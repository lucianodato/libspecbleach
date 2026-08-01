import os
import glob
import sys
import time
import argparse
import subprocess
import numpy as np
import soundfile as sf
import fast_bss_eval
import librosa
import scipy.signal

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR) if os.path.basename(SCRIPT_DIR) == "autoresearch" else SCRIPT_DIR
os.chdir(REPO_ROOT)

BENCHMARK_DIR = "autoresearch/research_data/benchmark_suite"
TIMEOUT_SEC = 15
NOISE_PREFIX_SEC = 2.5

# Standard Presets
PRESETS_ALL = [
    {"name": "light",    "args": ["--learn-frames", "200", "--reduction", "12.0", "--whitening", "10.0", "--smoothing", "0.05" ]},
    {"name": "moderate", "args": ["--learn-frames", "200", "--reduction", "18.0", "--whitening", "20.0", "--smoothing", "0.15" ]},
    {"name": "heavy",    "args": ["--learn-frames", "200", "--reduction", "24.0", "--whitening", "30.0", "--smoothing", "0.25" ]}
]

PRESET_FAST = [
    {"name": "moderate", "args": ["--learn-frames", "200", "--reduction", "18.0", "--whitening", "20.0", "--smoothing", "0.15" ]}
]

def compile_libspecbleach() -> bool:
    if not os.path.exists("build"):
        setup_cmd = [
            "cmake", "-B", "build",
            "-DCMAKE_BUILD_TYPE=Release",
            "-DENABLE_EXAMPLES=ON"
        ]
        if subprocess.run(setup_cmd, capture_output=True).returncode != 0:
            return False

    compile_cmd = ["cmake", "--build", "build", "--config", "Release", "-j4"]
    return subprocess.run(compile_cmd, capture_output=True).returncode == 0

def find_denoiser_executable() -> str:
    possible_paths = [
        "./build/examples/denoiser_demo",
        "./build/example/denoiser_demo",
        "./build/denoiser_demo",
        "./build/Release/denoiser_demo",
        "./build/examples/Release/denoiser_demo"
    ]
    for path in possible_paths:
        if os.path.exists(path):
            return path
    return ""

def compute_temporal_echo_penalty(s_clean: np.ndarray, s_proc: np.ndarray, sr: int = 48000) -> float:
    """
    Measures post-transient energy leakage (echoes/ghosting).
    Returns a penalty score >= 0.0 (higher = worse echo/smearing).
    """
    # 1. Compute frame energy envelopes
    hop = 256
    clean_env = librosa.feature.rms(y=s_clean, frame_length=512, hop_length=hop)[0]
    proc_env  = librosa.feature.rms(y=s_proc,  frame_length=512, hop_length=hop)[0]
    
    # 2. Find frames where clean energy drops rapidly (>15dB decay in 2-3 frames)
    clean_db = 20 * np.log10(clean_env + 1e-6)
    proc_db  = 20 * np.log10(proc_env + 1e-6)
    
    diff_clean = np.diff(clean_db)
    
    # Identify sharp drop-offs (phrase ends / post-transients)
    drop_indices = np.where(diff_clean < -15.0)[0]
    
    if len(drop_indices) == 0:
        return 0.0

    # 3. Measure if processed audio decays much slower than clean audio
    echo_leakage = 0.0
    for idx in drop_indices:
        if idx + 2 < len(proc_db):
            clean_decay = clean_db[idx] - clean_db[idx + 2]
            proc_decay  = proc_db[idx] - proc_db[idx + 2]
            
            # If clean decayed 20dB but processed only decayed 5dB -> Echo detected!
            decay_error = max(0.0, clean_decay - proc_decay)
            echo_leakage += decay_error

    mean_echo_leakage = echo_leakage / len(drop_indices)
    return float(mean_echo_leakage)

def compute_lsd(s_clean: np.ndarray, s_proc: np.ndarray, sr: int = 48000) -> float:
    S_clean = np.abs(librosa.stft(s_clean, n_fft=2048, hop_length=512))
    S_proc  = np.abs(librosa.stft(s_proc,  n_fft=2048, hop_length=512))
    log_clean = 20 * np.log10(S_clean + 1e-7)
    log_proc  = 20 * np.log10(S_proc  + 1e-7)
    return float(np.mean(np.sqrt(np.mean((log_clean - log_proc) ** 2, axis=0))))

def compute_mr_stft_loss(s_clean: np.ndarray, s_proc: np.ndarray) -> float:
    losses = []
    for n_fft in [512, 1024, 2048]:
        S_clean = np.abs(librosa.stft(s_clean, n_fft=n_fft, hop_length=n_fft // 4))
        S_proc  = np.abs(librosa.stft(s_proc,  n_fft=n_fft, hop_length=n_fft // 4))
        losses.append(np.linalg.norm(S_clean - S_proc, ord='fro') / (np.linalg.norm(S_clean, ord='fro') + 1e-7))
    return float(np.mean(losses))

def measure_empirical_latency_fast(mix_signal: np.ndarray, proc_signal: np.ndarray, sr: int = 48000) -> float:
    """Optimized fast cross-correlation using a 0.25-second window."""
    n_samples = min(len(mix_signal), len(proc_signal), int(sr * 0.25))
    mix_seg  = mix_signal[:n_samples]
    proc_seg = proc_signal[:n_samples]
    
    correlation = scipy.signal.correlate(proc_seg, mix_seg, mode='full')
    lags = scipy.signal.correlation_lags(len(proc_seg), len(mix_seg), mode='full')
    
    delay_samples = lags[np.argmax(correlation)]
    return max(0.0, (delay_samples / sr) * 1000.0)

def process_case(case_dir: str, presets: list, skip_perf: bool = False, fast_metrics: bool = False) -> dict:
    mix_path   = os.path.join(case_dir, "mix.wav")
    clean_path = os.path.join(case_dir, "clean.wav")
    noise_path = os.path.join(case_dir, "noise.wav")

    executable = find_denoiser_executable()
    if not executable:
        return {"score": -999.0}

    s_clean, sr = sf.read(clean_path)
    s_noise, _  = sf.read(noise_path)
    s_mix, _    = sf.read(mix_path)

    scores, qualities, rtfs, latencies = [], [], [], []

    for preset in presets:
        out_path = os.path.join(case_dir, f"processed_{preset['name']}.wav")
        cmd = [executable] + preset["args"] + [mix_path, out_path]
        
        try:
            t0 = time.perf_counter()
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT_SEC)
            exec_time_sec = time.perf_counter() - t0

            if res.returncode != 0 or not os.path.exists(out_path):
                return {"score": -999.0}
        except subprocess.TimeoutExpired:
            return {"score": -999.0}

        s_proc, _ = sf.read(out_path)

        # Performance evaluation
        if not skip_perf:
            rtf = exec_time_sec / (len(s_mix) / sr)
            latency_ms = measure_empirical_latency_fast(s_mix, s_proc, sr)
        else:
            rtf = 0.0
            latency_ms = 0.0

        # Quality evaluation window
        start_sample = int(NOISE_PREFIX_SEC * sr)
        s_clean_m = s_clean[start_sample:]
        s_noise_m = s_noise[start_sample:]
        s_proc_m  = s_proc[start_sample:]

        min_len = min(len(s_clean_m), len(s_noise_m), len(s_proc_m))
        if min_len <= 0:
            return {"score": -999.0}

        s_clean_m = s_clean_m[:min_len]
        s_noise_m = s_noise_m[:min_len]
        s_proc_m  = s_proc_m[:min_len]

        # Fast BSS Eval
        ref = np.vstack([s_clean_m, s_noise_m])
        est_noise = (s_clean_m + s_noise_m) - s_proc_m
        est = np.vstack([s_proc_m, est_noise])

        sdr, sir, sar, _ = fast_bss_eval.bss_eval_sources(ref, est)
        sar_val, sir_val = float(sar[0]), float(sir[0])

        clean_onset = librosa.onset.onset_strength(y=s_clean_m, sr=sr)
        proc_onset  = librosa.onset.onset_strength(y=s_proc_m, sr=sr)
        onset_corr  = float(np.corrcoef(clean_onset, proc_onset)[0, 1])

        # Skip heavy STFT metrics in extreme fast mode
        if not fast_metrics:
            lsd_val = compute_lsd(s_clean_m, s_proc_m, sr)
            mr_stft = compute_mr_stft_loss(s_clean_m, s_proc_m)
        else:
            lsd_val = 0.0
            mr_stft = 0.0

        quality_score = (
            (0.35 * sar_val) +
            (0.30 * sir_val) +
            (0.15 * (onset_corr * 10.0)) -
            (0.10 * lsd_val) -
            (0.10 * (mr_stft * 10.0))
        )

        echo_penalty = compute_temporal_echo_penalty(s_clean_m, s_proc_m, sr)

        # Subtract echo penalty from the quality score
        quality_score = quality_score - (0.05 * echo_penalty)

        cpu_penalty = 10.0 * rtf
        latency_penalty = 0.05 * max(0.0, latency_ms - 10.0)

        preset_score = quality_score - cpu_penalty - latency_penalty

        scores.append(preset_score)
        qualities.append(quality_score)
        rtfs.append(rtf)
        latencies.append(latency_ms)

    return {
        "score": float(np.mean(scores)),
        "quality": float(np.mean(qualities)),
        "rtf": float(np.mean(rtfs)),
        "latency_ms": float(np.mean(latencies))
    }

def evaluate_suite(fast_mode: bool = False, max_cases: int = 0, skip_perf: bool = False) -> float:
    if not compile_libspecbleach():
        print("[FAIL] CMake build failed")
        return -999.0

    case_dirs = sorted(glob.glob(f"{BENCHMARK_DIR}/case_*"))
    if not case_dirs:
        print("[FAIL] Benchmark suite empty!")
        return -999.0

    # Apply fast mode defaults
    presets = PRESET_FAST if fast_mode else PRESETS_ALL
    if max_cases > 0:
        case_dirs = case_dirs[:max_cases]

    scores, rtfs, latencies = [], [], []
    for c_dir in case_dirs:
        res = process_case(c_dir, presets=presets, skip_perf=skip_perf, fast_metrics=fast_mode)
        if res["score"] == -999.0:
            print(f"[FAIL] Execution error in case: {os.path.basename(c_dir)}")
            return -999.0
        scores.append(res["score"])
        rtfs.append(res["rtf"])
        latencies.append(res["latency_ms"])

    mean_score = float(np.mean(scores))
    mean_rtf = float(np.mean(rtfs))
    mean_lat = float(np.mean(latencies))

    print(f"[SUCCESS] Evaluated {len(case_dirs)} cases ({'FAST' if fast_mode else 'FULL'} mode).")
    if not skip_perf:
        print(f" -> Avg RTF: {mean_rtf:.4f} | Avg Latency: {mean_lat:.2f} ms")
    print(f" -> Final Score: {mean_score:.4f}")
    return mean_score

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate libspecbleach performance and quality.")
    parser.add_argument("--fast", action="store_true", help="Run 1 preset, fast cross-correlation, and skip heavy STFTs.")
    parser.add_argument("--max-cases", type=int, default=0, help="Cap evaluation to first N test cases (e.g., 6).")
    parser.add_argument("--skip-perf", action="store_true", help="Disable CPU RTF and Latency penalties completely.")
    
    args = parser.parse_args()
    
    score = evaluate_suite(fast_mode=args.fast, max_cases=args.max_cases, skip_perf=args.skip_perf)
    print(f"METRIC_SCORE: {score}")
