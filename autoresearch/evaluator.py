import os
import glob
import sys
import subprocess
import numpy as np
import soundfile as sf
import fast_bss_eval
import librosa

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR) if os.path.basename(SCRIPT_DIR) == "autoresearch" else SCRIPT_DIR
os.chdir(REPO_ROOT)

BENCHMARK_DIR = "autoresearch/research_data/benchmark_suite"
TIMEOUT_SEC = 15
NOISE_PREFIX_SEC = 2.5  # Must match build_dataset.py

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

def compute_lsd(s_clean: np.ndarray, s_proc: np.ndarray, sr: int = 48000) -> float:
    """Computes Log-Spectral Distance (lower score = higher fidelity)."""
    S_clean = np.abs(librosa.stft(s_clean, n_fft=2048, hop_length=512))
    S_proc  = np.abs(librosa.stft(s_proc,  n_fft=2048, hop_length=512))
    
    # Avoid log(0) with small epsilon
    log_clean = 20 * np.log10(S_clean + 1e-7)
    log_proc  = 20 * np.log10(S_proc  + 1e-7)
    
    # Frequency-wise RMS distance averaged over time frames
    frame_lsd = np.sqrt(np.mean((log_clean - log_proc) ** 2, axis=0))
    return float(np.mean(frame_lsd))

def compute_mr_stft_loss(s_clean: np.ndarray, s_proc: np.ndarray) -> float:
    """Computes Multi-Resolution STFT Spectral Convergence Error across 3 window sizes."""
    fft_sizes = [512, 1024, 2048]
    losses = []
    
    for n_fft in fft_sizes:
        hop_length = n_fft // 4
        S_clean = np.abs(librosa.stft(s_clean, n_fft=n_fft, hop_length=hop_length))
        S_proc  = np.abs(librosa.stft(s_proc,  n_fft=n_fft, hop_length=hop_length))
        
        # Spectral Convergence Loss
        convergence = np.linalg.norm(S_clean - S_proc, ord='fro') / (np.linalg.norm(S_clean, ord='fro') + 1e-7)
        losses.append(convergence)
        
    return float(np.mean(losses))

def process_case(case_dir: str) -> float:
    mix_path   = os.path.join(case_dir, "mix.wav")
    clean_path = os.path.join(case_dir, "clean.wav")
    noise_path = os.path.join(case_dir, "noise.wav")
    out_path   = os.path.join(case_dir, "processed_out.wav")

    executable = find_denoiser_executable()
    if not executable:
        return -999.0

    # Execute denoiser using --learn-frames 200 (~2.0-2.5s of STFT frames)
    cmd = [
        executable,
        "--learn-frames", "200",
        "--reduction", "18.0",
        mix_path,
        out_path
    ]
    
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT_SEC)
        if res.returncode != 0 or not os.path.exists(out_path):
            return -999.0
    except subprocess.TimeoutExpired:
        return -999.0

    s_clean, sr = sf.read(clean_path)
    s_noise, _  = sf.read(noise_path)
    s_proc, _   = sf.read(out_path)

    # Cut off initial 2.5s noise prefix so metrics only score music quality
    start_sample = int(NOISE_PREFIX_SEC * sr)
    s_clean = s_clean[start_sample:]
    s_noise = s_noise[start_sample:]
    s_proc  = s_proc[start_sample:]

    min_len = min(len(s_clean), len(s_noise), len(s_proc))
    if min_len <= 0:
        return -999.0

    s_clean, s_noise, s_proc = s_clean[:min_len], s_noise[:min_len], s_proc[:min_len]

    # Calculate BSSEval metrics on the music section
    ref = np.vstack([s_clean, s_noise])
    est_noise = (s_clean + s_noise) - s_proc
    est = np.vstack([s_proc, est_noise])

    sdr, sir, sar, _ = fast_bss_eval.bss_eval_sources(ref, est)
    sar_val, sir_val, sdr_val = float(sar[0]), float(sir[0]), float(sdr[0])

    # 1. BSSEval (SAR & SIR)
    ref = np.vstack([s_clean, s_noise])
    est_noise = (s_clean + s_noise) - s_proc
    est = np.vstack([s_proc, est_noise])

    sdr, sir, sar, _ = fast_bss_eval.bss_eval_sources(ref, est)
    sar_val, sir_val = float(sar[0]), float(sir[0])

    # 2. Transient Onset Correlation
    clean_onset = librosa.onset.onset_strength(y=s_clean, sr=sr)
    proc_onset  = librosa.onset.onset_strength(y=s_proc, sr=sr)
    onset_corr  = float(np.corrcoef(clean_onset, proc_onset)[0, 1])

    # 3. Log-Spectral Distance & Multi-Resolution STFT Loss
    lsd_val = compute_lsd(s_clean, s_proc, sr)
    mr_stft = compute_mr_stft_loss(s_clean, s_proc)

    # --- Balanced Composite Score Calculation ---
    # + SAR (35%): Freedom from musical noise / burbling
    # + SIR (30%): Deep background noise suppression
    # + Onset (15%): Drum / pick attack preservation
    # - LSD (10%): Overall spectral envelope deviation
    # - MR-STFT (10%): Multi-resolution time-frequency smearing
    case_score = (
        (0.35 * sar_val) +
        (0.30 * sir_val) +
        (0.15 * (onset_corr * 10.0)) -
        (0.10 * lsd_val) -
        (0.10 * (mr_stft * 10.0))
    )


    return case_score

def evaluate_suite() -> float:
    if not compile_libspecbleach():
        print("[FAIL] CMake build failed")
        return -999.0

    case_dirs = sorted(glob.glob(f"{BENCHMARK_DIR}/case_*"))
    if not case_dirs:
        print("[FAIL] Benchmark suite empty. Run build_dataset.py first!")
        return -999.0

    scores = []
    for c_dir in case_dirs:
        score = process_case(c_dir)
        if score == -999.0:
            print(f"[FAIL] Execution error in case: {os.path.basename(c_dir)}")
            return -999.0
        scores.append(score)

    mean_score = float(np.mean(scores))
    print(f"[SUCCESS] Evaluated {len(scores)} cases using profile learning. Average Score: {mean_score:.4f}")
    return mean_score

if __name__ == "__main__":
    score = evaluate_suite()
    print(f"METRIC_SCORE: {score}")
