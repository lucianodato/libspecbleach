import os
import glob
import numpy as np
import soundfile as sf
import librosa

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR) if os.path.basename(SCRIPT_DIR) == "autoresearch" else SCRIPT_DIR
os.chdir(REPO_ROOT)

CLEAN_DIR = "autoresearch/audio_dataset/clean"
NOISE_DIR = "autoresearch/audio_dataset/noise"
BENCHMARK_DIR = "autoresearch/research_data/benchmark_suite"

TARGET_SR = 48000
TARGET_SNRS = [6.0, 12.0, 18.0]
NOISE_PREFIX_SEC = 2.5  # Pure noise lead-in duration for profile capture
SUPPORTED_EXTS = ["*.wav", "*.mp3", "*.flac", "*.ogg", "*.m4a"]

def find_audio_files(directory: str) -> list[str]:
    files = []
    for ext in SUPPORTED_EXTS:
        files.extend(glob.glob(os.path.join(directory, ext)))
        files.extend(glob.glob(os.path.join(directory, ext.upper())))
    return sorted(files)

def calculate_rms(audio: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(audio)) + 1e-12))

def match_noise_length(s_noise: np.ndarray, target_length: int) -> np.ndarray:
    noise_length = len(s_noise)
    if noise_length < target_length:
        repeats = int(np.ceil(target_length / noise_length))
        s_noise = np.tile(s_noise, repeats)
    return s_noise[:target_length]

def mix_at_snr(clean: np.ndarray, noise: np.ndarray, snr_db: float) -> tuple[np.ndarray, np.ndarray]:
    clean_rms = calculate_rms(clean)
    noise_rms = calculate_rms(noise)
    target_noise_rms = clean_rms / (10 ** (snr_db / 20.0))
    scale_factor = target_noise_rms / noise_rms
    
    scaled_noise = noise * scale_factor
    mix = clean + scaled_noise
    
    max_amp = np.max(np.abs(mix))
    if max_amp > 0.99:
        mix /= max_amp
        clean /= max_amp
        scaled_noise /= max_amp
        
    return mix, scaled_noise

def prepare_benchmark_suite():
    os.makedirs(BENCHMARK_DIR, exist_ok=True)
    
    clean_files = find_audio_files(CLEAN_DIR)
    noise_files = find_audio_files(NOISE_DIR)
    
    if not clean_files or not noise_files:
        raise FileNotFoundError("Populate clean and noise folders first!")

    prefix_samples = int(NOISE_PREFIX_SEC * TARGET_SR)
    print(f"Generating test cases with {NOISE_PREFIX_SEC}s pure noise prefix for profile learning...\n")

    pair_id = 0
    for clean_path in clean_files:
        for noise_path in noise_files:
            snr = TARGET_SNRS[pair_id % len(TARGET_SNRS)]
            pair_id += 1
            
            clean_name = os.path.splitext(os.path.basename(clean_path))[0]
            noise_name = os.path.splitext(os.path.basename(noise_path))[0]
            
            s_clean_body, _ = librosa.load(clean_path, sr=TARGET_SR, mono=True)
            s_noise_raw, _  = librosa.load(noise_path, sr=TARGET_SR, mono=True)
            
            # Prepare clean track with initial silence
            silence_prefix = np.zeros(prefix_samples, dtype=np.float32)
            s_clean_full = np.concatenate([silence_prefix, s_clean_body])
            
            # Prepare noise track to match full length
            s_noise_full = match_noise_length(s_noise_raw, target_length=len(s_clean_full))
            
            # Mix music body at target SNR
            mix_body, scaled_noise_body = mix_at_snr(s_clean_body, s_noise_full[prefix_samples:], snr)
            
            # Combine prefix noise + mixed music
            noise_prefix = s_noise_full[:prefix_samples]
            
            full_mix = np.concatenate([noise_prefix, mix_body])
            full_clean = s_clean_full
            full_noise = np.concatenate([noise_prefix, scaled_noise_body])
            
            case_dir = os.path.join(BENCHMARK_DIR, f"case_{pair_id:02d}_{clean_name}_vs_{noise_name}")
            os.makedirs(case_dir, exist_ok=True)
            
            sf.write(os.path.join(case_dir, "clean.wav"), full_clean, TARGET_SR)
            sf.write(os.path.join(case_dir, "noise.wav"), full_noise, TARGET_SR)
            sf.write(os.path.join(case_dir, "mix.wav"), full_mix, TARGET_SR)
            
            print(f"Case {pair_id:02d}: [{clean_name}] + [{noise_name}] ({NOISE_PREFIX_SEC}s noise lead-in @ {snr}dB SNR)")

    print(f"\nGenerated {pair_id} test cases with noise profile prefixes!")

if __name__ == "__main__":
    prepare_benchmark_suite()
