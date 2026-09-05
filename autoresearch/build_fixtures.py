"""Generates the committed real-world fixtures for the quality-metrics test.

Decodes/resamples/downmixes sources with macOS-native `afconvert` (48 kHz,
mono, 16-bit), then trims/tiles with numpy and writes raw `clean.wav` (signal
body, no lead-in) and `noise.wav` (tiled, full length). The C test
(test_integration_realworld_quality.c) synthesizes the mix itself from these
two files, which lets it sweep SNRs without extra committed data.

Requires: soundfile (fixture writing) + afconvert (source decoding).
Usage: python3 autoresearch/build_fixtures.py
"""

import os
import subprocess
import tempfile

import numpy as np
import soundfile as sf

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
DATASET_DIR = os.path.join(REPO_ROOT, "autoresearch", "audio_dataset")
FIXTURE_DIR = os.path.join(REPO_ROOT, "tests", "test_data", "realworld")

TARGET_SR = 48000
NOISE_LEAD_IN_SEC = 1.5  # Must match NOISE_LEAD_IN_SEC in the C test
BODY_SEC = 8.0

# (case_name, clean_file, clean_offset_sec, noise_file, noise_offset_sec)
# SNR designation lives in the C test table.
CASES = [
    ("case_voice_fan", "578075__bigspider95__italian_voice.wav", 1.0,
     "freesound_community-noisy-bathroom-exhaust-fan-31078.mp3", 0.5),
    ("case_voice_city", "578075__bigspider95__italian_voice.wav", 10.0,
     "freesound_community-city-skyline-rooftops-noisy-london-19916.mp3", 2.0),
    ("case_guitar_electric",
     "320800__ajubamusic__acoustic-picking-guitar-100bpm-key-of-c.mp3", 2.0,
     "Electric noise.mp3", 0.0),
    ("case_drums_fridge",
     "710275__lucianodato__headphone-test-drums-2.wav", 5.0,
     "freesound_community-fridge-buzz-loop-39612.mp3", 0.0),
]


def decode_48k_mono(path: str) -> np.ndarray:
    with tempfile.TemporaryDirectory() as tmp:
        tmp_wav = os.path.join(tmp, "decoded.wav")
        subprocess.run(["afconvert", "-f", "WAVE", "-d", "LEI16@48000",
                        "-c", "1", os.path.abspath(path), tmp_wav],
                       check=True, capture_output=True)
        audio, sr = sf.read(tmp_wav, dtype="float32", always_2d=True)
    assert sr == TARGET_SR, f"{path}: expected {TARGET_SR} Hz, got {sr}"
    return audio.mean(axis=1)


def match_noise_length(s_noise: np.ndarray, target_length: int) -> np.ndarray:
    if len(s_noise) < target_length:
        repeats = int(np.ceil(target_length / len(s_noise)))
        s_noise = np.tile(s_noise, repeats)
    return s_noise[:target_length]


def build_case(case_name: str, clean_file: str, clean_offset: float,
               noise_file: str, noise_offset: float) -> None:
    s_clean = decode_48k_mono(os.path.join(DATASET_DIR, "clean", clean_file))
    s_noise = decode_48k_mono(os.path.join(DATASET_DIR, "noise", noise_file))

    start = int(clean_offset * TARGET_SR)
    end = start + int(BODY_SEC * TARGET_SR)
    assert end <= len(s_clean), f"{case_name}: clean source too short"
    s_clean = s_clean[start:end]

    n_start = int(noise_offset * TARGET_SR)
    assert n_start < len(s_noise), f"{case_name}: noise offset out of range"
    # Noise must cover lead-in + body (the C test mixes from these buffers)
    target = int((BODY_SEC + NOISE_LEAD_IN_SEC) * TARGET_SR)
    s_noise = match_noise_length(s_noise[n_start:], target)

    clean_rms = float(np.sqrt(np.mean(np.square(s_clean)) + 1e-12))
    noise_rms = float(np.sqrt(np.mean(np.square(s_noise)) + 1e-12))
    assert clean_rms > 1e-4, f"{case_name}: clean trim is (near) silent"

    case_dir = os.path.join(FIXTURE_DIR, case_name)
    os.makedirs(case_dir, exist_ok=True)
    sf.write(os.path.join(case_dir, "clean.wav"), s_clean, TARGET_SR,
             subtype="PCM_16")
    sf.write(os.path.join(case_dir, "noise.wav"), s_noise, TARGET_SR,
             subtype="PCM_16")
    print(f"{case_name}: clean {len(s_clean) / TARGET_SR:.1f}s "
          f"(rms {clean_rms:.4f}), noise {len(s_noise) / TARGET_SR:.1f}s "
          f"(rms {noise_rms:.4f})")


def main() -> None:
    for case_name, clean_file, clean_offset, noise_file, noise_offset in CASES:
        build_case(case_name, clean_file, clean_offset, noise_file,
                   noise_offset)
    print(f"\nFixtures written to {FIXTURE_DIR}")


if __name__ == "__main__":
    main()
