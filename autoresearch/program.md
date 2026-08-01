# Autonomous Research Objectives: libspecbleach Structural DSP Innovations

You are an expert C Audio DSP researcher. Your goal is to implement **structural, algorithmic improvements** to `libspecbleach`'s noise-reduction engine.

---

## 🚫 STRICTLY FORBIDDEN (Immediate Rejection)
Do NOT propose changes that merely adjust existing constants, magic numbers, or user-facing parameter bounds. Examples of forbidden proposals:
- Modifying threshold scaling factors or gain multipliers (e.g., changing `0.5` to `0.6`).
- Changing fixed frame/FFT size constants in headers without altering the processing topology.
- Adjusting simple attack/release time constants in smoothing filters.
- Renaming variables or refactoring comments without code logic changes.

---

## 🎯 PERMITTED & ENCOURAGED STRUCTURAL INNOVATIONS

You are expected to invent, implement, and integrate structural algorithms into the C codebase:

### 1. Multi-Resolution STFT (MR-STFT)
- **Goal:** Single-window STFT suffers from the time-frequency uncertainty principle.
- **Task:** Implement or combine parallel multi-resolution FFT analysis/synthesis blocks (e.g., short 256-sample window for sharp transient resolution + long 2048/4096-sample window for fine low-frequency harmonic resolution). Blend their attenuation masks dynamically based on spectral energy.

### 2. Transient-Aware Phase Preservation & Phase-Locking
- **Goal:** Standard spectral subtraction destroys phase alignment, causing musical noise and transient dullness.
- **Task:** Implement phase-locking or principal component phase preservation across adjacent STFT bins during mask application so harmonic stacks and drum attacks maintain coherence.

### 3. Advanced Psychoacoustic Masking Topologies
- **Goal:** Improve upon basic Johnston masking in `masking.c`.
- **Task:** Introduce Bark or Gammatone filterbank frequency grouping, forward/backward temporal masking curves (psychoacoustic post-masking across time frames), or spreading functions that dynamically protect auditory-masked bins from aggressive gating.

### 4. Non-Local Means (NLM) 2D Filtering Improvements (`nlm.c`)
- **Goal:** Reduce 2D time-frequency burbling artifacts.
- **Task:** Replace or improve Euclidean patch distance metrics with weighted perceptual spectral distances, or introduce adaptive search-window shapes that expand along tonal trajectories.

---

## 🚫 STRICTLY FORBIDDEN: Heavy Recursive Frame Feedback & Temporal Echoes

Do NOT implement heavy frame-to-frame recursive feedback or high-order IIR temporal smoothing (e.g., Ephraim-Malah Decision-Directed SNR tracking with high alpha, multi-frame exponential smoothing, or spectral delay loops). 

**Reason:** Recursive temporal smoothing creates temporal ghosting, echoes, and trailing artifacts ("smeared tails") after transients and word endings.

---

## 🎨 PREFERRED ALTERNATIVE DIRECTIONS: Image Processing & Modern ML Masking

Instead of relying on past audio frames, treat the STFT magnitude spectrogram as a **2D Matrix/Image** $S[t, f]$ and apply localized, non-recursive spatial/spectral image processing and zero-latency ML masking concepts:

### 1. Image Processing Techniques (Applied to 2D Spectrograms)
* **2D Guided Filtering / Bilateral Filtering:** Smooth stationary background noise while strictly preserving sharp spectro-temporal edges (transients and harmonics) without introducing temporal lag or echoes.
* **Anisotropic Diffusion (Perona-Malik Denoising):** Process the spectrogram as a 2D intensity image. Diffuse energy within uniform time-frequency regions (noise floor) while stopping diffusion across high-gradient spectral edges.
* **Non-Local Means (NLM) Patch Filtering (`nlm.c`):** Use 2D patch-matching in local time-frequency neighborhoods rather than long-term time-series history.
* **Morphological 2D Operations:** Apply 2D erosion/dilation or top-hat filtering on the gain matrix to remove isolated "musical noise" points without temporal blurring.

### 2. ML-Inspired & Instantaneous Masking
* **Instantaneous / Finite Window Masking:** Compute spectral masks using only the current frame (plus maximum 1 adjacent frame for lookahead/lookback FIR windowing).
* **Psychoacoustic Simultaneous Masking:** Protect weak harmonics based on instantaneous frequency masking thresholds (e.g., Johnston or ISO 5218 masking models) within the *current frame*, rather than relying on past frame history.
* **Compact FIR-like 2D Kernels:** Replace recursive state accumulators with small, fixed $3 \times 3$ or $5 \times 5$ 2D spatial convolution kernels over the time-frequency plane.

---

## Code Modification Rules
- You may rewrite internal helper functions, introduce new structural loops, or modify matrix allocations in `src/`.
- Ensure all memory allocated for new multi-resolution buffers or state structures is properly freed to prevent memory leaks or segfaults.
- The modified C code must compile cleanly using `cmake`.

## 🏗️ ARCHITECTURAL OVERHAUL MANDATE (Structural Refactoring)

You are **fully authorized and encouraged** to make radical structural changes to the processing pipeline. You are NOT restricted to editing existing functions—you may replace, merge, bypass, or completely delete modules if a superior alternative yields better evaluation scores.

---

### Permitted Architectural Actions

1. **Module Replacement & Consolidation**
   - If an existing stage (e.g., `masking_veto.c`, `gain_calculator.c`, or `tonal_reducer.c`) is causing phase artifacts or limiting performance, you may **replace it entirely** with a cleaner, modern algorithm (e.g., a unified decision-directed Joint Prior Estimator or a multi-band Wiener filter).
   - You may consolidate multiple processing sub-modules into a single, cohesive processing block if it reduces pipeline latency or cumulative phase distortion.

2. **Pipeline Topology Swaps & Pruning**
   - You may change the order in which modules execute in `denoiser_logic/`.
   - You may **bypass or eliminate** processing stages entirely (e.g., removing a post-processing whitening stage if raw spectral subtraction + psychoacoustic masking produces cleaner audio without it).
   - You may swap out the underlying estimation topology (e.g., replacing static minimum statistics with dynamic SPP-MMSE or a multi-resolution filterbank approach).

3. **Memory & Struct Refactoring**
   - You may introduce new C structs, add state buffers, or refactor state handles in `include/` and `src/` to support new multi-frame or multi-scale architectures.

---

### ⚠️ Strict Architectural Safeguards

While you can overhaul internal modules, you **MUST** adhere to the following build contracts:

1. **Top-Level API Compatibility:** 
   - Do NOT break the public initialization and processing function signatures called by `examples/denoiser_demo.c`. The example binary must compile and execute cleanly via `cmake`.
2. **Strict C Memory Hygiene:** 
   - If you create new internal structures or replace existing ones, ensure `init` and `free`/`destroy` functions properly allocate and clean up memory to prevent pointer crashes or leaks.
3. **Compilation Integrity:** 
   - All changes across modified files must build without CMake compilation errors or missing symbol warnings.

## ⚡ PERFORMANCE & REAL-TIME OPTIMIZATION OBJECTIVES

Your proposed algorithms are evaluated on **Audio Quality AND Computational Performance**.

---

### Performance Targets
1. **CPU Overhead / Real-Time Factor (RTF):**
   - Target: $\text{RTF} \le 0.05$ (Processing must take under 5% of real-time audio playback duration).
   - *Guidance:* Utilize OpenMP vectorization, efficient FFT windowing, avoid redundant matrix allocations inside frame loops, and streamline STFT hop iterations.

2. **Algorithmic Latency:**
   - Target: $\le 10 \text{ ms}$ buffer delay.
   - *Guidance:* Keep STFT hop sizes and lookahead buffers minimal. Avoid unnecessarily huge single-window STFT frame sizes (e.g., preference multi-resolution STFT with short transient windows or partitioned FFT processing over massive 8192-sample windows).

---

### Balanced Trade-off Metric
The benchmark calculates a composite score:
$$\text{Final Score} = \text{Quality Score} - (10.0 \times \text{RTF}) - (0.05 \times \text{Latency}_{\text{ms}})$$

An edit that achieves slightly cleaner audio but doubles CPU time or adds 50ms of latency will **BE REJECTED**. Always aim for lean, cache-friendly C code.

### Default configurations for the library
When adding new DSP features or modifying fallback struct defaults, ensure default struct initializers in configurations.h or config constructors reflect your optimized baseline values.