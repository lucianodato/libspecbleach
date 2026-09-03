# libspecbleach Examples

Executable integration references. Read the one matching your use case,
then steal its structure — they are written to be copied.

## Which example should I read?

| Your situation | Read | Uses |
| :--- | :--- | :--- |
| Embedding the library anywhere, no file I/O | [`simple_embed.c`](simple_embed.c) | Core API, zero dependencies |
| Mono processing, offline or simple pipeline | [`denoiser_demo.c`](denoiser_demo.c) | Core API only (`specbleach_denoiser.h`) |
| Stereo/surround, smoothing-mode switching, "what the plugin does" | [`stereo_denoiser_demo.c`](stereo_denoiser_demo.c) | Extras layer (`specbleach_stereo.h`) |
| Real-time application with GUI-driven learn/switch | `stereo_denoiser_demo.c` **plus** the [Noise Repellent plugin source](https://github.com/lucianodato/noise-repellent/blob/master/Source/PluginProcessor.cpp) | Extras + threading patterns |

Both file demos process audio via libsndfile so they are runnable and
verifiable; everything except file I/O maps 1:1 to real-time callback code,
and the demos call out exactly where the real-time differences are.
`simple_embed.c` has no external dependencies at all.

## The 60-second version

```c
#include <math.h>
#include <specbleach_denoiser.h>

// 1. CREATE — one instance per channel; frame_size_ms = STFT window (20-100)
specbleach_denoiser* denoiser =
    specbleach_denoiser_initialize(sample_rate, 46.0f);

// 2. CONFIGURE — start from documented-safe defaults, override what you need.
//    Never "= {0}": zero-initialized reduction_gain means MAXIMUM reduction.
//    The library copies everything on load, including curve bias arrays.
SpecbleachDenoiserParameters p = specbleach_denoiser_get_default_parameters();
p.learn_noise = SPECBLEACH_LEARN_ALL;
p.reduction_gain = powf(10.0f, -20.0f / 20.0f);   // -20 dB gain floor
specbleach_denoiser_load_parameters(denoiser, &p, sizeof(p));

// 3. LEARN — feed blocks containing only noise; each channel learns its own
specbleach_denoiser_process(denoiser, nframes, in, out);

// 4. FINALIZE — profiles are only usable after learning turns OFF
p.learn_noise = SPECBLEACH_LEARN_OFF;
specbleach_denoiser_load_parameters(denoiser, &p, sizeof(p));

// 5. REDUCE — same call shape forever after; any block size
specbleach_denoiser_process(denoiser, nframes, in, out);

// 6. CLEANUP
specbleach_denoiser_free(denoiser);
```

For multi-channel, wrap steps in a group instead:
`specbleach_stereo_initialize(sr, ms, channels)` then load/process
once for all channels (deinterleaved pointer arrays).

## Building the examples

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release \
      -DENABLE_EXAMPLES=ON -DSPECBLEACH_BUILD_EXTRAS=ON
cmake --build build --config Release

./build/stereo_denoiser_demo --switch-smoothing noisy.wav clean.wav
```

Linking from your own project (installed package):

```cmake
find_package(libspecbleach 0.4 REQUIRED)                   # core only
find_package(libspecbleach 0.4 REQUIRED COMPONENTS extras) # + stereo layer
target_link_libraries(myapp PRIVATE libspecbleach::libspecbleach)
target_link_libraries(myapp PRIVATE libspecbleach::extras)
```

## Switching smoothing modes

There is only ONE denoiser family. The `smoothing_mode` field of
`SpecbleachDenoiserParameters` selects the temporal (1D) or NLM 2D smoothing strategy, and the library owns
the runtime transition: it crossfades internally over
`SMOOTHING_TRANSITION_SECONDS` (30 ms) with zero allocations on the audio
thread. Both modes report the same constant latency, so hosts never need to
re-anchor delay compensation. Switching is just a parameter reload.

C++ integrators: `<specbleach.hpp>` provides header-only RAII ownership
for the core handle (`make_denoiser`); `<specbleach_stereo.hpp>` adds the
extras group (`make_stereo_group`).
It wraps lifetime management only and has no framework dependencies, so it
suits JUCE, raw VST3/CLAP/LV2, DAW codebases, and standalone apps alike.

## Pitfalls these examples exist to prevent

1. **Profiles are not ready until learning turns OFF.** Capture modes are
   finalized on the learn → off transition, not while learning — and only
   fully usable after at least one block is processed afterwards. Both demos
   reload parameters with `SPECBLEACH_LEARN_OFF` before reducing.
2. **Each channel learns independently.** Do not average stereo inputs into
   one engine; use `specbleach_stereo`, which keeps per-channel profiles and
   can fallback-fill gaps with `specbleach_stereo_sync_profiles()`.
3. **Latency reports belong on a message thread.** Hosts compensate delay
   based on what you tell them, and several hosts restart plugin processing
   when the reported value changes synchronously mid-callback. (Smoothing
   mode switches never change latency, so this only applies at prepare time.)
4. **Keep control work off the audio thread.** `process()` is real-time
   safe; parameter loads and profile sync/serialization are not.
5. **Wrong-size parameter loads fail cleanly by design.** Always pass
   `sizeof(the_exact_struct)` (or `SPECBLEACH_PARAMETERS_SIZE`); this protects
   you across library upgrades. Prefer `specbleach_denoiser_get_default_parameters()`
   over `= {0}` for the same reason.
6. **Mixing handle types does not compile** — that is intentional. Each
   handle type has its own opaque type; if it compiles, it is correct.

## Further reference

- Public API contracts: `include/*.h` — every function documents threading
  and ownership behavior.
- Module-level behavior: `tests/test_*.c` double as executable
  specification; `test_specbleach_stereo.c` shows edge cases (wrong sizes,
  NULL handling).
- Real-time orchestration reference:
  [noise-repellent `PluginProcessor.cpp`](https://github.com/lucianodato/noise-repellent/blob/master/Source/PluginProcessor.cpp).
