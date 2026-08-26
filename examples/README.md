# libspecbleach Examples

Executable integration references. Read the one matching your use case,
then steal its structure — they are written to be copied.

## Which example should I read?

| Your situation | Read | Uses |
| :--- | :--- | :--- |
| Mono processing, offline or simple pipeline | [`denoiser_demo.c`](denoiser_demo.c) | Core API only (`specbleach_denoiser.h`) |
| Stereo/surround, engine switching, "what the plugin does" | [`stereo_denoiser_demo.c`](stereo_denoiser_demo.c) | Extras layer (`specbleach_stereo.h`, `specbleach_delay_line.h`, `specbleach_profile_migration.h`) |
| Real-time application with GUI-driven learn/switch | `stereo_denoiser_demo.c` **plus** the [Noise Repellent plugin source](https://github.com/lucianodato/noise-repellent/blob/master/Source/PluginProcessor.cpp) | Extras + threading patterns |

Both demos process audio files via libsndfile so they are runnable and
verifiable; everything except file I/O maps 1:1 to real-time callback code,
and the demos call out exactly where the real-time differences are.

## The 60-second version

```c
#include <specbleach_denoiser.h>

// 1. CREATE — one instance per channel; frame_size_ms = STFT window (20-100)
specbleach_denoiser* denoiser =
    specbleach_denoiser_initialize(sample_rate, 46.0f);

// 2. CONFIGURE — library copies everything, including curve bias arrays
SpecbleachDenoiserParameters p = {0};
p.learn_noise = SPECBLEACH_LEARN_ALL;
p.reduction_gain = powf(10.0f, -20.0f / 20.0f);   // -20 dB
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
`specbleach_stereo_initialize(sr, ms, channels, ENGINE)` then load/process
once for all channels (deinterleaved pointer arrays).

## Building the examples

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release \
      -DENABLE_EXAMPLES=ON -DSPECBLEACH_BUILD_EXTRAS=ON
cmake --build build --config Release

./build/stereo_denoiser_demo --switch-engine noisy.wav clean.wav
```

Linking from your own project:

```bash
pkg-config --cflags --libs specbleach          # core only
pkg-config --cflags --libs specbleach-extras   # core + extras layer
```

```cmake
find_package(libspecbleach 0.4 REQUIRED COMPONENTS extras)
target_link_libraries(myapp PRIVATE libspecbleach::extras)
```

## Switching between processors

The engine families expose different opaque handles and report different
algorithmic latencies (e.g., spectral vs. 2D NLM). Switching between them at
runtime is YOUR policy as the host application; extras deliberately provides
no switch widget because no in-library blend can hide the moment your host
re-anchors its delay compensation. The two recipes below are both proven in
production use.

**Recipe A — structural alignment (seamless, costs latency):**

1. Wrap the shorter-latency family in a permanent delay ring so BOTH
   families always emit at `max(latency_a, latency_b)` samples of delay.
2. Report that constant maximum to your host and never change it.
3. Switch by rendering both families and applying an equal-power crossfade
   (`w_old = cos(t*pi/2)`, `w_new = sin(t*pi/2)`, t = fade progress) over
   40-100 ms. Streams are time-aligned, so the blend is artifact-free.

Costs: every session runs at the slower family's latency even when only the
fast one is active.

**Recipe B — mute-and-warm (native latency, brief silence):** this is what
the Noise Repellent plugin ships.

1. Ramp output down (~100 ms), then mute completely.
2. While muted, render ONLY the target family silently for ~700 ms so its
   internal history buffers warm up.
3. Move your reported latency to the target's native value while silent —
   and deliver the notification from a message thread, never mid-callback
   (hosts like Reaper suspend/resume plugins on synchronous latency-change
   notifications).
4. Wait for old buffered tails to drain through the host before reporting a
   latency DECREASE; then ramp back up over ~100 ms.

Costs: ~1 second of muted audio per switch; rewards: native latencies and
single-engine CPU in steady state.

In both recipes run profile migration (`specbleach_stereo_migrate_profiles_from`)
before the target renders, or it starts deaf.

For Recipe A's alignment stage there is a ready-made building block:
`extras/specbleach_delay_line.h` is a policy-free multi-channel single-tap
delay line (real-time safe process, any block interleaving). Size it to
the latency difference, park the shorter-latency family in it, report
max(latency) once, and crossfade.

C++ integrators: `extras/specbleach.hpp` provides header-only RAII
ownership for every handle type (`make_denoiser`, `make_stereo_group`,
`make_delay_line`, ...). It wraps lifetime management only and has no
framework dependencies, so it suits JUCE, raw VST3/CLAP/LV2, DAW
codebases, and standalone apps alike.

## Pitfalls these examples exist to prevent

1. **Profiles are not ready until learning turns OFF.** Capture modes are
   finalized on the learn → off transition, not while learning. Both demos
   reload parameters with `SPECBLEACH_LEARN_OFF` before reducing.
2. **Each channel learns independently.** Do not average stereo inputs into
   one engine; use `specbleach_stereo`, which keeps per-channel profiles and
   can fallback-fill gaps with `specbleach_stereo_sync_profiles()`.
3. **Engine switches need profile migration**, otherwise the new family
   starts deaf. Call `specbleach_stereo_migrate_profiles_from()` before the
   target renders (see "Switching between processors" above).
4. **Latency reports belong on a message thread.** Hosts compensate delay
   based on what you tell them, and several hosts restart plugin processing
   when the reported value changes synchronously mid-callback.
5. **Keep control work off the audio thread.** `process()` is real-time
   safe; parameter loads, profile migration/sync, and serialization are
   not.
6. **Wrong-size parameter loads fail cleanly by design.** Always pass
   `sizeof(the_exact_struct)`; this protects you across library upgrades.
7. **Mixing handle types does not compile** — that is intentional. Each
   engine family has its own opaque type; if it compiles, it is correct.

## Further reference

- Public API contracts: `include/*.h` — every function documents threading
  and ownership behavior.
- Module-level behavior: `tests/test_*.c` double as executable
  specification; `test_specbleach_stereo.c` and `test_specbleach_transition.c`
  show edge cases (wrong sizes, mid-fade re-begin, NULL handling).
- Real-time orchestration reference:
  [noise-repellent `PluginProcessor.cpp`](https://github.com/lucianodato/noise-repellent/blob/master/Source/PluginProcessor.cpp).
