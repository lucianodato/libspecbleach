# libspecbleach Examples

Executable integration references. Read the one matching your use case,
then steal its structure — they are written to be copied.

## Which example should I read?

| Your situation | Read | Uses |
| :--- | :--- | :--- |
| Mono processing, offline or simple pipeline | [`denoiser_demo.c`](denoiser_demo.c) | Core API only (`specbleach_denoiser.h`) |
| Stereo/surround, engine switching, "what the plugin does" | [`stereo_denoiser_demo.c`](stereo_denoiser_demo.c) | Extras layer (`specbleach_stereo.h`, `specbleach_transition.h`, `specbleach_profile_migration.h`) |
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

## Pitfalls these examples exist to prevent

- **Switching engines with different algorithmic latencies**: no amount of
  in-plugin blending hides the moment the host re-anchors delay
  compensation. The proven pattern is structural: run the shorter-latency
  family through a permanent delay ring, report max(latency) constantly,
  and blend aligned streams. See the KNOWN LIMITATION note in
  `extras/include/specbleach_transition.h`.

1. **Profiles are not ready until learning turns OFF.** Capture modes are
   finalized on the learn → off transition, not while learning. Both demos
   reload parameters with `SPECBLEACH_LEARN_OFF` before reducing.
2. **Each channel learns independently.** Do not average stereo inputs into
   one engine; use `specbleach_stereo`, which keeps per-channel profiles and
   can fallback-fill gaps with `specbleach_stereo_sync_profiles()`.
3. **Engine switches need profile migration**, otherwise the new family
   starts deaf. `specbleach_stereo_migrate_profiles_from()` before
   `specbleach_transition_begin()`.
4. **Report latency when you call `transition_begin()`.** Hosts compensate
   delay based on what you tell them; the transition module aligns audio
   internally but cannot inform your host for you.
5. **Keep control work off the audio thread.** `process()` is real-time
   safe; `transition_begin()`, profile migration/sync, and serialization are
   not. The Noise Repellent plugin shows the pause-gate handshake pattern.
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
