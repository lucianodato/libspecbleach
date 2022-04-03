# libspecbleach

C library for audio noise reduction

## Background

This library is based on the extraction of the noise-repellent algorithm into a standalone library. It was re-written to be extensible and in a modular fashion. It uses the concept of a spectral processor which uses an short time Fourier transform (STFT) to process the audio. There are two initial processors in place, one which uses the adaptive part of noise repellent and one that uses the manual capturing profile based denoising. The library could be extended with more spectral processors using any STFT based algorithm such as de-crackle, de-click and other audio restoration algorithms.

## De-noise algorithms

There several techniques implemented in the library that are being used in the denoisers, such as masking thresholds estimation, onset detectors, etc. All these are being used in conjuntion to improve the very basic spectral substraction algorithm. Most of the papers used are listed in the wiki of the project. Also a block diagram is provided to explain the reduction architecture.
