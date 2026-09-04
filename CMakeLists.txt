cmake_minimum_required(VERSION 3.22)
project(libspecbleach VERSION 0.4.0 LANGUAGES C)

# Guard against drift between the project version and the static
# include/specbleach_version.h macros (single source of truth: this file).
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/include/specbleach_version.h"
     SPECBLEACH_VERSION_HEADER)
string(REGEX MATCH "#define SPECBLEACH_VERSION_STRING \"([^\"]+)\""
       _ "${SPECBLEACH_VERSION_HEADER}")
set(SPECBLEACH_HEADER_VERSION "${CMAKE_MATCH_1}")
if(NOT SPECBLEACH_HEADER_VERSION VERSION_EQUAL PROJECT_VERSION)
    message(FATAL_ERROR
        "specbleach_version.h reports ${SPECBLEACH_HEADER_VERSION} but "
        "CMake project() is ${PROJECT_VERSION}; keep them in sync.")
endif()

# Optional C++ support (used only by the header-wrapper smoke test; the
# library itself is pure C and builds without a C++ compiler).
include(CheckLanguage)
check_language(CXX)
if(CMAKE_CXX_COMPILER)
    enable_language(CXX)
endif()

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ==============================================================================
# Build Options
# ==============================================================================
option(BUILD_SHARED_LIBS "Build shared libraries" ON)
option(ENABLE_AVX "Enable AVX optimizations on x86_64" ON)
option(ENABLE_TESTS "Build test suite" OFF)
option(ENABLE_EXAMPLES "Build example executables" OFF)
option(ENABLE_SANITIZERS "Enable Address and Undefined Behavior Sanitizers" OFF)
option(ENABLE_COVERAGE "Enable code coverage instrumentation" OFF)
# ==============================================================================
# 1. Vendored Third-Party Dependencies
# ==============================================================================
set(PFFFT_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/pffft/src/pffft.c"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/pffft/src/pffft_common.c"
)
set(PFFFT_INCLUDE_DIRS
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/pffft/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/pffft/src"
)
set(PKG_CONFIG_REQUIRES_PRIVATE "")

# ==============================================================================
# 2. Math & Threading Dependencies
# ==============================================================================
find_library(MATH_LIB m)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
message(STATUS "Threading: internal worker pool via ${Threads_LIBRARIES}")

# ==============================================================================
# 3. Optional Dependencies (Examples & Tests)
# ==============================================================================
if(ENABLE_EXAMPLES OR ENABLE_TESTS)
    find_package(PkgConfig QUIET)
    if(PkgConfig_FOUND)
        pkg_check_modules(SNDFILE IMPORTED_TARGET sndfile)
    endif()
endif()

# ==============================================================================
# 4. Compiler Flags & Sanitizers
# ==============================================================================
if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Wno-unused-parameter -Wno-unused-variable)
endif()

if(ENABLE_SANITIZERS AND NOT MSVC)
    add_compile_options(-fsanitize=address -fsanitize=undefined)
    add_link_options(-fsanitize=address -fsanitize=undefined)
endif()

if(ENABLE_COVERAGE AND NOT MSVC)
    add_compile_options(--coverage -O0)
    add_link_options(--coverage)
endif()

# ==============================================================================
# 5. Source Files Collection
# ==============================================================================
file(GLOB_RECURSE SPECLEACH_SHARED_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/shared/*.c"
)
list(FILTER SPECLEACH_SHARED_SOURCES EXCLUDE REGEX ".*nlm_filter_avx\\.c$")

file(GLOB_RECURSE SPECLEACH_PROCESSOR_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/processors/*.c"
)

# ==============================================================================
# 6. Target Definition
# ==============================================================================
add_library(libspecbleach
    ${SPECLEACH_SHARED_SOURCES}
    ${SPECLEACH_PROCESSOR_SOURCES}
    ${PFFFT_SOURCES}
)

set_target_properties(libspecbleach PROPERTIES
    OUTPUT_NAME "specbleach"
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    C_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN ON
)

# Symbol visibility contract (see include/specbleach_export.h):
# - Building the shared library exports SPECBLEACH_API symbols only.
# - Consuming a static archive must not use dllimport on Windows.
if(BUILD_SHARED_LIBS)
    target_compile_definitions(libspecbleach PRIVATE SPECBLEACH_EXPORTS)
else()
    target_compile_definitions(libspecbleach PUBLIC SPECBLEACH_STATIC)
endif()
target_compile_definitions(libspecbleach PRIVATE PFFFT_STATIC_DEFINE)

target_include_directories(libspecbleach PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_include_directories(libspecbleach PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(libspecbleach SYSTEM PRIVATE
    ${PFFFT_INCLUDE_DIRS}
)

if(MATH_LIB)
    target_link_libraries(libspecbleach PRIVATE ${MATH_LIB})
endif()

# Threading is an internal implementation detail -> PRIVATE
target_link_libraries(libspecbleach PRIVATE Threads::Threads)

# AVX SIMD Target
if(ENABLE_AVX)
    set(AVX_C_FILE "${CMAKE_CURRENT_SOURCE_DIR}/src/shared/denoiser_logic/processing/nlm_filter_avx.c")
    if(EXISTS "${AVX_C_FILE}" AND (CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64|amd64|AMD64|x86-64)" OR "x86_64" IN_LIST CMAKE_OSX_ARCHITECTURES))
        add_library(specbleach_avx STATIC "${AVX_C_FILE}")
        if(MSVC)
            target_compile_options(specbleach_avx PRIVATE /arch:AVX)
        else()
            target_compile_options(specbleach_avx PRIVATE -mavx)
        endif()
        target_include_directories(specbleach_avx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/src)
        target_link_libraries(libspecbleach PRIVATE specbleach_avx)
    endif()
endif()

# Alias Target for Subdirectories/FetchContent
add_library(libspecbleach::libspecbleach ALIAS libspecbleach)
add_library(specbleach ALIAS libspecbleach)

# ==============================================================================
# 6b. Optional Extras Layer (multi-channel orchestration)
# ==============================================================================
option(SPECBLEACH_BUILD_EXTRAS
       "Build optional orchestration helpers (stereo groups, profile migration)"
       OFF)

if(SPECBLEACH_BUILD_EXTRAS)
    add_library(libspecbleach_extras
        extras/src/specbleach_stereo.c
    )

    set_target_properties(libspecbleach_extras PROPERTIES
        OUTPUT_NAME "specbleach_extras"
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
        C_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN ON
    )

    target_include_directories(libspecbleach_extras PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/extras/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )

    target_link_libraries(libspecbleach_extras PUBLIC libspecbleach)

    if(BUILD_SHARED_LIBS)
        target_compile_definitions(libspecbleach_extras PRIVATE SPECBLEACH_EXPORTS)
    else()
        target_compile_definitions(libspecbleach_extras PUBLIC SPECBLEACH_STATIC)
    endif()

    add_library(libspecbleach::extras ALIAS libspecbleach_extras)
endif()

# ==============================================================================
# 7. Examples & Tests
# ==============================================================================
if(ENABLE_EXAMPLES)
    # Zero-dependency embed example: no libsndfile required
    add_executable(simple_embed "${CMAKE_CURRENT_SOURCE_DIR}/examples/simple_embed.c")
    target_include_directories(simple_embed PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(simple_embed PRIVATE libspecbleach)
    if(MATH_LIB)
        target_link_libraries(simple_embed PRIVATE ${MATH_LIB})
    endif()
endif()

if(ENABLE_EXAMPLES AND TARGET PkgConfig::SNDFILE)
    add_executable(denoiser_demo "${CMAKE_CURRENT_SOURCE_DIR}/examples/denoiser_demo.c")
    target_include_directories(denoiser_demo PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(denoiser_demo PRIVATE libspecbleach PkgConfig::SNDFILE)
    if(MATH_LIB)
        target_link_libraries(denoiser_demo PRIVATE ${MATH_LIB})
    endif()
    if(SPECBLEACH_BUILD_EXTRAS)
        add_executable(stereo_denoiser_demo "${CMAKE_CURRENT_SOURCE_DIR}/examples/stereo_denoiser_demo.c")
        target_include_directories(stereo_denoiser_demo PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
        target_link_libraries(stereo_denoiser_demo PRIVATE libspecbleach_extras PkgConfig::SNDFILE)
        if(MATH_LIB)
            target_link_libraries(stereo_denoiser_demo PRIVATE ${MATH_LIB})
        endif()
    endif()
    if(NOT TARGET PkgConfig::SNDFILE)
        message(WARNING
            "libsndfile not found. denoiser_demo and test_audio_file_regression will not be built.")
    endif()
endif()

if(ENABLE_TESTS)
    enable_testing()
    set(TEST_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_data/")

    # Unit tests exercise internal APIs directly, which are hidden in the
    # shared library by design. Link them against a static full-symbol copy of
    # the library sources instead. The public shared artifact is untouched.
    if(BUILD_SHARED_LIBS)
        add_library(specbleach_test_internal STATIC EXCLUDE_FROM_ALL
            ${SPECLEACH_SHARED_SOURCES}
            ${SPECLEACH_PROCESSOR_SOURCES}
            ${PFFFT_SOURCES}
        )
        target_compile_definitions(specbleach_test_internal PUBLIC
            SPECBLEACH_STATIC)
        target_compile_definitions(specbleach_test_internal PRIVATE
            PFFFT_STATIC_DEFINE)
        target_include_directories(specbleach_test_internal PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/src
            ${CMAKE_CURRENT_SOURCE_DIR}/include
        )
        target_include_directories(specbleach_test_internal SYSTEM PRIVATE
            ${PFFFT_INCLUDE_DIRS}
        )
        target_link_libraries(specbleach_test_internal PRIVATE
            Threads::Threads)
        if(TARGET specbleach_avx)
            target_link_libraries(specbleach_test_internal PRIVATE specbleach_avx)
        endif()
        set(SPECBLEACH_TEST_LIB specbleach_test_internal)
    else()
        set(SPECBLEACH_TEST_LIB libspecbleach)
    endif()

    macro(add_specbleach_test test_name src_file)
        add_executable(${test_name} "${CMAKE_CURRENT_SOURCE_DIR}/tests/${src_file}")
        target_compile_definitions(${test_name} PRIVATE "TEST_DATA_DIR=\"${TEST_DATA_DIR}\"")
        target_include_directories(${test_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/src)
        target_link_libraries(${test_name} PRIVATE ${SPECBLEACH_TEST_LIB} Threads::Threads)
        if(MATH_LIB)
            target_link_libraries(${test_name} PRIVATE ${MATH_LIB})
        endif()
        add_test(NAME ${test_name} COMMAND ${test_name})
    endmacro()

    add_specbleach_test(test_utils "test_utils.c")
    add_specbleach_test(test_integration "test_integration.c")
    add_specbleach_test(test_audio_regression "test_audio_regression.c")
    add_specbleach_test(test_audio_files "test_audio_files.c")
    add_specbleach_test(test_fft "test_fft.c")
    add_specbleach_test(test_spectral_whitening "test_spectral_whitening.c")
    add_specbleach_test(test_stft_latency "test_stft_latency.c")
    add_specbleach_test(test_gain_calculator "test_gain_calculator.c")
    add_specbleach_test(test_shared_utils "test_shared_utils.c")
    add_specbleach_test(test_frame_rate_norm "test_frame_rate_norm.c")
    add_specbleach_test(test_spectral_smoother "test_spectral_smoother.c")
    add_specbleach_test(test_spectral_utils "test_spectral_utils.c")
    add_specbleach_test(test_noise_profile "test_noise_profile.c")
    add_specbleach_test(test_specbleach_denoiser "test_specbleach_denoiser.c")
    add_specbleach_test(test_spectral_denoiser "test_spectral_denoiser.c")
    add_specbleach_test(test_noise_floor_manager "test_noise_floor_manager.c")
    add_specbleach_test(test_adaptive_noise_estimator "test_adaptive_noise_estimator.c")
    add_specbleach_test(test_nlm_filter "test_nlm_filter.c")
    add_specbleach_test(test_dftt_filter "test_dftt_filter.c")
    add_specbleach_test(test_thread_pool "test_thread_pool.c")
    add_specbleach_test(test_suppression_engine "test_suppression_engine.c")
    add_specbleach_test(test_masking_veto "test_masking_veto.c")
    add_specbleach_test(test_transient_detector "test_transient_detector.c")
    add_specbleach_test(test_spectral_circular_buffer "test_spectral_circular_buffer.c")
    add_specbleach_test(test_martin "test_martin_noise_estimator.c")
    add_specbleach_test(test_brandt "test_brandt_noise_estimator.c")
    add_specbleach_test(test_spp_mmse "test_spp_mmse_noise_estimator.c")
    add_specbleach_test(test_tonal_detector "test_tonal_detector.c")
    add_specbleach_test(test_tonal_reducer "test_tonal_reducer.c")
    add_specbleach_test(test_hpss_filter "test_hpss_filter.c")
    add_specbleach_test(test_denoiser_post_process "test_denoiser_post_process.c")
    add_specbleach_test(test_specbleach_version "test_specbleach_version.c")

    if(TARGET PkgConfig::SNDFILE)
        add_specbleach_test(test_audio_file_regression "test_audio_file_regression.c")
        target_link_libraries(test_audio_file_regression PRIVATE PkgConfig::SNDFILE)
    endif()

    if(SPECBLEACH_BUILD_EXTRAS)
        add_specbleach_test(test_specbleach_stereo "test_specbleach_stereo.c")
        target_link_libraries(test_specbleach_stereo PRIVATE libspecbleach_extras)
        if(CMAKE_CXX_COMPILER)
            add_executable(test_specbleach_cpp_smoke
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_specbleach_cpp_smoke.cpp")
            target_include_directories(test_specbleach_cpp_smoke PRIVATE
                ${CMAKE_CURRENT_SOURCE_DIR}/extras/include
                ${CMAKE_CURRENT_SOURCE_DIR}/include
            )
            target_link_libraries(test_specbleach_cpp_smoke PRIVATE
                libspecbleach_extras)
            add_test(NAME test_specbleach_cpp_smoke
                     COMMAND test_specbleach_cpp_smoke)
        endif()
    endif()
endif()

# ==============================================================================
# 8. Installation & Packaging Targets
# ==============================================================================
option(SPECBLEACH_INSTALL "Install libspecbleach binaries, headers, and package config" ON)

if(SPECBLEACH_INSTALL)
    install(TARGETS libspecbleach
        EXPORT libspecbleachTargets
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

    install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

    # Generate and Install CMake Package Config
    configure_package_config_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/libspecbleachConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/libspecbleachConfig.cmake"
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libspecbleach
    )

    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/libspecbleachConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
    )

    install(EXPORT libspecbleachTargets
        FILE libspecbleachTargets.cmake
        NAMESPACE libspecbleach::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libspecbleach
    )

    install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/libspecbleachConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/libspecbleachConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libspecbleach
    )

    # Optional extras layer install (own component so packagers can split it)
    if(SPECBLEACH_BUILD_EXTRAS)
        install(TARGETS libspecbleach_extras
            EXPORT libspecbleachExtrasTargets
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
            INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        )

        install(DIRECTORY extras/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

        install(EXPORT libspecbleachExtrasTargets
            FILE libspecbleachExtrasTargets.cmake
            NAMESPACE libspecbleach::
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libspecbleach
        )

        set(PKG_CONFIG_EXTRAS_FILE "${CMAKE_CURRENT_BINARY_DIR}/specbleach-extras.pc")
        file(WRITE "${PKG_CONFIG_EXTRAS_FILE}"
"prefix=${CMAKE_INSTALL_PREFIX}\n"
"exec_prefix=\${prefix}\n"
"libdir=\${exec_prefix}/${CMAKE_INSTALL_LIBDIR}\n"
"includedir=\${prefix}/${CMAKE_INSTALL_INCLUDEDIR}\n"
"\n"
"Name: specbleach-extras\n"
"Description: Optional orchestration helpers for libspecbleach\n"
"Version: ${PROJECT_VERSION}\n"
"Requires: specbleach = ${PROJECT_VERSION}\n"
"Libs: -L\${libdir} -lspecbleach_extras\n"
"Cflags: -I\${includedir}\n"
        )
        install(FILES "${PKG_CONFIG_EXTRAS_FILE}"
            DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
        )
    endif()

    # Generate pkg-config file (specbleach.pc)
    set(PKG_CONFIG_PREFIX "\${pcfiledir}/../..")
    set(PKG_CONFIG_EXEC_PREFIX "\${prefix}")
    set(PKG_CONFIG_VERSION "${PROJECT_VERSION}")

    if(IS_ABSOLUTE "${CMAKE_INSTALL_LIBDIR}")
        set(PKG_CONFIG_LIBDIR "${CMAKE_INSTALL_LIBDIR}")
    else()
        set(PKG_CONFIG_LIBDIR "\${exec_prefix}/${CMAKE_INSTALL_LIBDIR}")
    endif()
    if(IS_ABSOLUTE "${CMAKE_INSTALL_INCLUDEDIR}")
        set(PKG_CONFIG_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}")
    else()
        set(PKG_CONFIG_INCLUDEDIR "\${prefix}/${CMAKE_INSTALL_INCLUDEDIR}")
    endif()

    set(PKG_CONFIG_LIBS_PRIVATE "-lm")
    if(CMAKE_THREAD_LIBS_INIT)
        string(APPEND PKG_CONFIG_LIBS_PRIVATE " ${CMAKE_THREAD_LIBS_INIT}")
    endif()

    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/specbleach.pc.in"
        "${CMAKE_CURRENT_BINARY_DIR}/specbleach.pc"
        @ONLY
    )

    install(FILES "${CMAKE_CURRENT_BINARY_DIR}/specbleach.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
    )
endif()

