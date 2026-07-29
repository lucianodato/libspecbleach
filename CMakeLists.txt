cmake_minimum_required(VERSION 3.16)
project(libspecbleach VERSION 0.3.0 LANGUAGES C)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)

option(BUILD_SHARED_LIBS "Build shared libraries" OFF)
option(ENABLE_AVX "Enable AVX optimizations on x86_64" ON)
option(ENABLE_TESTS "Build test suite" ON)
option(ENABLE_EXAMPLES "Build example executables" ON)
option(ENABLE_SANITIZERS "Enable Address and Undefined Behavior Sanitizers" OFF)

# Dependencies
find_library(MATH_LIB m)
find_path(FFTW3F_INCLUDE_DIR fftw3.h)
find_library(FFTW3F_LIB fftw3f REQUIRED)

find_package(OpenMP QUIET)

if(NOT OpenMP_C_FOUND)
    find_library(OMP_LIB omp HINTS /opt/homebrew/opt/libomp/lib /usr/local/lib)
    find_path(OMP_INCLUDE_DIR omp.h HINTS /opt/homebrew/opt/libomp/include /usr/local/include)
    if(OMP_LIB AND OMP_INCLUDE_DIR)
        add_library(OpenMP::OpenMP_C INTERFACE IMPORTED)
        set_target_properties(OpenMP::OpenMP_C PROPERTIES
            INTERFACE_COMPILE_OPTIONS "-Xpreprocessor;-fopenmp"
            INTERFACE_INCLUDE_DIRECTORIES "${OMP_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "${OMP_LIB}"
        )
        set(OpenMP_C_FOUND TRUE)
    endif()
endif()
find_package(PkgConfig QUIET)

find_path(SNDFILE_INCLUDE_DIR sndfile.h HINTS /opt/homebrew/include /usr/local/include)
find_library(SNDFILE_LIB sndfile HINTS /opt/homebrew/lib /usr/local/lib)
if(SNDFILE_LIB)
    set(SNDFILE_FOUND TRUE)
    set(SNDFILE_LIBRARIES ${SNDFILE_LIB})
endif()

# Compiler Warning Flags
if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic -Wno-unused-parameter -Wno-unused-variable)
endif()

if(ENABLE_SANITIZERS AND NOT MSVC)
    add_compile_options(-fsanitize=address -fsanitize=undefined)
    add_link_options(-fsanitize=address -fsanitize=undefined)
endif()

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

if(FFTW3F_INCLUDE_DIR)
    include_directories(${FFTW3F_INCLUDE_DIR})
endif()

# Collect all shared and processor sources
file(GLOB_RECURSE SPECLEACH_SHARED_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/shared/*.c"
)
list(FILTER SPECLEACH_SHARED_SOURCES EXCLUDE REGEX ".*nlm_filter_avx\\.c$")

file(GLOB_RECURSE SPECLEACH_PROCESSOR_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/processors/*.c"
)

# Main libspecbleach target
add_library(libspecbleach
    ${SPECLEACH_SHARED_SOURCES}
    ${SPECLEACH_PROCESSOR_SOURCES}
)

target_include_directories(libspecbleach PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

target_include_directories(libspecbleach PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

if(MATH_LIB)
    target_link_libraries(libspecbleach PUBLIC ${MATH_LIB})
endif()

target_link_libraries(libspecbleach PUBLIC ${FFTW3F_LIB})

if(OpenMP_C_FOUND)
    target_link_libraries(libspecbleach PRIVATE OpenMP::OpenMP_C)
endif()

# AVX compilation target for x86_64 architecture
if(ENABLE_AVX AND (CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64" OR CMAKE_SYSTEM_PROCESSOR MATCHES "AMD64"))
    add_library(specbleach_avx STATIC "${CMAKE_CURRENT_SOURCE_DIR}/src/shared/denoiser_logic/processing/nlm_filter_avx.c")
    target_compile_options(specbleach_avx PRIVATE -mavx)
    target_include_directories(specbleach_avx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/src)
    if(FFTW3F_INCLUDE_DIR)
        target_include_directories(specbleach_avx PRIVATE ${FFTW3F_INCLUDE_DIR})
    endif()
    target_link_libraries(specbleach_avx PRIVATE ${FFTW3F_LIB})
    target_link_libraries(libspecbleach PRIVATE specbleach_avx)
endif()

# Alias target
add_library(libspecbleach::libspecbleach ALIAS libspecbleach)

# Examples
if(ENABLE_EXAMPLES AND SNDFILE_FOUND)
    add_executable(denoiser_demo "${CMAKE_CURRENT_SOURCE_DIR}/examples/denoiser_demo.c")
    target_include_directories(denoiser_demo PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/src)
    if(SNDFILE_INCLUDE_DIRS)
        target_include_directories(denoiser_demo PRIVATE ${SNDFILE_INCLUDE_DIRS})
    endif()
    target_link_libraries(denoiser_demo PRIVATE libspecbleach ${SNDFILE_LIBRARIES})
endif()

# CTest Suite
if(ENABLE_TESTS)
    enable_testing()
    set(TEST_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_data/")

    macro(add_specbleach_test test_name src_file)
        add_executable(${test_name} "${CMAKE_CURRENT_SOURCE_DIR}/tests/${src_file}")
        target_compile_definitions(${test_name} PRIVATE "TEST_DATA_DIR=\"${TEST_DATA_DIR}\"")
        target_include_directories(${test_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/src)
        target_link_libraries(${test_name} PRIVATE libspecbleach ${FFTW3F_LIB})
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
    add_specbleach_test(test_spectral_utils "test_spectral_utils.c")
    add_specbleach_test(test_noise_profile "test_noise_profile.c")
    add_specbleach_test(test_specbleach_processor_core "test_specbleach_processor_core.c")
    add_specbleach_test(test_specbleach_denoiser "test_specbleach_denoiser.c")
    add_specbleach_test(test_noise_floor_manager "test_noise_floor_manager.c")
    add_specbleach_test(test_adaptive_noise_estimator "test_adaptive_noise_estimator.c")
    add_specbleach_test(test_nlm_filter "test_nlm_filter.c")
    add_specbleach_test(test_specbleach_2d_denoiser "test_specbleach_2d_denoiser.c")
    add_specbleach_test(test_suppression_engine "test_suppression_engine.c")
    add_specbleach_test(test_masking_veto "test_masking_veto.c")
    add_specbleach_test(test_transient_detector "test_transient_detector.c")
    add_specbleach_test(test_spectral_circular_buffer "test_spectral_circular_buffer.c")
    add_specbleach_test(test_martin "test_martin_noise_estimator.c")
    add_specbleach_test(test_brandt "test_brandt_noise_estimator.c")
    add_specbleach_test(test_spp_mmse "test_spp_mmse_noise_estimator.c")
    add_specbleach_test(test_tonal_detector "test_tonal_detector.c")
    add_specbleach_test(test_tonal_reducer "test_tonal_reducer.c")
    add_specbleach_test(test_denoiser_post_process "test_denoiser_post_process.c")

    if(SNDFILE_FOUND)
        add_specbleach_test(test_audio_file_regression "test_audio_file_regression.c")
        if(SNDFILE_INCLUDE_DIRS)
            target_include_directories(test_audio_file_regression PRIVATE ${SNDFILE_INCLUDE_DIRS})
        endif()
        target_link_libraries(test_audio_file_regression PRIVATE ${SNDFILE_LIBRARIES})
    endif()
endif()
