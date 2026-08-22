cmake_minimum_required(VERSION 3.22)
project(libspecbleach VERSION 0.3.2 LANGUAGES C)

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
option(USE_SYSTEM_FFTW "Use system-installed FFTW3 library" ON)

# ==============================================================================
# 1. FFTW3 Dependency
# ==============================================================================
if(USE_SYSTEM_FFTW)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(FFTW3F REQUIRED IMPORTED_TARGET fftw3f)
    set(FFTW_TARGET PkgConfig::FFTW3F)
    set(PKG_CONFIG_REQUIRES_PRIVATE "fftw3f")
else()
    include(FetchContent)
    set(BUILD_SHARED_LIBS_SAVED ${BUILD_SHARED_LIBS})
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
    set(ENABLE_FLOAT ON CACHE BOOL "" FORCE)
    set(ENABLE_THREADS ON CACHE BOOL "" FORCE)

    FetchContent_Declare(
        fftw3
        URL https://www.fftw.org/fftw-3.3.10.tar.gz
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_MakeAvailable(fftw3)
    set(FFTW_TARGET fftw3f)
    set(BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS_SAVED} CACHE BOOL "Build shared libraries" FORCE)
    set(PKG_CONFIG_REQUIRES_PRIVATE "")
endif()

# ==============================================================================
# 2. Math & OpenMP Dependencies
# ==============================================================================
find_library(MATH_LIB m)

if(APPLE)
    list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/opt/libomp" "/usr/local/opt/libomp")
endif()

find_package(OpenMP QUIET)
if(MSVC AND OpenMP_C_FOUND)
    set_target_properties(OpenMP::OpenMP_C PROPERTIES
        INTERFACE_COMPILE_OPTIONS "/openmp:llvm"
    )
elseif(NOT OpenMP_C_FOUND)
    message(STATUS "OpenMP not found. Searching Homebrew paths...")
    find_path(OMP_INCLUDE_DIR omp.h HINTS /opt/homebrew/opt/libomp/include /usr/local/opt/libomp/include)
    find_library(OMP_LIB NAMES omp gomp HINTS /opt/homebrew/opt/libomp/lib /usr/local/opt/libomp/lib)
    if(OMP_LIB AND OMP_INCLUDE_DIR)
        add_library(OpenMP::OpenMP_C INTERFACE IMPORTED)
        if(CMAKE_C_COMPILER_ID MATCHES "Clang" AND APPLE)
            message(STATUS "OpenMP found: ${OMP_LIB}")
            set_target_properties(OpenMP::OpenMP_C PROPERTIES
                INTERFACE_COMPILE_OPTIONS "-Xpreprocessor;-fopenmp"
                INTERFACE_INCLUDE_DIRECTORIES "${OMP_INCLUDE_DIR}"
                INTERFACE_LINK_LIBRARIES "${OMP_LIB}"
            )
        else()

            set_target_properties(OpenMP::OpenMP_C PROPERTIES
                INTERFACE_COMPILE_OPTIONS "-fopenmp"
                INTERFACE_INCLUDE_DIRECTORIES "${OMP_INCLUDE_DIR}"
                INTERFACE_LINK_LIBRARIES "${OMP_LIB}"
            )
        endif()
        set(OpenMP_C_FOUND TRUE)
    endif()
endif()

if(OpenMP_C_FOUND OR TARGET OpenMP::OpenMP_C)
    message(STATUS "OpenMP: ENABLED (Multi-threaded NLM 2D filter enabled)")
else()
    message(WARNING "OpenMP: NOT FOUND (NLM 2D filter will run single-threaded)")
endif()

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
)

set_target_properties(libspecbleach PROPERTIES
    OUTPUT_NAME "specbleach"
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    C_VISIBILITY_PRESET default
)

target_include_directories(libspecbleach PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_include_directories(libspecbleach PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

if(MATH_LIB)
    target_link_libraries(libspecbleach PRIVATE ${MATH_LIB})
endif()

# FFTW is an internal implementation detail (not exposed in public headers) -> PRIVATE
target_link_libraries(libspecbleach PRIVATE ${FFTW_TARGET})

if(OpenMP_C_FOUND OR TARGET OpenMP::OpenMP_C)
    target_link_libraries(libspecbleach PRIVATE OpenMP::OpenMP_C)
endif()

# AVX SIMD Target
if(ENABLE_AVX)
    set(AVX_C_FILE "${CMAKE_CURRENT_SOURCE_DIR}/src/shared/denoiser_logic/processing/nlm_filter_avx.c")
    if(EXISTS "${AVX_C_FILE}" AND CMAKE_SYSTEM_PROCESSOR MATCHES "(x86_64|amd64|AMD64|x86-64)")
        add_library(specbleach_avx STATIC "${AVX_C_FILE}")
        if(MSVC)
            target_compile_options(specbleach_avx PRIVATE /arch:AVX)
        else()
            target_compile_options(specbleach_avx PRIVATE -mavx)
        endif()
        target_include_directories(specbleach_avx PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/src)
        target_link_libraries(specbleach_avx PRIVATE ${FFTW_TARGET})
        if(OpenMP_C_FOUND OR TARGET OpenMP::OpenMP_C)
            target_link_libraries(specbleach_avx PRIVATE OpenMP::OpenMP_C)
        endif()
        target_link_libraries(libspecbleach PRIVATE specbleach_avx)
    endif()
endif()

# Alias Target for Subdirectories/FetchContent
add_library(libspecbleach::libspecbleach ALIAS libspecbleach)
add_library(specbleach ALIAS libspecbleach)

# ==============================================================================
# 7. Examples & Tests
# ==============================================================================
if(ENABLE_EXAMPLES AND TARGET PkgConfig::SNDFILE)
    add_executable(denoiser_demo "${CMAKE_CURRENT_SOURCE_DIR}/examples/denoiser_demo.c")
    target_include_directories(denoiser_demo PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/src)
    target_link_libraries(denoiser_demo PRIVATE libspecbleach PkgConfig::SNDFILE)
    if(MATH_LIB)
        target_link_libraries(denoiser_demo PRIVATE ${MATH_LIB})
    endif()
    if(NOT TARGET PkgConfig::SNDFILE)
        message(WARNING
            "libsndfile not found. denoiser_demo and test_audio_file_regression will not be built.")
    endif()
endif()

if(ENABLE_TESTS)
    enable_testing()
    set(TEST_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_data/")

    macro(add_specbleach_test test_name src_file)
        add_executable(${test_name} "${CMAKE_CURRENT_SOURCE_DIR}/tests/${src_file}")
        target_compile_definitions(${test_name} PRIVATE "TEST_DATA_DIR=\"${TEST_DATA_DIR}\"")
        target_include_directories(${test_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/src)
        target_link_libraries(${test_name} PRIVATE libspecbleach ${FFTW_TARGET})
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
    add_specbleach_test(test_spectral_smoother "test_spectral_smoother.c")
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
    add_specbleach_test(test_hpss_filter "test_hpss_filter.c")
    add_specbleach_test(test_denoiser_post_process "test_denoiser_post_process.c")

    if(TARGET PkgConfig::SNDFILE)
        add_specbleach_test(test_audio_file_regression "test_audio_file_regression.c")
        target_link_libraries(test_audio_file_regression PRIVATE PkgConfig::SNDFILE)
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

    set(PKG_CONFIG_OPENMP_FLAGS "")
    if(OpenMP_C_FOUND)
        set(PKG_CONFIG_OPENMP_FLAGS "${OpenMP_C_FLAGS}")
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

