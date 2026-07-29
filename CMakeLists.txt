cmake_minimum_required(VERSION 3.16)
project(libspecbleach VERSION 0.3.0 LANGUAGES C)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)

option(BUILD_SHARED_LIBS "Build shared libraries instead of static" OFF)
option(ENABLE_AVX "Enable AVX optimizations on x86_64" ON)
option(ENABLE_TESTS "Build test suite" OFF)

find_library(MATH_LIB m)
find_path(FFTW3F_INCLUDE_DIR fftw3.h)
find_library(FFTW3F_LIB fftw3f REQUIRED)

find_package(OpenMP)

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

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

if(FFTW3F_INCLUDE_DIR)
    target_include_directories(libspecbleach PRIVATE ${FFTW3F_INCLUDE_DIR})
endif()

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

# Create alias target for clean CMake usage
add_library(libspecbleach::libspecbleach ALIAS libspecbleach)
