# CMakeLists.txt
cmake_minimum_required(VERSION 3.20)
project(lasm)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Use the mold linker when available (much faster links). mold only links ELF,
# so verify the compiler actually accepts it before selecting it -- on macOS
# (Mach-O) the flag is present but unusable.
find_program(MOLD_LINKER mold)
if(MOLD_LINKER)
    include(CheckLinkerFlag)
    check_linker_flag(CXX "-fuse-ld=mold" LASM_HAVE_MOLD)
    if(LASM_HAVE_MOLD)
        if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.29)
            set(CMAKE_LINKER_TYPE MOLD)
        else()
            add_link_options(-fuse-ld=mold)
        endif()
        message(STATUS "Using mold linker: ${MOLD_LINKER}")
    endif()
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
# Build options
option(LASM_NATIVE "Optimize for the build host CPU (-march=native); produces non-portable binaries" OFF)

# Compiler flags
# Project defaults come first; a user's -DCMAKE_<LANG>_FLAGS_<CONFIG> is
# appended so it reaches the compile line and wins over the defaults.
# Debug flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG -Wall -Wextra -Wpedantic -fsanitize=undefined,null,bounds -fsanitize-trap=all ${CMAKE_CXX_FLAGS_DEBUG}")
# Release flags. -g1 keeps line tables (backtraces, profiles) and drops the
# variable locations, which cost about a quarter of the compile time.
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g1 -ftree-vectorize -DNDEBUG -Wall -Wextra -Wpedantic ${CMAKE_CXX_FLAGS_RELEASE}")
# RelWithDebInfo flags — includes gprof instrumentation (-pg, -no-pie)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -pg -DNDEBUG -Wall -Wextra -ftree-vectorize ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
set(CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO "-pg -no-pie ${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO}")
# MinSizeRel flags
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG ${CMAKE_CXX_FLAGS_MINSIZEREL}")

#We are not providing a public ABI, no need to fix this
add_compile_options($<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wno-interference-size>)

# Optimize for the host CPU only when explicitly requested (non-portable;
# default off so distributed binaries run on other machines). Portable x86-64
# builds target x86-64-v2 (SSSE3/SSE4.2/POPCNT; every x86-64 CPU since ~2009)
# so the SSSE3 SIMD paths are used instead of the scalar fallbacks. ARM64
# always has NEON and needs no flag.
option(LASM_X86_64_V2 "Portable x86-64 builds target x86-64-v2 (SSSE3/SSE4.2)" ON)
if(LASM_NATIVE)
    add_compile_options($<$<CONFIG:Release,RelWithDebInfo>:-march=native>)
    add_compile_options($<$<CONFIG:RelWithDebInfo>:-mtune=native>)
elseif(LASM_X86_64_V2 AND CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|AMD64|amd64)$")
    add_compile_options(-march=x86-64-v2)
endif()

if(PGO_GENERATE)
    message(STATUS "Building with PGO instrumentation")
    add_compile_options(-fprofile-generate)
    add_link_options(-fprofile-generate)
endif()
if(PGO_USE)
    add_compile_options(-fprofile-use=${PROFILE_DATA_FILE} -fprofile-correction)
    add_link_options(-fprofile-use=${PROFILE_DATA_FILE})
endif()
# # Remove problematic flags in debug mode
# if(CMAKE_BUILD_TYPE STREQUAL "Debug")
#     # Remove optimization that can interfere with debugging
#     string(REPLACE "-Wno-interference-size" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
# endif()
# libdeflate inflates BGZF blocks; zlib-ng (fetched below) streams plain gzip.
find_package(PkgConfig REQUIRED)
pkg_check_modules(DEFLATE REQUIRED IMPORTED_TARGET libdeflate)

set(SOURCE_FILES lasm.cpp asm2.cpp dump_cores.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_include_directories(${PROJECT_NAME} PRIVATE ./ext)

# External dependencies
include(FetchContent)


FetchContent_Declare(
    cli11_proj
    QUIET
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG v2.5.0
)
FetchContent_MakeAvailable(cli11_proj)
# CLI11's implementation is compiled once, in its own source file, instead of
# inline in every source file that includes CLI11. It is part of lasm rather
# than CLI11's precompiled library, which the Makefile generator would build
# before starting on lasm's sources.
set(CLI11_IMPL ${cli11_proj_SOURCE_DIR}/src/Precompile.cpp)
target_sources(${PROJECT_NAME} PRIVATE ${CLI11_IMPL})
set_source_files_properties(lasm.cpp dump_cores.cpp ${CLI11_IMPL}
                            PROPERTIES COMPILE_DEFINITIONS CLI11_COMPILE)

FetchContent_Declare(
    gtl
    GIT_REPOSITORY https://github.com/greg7mdp/gtl.git
    GIT_TAG v1.2.0
)
FetchContent_MakeAvailable(gtl)

# v1.0.5 orders an implicit-producer block's reuse after every dequeue from
# it (upstream 0d54c67, e74b07d); v1.0.4 left that a data race.
FetchContent_Declare(
    concurrentqueue
    GIT_REPOSITORY https://github.com/cameron314/concurrentqueue.git
    GIT_TAG v1.0.5
)
FetchContent_MakeAvailable(concurrentqueue)

# ksw2 (lh3/ksw2, minimap2's aligner): SIMD dual-affine z-drop extension for
# the end joins (align.h). It has no CMake build, so only its dual-gap SSE
# kernel is compiled here.
FetchContent_Declare(
    ksw2
    GIT_REPOSITORY https://github.com/lh3/ksw2.git
    GIT_TAG 289609bd9e5381a13b16239d0a7703f1ff03f9ca
)
FetchContent_MakeAvailable(ksw2)
add_library(ksw2_extd2 STATIC ${ksw2_SOURCE_DIR}/ksw2_extd2_sse.c)
target_include_directories(ksw2_extd2 SYSTEM PUBLIC ${ksw2_SOURCE_DIR})
target_compile_options(ksw2_extd2 PRIVATE -w)
if(NOT CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|i[3-6]86")
    # No SSE: translate the intrinsics with sse2neon, as minimap2 does on ARM.
    FetchContent_Declare(
        sse2neon
        GIT_REPOSITORY https://github.com/DLTcollab/sse2neon.git
        GIT_TAG v1.8.0
    )
    FetchContent_MakeAvailable(sse2neon)
    set(SSE2NEON_SHIM ${CMAKE_BINARY_DIR}/sse2neon_shim)
    file(WRITE ${SSE2NEON_SHIM}/emmintrin.h "#include \"sse2neon.h\"\n")
    file(WRITE ${SSE2NEON_SHIM}/smmintrin.h "#include \"sse2neon.h\"\n")
    target_include_directories(ksw2_extd2 BEFORE PRIVATE
                               ${SSE2NEON_SHIM} ${sse2neon_SOURCE_DIR})
    target_compile_definitions(ksw2_extd2 PRIVATE __SSE2__ KSW_SSE2_ONLY)
endif()

# zlib-ng in zlib-compatible mode: the zlib API behind <zlib.h>, with a faster
# inflate for plain gzip input. Built as a static library; it picks its SIMD
# code paths at run time, so the binary stays portable.
set(ZLIB_COMPAT ON)
set(BUILD_SHARED_LIBS OFF)
set(BUILD_TESTING OFF)
set(SKIP_INSTALL_ALL ON)
FetchContent_Declare(
    zlib_ng
    GIT_REPOSITORY https://github.com/zlib-ng/zlib-ng.git
    GIT_TAG 12731092979c6d07f42da27da673a9f6c7b13586 # 2.3.3
)
FetchContent_MakeAvailable(zlib_ng)

# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE
                      zlib-ng::zlibstatic PkgConfig::DEFLATE gtl CLI11::CLI11 ksw2_extd2
                      concurrentqueue pthread
                      $<$<PLATFORM_ID:Linux>:rt>)

# Link-time optimisation in Release builds with GCC: code generation moves to
# the link, where GCC splits it into partitions compiled in parallel, instead
# of running serially per source file (asm2.cpp alone dominated the build).
# Clang's ThinLTO only parallelises across source files, so it is left off.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set_property(TARGET ${PROJECT_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
endif()

# Convenience target: run the assembler on a local sample (data not in repo).
add_custom_target(rtest COMMAND lasm hifi test/human2.mapped.sorted_small_set.bam test.gfa -l3 -t8)
add_dependencies(rtest ${PROJECT_NAME})

# Debug targets
add_custom_target(debug
    COMMENT "Building in debug mode"
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${PROJECT_NAME}
)


# Compile definitions
target_compile_definitions(${PROJECT_NAME} PRIVATE
    __STDC_FORMAT_MACROS=1
)



# Debug-specific definitions
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(${PROJECT_NAME} PRIVATE DEBUG_MODE=1)

    # Enable additional debug features
    target_compile_definitions(${PROJECT_NAME} PRIVATE ENABLE_DEBUG_LOGGING=1)

endif()
