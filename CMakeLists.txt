# CMakeLists.txt
cmake_minimum_required(VERSION 3.20)
project(lasm)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_USING_LINKER_mold)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
enable_testing()
# Compiler flags
# Debug flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg -O0 -DDEBUG -Wall -Wextra -Wpedantic -fsanitize=undefined,null,bounds -fsanitize-trap=all")
# Release flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g -pg -march=native -ftree-vectorize -DNDEBUG -Wall -Wextra -Wpedantic")
# RelWithDebInfo flags
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG -Wall -Wextra")
# MinSizeRel flags
set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG")

if(PGO_GENERATE)
    message(STATUS "Building with PGO instrumentation")
    add_compile_options(-fprofile-generate)
    add_link_options(-fprofile-generate)
endif()
if(PGO_USE)
    add_compile_options(-fprofile-use=${PROFILE_DATA_FILE} -fprofile-correction)
    add_link_options(-fprofile-use=${PROFILE_DATA_FILE})
endif()
# Set linker flags for release builds
if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -no-pie")
endif()

# # Remove problematic flags in debug mode
# if(CMAKE_BUILD_TYPE STREQUAL "Debug")
#     # Remove optimization that can interfere with debugging
#     string(REPLACE "-Wno-interference-size" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
# endif()
find_package(PkgConfig REQUIRED)
pkg_check_modules(HTS REQUIRED htslib)
include_directories(${HTS_INCLUDE_DIRS})
link_directories(${HTS_LIBRARY_DIRS})
add_definitions(${HTS_CFLAGS_OTHER})

set(SOURCE_FILES lasm.cpp asm2.cpp stats.cpp lcmer.cpp sthist.cpp)
include_directories(./ext)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

# External dependencies
include(FetchContent)

FetchContent_Declare(
    cli11_proj
    QUIET
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG v2.5.0
)
FetchContent_MakeAvailable(cli11_proj)

FetchContent_Declare(
    gtl
    GIT_REPOSITORY https://github.com/greg7mdp/gtl.git
    GIT_TAG v1.2.0
)
FetchContent_MakeAvailable(gtl)

FetchContent_Declare(
    concurrentqueue
    GIT_REPOSITORY https://github.com/cameron314/concurrentqueue.git
    GIT_TAG v1.0.4
)
FetchContent_MakeAvailable(concurrentqueue)

FetchContent_Declare(
    libbloom
    GIT_REPOSITORY https://github.com/jvirkki/libbloom.git
    GIT_TAG v2.0
)
FetchContent_MakeAvailable(libbloom)


# Link libraries
target_link_libraries(${PROJECT_NAME} ${HTS_LIBRARIES} gtl CLI11::CLI11 concurrentqueue pthread deflate z)


# Tests and custom targets
add_test(NAME test COMMAND lasm hifi test/human2.mapped.sorted_small_set.bam test.gfa -l3 -t8)
add_custom_target(rt COMMAND lasm hifi test/human2.mapped.sorted_small_set.bam test.gfa -l3 -t8)
add_dependencies(rt ${PROJECT_NAME})

# Debug targets
add_custom_target(debug
    COMMENT "Building in debug mode"
    COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${PROJECT_NAME}
)


# Compile definitions
target_compile_definitions(${PROJECT_NAME} PUBLIC
    SYNCMER_RULE=DensitySyncmerRule
    LCMER_L=8
    SYNCMER_K=8
)


# Debug-specific definitions
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(${PROJECT_NAME} PUBLIC DEBUG_MODE=1)

    # Enable additional debug features
    target_compile_definitions(${PROJECT_NAME} PUBLIC ENABLE_DEBUG_LOGGING=1)

endif()
