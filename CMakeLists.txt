cmake_minimum_required(VERSION 3.22)
project(advanced_data_structures)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-std=c++20 -Wall -Wextra -mcx16 -lpthread)
# add_compile_options(-Wconversion -Wcast-align -Wunused -Wshadow -Wold-style-cast -Wpointer-arith -Wcast-qual -Wno-missing-braces)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -finline-functions -m64 -funroll-loops")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")

include_directories(${PROJECT_SOURCE_DIR}/src)


include(FetchContent)

# See: https://google.github.io/googletest/quickstart-cmake.html
#FetchContent_Declare(googletest URL https://github.com/google/googletest/archive/e2239ee6043f73722e7aa812a459f54a28552929.zip)
FetchContent_Declare(googletest GIT_REPOSITORY https://github.com/google/googletest.git GIT_TAG release-1.11.0)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

#FetchContent_Declare(benchmark URL https://github.com/google/benchmark/archive/refs/tags/v1.6.1.zip)
FetchContent_Declare(benchmark GIT_REPOSITORY https://github.com/google/benchmark.git GIT_TAG v1.6.1)

FetchContent_MakeAvailable(googletest benchmark)


## Targets

add_executable(ads_program_a src/main.cpp)

add_executable(generate_bv_data app/generate_bv_data.cpp)

## Tests

enable_testing()
include(GoogleTest)

add_executable(bv_test tests/bv_test.cpp)
target_link_libraries(bv_test gtest gtest_main)
gtest_discover_tests(bv_test)

add_executable(bp_test tests/bp_test.cpp)
target_link_libraries(bp_test gtest gtest_main)
gtest_discover_tests(bp_test)

add_executable(utils_test tests/utils_test.cpp)
target_link_libraries(utils_test gtest gtest_main)
gtest_discover_tests(utils_test)


## Benchmarks

add_executable(utils_benchmark benchmarks/utils_benchmark.cpp)
target_link_libraries(utils_benchmark benchmark::benchmark)

add_executable(bv_benchmark benchmarks/bv_benchmark.cpp)
target_link_libraries(bv_benchmark benchmark::benchmark)

add_executable(bp_benchmark benchmarks/bp_benchmark.cpp)
target_link_libraries(bp_benchmark benchmark::benchmark)
