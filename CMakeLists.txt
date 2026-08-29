cmake_minimum_required(VERSION 3.16)

project(SupersonicDiamondAirfoil LANGUAGES CXX)

add_executable(diamond_solver
    src/euler_solver/main.cpp
)

target_compile_features(diamond_solver PRIVATE cxx_std_17)

if(MSVC)
    target_compile_options(diamond_solver PRIVATE /W4)
else()
    target_compile_options(diamond_solver PRIVATE -Wall -Wextra -Wpedantic)
endif()
