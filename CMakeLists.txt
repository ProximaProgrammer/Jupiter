cmake_minimum_required(VERSION 3.16)
project(JupiterRT VERSION 0.4.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_executable(jupiter
  src/main.cpp
  src/config/SimConfig.cpp
  src/core/Grid.cpp
  src/physics/EOS.cpp
  src/physics/Dynamics.cpp
  src/physics/RadiativeTransfer.cpp
  src/physics/TimeStep.cpp
  src/io/Diagnostics.cpp
  src/io/Visualizer.cpp
)

target_include_directories(jupiter PRIVATE src)

target_compile_options(jupiter PRIVATE
  $<$<CXX_COMPILER_ID:AppleClang,Clang,GNU>:-Wall -Wextra -Wpedantic -Wno-unused-parameter>
)
