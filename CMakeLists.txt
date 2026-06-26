cmake_minimum_required(VERSION 3.20)

project(JupiterRT VERSION 0.3.0 LANGUAGES CXX)

option(JUPITER_ENABLE_CUDA "Build CUDA accelerated kernels for NVIDIA GPUs" OFF)
option(JUPITER_WARNINGS_AS_ERRORS "Treat warnings as errors" OFF)
option(JUPITER_NATIVE_OPT "Use native CPU optimization flags on local machine" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(APPLE)
  execute_process(COMMAND xcrun --show-sdk-path
                  OUTPUT_VARIABLE _macosx_sdk
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(_macosx_sdk)
    set(CMAKE_OSX_SYSROOT "${_macosx_sdk}" CACHE PATH "macOS SDK" FORCE)
  endif()
endif()

file(GLOB_RECURSE JUPITER_CPP CONFIGURE_DEPENDS "src/*.cpp")
list(FILTER JUPITER_CPP EXCLUDE REGEX ".*/src/legacy/.*")

add_executable(jupiter ${JUPITER_CPP})
target_include_directories(jupiter PRIVATE src)

target_compile_options(jupiter PRIVATE
  $<$<CXX_COMPILER_ID:Clang,AppleClang,GNU>:-Wall -Wextra -Wpedantic>
  $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

if(JUPITER_NATIVE_OPT)
  target_compile_options(jupiter PRIVATE
    $<$<CXX_COMPILER_ID:Clang,AppleClang,GNU>:-O3>
  )
endif()

if(JUPITER_WARNINGS_AS_ERRORS)
  target_compile_options(jupiter PRIVATE
    $<$<CXX_COMPILER_ID:Clang,AppleClang,GNU>:-Werror>
    $<$<CXX_COMPILER_ID:MSVC>:/WX>
  )
endif()

if(JUPITER_ENABLE_CUDA)
  enable_language(CUDA)
  file(GLOB_RECURSE JUPITER_CUDA CONFIGURE_DEPENDS "src/*.cu")
  target_sources(jupiter PRIVATE ${JUPITER_CUDA})
  target_compile_definitions(jupiter PRIVATE JUPITER_ENABLE_CUDA=1)
  set_property(TARGET jupiter PROPERTY CUDA_STANDARD 17)
  set_property(TARGET jupiter PROPERTY CUDA_STANDARD_REQUIRED ON)
  # RTX 3060 desktop/laptop = Ampere compute capability 8.6.
  set_property(TARGET jupiter PROPERTY CUDA_ARCHITECTURES 86)
  target_compile_options(jupiter PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--use_fast_math>)
endif()

install(TARGETS jupiter RUNTIME DESTINATION bin)
