# cmake_minimum_required(VERSION 3.10.0)
# project(weather VERSION 0.1.0 LANGUAGES C CXX)

# add_executable(weather main.cpp Log.cpp)

# include(CTest)
# enable_testing()

# set(CPACK_PROJECT_NAME ${PROJECT_NAME})
# set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
# include(CPack)



# cmake_minimum_required(VERSION 3.10)
# project(WeatherApp)

# # Automatically finds all .cpp files in your current folder
# file(GLOB SOURCES "*.cpp")

# add_executable(WeatherApp ${SOURCES})


cmake_minimum_required(VERSION 3.10)
project(playground)

# Force CMake to find the Mac System Headers
if(APPLE)
    execute_process(COMMAND xcrun --show-sdk-path OUTPUT_VARIABLE "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk" OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(CMAKE_OSX_SYSROOT ${/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk})
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

file(GLOB_RECURSE SOURCES
    "src/*.cpp"
)
add_executable(playground ${SOURCES})