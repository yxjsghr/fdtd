cmake_minimum_required(VERSION 3.16)
project(fdtd3d LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 让 <cmath> 中的 M_PI 可用
add_compile_definitions(_USE_MATH_DEFINES)

# 编译优化（可按需改为 RelWithDebInfo）
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

include_directories(include)

add_executable(fdtd3d
    src/main_3d.cpp
    src/FDTD3D.cpp
    src/OutputManager.cpp
)

