cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(custom_ops)

find_package(Torch REQUIRED)

add_executable(minimal-app minimal.cpp)
target_link_libraries(minimal-app "${TORCH_LIBRARIES}")
set_property(TARGET minimal-app PROPERTY CXX_STANDARD 11)
