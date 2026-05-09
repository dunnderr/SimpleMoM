cmake_minimum_required(VERSION 4.0)
project(SimpleMoM)

set(CMAKE_CXX_STANDARD 20)

add_executable(SimpleMoM main.cpp)

find_package(Threads REQUIRED)
target_link_libraries(SimpleMoM Threads::Threads m)
