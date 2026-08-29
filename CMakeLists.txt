cmake_minimum_required(VERSION 4.3)
project(NE684a2)

set(CMAKE_CXX_STANDARD 26)

add_executable(NE684a2 main.cpp
        neutron.h
        RandomManager.h
        CrossSections.h
)
