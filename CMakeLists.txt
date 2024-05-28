cmake_minimum_required(VERSION 3.22)
project(complex_predict)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY  ${CMAKE_CURRENT_SOURCE_DIR}/lib)

# 设置动态库路径
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/lib)
# 设置头文件路径
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# 源文件
set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
# 库 源文件
file(GLOB LIB_SOURCE ${SOURCE_DIR}/lib/*.cpp)

# 编译动态库
add_library(cpdp SHARED ${LIB_SOURCE})

# test
add_executable(test ${SOURCE_DIR}/test.cpp)
target_link_libraries(test cpdp)

# # cns
add_executable(cns ${SOURCE_DIR}/cns.cpp)
target_link_libraries(cns cpdp)

# bops
add_executable(bops ${SOURCE_DIR}/bops.cpp)
target_link_libraries(bops cpdp)

# dpin
add_executable(dpin ${SOURCE_DIR}/dpin.cpp)
target_link_libraries(dpin cpdp)

# spin
add_executable(spin ${SOURCE_DIR}/spin.cpp)
target_link_libraries(spin cpdp)

add_executable(kpin ${SOURCE_DIR}/kpin.cpp)
target_link_libraries(kpin cpdp)