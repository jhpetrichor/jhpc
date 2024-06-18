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

# # cns
add_executable(cns ${SOURCE_DIR}/cns.cpp)
target_link_libraries(cns cpdp)

# # cns
add_executable(tcns ${SOURCE_DIR}/two_cns.cpp)
target_link_libraries(tcns cpdp)