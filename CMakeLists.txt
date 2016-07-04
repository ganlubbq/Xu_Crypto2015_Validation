cmake_minimum_required(VERSION 3.3)
project(Xu_Crypto2015_Validation)

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS "-std=c++1z")
    set(NTL_DEPENDENT_LIBS gmp)
endif()

find_path(BM_INCLUDE BerlekampMassey.h /cygdrive/d/cxx/include d:/cxx/include)
message(STATUS "Find BerlekampMassey.h at ${BM_INCLUDE}")
find_path(NTL_INCLUDE NTL/version.h d:/cxx/wlib/include)
message(STATUS "Find NTL header at ${NTL_INCLUDE}")
find_package(Boost)

include_directories(${BM_INCLUDE} ${NTL_INCLUDE} ${Boost_INCLUDE_DIRS})


add_executable(validate main.cpp)


find_library(NTL_LIB ntl d:/cxx/wlib/lib)
message(STATUS "Find NTL lib at ${NTL_LIB}")

target_link_libraries(validate ${NTL_LIB} ${NTL_DEPENDENT_LIBS})
