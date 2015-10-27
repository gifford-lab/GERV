cmake_minimum_required(VERSION 2.8)

project(mpi_motif CXX)

set(LIB_FLAGS "-static-libstdc++ -fopenmp")
set(W_FLAGS "-Wall -Wextra -Wshadow -pedantic -Werror -Weffc++")
set(FLOAT_FLAGS "-mssse3 -mfpmath=sse -fexcess-precision=fast")
set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -fno-operator-names ${LIB_FLAGS} ${W_FLAGS} ${FLOAT_FLAGS}")


set(Boost_USE_STATIC_LIBS ON)
find_package(Boost COMPONENTS system program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
include_directories(${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})

set(SUBDIRS classes fmin util)
foreach(subdir ${SUBDIRS})
  file(GLOB LIB_SOURCES_${subdir} ${subdir}/*.cpp)
  file(GLOB LIB_HEADERS_${subdir} ${subdir}/*.h)
  set(LIB_SOURCES ${LIB_SOURCES} ${LIB_SOURCES_${subdir}})
  set(LIB_HEADERS ${LIB_HEADERS} ${LIB_HEADERS_${subdir}})
endforeach()

add_library(mpi_motif_lib ${LIB_SOURCES} ${LIB_HEADERS})

macro (add_executable _name)
  # invoke built-in add_executable
  _add_executable(${ARGV} consts.h)
  if (TARGET ${_name})
    target_link_libraries(${_name} mpi_motif_lib ${Boost_LIBRARIES} -lgsl -lgslcblas)
  endif()
endmacro()

add_executable(mpi_motif mpi_motif.cpp)
add_executable(genome_ignore_n preprocess/genome_ignore_n)
add_executable(reads preprocess/reads.cpp)
add_executable(ascending postprocess/ascending.cpp)
add_executable(validation postprocess/validation.cpp)
