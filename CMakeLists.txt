cmake_minimum_required(VERSION 3.1)
project(inv)

find_package(ALPSCore REQUIRED COMPONENTS hdf5 gf params)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${DMFT_DIR}/cmake")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH true)

find_package(LAPACK REQUIRED)
find_package(MPI REQUIRED)
find_package(NFFT3 REQUIRED)
find_package(BlasClasses REQUIRED)
find_package(Cluster REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(BLAS REQUIRED)

find_package(Boost COMPONENTS program_options REQUIRED)

include(UseGtest)

include_directories(${BLAS_CLASSES_INCLUDE_DIR})
include_directories(${CLUSTER_INCLUDE_DIR})
include_directories(${GTEST_INCLUDE_DIR})
include_directories(${NFFT3_INCLUDE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g -DNDEBUG -DBOOST_DISABLE_ASSERTS")

add_executable(inverter inverter.cpp mc_vertex.cpp vertex.cpp k_space_structure.cpp vertexaux.cpp
chi0_nambu.cpp  chi_nambu.cpp  k_space_structure.cpp  mc_vertex.cpp  phase_boundary_task.cpp  single_freq_gf.cpp  task.cpp  vertex.cpp F.cpp)
target_link_libraries(inverter 
${CLUSTER_LIBRARIES}
${ALPSCore_LIBRARIES}
${MPI_CXX_LIBRARIES}
${LAPACK_LIBRARIES}
${NFFT3_LIBRARIES} ${GTEST_LIBRARY} ${Boost_LIBRARIES})


install(TARGETS inverter RUNTIME DESTINATION "bin/")
