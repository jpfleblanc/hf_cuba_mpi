cmake_minimum_required(VERSION 3.1)
project(diagami CXX)

find_package(Boost 1.71.0)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

find_package( Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

#add_executable(${PROJECT_NAME} main.cpp ising.cpp)
link_directories(${Boost_LIBRARY_DIRS})


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O3 -std=c++11 -lstdc++fs -fopenmp")
add_executable(hf_cuba_mpi 
	./src/hf_cuba_mpi_rewrite.cpp 
	# ./src/hf_cuba_all.cpp 
)

# add_executable(${PROJECT_NAME}_mpi 
# ./src/main_diagami.cpp 
# ./src/diagami_ami.cpp 
# ./src/diagami_bc.cpp 
# ./src/sim_diagami.cpp 
# ./src/parameters_diagami.cpp
# ./src/diagami_V.cpp
# )


#set(EXEC_NAME ${PROJECT_NAME}_mpi)


# Request the ALPSCore package (with all components)
# The ALPSCore package must be installed in some
# standard place (like /usr/local),
# or somewhere in your PATH. 
# Otherwise, point ALPSCore_DIR environment variable
# to the ALPScore installation tree.
# find_package(ALPSCore 2.0 REQUIRED)


SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -std=c++11 -DUSE_MPI -lstdc++fs")# -ffast-math")

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../libami/install/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../libami/install/include)

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../libamicalc/install/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../libamicalc/install/include)
link_directories(${CMAKE_CURRENT_SOURCE_DIR}/../libamigraph/install/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../libamigraph/install/include)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../Cuba-4.2.2/install/include)

#set_target_properties(${EXEC_NAME} PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
# Use ALPSCore_LIBRARIES variable to link to ALPSCore 
#target_link_libraries(${PROJECT_NAME} ${ALPSCore_LIBRARIES})

#NOTE: The relative paths in libami and libamigraph compilation don't allow cmake to install. fix later

	
target_link_libraries(hf_cuba_mpi ${CMAKE_CURRENT_SOURCE_DIR}/../../../libami/install/lib/libami.so ${CMAKE_CURRENT_SOURCE_DIR}/../libamigraph/install/lib/libamigraph.so ${CMAKE_CURRENT_SOURCE_DIR}/../libamicalc/install/lib/libamicalc.so  m stdc++fs Eigen3::Eigen ${CMAKE_CURRENT_SOURCE_DIR}/../../../Cuba-4.2.2/install/lib/libcuba.a boost_system boost_filesystem)
install(TARGETS ${EXEC_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin )
#install (TARGETS ${EXEC_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin )
#install(TARGETS dmft RUNTIME DESTINATION "bin/")
