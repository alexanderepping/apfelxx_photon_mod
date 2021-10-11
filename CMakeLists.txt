cmake_minimum_required(VERSION 3.10)
project(apfelxx_photon_mod)

# use this for normal compiling
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 -O2")
enable_language(Fortran)
# use this for debugging
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 -g"))
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g")


#prepare LHAPDF for linking
set(LHAPDF_CONFIG "lhapdf-config")
exec_program(${LHAPDF_CONFIG} ARGS --libs OUTPUT_VARIABLE LHAPDF_LIB)
exec_program(${LHAPDF_CONFIG} ARGS --incdir OUTPUT_VARIABLE LHAPDF_INC)

include_directories(
    ./myApfelModified/apfelxx/inc
    ${LHAPDF_INC}
)

file(GLOB_RECURSE LIB_APFELXX_FILES "myApfelModified/apfelxx/src/*")

#compile APFELXX and store the lib into ./lib
add_library(apfelxx SHARED ${LIB_APFELXX_FILES})
set_target_properties(apfelxx PROPERTIES LIBRARY_OUTPUT_DIRECTORY ./lib)

#example for compiling executable
add_executable(evolution "${CMAKE_CURRENT_SOURCE_DIR}/testModifiedApfel/Evolution.cc")
target_link_libraries(evolution apfelxx ${LHAPDF_LIB})
set_target_properties(evolution PROPERTIES RUNTIME_OUTPUT_DIRECTORY ./bin)