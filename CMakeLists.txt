cmake_minimum_required(VERSION 3.23.1)

# add these options before PROJECT keyword
# disallow in-source build:
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(lossfraction)

# compiler requirements
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3 -Wall")

######################
# required libraries #
######################

# gyronimo requirement
set(gyronimo_dir /home/master/boris/gyrofork/gyronimo)
set(gyronimo_build /home/master/boris/gyrofork/build)
find_library(found gyronimo ${gyronimo_build})
  if(found)
    message(STATUS "Found ::gyronimo:: library")
  else()
    message(FATAL_ERROR "::gyronimo:: library not found, update gyronimo_dir path. See readme.md for installation details.")
  endif()
set(lossfraction_includes ${lossfraction_includes} ${gyronimo_dir} ${gyronimo_dir}/misc/apps/include)
set(lossfraction_depending_libs ${lossfraction_depending_libs} gyronimo)
link_directories(${gyronimo_build})

# mpi requirement
find_package(MPI REQUIRED)
  set(lossfraction_includes ${lossfraction_includes} ${MPI_INCLUDE_PATH})
  set(lossfraction_depending_libs ${lossfraction_depending_libs} ${MPI_C_LIBRARIES})

# gsl requirement
find_package(GSL REQUIRED)
  set(lossfraction_includes ${lossfraction_includes} ${GSL_INCLUDE_DIRS})
  set(lossfraction_depending_libs ${lossfraction_depending_libs} ${GSL_LIBRARIES})
  
  # boost requirement
find_package(Boost 1.73.0 REQUIRED)
  set(lossfraction_includes ${lossfraction_includes} ${Boost_INCLUDE_DIRS})

# netcdf requirement
find_program(NCXX4_CONFIG "ncxx4-config")
  if(NCXX4_CONFIG)
    message(STATUS "Found ncxx4-config")
    execute_process(COMMAND "ncxx4-config" --libs OUTPUT_VARIABLE netcdf_cxx_libs)
    execute_process(COMMAND "ncxx4-config" --cflags OUTPUT_VARIABLE netcdf_cxx_flags)
  else()
    message(SEND_ERROR "ncxx4-config not found, install Unidata's netCDF-CXX4")
  endif()
  string(STRIP ${netcdf_cxx_libs} netcdf_cxx_libs)
  string(STRIP ${netcdf_cxx_flags} netcdf_cxx_flags)
  set(lossfraction_includes ${lossfraction_includes} ${netcdf_cxx_flags})
  set(lossfraction_depending_libs ${lossfraction_depending_libs} ${netcdf_cxx_libs})

add_executable(loss-fraction ${CMAKE_SOURCE_DIR}/app/loss-fraction.cc)
target_include_directories(loss-fraction PUBLIC ${lossfraction_includes})
target_link_libraries(loss-fraction PRIVATE ${lossfraction_depending_libs} MPI::MPI_CXX)