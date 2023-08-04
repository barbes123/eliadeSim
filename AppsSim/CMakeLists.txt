# $Id: CMakeLists.txt 86065 2014-11-07 08:51:15Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(Apps)

set (CMAKE_CXX_STANDARD 11)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#

set(Geant4_DIR "/home/eliade/geant4/geant4-v11.1.1-install/lib/cmake/Geant4/" CACHE PATH "/home/eliade/geant4/geant4-v11.1.1-install/lib/cmake/Geant4/")

option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)


# Find ROOT (required package)
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

include_directories(${ROOT_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/include)
#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(exampleApps exampleApps.cc ${sources} ${headers})
target_link_libraries(exampleApps ${Geant4_LIBRARIES})
target_link_libraries(exampleApps ${ROOT_LIBRARIES})
target_link_libraries(exampleApps tbb)



#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build Apps. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEApps_SCRIPTS
  exampleApps.in
  exampleApps.out
  init_vis.mac
  run1.mac
  run2.mac
  vis.mac
  )

foreach(_script ${EXAMPLEApps_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(Apps DEPENDS exampleApps)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS exampleApps DESTINATION bin)


