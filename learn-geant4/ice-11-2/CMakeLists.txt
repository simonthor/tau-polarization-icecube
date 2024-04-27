#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(B1)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)


#----------------------------------------------------------------------------
# Find Tauola (required package)
#
# Ideally, should be done via find_package, e.g.. TODO make this work by writing a FindTauola++.cmake file
# find_package(Tauola++ REQUIRED)
# For now, set the parameters manually. This can also be changed to use the 
# - Try to find Tauola++
# Defines:
#
#  Tauola++_FOUND
#  Tauola++_INCLUDE_DIR
#  Tauola++_INCLUDE_DIRS (not cached)
#  Tauola++_<component>_LIBRARY
#  Tauola++_<component>_FOUND
#  Tauola++_LIBRARIES (not cached)
#  Tauola++_LIBRARY_DIRS (not cached)

# Check all if none is explicitly requested
if(NOT Tauola++_FIND_COMPONENTS)
  set(Tauola++_FIND_COMPONENTS Fortran CxxInterface HEPEVT)
endif()
set(Tauola++_FOUND TRUE)
foreach(component ${Tauola++_FIND_COMPONENTS})
  find_library(Tauola++_${component}_LIBRARY NAMES Tauola${component}
               HINTS ${Tauola++_DIR}/lib $ENV{TAUOLAPP_ROOT_DIR}/lib ${TAUOLAPP_ROOT_DIR}/lib
                     ${Tauola++_DIR}/lib64 $ENV{TAUOLAPP_ROOT_DIR}/lib64 ${TAUOLAPP_ROOT_DIR}/lib64
                      )
  if (Tauola++_${component}_LIBRARY)
    set(Tauola++_${component}_FOUND TRUE)
    list(APPEND Tauola++_LIBRARIES ${Tauola++_${component}_LIBRARY})
    get_filename_component(libdir ${Tauola++_${component}_LIBRARY} PATH)
    list(APPEND Tauola++_LIBRARY_DIRS ${libdir})
  else()
    set(Tauola++_${component}_FOUND FALSE)
    if (Tauola++_FIND_REQUIRED_${component})
    set(Tauola++_FOUND FALSE)
    endif()
  endif()
  message(STATUS "Tauola++_${component}_FOUND=${Tauola++_${component}_FOUND}, Tauola++_FIND_REQUIRED_${component}=${Tauola++_FIND_REQUIRED_${component}}  in ${Tauola++_${component}_LIBRARY}")
  mark_as_advanced(Tauola++_${component}_LIBRARY)
endforeach()

if(Tauola++_LIBRARY_DIRS)
  list(REMOVE_DUPLICATES Tauola++_LIBRARY_DIRS)
endif()

find_path(Tauola++_INCLUDE_DIR Tauola/Tauola.h
          HINTS ${Tauola++_DIR}/include
                $ENV{TAUOLAPP_ROOT_DIR}/include ${TAUOLAPP_ROOT_DIR}/include)
set(Tauola++_INCLUDE_DIRS ${Tauola++_INCLUDE_DIR})
mark_as_advanced(Tauola++_INCLUDE_DIR)

# handle the QUIETLY and REQUIRED arguments and set Tauola++_FOUND to TRUE if
# all listed variables are TRUE
if (Tauola++_FOUND)
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Tauola++ DEFAULT_MSG Tauola++_INCLUDE_DIR Tauola++_LIBRARIES)
endif()
mark_as_advanced(Tauola++_FOUND)

set(Tauola++_LIBRARIES ${Tauola++_CxxInterface_LIBRARY} ${Tauola++_Fortran_LIBRARY} ${Tauola++_HEPEVT_LIBRARY})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(exampleB1 exampleB1.cc ${sources} ${headers})
target_link_libraries(exampleB1 ${Geant4_LIBRARIES} ${Tauola++_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEB1_SCRIPTS
  exampleB1.in
  exampleB1.out
  init_vis.mac
  run1.mac
  run2.mac
  vis.mac
  tsg_offscreen.mac
  )

foreach(_script ${EXAMPLEB1_SCRIPTS})
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
add_custom_target(B1 DEPENDS exampleB1)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS exampleB1 DESTINATION bin)
