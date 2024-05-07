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
  set(Tauola++_FIND_COMPONENTS Fortran CxxInterface)
endif()
set(Tauola++_FOUND TRUE)
foreach(component ${Tauola++_FIND_COMPONENTS})
  find_library(Tauola++_${component}_LIBRARY NAMES Tauola${component}
               HINTS ${Tauola++_ROOT_DIR}/lib $ENV{TAUOLAPP_ROOT_DIR}/lib ${TAUOLAPP_ROOT_DIR}/lib
                     ${Tauola++_ROOT_DIR}/lib64 $ENV{TAUOLAPP_ROOT_DIR}/lib64 ${TAUOLAPP_ROOT_DIR}/lib64
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
  message(STATUS "EvtGen: Tauola++_${component}_FOUND=${Tauola++_${component}_FOUND}, Tauola++_FIND_REQUIRED_${component}=${Tauola++_FIND_REQUIRED_${component}}  in ${Tauola++_${component}_LIBRARY}")
  mark_as_advanced(Tauola++_${component}_LIBRARY)
endforeach()

if(Tauola++_LIBRARY_DIRS)
  list(REMOVE_DUPLICATES Tauola++_LIBRARY_DIRS)
endif()

find_path(Tauola++_INCLUDE_DIR Tauola/Tauola.h
          HINTS ${Tauola++_ROOT_DIR}/include
                $ENV{TAUOLAPP_ROOT_DIR}/include ${TAUOLAPP_ROOT_DIR}/include)
set(Tauola++_INCLUDE_DIRS ${Tauola++_INCLUDE_DIR})
mark_as_advanced(Tauola++_INCLUDE_DIR)

# Make a variable that contains all compiler flags related to Tauola++
set(Tauola++_LIBRARIES ${Tauola++_CxxInterface_LIBRARY} ${Tauola++_Fortran_LIBRARY} ${Tauola++_HEPEVT_LIBRARY})
mark_as_advanced(Tauola++_LIBRARIES)

# handle the QUIETLY and REQUIRED arguments and set Tauola++_FOUND to TRUE if
# all listed variables are TRUE
if (Tauola++_FOUND)
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Tauola++ DEFAULT_MSG Tauola++_INCLUDE_DIR Tauola++_LIBRARIES)
endif()
mark_as_advanced(Tauola++_FOUND)