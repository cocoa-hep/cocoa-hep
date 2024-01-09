# - Locate FastJet library
# Defines:
#
#  FASTJET_FOUND
#  FASTJET_INCLUDE_DIR
#  FASTJET_INCLUDE_DIRS (not cached)
#  FASTJET_LIBRARY
#  FASTJET_LIBRARIES (not cached)
#  FASTJET_LIBRARY_DIRS (not cached)

find_path(FASTJET_INCLUDE_DIR fastjet/version.hh
          HINTS $ENV{FASTJET_HOME}/include ${FASTJET_HOME}/include /usr/include /usr/local/include $ENV{FASTJET_ROOT_DIR}/include ${FASTJET_ROOT_DIR}/include)

find_library(FASTJET_LIBRARY NAMES fastjet
             HINTS $ENV{FASTJET_HOME}/lib ${FASTJET_HOME}/lib /usr/local/lib /usr/lib /usr/lib64 $ENV{FASTJET_ROOT_DIR}/lib ${FASTJET_ROOT_DIR}/lib)

find_library(FASTJETPLUGINS_LIBRARY NAMES fastjetplugins
             HINTS $ENV{FASTJET_HOME}/lib ${FASTJET_HOME}/lib /usr/local/lib /usr/lib /usr/lib64 $ENV{FASTJET_ROOT_DIR}/lib ${FASTJET_ROOT_DIR}/lib)


find_library(FASTJETTOOLS_LIBRARY NAMES fastjettools
             HINTS $ENV{FASTJET_HOME}/lib ${FASTJET_HOME}/lib /usr/local/lib /usr/lib /usr/lib64 $ENV{FASTJET_ROOT_DIR}/lib ${FASTJET_ROOT_DIR}/lib)

# handle the QUIETLY and REQUIRED arguments and set FASTJET_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FastJet DEFAULT_MSG FASTJET_INCLUDE_DIR FASTJET_LIBRARY)

mark_as_advanced(FASTJET_FOUND FASTJET_INCLUDE_DIR FASTJET_LIBRARY FASTJETPLUGINS_LIBRARY FASTJETTOOLS_LIBRARY )

set(FASTJET_INCLUDE_DIRS ${FASTJET_INCLUDE_DIR})
set(FASTJET_LIBRARIES ${FASTJET_LIBRARY} ${FASTJETPLUGINS_LIBRARY} ${FASTJETTOOLS_LIBRARY} )
get_filename_component(FASTJET_LIBRARY_DIRS ${FASTJET_LIBRARY} PATH)

if(FASTJET_FOUND AND FASTJET_LIBRARY AND NOT TARGET fastjet::fastjet)
    add_library(fastjet::fastjet UNKNOWN IMPORTED)
    set_target_properties(fastjet::fastjet PROPERTIES
        IMPORTED_LOCATION "${FASTJET_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${FASTJET_INCLUDE_DIRS}"
    )
endif()

if(FASTJET_FOUND AND FASTJETPLUGINS_LIBRARY AND NOT TARGET fastjet::fastjetplugins)
    add_library(fastjet::fastjetplugins UNKNOWN IMPORTED)
    set_target_properties(fastjet::fastjetplugins PROPERTIES
        IMPORTED_LOCATION "${FASTJETPLUGINS_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${FASTJET_INCLUDE_DIRS}"
    )
endif()

if(FASTJET_FOUND AND FASTJETTOOLS_LIBRARY AND NOT TARGET fastjet::fastjettools)
    add_library(fastjet::fastjettools UNKNOWN IMPORTED)
    set_target_properties(fastjet::fastjettools PROPERTIES
        IMPORTED_LOCATION "${FASTJETTOOLS_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${FASTJET_INCLUDE_DIRS}"
    )
endif()


