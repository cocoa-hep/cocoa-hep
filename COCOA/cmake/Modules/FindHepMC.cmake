# - Locate HepMC library
# in a directory defined via  HEPMC_HOME or HEPMC_DIR environment variable
# Defines:
#
#  HEPMC_FOUND
#  HEPMC_INCLUDE_DIR
#  HEPMC_PYTHIA8INTERFACE_INCLUDE_DIR
#  HEPMC_INCLUDE_DIRS (not cached)
#  HEPMC_LIBRARIES
#  HEPMC_FIO_LIBRARIES

find_path(HEPMC_INCLUDE_DIR HepMC3/GenEvent.h
          HINTS $ENV{HEPMC_HOME}/include ${HEPMC_HOME}/include
          $ENV{HEPMC_DIR}/include ${HEPMC_DIR}/include)

find_path(HEPMC_PYTHIA8INTERFACE_INCLUDE_DIR Pythia8/Pythia8ToHepMC3.h
          HINTS $ENV{HEPMC_HOME}/share/HepMC3/interfaces/pythia8/include/ ${HEPMC_HOME}/share/HepMC3/interfaces/pythia8/include/
          $ENV{HEPMC_DIR}/share/HepMC3/interfaces/pythia8/include/ ${HEPMC_DIR}/share/HepMC3/interfaces/pythia8/include/)

find_library(HEPMC_LIBRARIES NAMES HepMC3
             HINTS $ENV{HEPMC_HOME}/lib ${HEPMC_HOME}/lib
             HINTS $ENV{HEPMC_DIR}/lib ${HEPMC_DIR}/lib)

get_filename_component(HEPMC_LIBRARY_DIR ${HEPMC_LIBRARIES} PATH)
set(HEPMC_FIO_LIBRARIES "-L${HEPMC_LIBRARY_DIR}" -lHepMC3 -lHepMC3rootIO)

set(HEPMC_INCLUDE_DIRS ${HEPMC_INCLUDE_DIR} ${HEPMC_PYTHIA8INTERFACE_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set HEPMC_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(HepMC DEFAULT_MSG HEPMC_INCLUDE_DIR HEPMC_PYTHIA8INTERFACE_INCLUDE_DIR HEPMC_LIBRARIES)

mark_as_advanced(HEPMC_FOUND HEPMC_INCLUDE_DIR HEPMC_PYTHIA8INTERFACE_INCLUDE_DIR HEPMC_LIBRARIES)
