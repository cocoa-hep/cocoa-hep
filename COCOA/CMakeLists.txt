cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(COCOA LANGUAGES CXX C)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
set(CMAKE_BUILD_TYPE RelWithDebInfo)

include("GNUInstallDirs")
#----------------------------------------------------------------------------
# Options
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
option(COCOA_ENABLE_TESTS "Enables testing" ON)
message(STATUS "COCOA: WITH_GEANT4_UIVIS=${WITH_GEANT4_UIVIS}")
message(STATUS "COCOA: COCOA_ENABLE_TESTS=${COCOA_ENABLE_TESTS}")
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Compiler flags
add_compile_options("$<$<COMPILE_LANG_AND_ID:CXX,AppleClang,Clang>:-Wc++17-extensions;-Wno-shadow>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wall;-Wno-sign-compare>")
message(STATUS "COCOA: CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID}" )
message(STATUS "COCOA: CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS} CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}=${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}" )
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()
message(STATUS "COCOA: Geant4_FOUND=${Geant4_FOUND} Geant4_INCLUDE_DIR=${Geant4_INCLUDE_DiR} Geant4_LIBRARIES=${Geant4_LIBRARIES} Geant4_USE_FILE=${Geant4_USE_FILE}  Geant4_DATASETS=${Geant4_DATASETS}")
include(${Geant4_USE_FILE})
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Find FastJet (required package)
#
find_package(FastJet REQUIRED)
message(STATUS "COCOA: FASTJET_FOUND=${FASTJET_FOUND} FASTJET_INCLUDE_DIR=${FASTJET_INCLUDE_DIR} FASTJET_LIBRARIES=${FASTJET_LIBRARIES}")
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Find JsonCPP (required package)
#
find_package(JsonCpp REQUIRED)
message(STATUS "COCOA: JsonCpp_FOUND=${JsonCpp_FOUND} JsonCpp_INCLUDE_DIR=${JsonCpp_INCLUDE_DIR} JsonCpp_LIBRARIES=${JsonCpp_LIBRARIES}")
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Find HepMC2 (required package)
#
find_package(HepMC2 REQUIRED QUIET)
message(STATUS "COCOA: HepMC2_FOUND=${HepMC2_FOUND} HEPMC2_INCLUDE_DIR=${HEPMC2_INCLUDE_DIR} HEPMC2_LIBRARIES=${HEPMC2_LIBRARIES}")
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS Core Tree RIO Net Hist Graf Physics)
message(STATUS "COCOA: ROOT_FOUND=${ROOT_FOUND} ROOT_INCLUDE_DIRS=${ROOT_INCLUDE_DIRS} ROOT_LIBRARIES=${ROOT_LIBRARIES} ROOT_USE_FILE=${ROOT_USE_FILE}")
include("${ROOT_USE_FILE}")
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Find Pythia8 (required package)
#
find_package(Pythia8 REQUIRED)
message(STATUS "COCOA: PYTHIA8_VERSION=${PYTHIA8_VERSION} PYTHIA8_LIBRARIES=${PYTHIA8_LIBRARIES} PYTHIA8_INCLUDE_DIRS=${PYTHIA8_INCLUDE_DIRS}")
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Sources and headers for this project
#
set(sources 
Cells_data.cc
CalorimeterConstruction.cc
Cell_var.cc
Config_reader_func.cc
Config_reader_var.cc
CSVReader.cc
Debug_Particle_Flow_data.cc
Debug_Particle_Flow_func.cc
Detector_analysis_var.cc
DetectorConstruction.cc
EventAction.cc
Full_trajectory_info_data.cc
Graph_construction_data.cc
GraphConstructor.cc
HEPEvtcom.cc
HepMCG4AsciiReader.cc
HepMCG4AsciiReaderMessenger.cc
HepMCG4Interface.cc
HepMCG4Pythia8Interface.cc
HepMCG4Pythia8Messenger.cc
InnerConstruction.cc
Jet_Builder_data.cc
Jet_Builder_func.cc
OutputRunAction.cc
Particle_flow_data.cc
Particle_flow_func.cc
Particle_flow_var.cc
Pion_for_pflow_var.cc
PrimaryGeneratorAction.cc
PrimaryGeneratorMessenger.cc
ReduceResolution.cc
RunTest.cc
SteppingAction.cc
Supercluster.cc
Superclustering.cc
Superclustering_data.cc
Topo_clust_func.cc
Topo_clusts_data.cc
Topo_clust_var.cc
TrackingAction.cc
Tracking_func.cc
Tracks_data.cc
Track_var.cc
TruthRecordGraph.cc
)
list(TRANSFORM sources PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/")
set( headers CalorimeterConstruction.hh
Cells_data.hh
Cell_var.hh
Config_reader_func.hh
Config_reader_var.hh
CSVReader.hh
Debug_Particle_Flow_data.hh
Debug_Particle_Flow_func.hh
Detector_analysis_var.hh
DetectorConstruction.hh
DetectorGeometryDefinitions.hh
EventAction.hh
Full_trajectory_info_data.hh
Graph_construction_data.hh
GraphConstructor.hh
HepMCG4AsciiReader.hh
HepMCG4AsciiReaderMessenger.hh
HepMCG4Interface.hh
HepMCG4Pythia8Interface.hh
HepMCG4Pythia8Messenger.hh
InnerConstruction.hh
Jet_Builder_data.hh
Jet_Builder_func.hh
OutputRunAction.hh
Particle_flow_data.hh
Particle_flow_func.hh
Particle_flow_var.hh
Pion_for_pflow_var.hh
PrimaryGeneratorAction.hh
PrimaryGeneratorMessenger.hh
ReduceResolution.hh
RunTest.hh
SteppingAction.hh
Supercluster.hh
Superclustering_data.hh
Superclustering.hh
Topo_clust_func.hh
Topo_clusts_data.hh
Topo_clust_var.hh
TrackingAction.hh
Tracking_func.hh
Tracks_data.hh
Track_var.hh
TruthRecordGraph.hh)
list(TRANSFORM headers PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/include/")

#----------------------------------------------------------------------------
#---Create a shared library with geneated dictionary
ROOT_GENERATE_DICTIONARY(G__Event LINKDEF LinkDef.h )
add_library(Event SHARED G__Event.cxx)
target_link_libraries(Event PRIVATE ROOT::Core ROOT::Tree ROOT::Physics)
if(CMAKE_SYSTEM_NAME MATCHES "Linux")
  target_link_options(Event PRIVATE "LINKER:-z,defs")
  target_link_options(Event PRIVATE "LINKER:-z,now")
endif()
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
#---Create COCOA executable
add_executable(COCOA COCOAMain.cc ${sources} ${headers})
target_include_directories(COCOA PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_compile_definitions(COCOA PRIVATE -DG4LIB_USE_PYTHIA8)
target_link_libraries(COCOA PRIVATE
                      ${Geant4_LIBRARIES}
                      HepMC2::HepMC2 HepMC2::HepMC2fio
                      Pythia8::Pythia8
                      fastjet::fastjet 
                      JsonCpp::JsonCpp
                      ROOT::Core ROOT::Tree ROOT::Hist ROOT::Graf ROOT::Physics
                      Event
                      ${CMAKE_DL_LIBS}
                     )
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build COCOA. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.

set(COCOA_SCRIPTS 
    # macro
    init_vis.mac
    vis.mac
    )


add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/tracking_configuration/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/tracking_configuration)

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        ${CMAKE_SOURCE_DIR}/pflow_configuration/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/pflow_configuration)

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            ${CMAKE_SOURCE_DIR}/config/ $<TARGET_FILE_DIR:${PROJECT_NAME}>/config)

foreach(_script ${COCOA_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Install the targets
#
install(TARGETS COCOA DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS Event DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(FILES ${PROJECT_BINARY_DIR}/libEvent_rdict.pcm ${PROJECT_BINARY_DIR}/libEvent.rootmap DESTINATION  ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/config/  DESTINATION ${CMAKE_INSTALL_DATADIR}/COCOA/config/ COMPONENT data)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/pflow_configuration/  DESTINATION ${CMAKE_INSTALL_DATADIR}/COCOA/pflow_configuration/ COMPONENT data)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tracking_configuration/  DESTINATION ${CMAKE_INSTALL_DATADIR}/COCOA/tracking_configuration/ COMPONENT data)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/macro/  DESTINATION ${CMAKE_INSTALL_DATADIR}/COCOA/macro/ COMPONENT data)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/init_vis.mac ${CMAKE_CURRENT_SOURCE_DIR}/vis.mac  DESTINATION ${CMAKE_INSTALL_DATADIR}/COCOA COMPONENT data)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Testing
#
if (COCOA_ENABLE_TESTS)
  ENABLE_TESTING()
  add_subdirectory(test)
endif()
#----------------------------------------------------------------------------

