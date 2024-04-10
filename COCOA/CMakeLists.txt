#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(COCOA)
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    add_compile_options(-Wno-shadow)
    add_compile_options(-Wc++17-extensions)
endif()
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)
set(CMAKE_BUILD_TYPE RelWithDebInfo)
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

find_package( FastJet REQUIRED)


#find_package (Acts REQUIRED)

find_package(JsonCpp REQUIRED)
include_directories(${JsonCpp_INCLUDE_DIR})
#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find HepMC (required package)
#
include_directories(${HEPMC_INCLUDE_DIR} ${HEPMC_PYTHIA8INTERFACE_INCLUDE_DIR})
find_package(HepMC REQUIRED)
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net)
include("${ROOT_USE_FILE}")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include_directories(${ROOT_INCLUDE_DIR})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
ROOT_GENERATE_DICTIONARY(G__Event LINKDEF LinkDef.h )# H02RunAction.hh
#---Create a shared library with geneated dictionary
add_library(Event SHARED G__Event.cxx)
target_link_libraries(Event ${ROOT_LIBRARIES})


#----------------------------------------------------------------------------
# Find Pythia8 (required package)
#
find_package(Pythia8 REQUIRED)
if(Pythia8_FOUND)
  message(STATUS "G4 Examples: Pythia8 found. --> COCOA example with Pythia8 enabled.") 
  add_definitions(-DG4LIB_USE_PYTHIA8)
else()
  set(PYTHIA8_LIBRARIES "")
  #set(PYTHIA8_HEPMC_LIBRARY "")
endif()


#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR}
                    ${HEPMC_INCLUDE_DIR}
                    ${HEPMC_PYTHIA8INTERFACE_INCLUDE_DIR}
                    ${PYTHIA8_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIR}
                    ${ROOTSYS}
                    ${FASTJET_INCLUDE_DIR}
                    ${JsonCpp_INCLUDE_DIR}

                    )
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)


#----------------------------------------------------------------------------
#find_package(UUID REQUIRED uuid)
# Add the executable, and link it to the Geant4 libraries
#
add_executable(COCOA COCOAMain.cc ${sources} ${headers})


target_include_directories(COCOA PUBLIC ${UUID_INCLUDE_DIRS})

target_link_libraries(COCOA ${Geant4_LIBRARIES}
                      ${HEPMC_LIBRARIES}
                      ${HEPMC_FIO_LIBRARIES}
                      ${PYTHIA6_LIBRARIES}
                      ${PYTHIA8_LIBRARIES}
                      ${FASTJET_LIBRARY}
                      ${JsonCpp_LIBRARIES}
                      #HepMCdict
                      ${ROOT_LIBRARIES}
                      # jsoncpp_object
                      # jsoncpp_static
                      #HepMCdict
                      Event
                      #clangBasic
                      #clangLex
                      #clangTooling
                      #clangSerialization
                      ${CMAKE_DL_LIBS}
                      ${UUID_LIBRARIES}
                      ${CMAKE_THREAD_LIBS_INIT}
                      #libuuid.so
                      #-llzma
                     )

# if pythia is compiled with g77, link with -lg2c instead.
#target_link_libraries(COCOA ${Geant4_LIBRARIES}
#                      ${HEPMC_LIBRARIES} ${HEPMC_FIO_LIBRARIES}
#                      ${PYTHIA6_LIBRARIES} g2c)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build COCOA. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
# set(COCOA_SCRIPTS 
# #    macro/Pythia8/hepmc_pythia8_qcd.in macro/Pythia8/QuarkMacro.in macro/Pythia8/GluonMacro.in
#     )
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
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS COCOA DESTINATION bin)

