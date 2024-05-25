
find_path(JsonCpp_INCLUDE_DIR json/json.h
NAMES jsoncpp
HINTS $ENV{JSONCPP_HOME}/include ${JSONCPP_HOME}/include
REQUIRED)
mark_as_advanced(JsonCpp_INCLUDE_DIR)
find_library(JsonCpp_LIBRARY
NAMES jsoncpp
HINTS $ENV{JSONCPP_HOME}/build/lib ${JSONCPP_HOME}/build/lib $ENV{JSONCPP_HOME}/lib64 ${JSONCPP_HOME}/lib64
REQUIRED)
mark_as_advanced(JsonCpp_LIBRARY)


get_filename_component(JsonCpp_LIBRARY_DIR ${JsonCpp_LIBRARY} PATH)

# handle the QUIETLY and REQUIRED arguments and set JsonCpp_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(JsonCpp JsonCpp_INCLUDE_DIR JsonCpp_LIBRARY)

mark_as_advanced(JsonCpp_FOUND JsonCpp_INCLUDE_DIR JsonCpp_LIBRARIES)

# set(JsonCpp_INCLUDE_DIRS ${JsonCpp_INCLUDE_DIR})
set(JsonCpp_LIBRARIES ${JsonCpp_LIBRARY})
# get_filename_component(JsonCpp_LIBRARY_DIRS ${JsonCpp_LIBRARY} PATH)

if(JsonCpp_FOUND AND NOT TARGET JsonCpp::JsonCpp)
    add_library(JsonCpp::JsonCpp UNKNOWN IMPORTED)
    set_target_properties(JsonCpp::JsonCpp PROPERTIES
        IMPORTED_LOCATION "${JsonCpp_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${JsonCpp_INCLUDE_DIR}"
    )
endif()