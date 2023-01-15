# BASIC DEFINITIONS ========================================================================

# define minimum version of cmake
cmake_minimum_required (VERSION 3.5)

if (NOT DEFINED CMAKE_MACOSX_RPATH)
   set(CMAKE_MACOSX_RPATH 0)
endif()

# define project name and its language
project(masterthesis)

# define c++ standard and issue all the warning demanded by this standard
# set(CMAKE_BUILD_TYPE Release)
set(masterthesis_VERSION 0.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# export version to file
configure_file(
  "${PROJECT_SOURCE_DIR}/inc/masterthesis/version.h.in"
  "${PROJECT_SOURCE_DIR}/inc/masterthesis/version.h"
)

# Configuration script
set(prefix ${CMAKE_INSTALL_PREFIX})
set(includedir "${prefix}/include")
set(libdir "${prefix}/lib")
configure_file("${PROJECT_SOURCE_DIR}/bin/masterthesis-config.in" "${PROJECT_BINARY_DIR}/bin/masterthesis-config")
include_directories(${PROJECT_SOURCE_DIR}/inc)

# APFEL++
find_program(APFELXX_CONFIG apfelxx-config REQUIRED)
if (APFELXX_CONFIG)
  exec_program(${APFELXX_CONFIG}
    ARGS --cppflags
    OUTPUT_VARIABLE APFELXX_CXX_FLAGS
  )
  set(APFELXX_CXX_FLAGS ${APFELXX_CXX_FLAGS} CACHE STRING INTERNAL)
  exec_program(${APFELXX_CONFIG}
    ARGS --ldflags
    OUTPUT_VARIABLE APFELXX_LIBRARIES
  )
  set(APFELXX_LIBRARIES ${APFELXX_LIBRARIES} CACHE STRING INTERNAL)
else(APFELXX_CONFIG)
  message(FATAL_ERROR "APFEL++ not found!")
endif(APFELXX_CONFIG)

# GSL
# find_program(GSL_CONFIG gsl-config)
# if (GSL_CONFIG)
#   exec_program(${GSL_CONFIG}
#     ARGS --cflags
#     OUTPUT_VARIABLE GSL_CXX_FLAGS
#   )
#   set(GSL_CXX_FLAGS ${GSL_CXX_FLAGS} CACHE STRING INTERNAL)
#   exec_program(${GSL_CONFIG}
#     ARGS --libs
#     OUTPUT_VARIABLE GSL_LIBRARIES
#   )
#   set(GSL_LIBRARIES ${GSL_LIBRARIES} CACHE STRING INTERNAL)
# else(GSL_CONFIG)
#   message ("GSL not found. Some parts of the code will not be compiled.")
# endif(GSL_CONFIG)
# message("${GSL_CXX_FLAGS}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fPIC -std=c++17 -stdlib=libc++ -lgsl -lgslcblas -lm `apfelxx-config --cppflags --ldflags --cxxflags`")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GSL_CXX_FLAGS} ${APFELXX_CXX_FLAGS} -Wall -Wextra -fPIC")

find_package(GSL REQUIRED)
find_library(APFELXX apfelxx)

# FINALIZE ==================================================================================

# generate list of source files
file(GLOB_RECURSE  source_files src/*)

# define target library
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

add_library(masterthesis SHARED ${source_files})

# define libraries to be linked
target_link_libraries(masterthesis GSL::gsl GSL::gslcblas "${APFELXX}")

# pybind11
find_package(pybind11 QUIET)
if (pybind11_FOUND)
   pybind11_add_module(masterthesispy pywrap/pywrap.cc)
   target_link_libraries(masterthesispy PRIVATE masterthesis)
   find_package (Python3 REQUIRED)
   execute_process(COMMAND "${Python3_EXECUTABLE}" -c "import distutils.sysconfig; print(distutils.sysconfig.get_python_lib(prefix='${CMAKE_INSTALL_PREFIX}', plat_specific=True), end = '')" OUTPUT_VARIABLE PYTHON_INSTALLATION_PATH)
   install(TARGETS masterthesispy LIBRARY DESTINATION "${PYTHON_INSTALLATION_PATH}")
else (pybind11_FOUND)
   message ("-- pybind11 not found. The python wrapper will not be built.")
endif (pybind11_FOUND)

# install
install(TARGETS masterthesis
           RUNTIME DESTINATION bin
           LIBRARY DESTINATION lib)
install(DIRECTORY inc/masterthesis DESTINATION include)
install(FILES ${PROJECT_BINARY_DIR}/bin/masterthesis-config DESTINATION bin
PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
