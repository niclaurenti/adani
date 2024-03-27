# BASIC DEFINITIONS ========================================================================

# define minimum version of cmake
cmake_minimum_required (VERSION 3.15)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
        "MinSizeRel" "RelWithDebInfo")
endif()

# define project name and its language
project(adani)

# define c++ standard and issue all the warning demanded by this standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fPIC")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fPIC -lgsl -lgslcblas -lm")
enable_language(Fortran)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")

if (NOT DEFINED CMAKE_MACOSX_RPATH)
    set(CMAKE_MACOSX_RPATH 0)
endif()

# Configuration script
set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix "${prefix}")
set(includedir "${prefix}/include")
set(libdir "${prefix}/lib")

# export version to file
configure_file(
    "${PROJECT_SOURCE_DIR}/inc/adani/version.h.in"
    "${PROJECT_SOURCE_DIR}/inc/adani/version.h"
)

configure_file(
    "${PROJECT_SOURCE_DIR}/bin/adani-config.in"
    "${PROJECT_SOURCE_DIR}/bin/adani-config"
)
configure_file(
    "${PROJECT_SOURCE_DIR}/bin/adani.pc.in"
    "${PROJECT_SOURCE_DIR}/bin/adani.pc"
)

include_directories(inc)

# find gsl
find_package(GSL REQUIRED)

# FINALIZE ==================================================================================

option(PYTHON_ONLY "Compiles only for python." OFF)

# pybind11
if(PYTHON_ONLY)
    if(SKBUILD)
    message("Found skbuild")
        execute_process(
            COMMAND "${PYTHON_EXECUTABLE}" -c
                    "import pybind11; print(pybind11.get_cmake_dir())"
        OUTPUT_VARIABLE _tmp_dir
        OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ECHO STDOUT)
        list(APPEND CMAKE_PREFIX_PATH "${_tmp_dir}")
    endif(SKBUILD)
    find_package(pybind11 CONFIG REQUIRED)
    pybind11_add_module(
        _core MODULE
        pywrap/pywrap.cc
        src/ApproximateCoefficientFunction.cc
        src/AsymptoticCoefficientFunction.cc
        src/CoefficientFunction.cc
        src/Convolution.cc
        src/ExactCoefficientFunction.cc
        src/HighEnergyCoefficientFunction.cc
        src/HighScaleCoefficientFunction.cc
        src/HighScaleSplitLogs.cc
        src/hplog.f
        src/hqcoef.f
        src/MasslessCoefficientFunction.cc
        src/MatchingCondition.cc
        src/SpecialFunctions.cc
        src/SplittingFunction.cc
        src/ThresholdCoefficientFunction.cc
        src/Value.cc
    )
    target_link_libraries(_core PRIVATE GSL::gsl GSL::gslcblas)
    # target_link_libraries(_core PRIVATE gsl gslcblas)
    target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})
    install(TARGETS _core DESTINATION .)

else(PYTHON_ONLY)
    add_library(
        adani SHARED
        src/ApproximateCoefficientFunction.cc
        src/AsymptoticCoefficientFunction.cc
        src/CoefficientFunction.cc
        src/Convolution.cc
        src/ExactCoefficientFunction.cc
        src/HighEnergyCoefficientFunction.cc
        src/HighScaleCoefficientFunction.cc
        src/HighScaleSplitLogs.cc
        src/hplog.f
        src/hqcoef.f
        src/MasslessCoefficientFunction.cc
        src/MatchingCondition.cc
        src/SpecialFunctions.cc
        src/SplittingFunction.cc
        src/ThresholdCoefficientFunction.cc
        src/Value.cc
    )

    # define libraries to be linked
    target_link_libraries(adani GSL::gsl GSL::gslcblas)
    # target_link_libraries(adani gsl gslcblas)

    # install
    install(FILES ${PROJECT_SOURCE_DIR}/bin/adani-config
    DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
    install(FILES ${PROJECT_SOURCE_DIR}/bin/adani.pc DESTINATION lib/pkgconfig)
    install(DIRECTORY inc/adani DESTINATION include)
    install(TARGETS adani DESTINATION lib)

endif (PYTHON_ONLY)
