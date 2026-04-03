cmake_minimum_required(VERSION 3.15)
project(adani LANGUAGES CXX Fortran)

# --- Build settings ----------------------------------------------------------

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Conan will inject include paths, so no manual include_directories()
# include_directories(inc)  <-- REMOVE THIS

# --- Dependencies ------------------------------------------------------------

find_package(gsl REQUIRED)   # Conan package name is lowercase "gsl"

# --- Library ----------------------------------------------------------------

add_library(adani SHARED
    src/ApproximateCoefficientFunction.cc
    src/AQg3mod.f90
    src/AQg3red.f
    src/AsymptoticCoefficientFunction.cc
    src/CoefficientFunction.cc
    src/CommonTypes.cc
    src/Convolution.cc
    src/ExactCoefficientFunction.cc
    src/Exceptions.cc
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

# PUBLIC so consumers inherit GSL include dirs + libs
target_link_libraries(adani PUBLIC gsl::gsl)

# PUBLIC include directory for consumers
target_include_directories(adani PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
    $<INSTALL_INTERFACE:include>
)

# --- Install headers ---------------------------------------------------------

install(DIRECTORY inc/adani DESTINATION include)

# --- Export CMake targets ----------------------------------------------------

install(
    TARGETS adani
    EXPORT adaniTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
)

install(
    EXPORT adaniTargets
    FILE adaniTargets.cmake
    NAMESPACE adani::
    DESTINATION lib/cmake/adani
)

# --- Config file -------------------------------------------------------------

include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/adaniConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/adaniConfig.cmake
    INSTALL_DESTINATION lib/cmake/adani
)

install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/adaniConfig.cmake
    DESTINATION lib/cmake/adani
)
