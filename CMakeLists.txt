cmake_minimum_required(VERSION 3.15)
project(adani LANGUAGES CXX Fortran)

# --- Build settings ----------------------------------------------------------

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

file(GLOB_RECURSE SOURCES src/*.cc src/*.f src/*.f90)

option(PYTHON_BUILD "Compiles only for python." OFF)

# pybind11
if(PYTHON_BUILD)
    include_directories(inc)
    find_package(GSL REQUIRED)
    # --- Python Package ----------------------------------------------------------------
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
        ${SOURCES}
    )
    target_link_libraries(_core PRIVATE GSL::gsl)
    target_compile_definitions(_core PRIVATE VERSION_INFO=${PROJECT_VERSION})
    install(TARGETS _core DESTINATION .)

else(PYTHON_BUILD)
    if(DEFINED CONAN_TOOLCHAIN_FILE) # Check if we're being built by Conan
        find_package(GSL REQUIRED)   # Conan package name is lowercase "gsl"
    else() # Assume we're being built by CMake directly, so use the system package
        option(BUILD_SHARED_LIBS "Build shared instead of static" ON)
        find_package(GSL REQUIRED)   # System package name is uppercase "GSL"
    endif()
    # --- C++ Library ----------------------------------------------------------------
    add_library(adani
        ${SOURCES}
    )

    # PUBLIC so consumers inherit GSL include dirs + libs
    target_link_libraries(adani PUBLIC GSL::gsl)

    # PUBLIC include directory for consumers
    target_include_directories(adani PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<INSTALL_INTERFACE:include>
    )

    # --- Install headers ---------------------------------------------------------

    install(DIRECTORY inc/adani/
        DESTINATION include/adani
        FILES_MATCHING PATTERN "*.h")


    # --- Export CMake targets ----------------------------------------------------

    install(
        TARGETS adani
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin
    )

endif (PYTHON_BUILD)
