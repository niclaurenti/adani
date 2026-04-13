from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, CMakeDeps, cmake_layout
from conan.tools.files import copy
import subprocess
import os

class AdaniConan(ConanFile):
    name = "adani"
    package_type = "library"

    # Version is set dynamically from git
    def set_version(self):
        try:
            version = subprocess.check_output(
                ["git", "describe", "--tags", "--dirty", "--always"],
                text=True
            ).strip()
            if version.startswith("v"):
                version = version[1:]
        except Exception:
            version = "0.0.0"
        self.version = "0.0.0"

    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fPIC": [True, False], "PYTHON_BUILD": [False]}
    default_options = {"shared": True, "fPIC": True, "PYTHON_BUILD": False}

    exports_sources = (
        "CMakeLists.txt",
        "cmake/*",
        "inc/adani/*",
        "src/*"
    )

    def requirements(self):
        self.requires("gsl/2.7.1")

    def layout(self):
        cmake_layout(self, src_folder=".")

    def generate(self):
        tc = CMakeToolchain(self)
        tc.cache_variables["CMAKE_EXPORT_COMPILE_COMMANDS"] = True
        tc.cache_variables["CMAKE_INSTALL_RPATH"] = "@loader_path"
        tc.cache_variables["CMAKE_BUILD_RPATH"] = "@loader_path"
        tc.generate()
        deps = CMakeDeps(self)
        deps.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.set_property("cmake_file_name", "adani")
        self.cpp_info.set_property("cmake_target_name", "adani::adani")
        self.cpp_info.libs = ["adani"]
        self.cpp_info.requires = ["gsl::gsl"]
