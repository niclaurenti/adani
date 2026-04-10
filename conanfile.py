from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, CMakeDeps, cmake_layout
import subprocess

class AdaniConan(ConanFile):
    name = "adani"
    package_type = "library"

    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "PYTHON_BUILD": [False]}
    default_options = {"shared": True, "PYTHON_BUILD": False}

    generators = "CMakeDeps"

    exports_sources = "CMakeLists.txt", "cmake/*", "inc/adani/*", "src/*"

    def requirements(self):
        self.requires("gsl/2.7.1", transitive_headers=True, transitive_libs=True)

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.cache_variables["CMAKE_EXPORT_COMPILE_COMMANDS"] = True
        tc.generate()

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
        self.version = version

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_id(self):
        del self.info.settings.compiler

    def package_info(self):
        self.cpp_info.set_property("cmake_file_name", "adani")
        self.cpp_info.set_property("cmake_target_name", "adani::adani")
        self.cpp_info.libs = ["adani"]
        self.cpp_info.include_dirs = ["include", "include/adani"]
        # self.cpp_info.requires = ["gsl::gsl"]
