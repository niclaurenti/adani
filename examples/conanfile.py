from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, CMakeDeps, cmake_layout
import subprocess

class TestConan(ConanFile):
    name = "test"
    version = "1.0"
    package_type = "application"
    settings = "os", "compiler", "build_type", "arch"

    exports_sources = "CMakeLists.txt", "test.cpp"

    def requirements(self):
        self.requires("adani/" + self.get_commit())

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.cache_variables["CMAKE_EXPORT_COMPILE_COMMANDS"] = True
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
        self.cpp_info.libs = ["test"]

    def get_commit(self):
        try:
            version = subprocess.check_output(
                ["git", "describe", "--tags", "--dirty", "--always"],
                text=True
            ).strip()
            if version.startswith("v"):
                version = version[1:]
            return version
        except Exception:
            return "0.0.0"
