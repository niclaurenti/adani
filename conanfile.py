from conan import ConanFile
from conan.tools.cmake import CMake, cmake_layout
import subprocess

def git_version():
    try:
        return subprocess.check_output(
            ["git", "describe", "--tags", "--dirty", "--always"],
            text=True
        ).strip()
    except Exception:
        return "0.0.0"

class AdaniConan(ConanFile):
    name = "adani"
    version = git_version()
    package_type = "library"

    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True]}
    default_options = {"shared": True}

    requires = "gsl/2.7.1"
    generators = "CMakeDeps", "CMakeToolchain"

    exports_sources = "*"

    def layout(self):
        cmake_layout(self)

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["adani"]
