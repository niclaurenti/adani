from conan import ConanFile
from conan.tools.cmake import CMake, cmake_layout

class TestConan(ConanFile):
    name = "test"
    version = "1.0"

    settings = "os", "compiler", "build_type", "arch"

    # Your test depends on the adani library
    requires = "adani/v1.0.6-7-g4ba624c-dirty"  # or a specific version

    generators = "CMakeDeps", "CMakeToolchain"

    exports_sources = "CMakeLists.txt", "test.cpp"

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
        self.cpp_info.libs = ["test"]
