from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, CMakeDeps, cmake_layout

class TestConan(ConanFile):
    name = "test"
    version = "1.0"
    package_type = "application"
    settings = "os", "compiler", "build_type", "arch"

    exports_sources = "CMakeLists.txt", "test.cpp"

    def requirements(self):
        # Your test depends on the adani library
        self.requires("adani/v1.0.6-8-gaafc7cb-dirty", transitive_headers=True, transitive_libs=True)

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
