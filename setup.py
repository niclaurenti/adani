import sys, os
import versioneer

try:
    from skbuild import setup
except ImportError:
    print(
        "Please update pip, you need pip 10 or greater,\n"
        " or you need to install the PEP 518 requirements in pyproject.toml yourself",
        file=sys.stderr,
    )
    raise

from setuptools import find_packages

# load long description from README
this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="adani",
    version=versioneer.get_version(),
    description="Code computing approximate DIS N3LO coefficients",
    author="Niccolò Laurenti",
    license="AGPLv3",
    packages=find_packages(where="inc"),
    package_dir={"": "inc"},
    cmake_install_dir="inc/adani",
    cmake_args=['-DPYTHON_ONLY:BOOL=ON'],
    python_requires=">=3.8",
    long_description=long_description,
    long_description_content_type="text/markdown",
)
