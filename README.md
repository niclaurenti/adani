# ADANI

ADANI (Approximate DIS At N3LO Implementation) is a C++ code that computes an approximation for the DIS coefficient functions at N3LO in heavy quark production, that are not fully known yet.

## Dependencies
The codes depend on ```apfelxx``` and ```gsl```. In order to build also the Python bindings, also ```pybind11``` and the python module ```scikit-build``` are required.

## Installation
```bash
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/your/installation/path/ ..
make && make install
```
This will install both the C++ libary and the Python bindings. In order to build just the Python bindings run
```bash
pip install .
```