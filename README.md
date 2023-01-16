# ADANI

ADANI (Approximate DIS At N3LO Implementation) is a C++ code that computes an approximation for the DIS coefficient functions at N3LO in heavy quark production, that are not fully known yet.

## Dependencies
The codes depend on ```apfelxx``` and ```gsl```. In order to build also the Python bindings, also ```pybind11``` and the python module ```scikit-build``` are required.

## Installation
In order to install both the C++ libary and the Python bindings run
```bash
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/your/installation/path/ ..
make && make install
```
In order to build just the Python bindings run
```bash
pip install .
```

## Words of our prophet

> La garra charrúa! L'ultima parola agli uruguagi, sempre loro! L'ultima parola nel calcio è la loro: hanno un cuore differente, lo capisci o no? L'artiglio che graffia,
> che lascia il segno nella storia dell'Inter: questa è la storia che si ripete! [...] Il graffio che aveva portato l'Inter in Champions serve per rimarcare il territorio: 
> questo è l'Uruguay quando va in campo con tutto se stesso...ecco chi è Vecino: stanco, si, ma lascia in campo Vecino che parla alla fine, lascialo in campo, che la dice 
> lui l'ultima cosa nel calcio! Stanno a insegnare cos'è il calcio agli uruguagi, ma vedi un po' te.