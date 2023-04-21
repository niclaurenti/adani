# ADANI

ADANI (Approximate DIS At N3LO Implementation) is a C++ code that computes an approximation for the DIS coefficient functions at N3LO in heavy quark production, that are not fully known yet.

## Dependencies

The code depends on the public library ```gsl```.

Optional dependencies are the library ```pybind11``` and the Python module ```scikit-build``` (both public), that are required for building the Pyhton bindings.

## Installation

In order to install the C++ library run
```bash
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/your/installation/path/ ..
make && make install
```

In order to install the Python package run
```bash
pip install adani
```

## Compile a program

In order to compile a simple program run
```bash
g++ -Wall -o test.exe test.cpp -ladani `adani-config --cppflags --ldflags --cxxflags`
```
or
```bash
g++ -Wall -I/your/installation/path/include -L/your/installation/path/lib/ -o test.exe test.cpp -ladani
```
In the first case remember to run
```bash
export PATH=$PATH:/your/installation/path/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/your/installation/path/lib
```
For MacOS users: add the flags ```-std=c++17 -stdlib=libc++```.

## Python module

In order to use the Python module do
```bash
import adani
```
remembering to run
```bash
export PYTHONPATH=$PYTHONPATH:/your/installation/path/lib/python3.X/site-packages/
```
where ```X``` is your Python version.

## Contacts

Contact me at niccolo.laurenti@mi.infn.it.

## Words of our prophet

> La garra charrúa! L'ultima parola agli uruguagi, sempre loro! L'ultima parola nel calcio è la loro: hanno un cuore differente, lo capisci o no? L'artiglio che graffia,
> che lascia il segno nella storia dell'Inter: questa è la storia che si ripete! [...] Il graffio che aveva portato l'Inter in Champions serve per rimarcare il territorio:
> questo è l'Uruguay quando va in campo con tutto se stesso...ecco chi è Vecino: stanco, si, ma lascia in campo Vecino che parla alla fine, lascialo in campo, che la dice
> lui l'ultima cosa nel calcio! Stanno a insegnare cos'è il calcio agli uruguagi, ma vedi un po' te.
>
> 18/09/2018

> Messi! Messi! Leo Messi! Il sinistro migliore del mondo. Da Di Maria a Messi, dalla Bajada alla Perdriel sempre Rosario, la città del calcio. Uno per l'altro. Si sblocca
> la partita. Football! [...] Tutti in piedi per il miglior giocatore del mondo. Rispetto per il numero uno. [...] Rosario, città del calcio, per questo calciatore che troppe
> volte è stato criticato. Anche quando era a Barcellona, quando è arrivato, sofferente, da giovane, non ha mai dimenticato l'argentino. Parla rosarino, sente l'argentino nel
> sangue, ha pianto per la Seleccion ed è lui che la tiene in vita. [...] La mistica che entra in campo. Abbiamo nominato Diego [Armando Maradona] 10 minuti fa. Con Diego dentro
> tutto è possibile ed è col pianto dell'Argentina e gli occhi spiritati del miglior giocatore del mondo. [...] VAMO!
>
> 26/11/2022
