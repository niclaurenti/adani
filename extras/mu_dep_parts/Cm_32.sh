clang++ -Wall -I/Users/niccololaurenti/Master-thesis/include -L/Users/niccololaurenti/Master-thesis/ -L/Users/niccololaurenti/.local/lib/ -o Cm_32.exe Cm_32.cpp -lMyLib -lapfelxx -std=c++17 -stdlib=libc++ -lgsl -lgslcblas -lm `apfelxx-config --cppflags --ldflags --cxxflags`

./Cm_32.exe $1 3 prova.dat