export PATH=$PATH:/home/niccolo/external/apfelxx/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/niccolo/external/apfelxx/lib


export PATH=$PATH:/home/niccolo/.local/bin/gsl
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/niccolo/.local/lib

g++ -Wall -I/home/niccolo/.local/include/gsl -L/home/niccolo/.local/lib -c src/Convolutions.cc src/HighEnergyCoefficientFunctions.cc src/MassiveCoefficientFunctions.cc src/HighScaleCoefficientFunctions.cc src/MasslessCoefficientFunctions.cc src/MatchingConditions.cc src/ThresholdCoefficientFunctions.cc src/ApproximateCoefficientFunctions.cc src/AsymptoticCoefficientFunctions.cc src/SpecialFunctions.cc -lgsl -lgslcblas -lm `apfelxx-config --cppflags --ldflags --cxxflags`

ar -r libMyLib.a Convolutions.o HighEnergyCoefficientFunctions.o MassiveCoefficientFunctions.o HighScaleCoefficientFunctions.o MasslessCoefficientFunctions.o MatchingConditions.o  ThresholdCoefficientFunctions.o ApproximateCoefficientFunctions.o AsymptoticCoefficientFunctions.o SpecialFunctions.o

ar tv libMyLib.a
