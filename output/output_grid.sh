clang++ -Wall -I/Users/niccololaurenti/Master-thesis/include -L/Users/niccololaurenti/Master-thesis/ -L/Users/niccololaurenti/.local/lib/ -o output_grid.exe output_grid.cpp -lMyLib -lapfelxx -std=c++17 -stdlib=libc++ `apfelxx-config --cppflags --ldflags --cxxflags`
for channel in "2g" "2q" "Lg" "Lq";
do
./output_grid.exe $1 $2 $channel
done
#./output_grid.exe mufrac= mu/Q m filename