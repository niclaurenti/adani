#source ../lib.sh
clang++ -Wall -I/Users/niccololaurenti/Master-thesis/include -L/Users/niccololaurenti/Master-thesis/ -L/Users/niccololaurenti/.local/lib/ -o damptest.exe damptest.cpp -lMyLib -lapfelxx -std=c++17 -stdlib=libc++ `apfelxx-config --cppflags --ldflags --cxxflags`

for i in 2 5 10 15 20 25 30 50 100 200 500 800 1000 10000 100000;
do ./damptest.exe $i fileeta.dat
gnuplot -e "xi=$i; filename='fileeta.dat' " damptest.gnp ;
done


gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=approximation_g_N2LO_FL.pdf damptest_g_N2LO_vs_eta_xi=2.pdf damptest_g_N2LO_vs_eta_xi=5.pdf damptest_g_N2LO_vs_eta_xi=10.pdf damptest_g_N2LO_vs_eta_xi=15.pdf damptest_g_N2LO_vs_eta_xi=20.pdf damptest_g_N2LO_vs_eta_xi=25.pdf damptest_g_N2LO_vs_eta_xi=30.pdf damptest_g_N2LO_vs_eta_xi=50.pdf damptest_g_N2LO_vs_eta_xi=100.pdf damptest_g_N2LO_vs_eta_xi=200.pdf damptest_g_N2LO_vs_eta_xi=500.pdf damptest_g_N2LO_vs_eta_xi=800.pdf damptest_g_N2LO_vs_eta_xi=1000.pdf damptest_g_N2LO_vs_eta_xi=10000.pdf damptest_g_N2LO_vs_eta_xi=100000.pdf
rm damptest_g_N2LO_vs_eta_xi=2.pdf
rm damptest_g_N2LO_vs_eta_xi=5.pdf
rm damptest_g_N2LO_vs_eta_xi=10.pdf
rm damptest_g_N2LO_vs_eta_xi=15.pdf
rm damptest_g_N2LO_vs_eta_xi=20.pdf
rm damptest_g_N2LO_vs_eta_xi=25.pdf
rm damptest_g_N2LO_vs_eta_xi=30.pdf
rm damptest_g_N2LO_vs_eta_xi=50.pdf
rm damptest_g_N2LO_vs_eta_xi=100.pdf
rm damptest_g_N2LO_vs_eta_xi=200.pdf
rm damptest_g_N2LO_vs_eta_xi=500.pdf
rm damptest_g_N2LO_vs_eta_xi=800.pdf
rm damptest_g_N2LO_vs_eta_xi=1000.pdf
rm damptest_g_N2LO_vs_eta_xi=10000.pdf
rm damptest_g_N2LO_vs_eta_xi=100000.pdf

gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=approximation_ps_N2LO_FL.pdf damptest_ps_N2LO_vs_eta_xi=2.pdf damptest_ps_N2LO_vs_eta_xi=5.pdf damptest_ps_N2LO_vs_eta_xi=10.pdf damptest_ps_N2LO_vs_eta_xi=15.pdf damptest_ps_N2LO_vs_eta_xi=20.pdf damptest_ps_N2LO_vs_eta_xi=25.pdf damptest_ps_N2LO_vs_eta_xi=30.pdf damptest_ps_N2LO_vs_eta_xi=50.pdf damptest_ps_N2LO_vs_eta_xi=100.pdf damptest_ps_N2LO_vs_eta_xi=200.pdf damptest_ps_N2LO_vs_eta_xi=500.pdf damptest_ps_N2LO_vs_eta_xi=800.pdf damptest_ps_N2LO_vs_eta_xi=1000.pdf damptest_ps_N2LO_vs_eta_xi=10000.pdf damptest_ps_N2LO_vs_eta_xi=100000.pdf
rm damptest_ps_N2LO_vs_eta_xi=2.pdf
rm damptest_ps_N2LO_vs_eta_xi=5.pdf
rm damptest_ps_N2LO_vs_eta_xi=10.pdf
rm damptest_ps_N2LO_vs_eta_xi=15.pdf
rm damptest_ps_N2LO_vs_eta_xi=20.pdf
rm damptest_ps_N2LO_vs_eta_xi=25.pdf
rm damptest_ps_N2LO_vs_eta_xi=30.pdf
rm damptest_ps_N2LO_vs_eta_xi=50.pdf
rm damptest_ps_N2LO_vs_eta_xi=100.pdf
rm damptest_ps_N2LO_vs_eta_xi=200.pdf
rm damptest_ps_N2LO_vs_eta_xi=500.pdf
rm damptest_ps_N2LO_vs_eta_xi=800.pdf
rm damptest_ps_N2LO_vs_eta_xi=1000.pdf
rm damptest_ps_N2LO_vs_eta_xi=10000.pdf
rm damptest_ps_N2LO_vs_eta_xi=100000.pdf
