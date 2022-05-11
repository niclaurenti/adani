clang++ -Wall -I/Users/niccololaurenti/Master-thesis/include -L/Users/niccololaurenti/Master-thesis/ -L/Users/niccololaurenti/.local/lib/ -o approxN3LO.exe approxN3LO.cpp -lMyLib -lapfelxx -std=c++17 -stdlib=libc++ `apfelxx-config --cppflags --ldflags --cxxflags`

for i in 2 5 10 15 20 25 30 50 100 200 500 800 1000 10000 100000;
do ./approxN3LO.exe $i fileeta.dat filex.dat
gnuplot -e "xi=$i; filename1='filex.dat'; filename2='fileeta.dat' " approxN3LO.gnp ;
done


gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=ApproxN3LO_g_FL.pdf Approx_N3LO_g_vs_eta_xi=2.pdf Approx_N3LO_g_vs_eta_xi=5.pdf Approx_N3LO_g_vs_eta_xi=10.pdf Approx_N3LO_g_vs_eta_xi=15.pdf Approx_N3LO_g_vs_eta_xi=20.pdf Approx_N3LO_g_vs_eta_xi=25.pdf Approx_N3LO_g_vs_eta_xi=30.pdf Approx_N3LO_g_vs_eta_xi=50.pdf Approx_N3LO_g_vs_eta_xi=100.pdf Approx_N3LO_g_vs_eta_xi=200.pdf Approx_N3LO_g_vs_eta_xi=500.pdf Approx_N3LO_g_vs_eta_xi=800.pdf Approx_N3LO_g_vs_eta_xi=1000.pdf Approx_N3LO_g_vs_eta_xi=10000.pdf Approx_N3LO_g_vs_eta_xi=100000.pdf
rm Approx_N3LO_g_vs_eta_xi=2.pdf 
rm Approx_N3LO_g_vs_eta_xi=5.pdf 
rm Approx_N3LO_g_vs_eta_xi=10.pdf 
rm Approx_N3LO_g_vs_eta_xi=15.pdf 
rm Approx_N3LO_g_vs_eta_xi=20.pdf 
rm Approx_N3LO_g_vs_eta_xi=25.pdf 
rm Approx_N3LO_g_vs_eta_xi=30.pdf 
rm Approx_N3LO_g_vs_eta_xi=50.pdf 
rm Approx_N3LO_g_vs_eta_xi=100.pdf 
rm Approx_N3LO_g_vs_eta_xi=200.pdf 
rm Approx_N3LO_g_vs_eta_xi=500.pdf 
rm Approx_N3LO_g_vs_eta_xi=800.pdf 
rm Approx_N3LO_g_vs_eta_xi=1000.pdf 
rm Approx_N3LO_g_vs_eta_xi=10000.pdf 
rm Approx_N3LO_g_vs_eta_xi=100000.pdf

gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=ApproxN3LO_ps_FL.pdf Approx_N3LO_ps_vs_eta_xi=2.pdf Approx_N3LO_ps_vs_eta_xi=5.pdf Approx_N3LO_ps_vs_eta_xi=10.pdf Approx_N3LO_ps_vs_eta_xi=15.pdf Approx_N3LO_ps_vs_eta_xi=20.pdf Approx_N3LO_ps_vs_eta_xi=25.pdf Approx_N3LO_ps_vs_eta_xi=30.pdf Approx_N3LO_ps_vs_eta_xi=50.pdf Approx_N3LO_ps_vs_eta_xi=100.pdf Approx_N3LO_ps_vs_eta_xi=200.pdf Approx_N3LO_ps_vs_eta_xi=500.pdf Approx_N3LO_ps_vs_eta_xi=800.pdf Approx_N3LO_ps_vs_eta_xi=1000.pdf Approx_N3LO_ps_vs_eta_xi=10000.pdf Approx_N3LO_ps_vs_eta_xi=100000.pdf
rm Approx_N3LO_ps_vs_eta_xi=2.pdf 
rm Approx_N3LO_ps_vs_eta_xi=5.pdf 
rm Approx_N3LO_ps_vs_eta_xi=10.pdf 
rm Approx_N3LO_ps_vs_eta_xi=15.pdf 
rm Approx_N3LO_ps_vs_eta_xi=20.pdf 
rm Approx_N3LO_ps_vs_eta_xi=25.pdf 
rm Approx_N3LO_ps_vs_eta_xi=30.pdf 
rm Approx_N3LO_ps_vs_eta_xi=50.pdf 
rm Approx_N3LO_ps_vs_eta_xi=100.pdf 
rm Approx_N3LO_ps_vs_eta_xi=200.pdf 
rm Approx_N3LO_ps_vs_eta_xi=500.pdf 
rm Approx_N3LO_ps_vs_eta_xi=800.pdf 
rm Approx_N3LO_ps_vs_eta_xi=1000.pdf 
rm Approx_N3LO_ps_vs_eta_xi=10000.pdf 
rm Approx_N3LO_ps_vs_eta_xi=100000.pdf
