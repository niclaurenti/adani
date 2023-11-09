// See README.md for compilation suggestions

#include "adani/adani.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <ctime>
#include <sstream>
#include <string.h>
#include <unistd.h>

using namespace std;

std::string print_time(time_t seconds);

int main(int argc, char** argv) {

    if(argc!=7) {
        cout<< "ERROR!\nUsage: ./output_grid.exe mufrac = mu/Q m nf v channel filename\nExiting..." <<endl;
        return -1;
    }

    double mufrac = atof(argv[1]);
    double m = atof(argv[2]);

    int nf = atoi(argv[3]) ;
    int v = atoi(argv[4]) ;
    if (v < -1 || v > 1) {
        cout << "Choose v = {-1, 0, 1}" << endl ;
        exit(-1);
    }
    string channel = argv[5];

    string filename = argv[6] ;
    ifstream inputQ;
    inputQ.open("Q.txt");

    ifstream inputx;
    inputx.open("x.txt");

    std::vector<double> Q, x ;

    double xtmp ;
    inputx >> xtmp ;
    while(!inputx.eof()) {
        x.push_back(xtmp) ;
        inputx >> xtmp ;
    }

    double Qtmp ;
    inputQ >> Qtmp ;
    while(!inputQ.eof()) {
        Q.push_back(Qtmp) ;
        inputQ >> Qtmp ;
    }

    cout    << "Computation of the grid for the coefficient function C"<< channel
            << " for m = " << m << " GeV, nf = "<< nf << " and µ/Q = "<< mufrac
            << endl ;

    cout << "Size of the grid (x,Q) = (" << x.size() <<","<< Q.size() << ")" <<endl ;


    ofstream output;
    output.open(filename);
    if (! output.is_open()) {
        cout<<"Problems in opening "<<filename<<endl ;
        exit(-1);
    } else {
        cout << "Saving grid in " << filename << endl;
    }

    double res ;

    int k  = 1;

    time_t total = 0, mean;

    time_t starting_time = time(NULL) ;

    for(double Q_ : Q) {
        time_t t1 = time(NULL) ;
        for(double x_ : x) {
            double m2Q2, mu, m2mu2 ;
            m2Q2 = pow(m/Q_, 2) ;
            mu = mufrac * Q_ ;
            m2mu2 = pow(m/mu, 2) ;
            if(channel == "2g") res = C2_g3_approximation(x_, m2Q2, m2mu2, nf, v, 1) ;
            else if(channel == "2q") res = C2_ps3_approximation(x_, m2Q2, m2mu2, nf, v) ;
            else if(channel == "Lg") res = CL_g3_approximation(x_, m2Q2, m2mu2, nf, v, 1)  ;
            else if(channel == "Lq") res = CL_ps3_approximation(x_, m2Q2, m2mu2, nf, v) ;
            else {
                cout<< "ERROR!\nUsage: channel = {2g, 2q, Lg, Lq}\nExiting..." <<endl;
                exit(-1);
            }
            output << res << "   ";
        }
        output << endl ;
    }

    time_t ending_time = time(NULL) ;

    cout << "Total running time is " << print_time(ending_time - starting_time) << endl ;

    output.close();
    inputQ.close() ;
    inputx.close() ;

    return 0;

}

std::string print_time(time_t seconds) {

    stringstream ss;

    int hour = (int) seconds / 3600 ;
    int minute = (int) (seconds - hour * 3600) / 60 ;
    int second = seconds - hour * 3600 - minute * 60 ;

    ss << hour <<"h:"<< minute <<"m:"<< second <<"s" ;

    return ss.str();

}
