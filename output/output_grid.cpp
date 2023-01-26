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
        cout<< "ERROR!\nUsage: ./output_grid.exe mufrac = mu/Q m nf channel calls filename\nExiting..." <<endl;
        return -1;
    }

    int nf = atoi(argv[3]) ;
    string channel = argv[4];
    int calls = atoi(argv[5]) ;

    string filename = argv[6] ;
    ifstream inputQ;
    inputQ.open("Q.txt");

    ifstream inputx;
    inputx.open("x.txt");

    double mufrac = atof(argv[1]);
    double m = atof(argv[2]);

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
    //for(std::vector<Datum>::iterator d = data.begin(); d != data.end(); d++) {
    for(double Q_ : Q) {
        time_t t1 = time(NULL) ;
        for(double x_ : x) {
            double mQ, mu, mMu ;
            mQ = pow(m/Q_, 2) ;
            mu = mufrac * Q_ ;
            mMu = pow(m/mu, 2) ;
            if(channel == "2g") res = C2m_g3_approximation(x_,mQ,mMu,nf,1,calls) ;
            else if(channel == "2q") res = C2m_ps3_approximation(x_,mQ,mMu,nf) ;
            else if(channel == "Lg") res = CLm_g3_approximation(x_,mQ,mMu,nf,1,calls)  ;
            else if(channel == "Lq") res = CLm_ps3_approximation(x_,mQ,mMu,nf) ;
            else {
                cout<< "ERROR!\nUsage: channel = {2g, 2q, Lg, Lq}\nExiting..." <<endl;
                exit(-1);
            }
            output << res << "   ";
        }
        output << endl ;
        time_t t2 = time(NULL) ;
        cout << "Time for evaluating Q = "<< Q_ <<" is " << print_time(t2 - t1) << endl ;
        total += t2 - t1 ;
        mean = total / k ;
        cout << "Expected remaining time is " << print_time((Q.size() - k) * mean) << endl ;
        k++ ;
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
