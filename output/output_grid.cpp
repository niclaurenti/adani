#include "/Users/niccololaurenti/Master-thesis/include/MyLib.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

int main(int argc, char** argv) {
	
	if(argc!=4) {
		cout<< "ERROR!\nUsage: ./output_grid.exe mufrac = mu/Q m channel\nExiting..." <<endl;
		return -1;
	}
	
	string channel = argv[3];
	string filename;
	if(channel == "2g") filename = "C2g.dat";
	else if(channel == "2q") filename = "C2q.dat";
	else if(channel == "Lg") filename = "CLg.dat";
	else if(channel == "Lq") filename = "CLq.dat";	
	else {
		cout<< "ERROR!\nUsage: channel = {2g, 2q, Lg, Lq}\nExiting..." <<endl;
		return -1;
	}


	ofstream output;
	output.open(filename);

	ifstream inputQ;
	inputQ.open("Q.txt");

	ifstream inputx;
	inputx.open("x.txt");
	
    double mufrac = atof(argv[1]);
    double m = atof(argv[2]);

	double res ;
	
	int nf=4;

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
	
    //for(std::vector<Datum>::iterator d = data.begin(); d != data.end(); d++) {
    for(double Q_ : Q) {
		for(double x_ : x) {
			double mQ, mu, mMu ;
			mQ = pow(m/Q_, 2) ;
			mu = mufrac * Q_ ;
			mMu = pow(m/mu, 2) ;
			if(channel == "2g") res = C2m_g3_approximation(x_,mQ,mMu,nf) ;
			if(channel == "2q") res = C2m_ps3_approximation(x_,mQ,mMu,nf) ;
			if(channel == "Lg") res = CLm_g3_approximation(x_,mQ,mMu,nf)  ;
			if(channel == "Lq") res = CLm_ps3_approximation(x_,mQ,mMu,nf) ;
			cout <<"Channel = "<< channel<< "   Q = " << Q_ << "   x = " << x_ << " res = "<< res << endl ;	
		}
		output << endl ;	
	}
	output.close();
	inputQ.close() ;
	inputx.close() ;
	
	return 0;

}