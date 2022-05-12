#include "/Users/niccololaurenti/Master-thesis/include/MyLib.h"

#include <iostream>
#include <fstream>
#include <cmath>

#include "apfel/massivecoefficientfunctionsunp_sl.h"

using namespace std;
using namespace apfel;

int main(int argc, char** argv) {

    if(argc!=4) {
		cout << "ERROR\nUsage: ./c2_21.exe xi nf filename\nExiting..." << endl;
		return -1;
	}

    double xi = atof(argv[1]);
    int nf = atoi(argv[2]);
 
	ofstream output;
	output.open(argv[3]);

    double mQ=1/xi;

    double x, dx = 0.001, xmax = 1./(1.+4./xi);

    for(x=dx; x<xmax; x+=dx) {
		//cout << x << "   " << gluon.Regular(x*(1+4*mQ))/ norm << "   " << quark.Regular(x*(1+4*mQ))/ norm << endl;
        cout << x << "   "
               << C2m_g32(x, mQ, nf) << "   " 
               //<< C2m_ps32(x, mQ, nf)<< "   " 
               //<< CLm_g32(x, mQ, nf) << "   " 
               //<< CLm_ps32(x, mQ, nf) 
               << endl;
	}

    output.close();

    return 0;

}