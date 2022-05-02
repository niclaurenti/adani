#include "/Users/niccololaurenti/Master-thesis/include/MyLib.h"

#include <iostream>
#include <fstream>
#include <cmath>

#include "apfel/massivecoefficientfunctionsunp_sl.h"

using namespace std;
using namespace apfel;

int main(int argc, char** argv) {

    if(argc!=3) {
		cout << "ERROR\nUsage: ./c2_21.exe xi filename\nExiting..." << endl;
		return -1;
	}

    double xi = atof(argv[1]);
    //double D= atof(argv[2]);
	ofstream output;
	output.open(argv[2]);

    double norm = 16 * M_PI * M_PI ;

    double mQ=1/xi;

    Cm22bargNC gluon(1./(1+4*mQ));
    CmL2barpsNC quarkL(1./(1+4*mQ));

    double x, dx = 0.001, xmax = 1./(1.+4./xi);

    for(x=dx; x<xmax; x+=dx) {
		cout << x << "   " << gluon.Regular(x*(1+4*mQ))/ norm << "   " << quarkL.Regular(x*(1+4*mQ))/ norm << endl;
        output << x <<  "   " << gluon.Regular(x*(1+4*mQ)) / norm << "   " << quarkL.Regular(x*(1+4*mQ))/ norm  << endl;		
	}

    output.close();

    return 0;

}