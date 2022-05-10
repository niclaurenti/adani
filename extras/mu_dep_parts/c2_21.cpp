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

    Cm22bargNC gluonF2(1./(1+4*mQ));
    Cm22barpsNC quarkF2(1./(1+4*mQ));
    CmL2bargNC gluonFL(1./(1+4*mQ));
    CmL2barpsNC quarkFL(1./(1+4*mQ));

    double x, dx = 0.001, xmax = 1./(1.+4./xi);

    for(x=dx; x<xmax; x+=dx) {
		//cout << x << "   " << gluon.Regular(x*(1+4*mQ))/ norm << "   " << quark.Regular(x*(1+4*mQ))/ norm << endl;
        output << x << "   " << gluonF2.Regular(x*(1+4*mQ))/ norm << "   " << quarkF2.Regular(x*(1+4*mQ))/ norm << "   " << gluonFL.Regular(x*(1+4*mQ))/ norm << "   " << quarkFL.Regular(x*(1+4*mQ))/ norm << endl;
	}

    output.close();

    return 0;

}