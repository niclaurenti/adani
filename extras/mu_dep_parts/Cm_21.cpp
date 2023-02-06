#include "adani/adani.h"

#include <iostream>
#include <fstream>
#include <cmath>

#include "apfel/massivecoefficientfunctionsunp_sl.h"

using namespace std;
using namespace apfel;

int main(int argc, char** argv) {

    if(argc!=2) {
        cout << "ERROR\nUsage: ./c2_21.exe xi\nExiting..." << endl;
        return -1;
    }

    double xi = atof(argv[1]);
    //double D= atof(argv[2]);

    double norm = 16 * M_PI * M_PI ;

    double mQ=1/xi;

    Cm22bargNC gluonF2(1./(1+4*mQ));
    Cm22barpsNC quarkF2(1./(1+4*mQ));
    CmL2bargNC gluonFL(1./(1+4*mQ));
    CmL2barpsNC quarkFL(1./(1+4*mQ));

    double x, dx = 0.001;

    for(x=dx; x<1; x+=dx) {
        //cout << x << "   " << gluon.Regular(x*(1+4*mQ))/ norm << "   " << quark.Regular(x*(1+4*mQ))/ norm << endl;
        cout << x << "   " << gluonF2.Regular(x*(1+4*mQ))/ norm << "   " << quarkF2.Regular(x*(1+4*mQ))/ norm << "   " << gluonFL.Regular(x*(1+4*mQ))/ norm << "   " << quarkFL.Regular(x*(1+4*mQ))/ norm << "   " << C2m_g21(x, mQ) << "   " <<C2m_ps21(x, mQ)<< "   " << CLm_g21(x, mQ) << "   " <<CLm_ps21(x, mQ) << endl;
    }

    return 0;

}
