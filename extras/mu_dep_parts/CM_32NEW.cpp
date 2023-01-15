#include "masterthesis/MassiveCoefficientFunctions.h"

#include <iostream>
#include <iomanip>
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

    double mQ=1./xi;

    double norm = 4. * xi ;

    double x;
    //double xmax = 1./(1.+4./xi);

    double eta, logeta, logeta_min=-4, logeta_max=4;	
	
	int N=500;
	double dlog=(logeta_max - logeta_min)/N;

    for(int i=0; i<N; i++) {
		
		logeta=logeta_min + i*dlog;
		
		eta=pow(10, logeta);
        //if (eta == 1) eta = (double)eta ;
		
		x=1/(1+4*mQ*(eta+1));

        output << eta << "   "
               << x * C2m_g32(x, mQ, nf) / norm << "   " 
               << x * C2m_ps32(x, mQ, nf) / norm << "   " 
               << x * CLm_g32(x, mQ, nf) / norm << "   " 
               << x * CLm_ps32(x, mQ, nf) / norm << "   "
               << endl;
	}

    output.close();

    return 0;

}