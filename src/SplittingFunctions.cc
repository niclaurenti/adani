#include "../include/SplittingFunctions.h"
#include "../include/ColorFactors.h"
#include "../include/SpecialFunctions.h"
#include "apfel/specialfunctions.h"
#include <cmath>

using namespace std;
using namespace apfel;

double pgq(double x){

    if (x<0 || x>1) return 0;

    return 2. / x - 2. + x ;
}

//_________________________________________________

double Pgq0(double x) {

    if (x<0 || x>1) return 0;

    return 2. * CF * pgq(x) / 4 / M_PI ;
}

//_________________________________________________

double Pgg0reg(double x) {

    if (x<0 || x>1) return 0;

    double tmp = 1. / x - 2. + x - x * x ;

    return CA * 4. * tmp / 4 / M_PI ;    

}

//_________________________________________________________

double Pgg0loc(int nf) {

    double tmp = 11. / 3 * CA - 2. / 3 * nf ;

    return tmp / 4. / M_PI ;

}

//_____________________________________________________________

double Pgg0sing (double x) {
    if (x<0 || x>=1) return 0;

    return CA / M_PI / (1 - x) ;

}

//_________________________________________________________

double Pgq1(double x, int nf) {
    if (x<0 || x>=1) return 0;

    double H0 = H(x, 0) ;
    double H1 = H(x, 1) ;

    double H00 = H(x, 0, 0) ;
    double H10 = H(x, 1, 0) ;
    double Hm10 = H(x, -1, 0) ;
    double H01 = H(x, 0, 1) ;
    double H11 = H(x, 1, 1) ;

    double norm = (16 * M_PI * M_PI) ;

    double tmp_CACF = ( 1. / x + 2. * pgq(x) * (H10 + H11 + H01 - 11. / 6 * H1) 
                        - x*x * (8. / 3 * H0 - 44./ 9) + 4. * zeta(2) - 2. - 7. * H0 
                        + 2. * H00 - 2. * H1 * x + (1 + x) * (2. * H00 - 5 * H0 + 37. / 9) 
                        - 2. * pgq(-x) * Hm10 ) ;
    
    double tmp_CFnf = 2. / 3 * x - pgq(x) * (2. / 3 * H1 - 10. / 9) ;

    double tmp_CFCF = (
        pgq(x) * (3. * H1 - 2. * H11) + (1 + x) * (H00 - 7. / 2 + 7. / 2 * H0) 
        - 3. * H00 + 1 - 3. / 2* H0 + 2. * H1 * x
    ) ; 

    return 4. * (CA * CF * tmp_CACF + CF * nf * tmp_CFnf + CF * CF * tmp_CFCF) / norm ;
}


//______________________________________________________________________

double Pqq0reg(double x) {
    if (x<0 || x>=1) return 0;

    return CF * 2. * (- 1 - x )/ 4 / M_PI ;

}

//___________________________________________________________________________

double Pqq0loc() {

    return CF * 3. / 4 / M_PI ;

}

//___________________________________________________________________________

double Pqq0sing(double x) {
    if (x<0 || x>=1) return 0;

    return CF * 2. * 2. / ( 1 - x )/ 4. / M_PI ;

}

//___________________________________________________________________________

