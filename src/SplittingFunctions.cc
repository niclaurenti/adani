#include "../include/SplittingFunctions.h"
#include "../include/ColorFactors.h"
#include "apfel/specialfunctions.h"
#include <cmath>

using namespace std;
using namespace apfel;

double pgq(double x){

    if (x<0 || x>1) return 0;

    return 2. / x - 2. + x ;
}

double Pgq0(double x) {

    if (x<0 || x>1) return 0;

    return 2. * CF * pgq(x) / 4 / M_PI ;
}

double Pgg0reg(double x) {

    if (x<0 || x>1) return 0;

    double tmp = 1. / x - 2. + x - x * x ;

    return CA * 4. * tmp / 4 / M_PI ;    

}

double Pgg0loc(int nf) {

    double tmp = 11. / 3 * CA - 2. / 3 * nf ;

    return tmp / 4. / M_PI ;

}

double Pgg0sing (double x) {
    if (x<0 || x>=1) return 0;

    return CA / M_PI / (1 - x) ;

}


