#include "adani/AsymptoticCoefficientFunctions.h"
#include "adani/HighEnergyCoefficientFunctions.h"
#include "adani/HighScaleCoefficientFunctions.h"
#include <cmath>

//==========================================================================================//
//  Asymptotic limit of the gluon coefficient function for F2 at O(alpha_s).
//------------------------------------------------------------------------------------------//

double C2_g1_asymptotic(double x, double m2Q2) {

    return C2_g1_highscale(x, m2Q2);
}

//==========================================================================================//
//  Asymptotic limit of the gluon coefficient function for FL at O(alpha_s).
//------------------------------------------------------------------------------------------//

double CL_g1_asymptotic(double x) { return CL_g1_highscale(x); }

//==========================================================================================//
//  Asymptotic limit of the gluon coefficient function for F2 at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double C2_g2_asymptotic(double x, double m2Q2, double m2mu2) {

    return C2_g2_highscale(x, m2Q2, m2mu2) + C2_g2_power_terms(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  Asymptotic limit of the quark coefficient function for F2 at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double C2_ps2_asymptotic(double x, double m2Q2, double m2mu2) {

    return C2_ps2_highscale(x, m2Q2, m2mu2)
           + C2_ps2_power_terms(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  Asymptotic limit of the gluon coefficient function for FL at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double CL_g2_asymptotic(double x, double m2Q2, double m2mu2) {

    return CL_g2_highscale(x, m2Q2, m2mu2) + CL_g2_power_terms(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  Asymptotic limit of the quark coefficient function for FL at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double CL_ps2_asymptotic(double x, double m2Q2, double m2mu2) {

    return CL_ps2_highscale(x, m2Q2, m2mu2)
           + CL_ps2_power_terms(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  Asymptotic limit of the gluon coefficient function for F2 at O(alpha_s^3) at
//  leading logarithm.
//------------------------------------------------------------------------------------------//

// double C2_g3_asymptoticLL(double x, double m2Q2, double m2mu2, int nf, int v)
// {

//     return C2_g3_highscale(x, m2Q2, m2mu2, nf, v) + C2_g3_power_termsLL(x,
//     m2Q2, m2mu2);

// }

//==========================================================================================//
//  Asymptotic limit of the gluon coefficient function for F2 at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double
C2_g3_asymptotic(double x, double m2Q2, double m2mu2, int nf, int v1, int v2) {

    return C2_g3_highscale(x, m2Q2, m2mu2, nf, v1)
           + C2_g3_power_terms(x, m2Q2, m2mu2, nf, v2);
}

//==========================================================================================//
//  Asymptotic limit of the quark coefficient function for F2 at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double C2_ps3_asymptotic(double x, double m2Q2, double m2mu2, int nf, int v) {

    return C2_ps3_highscale(x, m2Q2, m2mu2, nf)
           + C2_ps3_power_terms(x, m2Q2, m2mu2, nf, v);
}

//==========================================================================================//
//  Asymptotic limit of the gluon coefficient function for FL at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double CL_g3_asymptotic(double x, double m2Q2, double m2mu2, int nf, int v) {

    return CL_g3_highscale(x, m2Q2, m2mu2, nf)
           + CL_g3_power_terms(x, m2Q2, m2mu2, nf, v);
}

//==========================================================================================//
//  Asymptotic limit of the quark coefficient function for FL at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double CL_ps3_asymptotic(double x, double m2Q2, double m2mu2, int nf, int v) {

    return CL_ps3_highscale(x, m2Q2, m2mu2, nf)
           + CL_ps3_power_terms(x, m2Q2, m2mu2, nf, v);
}
