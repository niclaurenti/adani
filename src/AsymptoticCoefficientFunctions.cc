#include "adani/AsymptoticCoefficientFunctions.h"
#include "adani/HighScaleCoefficientFunctions.h"
#include "adani/HighEnergyCoefficientFunctions.h"
#include <cmath>

//==========================================================================================//
//  Asymptotic limit of the gluon coefficient function for F2 at O(alpha_s).
//------------------------------------------------------------------------------------------//

double C2_g1_asymptotic(double x, double mQ) {

    return C2_g1_highscale(x, mQ);

}

//==========================================================================================//
//  Asymptotic limit of the gluon coefficient function for FL at O(alpha_s).
//------------------------------------------------------------------------------------------//

double CL_g1_asymptotic(double x) {

    return CL_g1_highscale(x);

}

//==========================================================================================//
//  Asymptotic limit of the gluon coefficient function for F2 at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double C2_g2_asymptotic(double x, double mQ, double mMu) {

    return C2_g2_highscale(x, mQ, mMu) + C2_g2_power_terms(x, mQ, mMu);

}

//==========================================================================================//
//  Asymptotic limit of the quark coefficient function for F2 at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double C2_ps2_asymptotic(double x, double mQ, double mMu) {

    return C2_ps2_highscale(x, mQ, mMu) + C2_ps2_power_terms(x, mQ, mMu);

}

//==========================================================================================//
//  Asymptotic limit of the gluon coefficient function for FL at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double CL_g2_asymptotic(double x, double mQ, double mMu) {

    return CL_g2_highscale(x, mQ, mMu) + CL_g2_power_terms(x, mQ, mMu);

}

//==========================================================================================//
//  Asymptotic limit of the quark coefficient function for FL at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double CL_ps2_asymptotic(double x, double mQ, double mMu) {

    return CL_ps2_highscale(x, mQ, mMu) + CL_ps2_power_terms(x, mQ, mMu);

}

//==========================================================================================//
//  Asymptotic limit of the gluon coefficient function for F2 at O(alpha_s^3) at leading logarithm.
//------------------------------------------------------------------------------------------//

// double C2_g3_asymptoticLL(double x, double mQ, double mMu, int nf, int v) {

//     return C2_g3_highscale(x, mQ, mMu, nf, v) + C2_g3_power_termsLL(x, mQ, mMu);

// }

//==========================================================================================//
//  Asymptotic limit of the gluon coefficient function for F2 at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double C2_g3_asymptotic(double x, double mQ, double mMu, int nf, int v1, int v2) {

    return C2_g3_highscale(x, mQ, mMu, nf, v1) + C2_g3_power_terms(x, mQ, mMu, nf, v2);

}

//==========================================================================================//
//  Asymptotic limit of the quark coefficient function for F2 at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double C2_ps3_asymptotic(double x, double mQ, double mMu, int nf, int v) {

    return C2_ps3_highscale(x, mQ, mMu, nf) + C2_ps3_power_terms(x, mQ, mMu, nf, v);

}

//==========================================================================================//
//  Asymptotic limit of the gluon coefficient function for FL at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double CL_g3_asymptotic(double x, double mQ, double mMu, int nf, int v) {

    return CL_g3_highscale(x, mQ, mMu, nf) + CL_g3_power_terms(x, mQ, mMu, nf, v);

}

//==========================================================================================//
//  Asymptotic limit of the quark coefficient function for FL at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double CL_ps3_asymptotic(double x, double mQ, double mMu, int nf, int v) {

    return CL_ps3_highscale(x, mQ, mMu, nf) + CL_ps3_power_terms(x, mQ, mMu, nf, v);

}
