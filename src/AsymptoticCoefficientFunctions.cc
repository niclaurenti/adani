#include "adani/AsymptoticCoefficientFunctions.h"
#include "adani/HighScaleCoefficientFunctions.h"
#include "adani/HighEnergyCoefficientFunctions.h"
#include "apfel/massivecoefficientfunctionsunp_sl.h"
#include <cmath>

using namespace apfel;

//==========================================================================================//
//  Asymptotic limit of the gluon coefficient function for F2 at O(alpha_s).
//------------------------------------------------------------------------------------------//

double C2m_g1_asymptotic(double x, double mQ) {

    return C2m_g1_highscale(x,mQ);

}

//==========================================================================================//
//  Asymptotic limit of the gluon coefficient function for FL at O(alpha_s).
//------------------------------------------------------------------------------------------//

double CLm_g1_asymptotic(double x) {
    
    return CLm_g1_highscale(x);

}

//==========================================================================================//
//  Asymptotic limit of the gluon coefficient function for F2 at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double C2m_g2_asymptotic(double x, double mQ, double mMu) {

    return C2m_g2_highscale(x,mQ,mMu) + C2m_g2_power_terms(x,mQ,mMu);

}

//==========================================================================================//
//  Asymptotic limit of the quark coefficient function for F2 at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double C2m_ps2_asymptotic(double x, double mQ, double mMu) {

    return C2m_ps2_highscale(x,mQ,mMu) + C2m_ps2_power_terms(x,mQ,mMu);

}

//==========================================================================================//
//  Asymptotic limit of the gluon coefficient function for FL at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double CLm_g2_asymptotic(double x, double mQ, double mMu) {

    return CLm_g2_highscale(x,mQ,mMu) + CLm_g2_power_terms(x,mQ,mMu);

}

//==========================================================================================//
//  Asymptotic limit of the quark coefficient function for FL at O(alpha_s^2).
//------------------------------------------------------------------------------------------//

double CLm_ps2_asymptotic(double x, double mQ, double mMu) {

    return CLm_ps2_highscale(x,mQ,mMu) + CLm_ps2_power_terms(x,mQ,mMu);

}

//==========================================================================================//
//  Asymptotic limit of the gluon coefficient function for F2 at O(alpha_s^3) at leading logarithm.
//------------------------------------------------------------------------------------------//

double C2m_g3_asymptoticLL(double x, double mQ, double mMu, int nf, int v) {

    return C2m_g3_highscale(x,mQ,mMu,nf,v) + C2m_g3_power_termsLL(x,mQ,mMu);

}

//==========================================================================================//
//  Asymptotic limit of the gluon coefficient function for F2 at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double C2m_g3_asymptotic(double x, double mQ, double mMu, int nf, int v1, int v2) {

    return C2m_g3_highscale(x,mQ,mMu,nf,v1) + C2m_g3_power_terms(x,mQ,mMu,nf,v2);

}

//==========================================================================================//
//  Asymptotic limit of the quark coefficient function for F2 at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double C2m_ps3_asymptotic(double x, double mQ, double mMu, int nf) {

    return C2m_ps3_highscale(x,mQ,mMu,nf) + C2m_ps3_power_terms(x,mQ,mMu,nf);

}

//==========================================================================================//
//  Asymptotic limit of the gluon coefficient function for FL at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double CLm_g3_asymptotic(double x, double mQ, double mMu, int nf) {

    return CLm_g3_highscale(x,mQ,mMu,nf) + CLm_g3_power_terms(x,mQ,mMu,nf);

}

//==========================================================================================//
//  Asymptotic limit of the quark coefficient function for FL at O(alpha_s^3).
//------------------------------------------------------------------------------------------//

double CLm_ps3_asymptotic(double x, double mQ, double mMu, int nf) {

    return CLm_ps3_highscale(x,mQ,mMu,nf) + CLm_ps3_power_terms(x,mQ,mMu,nf);

}
