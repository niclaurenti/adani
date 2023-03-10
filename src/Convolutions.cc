#include "adani/Convolutions.h"
#include "adani/Constants.h"
#include "adani/ExactCoefficientFunctions.h"
#include "adani/MasslessCoefficientFunctions.h"
#include "adani/MatchingConditions.h"
#include "adani/SpecialFunctions.h"
#include "adani/SplittingFunctions.h"

#include <cmath>
#include <iostream>

// TODO : in all the numerical convolutions the gsl default error is first switched off and
// then swithced on again: see if it can be done once for all

//==========================================================================================//
//  Convolution between first order massless quark coefficient functions for F2 and first
//  order matching KQg1
//------------------------------------------------------------------------------------------//

double C2_b1_x_K_Qg1(double x, double m2mu2) {

    if (x<0 || x>=1) return 0;

    double x2 = x * x ;
    double x1 = 1 - x ;
    double L = log(x) ;
    double L1 = log(x1) ;

    double res = (
        -5. / 2. +  2. * (3. - 4. * x) * x + zeta2 * (-1. + 2. * x - 4. * x2)
        + L1 * L1 * (1. - 2. * x1 * x) - 0.5 * L1 * (7. + 4. * x * (-4. + 3. * x) + (4. - 8. * x1 * x) * L)
        + 0.5 * L * (-1. + 4. * x * (3. * x - 2.) + (1. - 2. * x + 4. * x2) * L) + (2. * x - 1.) * Li2(x)
    ) ;

    return 4. * CF * TR * res * log(1./m2mu2) ;

}

//==========================================================================================//
//  Convolution between first order massless quark coefficient functions for FL and first
//  order matching KQg1
//------------------------------------------------------------------------------------------//

double CL_b1_x_K_Qg1(double x, double m2mu2) {

    if (x<0 || x>=1) return 0;

    double x2 = x * x ;
    double L = log(x) ;

    return 8. * CF * TR * (1. + x - 2. * x2 + 2. * x * L) * log(1./m2mu2) ;

}

//==========================================================================================//
//  Convolution between first order matching KQg1 and the second order matching Kgg2
//------------------------------------------------------------------------------------------//

// double K_Qg1_x_K_gg2(double x, double m2mu2) {

//     if (x<0 || x>=1) return 0;

//     double Lmu=log(1./m2mu2);
//     double Lmu2=Lmu*Lmu;

//     double pi2=M_PI*M_PI;
//     double pi3=pi2*M_PI;

//     double x2=x*x;
//     double x3=x2*x;

//     double L=log(x);
//     double L2=L*L;
//     double L3=L2*L;
//     double L4=L3*L;

//     double Lm=log(1-x);
//     //double Lm2=Lm*Lm;

//     double K_log2_TR2 = 16./9.*(1-2*x+2*x*x);

//     double K_log2_CATR = 4./9./x*(4 + 3*x +24*x2 - 31*x3 + 6*x*(1 - 2*x +2*x2)*Lm + 6*x*(1+4*x)*L);

//     double K_log2_CFTR = 4./9.*(-81 + 8./x +135*x -62*x2 + 3*(-9+8*x2)*L + 9*(-1+2*x)*L2);

//     double K_log2= TR*TR*K_log2_TR2 + CA*TR*K_log2_CATR + CF*TR*K_log2_CFTR;

//     double K_log_CATR = 4./27./x*(92 - 33*x +528*x2 - 551*x3 + 60*x*(1-2*x+2*x2)*Lm + 6*x*(13+64*x+26*x2)*L + 18*x*(-1+2*x)*L2);

//     double K_log_CFTR = 4./9./x*(-8 - 423*x +792*x2 - 352*x3 + 6*x*(-42 + 9*x + 20*x2)*L + 9*x*(-5+12*x)*L2 + 6*x*(-1+2*x)*L3);

//     double K_log= K_log_CATR*CA*TR + K_log_CFTR*CF*TR;

//     double K_const_CFTR = 1./3/x*(-16 -660*x + 999*x2 - 368*x3 + 9*x*(-48 - 9*x +16*x2)*L + 3*x*(-28 +23*x)*L2 + 2*x*(-5+12*x)*L3+x*(-1+2*x)*L4);

//     double K_const_CATR = 2./81./x*(556 - 660*x + 5052*x2 - 18*pi2*x2 - 4903*x3 + 6*x*(47 - 121*x + 130*x2)*Lm + 210*x*L + 3036*x2*L + 1320*x3*L - 171*x*L2 + 450*x2*L2 - 18*x*L3 + 36*x2*L3 +108*x2*Li2(x));

//     double K_const= K_const_CFTR*CF*TR + K_const_CATR*CA*TR;

//     return 2*TR*Lmu*(K_log2*Lmu2 + K_log*(-Lmu) + K_const)/64/pi3;

// }

// Requires Li4 that is not implemented

// double C2_b2_x_K_bg1(double x, double m2mu2, int nf) {

//     if (x<0 || x>=1) return 0;

//     double Lmu = log(1./m2mu2) ;
//     double pi3 = M_PI * M_PI * M_PI ;

//     double x2 = x * x ;
//     //double x3=x2*x;
//     //double x4=x3*x;

//     double L = log(x) ;
//     double L2 = L * L ;
//     double L3 = L2 * L ;
//     double L4 = L3 * L ;
//     double L5 = L4 * L ;

//     double xm = 1 - x ;
//     double xm2 = xm * xm ;

//     double Lm = log(xm) ;
//     double Lm2 = Lm * Lm ;
//     double Lm3 = Lm2 * Lm ;
//     double Lm4 = Lm3 * Lm ;

//     double Li2x = Li2(x) ;
//     double Li2xm = Li2(xm) ;
//     double Li2x_xm = Li2(-x / xm) ;

//     double Li3x = Li3(x) ;
//     double Li3xm = Li3(xm) ;
//     double Li3x_xm = Li3(-x / xm) ;

//     double Li4x = Li4(x) ;
//     double Li4xm = Li4(xm) ;
//     double Li4x_xm = Li4(-x / xm) ;

//     double c_const = -222.268 - 67.072 * x + 357.041 * x2 + 288.96 * x2 * atanh(1 - 2 * x);

//     double c_Lm3 = -16.9267 + 62.2978 * x - 65.8156 * x2 ;

//     double c_Lm4 = 7.11111 * (-4.67148 - x + x2) ;

//     double c_L3 = 9.86296 - 16.5387 * x + 12.5833 * x2 ;

//     double c_L4 = 0.740741 + 0.719 * x ;

//     double c_L5 = 0.2876 * x ;

//     double c_Li2xm = x * (437.4 - 2.91 * x) ;

//     double c_Li2x = - 659.905 + 105.837 * x - 47.5571 * x2 ;

//     double c_Lm2 = 194.191 - 353.395 * x + 143.651 * x2 - 85.3333 * (- 1.8475 -  x + x2) * Li2xm - 147.1 * Li2x_xm ;

//     double c_L2 = 48.2253 - 104.622 * x + 165.975 * x2 - 37.75 * xm2 * Lm - 73.55 * Lm2 - 109.35 * Li2x - 147.1 * Li2x_xm ;

//     double c_Li3xm = 21.8133 - 385.28 * x - 100.693 * x2 ;

//     double c_Li3x = (437.4 - 218.7 * x) * x ;

//     double c_Lm = 716.421 - 1188.82 * x + 656.924 * x2 + (- 21.8133 + 385.28 * x + 15.36 * x2) * Li2xm - 256. / 3. * x2 * Li2x + 170.667 * (- 1.8475 - x + x2) * Li3xm - 294.2 * Li3x_xm ;

//     double c_L = 113.954 - 128.79 * x - 516.556 * x2 + (136.193 - 101.56 * x + 112.113 * x2) * Lm2 - 28.4444 * (- 5.29516 - x + x2) * Lm3 + 218.7 * (- 2. + x) * x * Li2x + Lm * (- 113.25 + 882.6 * x - 334.86 * x2 + 294.2 * Li2x_xm) + 218.7 * Li3x + 294.2 * Li3x_xm ;

//     double c_Li4xm = 315.307 + 170.667 * x - 170.667 * x2 ;

//     double c_Li4x = 218.7 ;

//     double c_Li4x_xm = 294.2 ;

//     double c = c_const + c_Lm3 * Lm3 + c_Lm4 * Lm4 + c_L3 * L3 + c_L4 * L4 + c_L5 * L5 + c_Li2xm * Li2xm + c_Li2x * Li2x + c_Lm2 * Lm2 + c_L2 * L2 + c_Li3xm * Li3xm + c_Li3x * Li3x + c_Lm * Lm + c_L * L + c_Li4xm * Li4xm + c_Li4x * Li4x + c_Li4x_xm * Li4x_xm ;

//     double c_nf_const = - 7.68063 - 144.821 * x + 179.77 * x2 - 6 * x2 * atanh(1 - 2 * x) ;

//     double c_nf_Lm2 = - 4.57407 + 12.7037 * x - 12.4259 * x2 ;

//     double c_nf_Lm3 = 16./27. * (1 - 2 * x + 2 * x2) ;

//     double c_nf_L2 = - 4.88889 - 9.778 * x + 4.0565 * x2 ;

//     double c_nf_L3 = - 0.740741 - 0.185 * x ;

//     double c_nf_L4 = 0.0925 * x ;

//     double c_nf_Lm = - 22.098 + 18.3423 * x + 13.1046 * x2 -7.11111 * (0.078125 - x + x2) * Li2xm ;

//     double c_nf_Li2x = 16.2774 + 34.5223 * x - 32.9649 * x2 ;

//     double c_nf_L = - 10.538 - 80.9461 * x + 12.8354 * x2 - 8.113 * xm2 * Lm - 3.55556 * (0.078125 - x + x2) * Lm2 + 8.113 * Li2x ;

//     double c_nf_Li3xm = 7.11111 * (0.078125 - x + x2);

//     double c_nf_Li3x = 8.113;

//     double c_nf = c_nf_const + c_nf_Lm2 * Lm2 + c_nf_Lm3 * Lm3 + c_nf_L2 * L2 + c_nf_L3 * L3 + c_nf_L4 * L4 + c_nf_Lm * Lm + c_nf_Li2x * Li2x + c_nf_L * L + c_nf_Li3xm * Li3xm + c_nf_Li3x * Li3x ;

//     return 2 * TR * Lmu * (c + nf * c_nf) / 64. / pi3 ;

// }

//==========================================================================================//
//  Convolution between first order massive gluon coefficient functions for F2 and splitting
//  function Pgq0
//------------------------------------------------------------------------------------------//

//  Integrand

double C2_g1_x_Pgq0_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    return C2_g1(z, m2Q2) * Pgq0(x / z) / z ;

}

// Result

double C2_g1_x_Pgq0(double x, double m2Q2) {

    double result, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, static_cast<int>(nan(""))};
    //It is not dependent on nf so nf is put to nan

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    gsl_function F;
    F.function = &C2_g1_x_Pgq0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &result, &error);

    gsl_set_error_handler (old_handler);

    gsl_integration_workspace_free (w);

    return result;

}

//==========================================================================================//
//  Convolution between first order massive gluon coefficient functions for FL and splitting
//  function Pgq0
//------------------------------------------------------------------------------------------//

// Integrand

double CL_g1_x_Pgq0_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    return CL_g1(z, m2Q2) * Pgq0(x / z) / z;

}

// Result

double CL_g1_x_Pgq0(double x, double m2Q2) {

    double result, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, static_cast<int>(nan(""))};
    //It is not dependent on nf so nf is put to nan

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    gsl_function F;
    F.function = &CL_g1_x_Pgq0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &result, &error);

    gsl_set_error_handler (old_handler);

    gsl_integration_workspace_free (w);

    return result;

}

//==========================================================================================//
//  Convolution between first order massive gluon coefficient functions for F2 and splitting
//  function Pgg0
//------------------------------------------------------------------------------------------//

// Integrand of the regular part

double C2_g1_x_Pgg0_reg_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    return C2_g1(z, m2Q2) * Pgg0reg(x / z) / z;

}

// Integrand of the singular part

double C2_g1_x_Pgg0_sing_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    return Pgg0sing(z) * (C2_g1(x / z, m2Q2) / z - C2_g1(x , m2Q2)) ;

}

// Needed to integrate singular part

double Pgg0sing_integrated(double x) {

    return - Pgg0sing(0.) * log(1. - x) ;

}

// Result

double C2_g1_x_Pgg0(double x, double m2Q2, int nf) {

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double regular, singular, local, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, nf};

    gsl_function F;
    F.function = &C2_g1_x_Pgg0_reg_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &regular, &error);

    F.function = &C2_g1_x_Pgg0_sing_integrand;
    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &singular, &error);

    gsl_set_error_handler (old_handler);

    local = C2_g1(x, m2Q2) * (Pgg0loc(nf) - Pgg0sing_integrated(x)) ;

    gsl_integration_workspace_free (w);

    return regular + singular + local ;

}

//==========================================================================================//
//  Convolution between first order massive gluon coefficient functions for F2 and splitting
//  function Pgg0
//------------------------------------------------------------------------------------------//

// Integrand of the regular part

double CL_g1_x_Pgg0_reg_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    return CL_g1(z, m2Q2) * Pgg0reg(x / z) / z;

}

// Integrand of the singular part

double CL_g1_x_Pgg0_sing_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    return Pgg0sing(z) * (CL_g1(x / z, m2Q2) / z - CL_g1(x , m2Q2)) ;

}

// Result

double CL_g1_x_Pgg0(double x, double m2Q2, int nf) {

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double regular, singular, local, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, nf};

    gsl_function F;
    F.function = &CL_g1_x_Pgg0_reg_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &regular, &error);

    F.function = &CL_g1_x_Pgg0_sing_integrand;
    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &singular, &error);

    gsl_set_error_handler (old_handler);

    local = CL_g1(x, m2Q2) * (Pgg0loc(nf) - Pgg0sing_integrated(x)) ;

    gsl_integration_workspace_free (w);

    return regular + singular + local ;

}

//==========================================================================================//
//  Convolution between first order massive gluon coefficient functions for F2 and splitting
//  function Pgq1
//------------------------------------------------------------------------------------------//

// Integrand

double C2_g1_x_Pgq1_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return C2_g1(z, m2Q2) * Pgq1(x / z, nf) / z ;

}

// Result

double C2_g1_x_Pgq1(double x, double m2Q2, int nf) {

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double result, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, nf};

    gsl_function F;
    F.function = &C2_g1_x_Pgq1_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &result, &error);

    gsl_set_error_handler (old_handler);

    gsl_integration_workspace_free (w);

    return result ;

}

//==========================================================================================//
//  Convolution between first order massive gluon coefficient functions for FL and splitting
//  function Pgq1
//------------------------------------------------------------------------------------------//

// Integrand

double CL_g1_x_Pgq1_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return CL_g1(z, m2Q2) * Pgq1(x / z, nf) / z;

}

// Result

double CL_g1_x_Pgq1(double x, double m2Q2, int nf) {

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double result, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, nf};

    gsl_function F;
    F.function = &CL_g1_x_Pgq1_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &result, &error);

    gsl_set_error_handler (old_handler);

    gsl_integration_workspace_free (w);

    return result ;

}

//==========================================================================================//
//  Convolution between the mu-independent term of the second order massive gluon coefficient
//  functions for F2 and splitting function Pgq0
//------------------------------------------------------------------------------------------//

// Integrand

double C2_g20_x_Pgq0_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    return C2_g20(z, m2Q2) * Pgq0(x / z) / z ;

}

// Result

double C2_g20_x_Pgq0(double x, double m2Q2) {

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double result, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, static_cast<int>(nan(""))};

    gsl_function F;
    F.function = &C2_g20_x_Pgq0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &result, &error);

    gsl_set_error_handler (old_handler);

    gsl_integration_workspace_free (w);

    return result ;

}

//==========================================================================================//
//  Convolution between the mu-independent term of the second order massive gluon coefficient
//  functions for FL and splitting function Pgq0
//------------------------------------------------------------------------------------------//

// Integrand

double CL_g20_x_Pgq0_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    return CL_g20(z, m2Q2) * Pgq0(x / z) / z;

}

// Result

double CL_g20_x_Pgq0(double x, double m2Q2) {

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double result, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, static_cast<int>(nan(""))};

    gsl_function F;
    F.function = &CL_g20_x_Pgq0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &result, &error);

    gsl_set_error_handler (old_handler);

    gsl_integration_workspace_free (w);

    return result ;

}

//==========================================================================================//
//  Convolution between the mu-independent term of the second order massive quark coefficient
//  functions for F2 and splitting function Pqq0
//------------------------------------------------------------------------------------------//

// Integrand of the regular part

double C2_ps20_x_Pqq0_reg_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    return C2_ps20(z, m2Q2) * Pqq0reg(x / z) / z;

}

// Integrand of the singular part

double C2_ps20_x_Pqq0_sing_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    return Pqq0sing(z) * (C2_ps20(x / z, m2Q2) / z - C2_ps20(x , m2Q2)) ;

}

// Needed to integrate singular part

double Pqq0sing_integrated(double x) {

    return - Pqq0sing(0.) * log(1. - x) ;

}

// Result

double C2_ps20_x_Pqq0(double x, double m2Q2, int nf) {

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double regular, singular, local, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, nf};

    gsl_function F;
    F.function = &C2_ps20_x_Pqq0_reg_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &regular, &error);

    F.function = &C2_ps20_x_Pqq0_sing_integrand;
    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &singular, &error);

    gsl_set_error_handler (old_handler);

    local = C2_ps20(x, m2Q2) * (Pqq0loc() - Pqq0sing_integrated(x));

    gsl_integration_workspace_free (w);

    return regular + singular + local ;

}

//==========================================================================================//
//  Convolution between the mu-independent term of the second order massive quark coefficient
//  functions for FL and splitting function Pqq0
//------------------------------------------------------------------------------------------//

// Integrand of the regular part

double CL_ps20_x_Pqq0_reg_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    return CL_ps20(z, m2Q2) * Pqq0reg(x / z) / z;

}

// Integrand of the singular part

double CL_ps20_x_Pqq0_sing_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    return Pqq0sing(z) * (CL_ps20(x / z, m2Q2) / z - CL_ps20(x , m2Q2)) ;

}

// Result

double CL_ps20_x_Pqq0(double x, double m2Q2, int nf) {

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double regular, singular, local, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, nf};

    gsl_function F;
    F.function = &CL_ps20_x_Pqq0_reg_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &regular, &error);

    F.function = &CL_ps20_x_Pqq0_sing_integrand;
    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &singular, &error);

    gsl_set_error_handler (old_handler);

    local = CL_ps20(x, m2Q2) * (Pqq0loc() - Pqq0sing_integrated(x));

    gsl_integration_workspace_free (w);

    return regular + singular + local ;

}

//==========================================================================================//
//  Analytical onvolution between the splitting functions Pgg0 and Pgq0
//------------------------------------------------------------------------------------------//

double Pgg0_x_Pgq0(double x, int nf) {

    return (
        - 4. * CF * nf * (2. + (- 2. + x) * x)
        + 2. * CA * CF * (
            - 40. + x * (26. + x * (17. + 8. * x))
            + 12. * (2. + (- 2. + x) * x) * log(1. - x)
            - 24. * (1. + x + x * x) * log(x)
        )
    ) / 3. / x ;

}

//==========================================================================================//
//  Analytical onvolution between the splitting functions Pqq0 and Pgq0
//------------------------------------------------------------------------------------------//

double Pqq0_x_Pgq0(double x) {

    return (
        2. * CF * CF * (
            4. * (2. + (- 2. + x) * x) * log(1. - x)
            - x * (- 4. + x + 2. * (- 2. + x) * log(x))
        )
    ) / x ;

}

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions for F2 and the
//  convolution between the splitting functions Pgg0 and Pgq0
//------------------------------------------------------------------------------------------//

// Integrand

double C2_g1_x_Pgg0_x_Pgq0_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return C2_g1(z, m2Q2) * Pgg0_x_Pgq0(x / z, nf) / z;

}

// Result

double C2_g1_x_Pgg0_x_Pgq0(double x, double m2Q2, int nf) {

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double result, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, nf};

    gsl_function F;
    F.function = &C2_g1_x_Pgg0_x_Pgq0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &result, &error);

    gsl_set_error_handler (old_handler);

    gsl_integration_workspace_free (w);

    return result ;

}

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions for FL and the
//  convolution between the splitting functions Pgg0 and Pgq0
//------------------------------------------------------------------------------------------//

// Integrand

double CL_g1_x_Pgg0_x_Pgq0_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return CL_g1(z, m2Q2) * Pgg0_x_Pgq0(x / z, nf) / z;

}

// Result

double CL_g1_x_Pgg0_x_Pgq0(double x, double m2Q2, int nf) {

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double result, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, nf};

    gsl_function F;
    F.function = &CL_g1_x_Pgg0_x_Pgq0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &result, &error);

    gsl_set_error_handler (old_handler);

    gsl_integration_workspace_free (w);

    return result ;

}

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions for F2 and the
//  convolution between the splitting functions Pqq0 and Pgq0
//------------------------------------------------------------------------------------------//

// Integrand

double C2_g1_x_Pqq0_x_Pgq0_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    return C2_g1(z, m2Q2) * Pqq0_x_Pgq0(x / z) / z;

}

// Result

double C2_g1_x_Pqq0_x_Pgq0(double x, double m2Q2, int nf) {

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double result, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, nf};

    gsl_function F;
    F.function = &C2_g1_x_Pqq0_x_Pgq0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &result, &error);

    gsl_set_error_handler (old_handler);

    gsl_integration_workspace_free (w);

    return result ;

}

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions for FL and the
//  convolution between the splitting functions Pqq0 and Pgq0
//------------------------------------------------------------------------------------------//

// Integrand

double CL_g1_x_Pqq0_x_Pgq0_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    return CL_g1(z, m2Q2) * Pqq0_x_Pgq0(x / z) / z;

}

// Result

double CL_g1_x_Pqq0_x_Pgq0(double x, double m2Q2, int nf) {

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double result, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, nf};

    gsl_function F;
    F.function = &CL_g1_x_Pqq0_x_Pgq0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &result, &error);

    gsl_set_error_handler (old_handler);

    gsl_integration_workspace_free (w);

    return result ;

}

//==========================================================================================//
//  Convolution between first order massive gluon coefficient functions for F2 and splitting
//  function Pgg1
//------------------------------------------------------------------------------------------//

// Integrand of the regular part

double C2_g1_x_Pgg1_reg_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return C2_g1(z, m2Q2) * Pgg1reg(x / z, nf) / z;

}

// Integrand of the singular part

double C2_g1_x_Pgg1_sing_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return Pgg1sing(z, nf) * (C2_g1(x / z, m2Q2) / z - C2_g1(x , m2Q2)) ;

}

// Needed to integrate singular part

double Pgg1sing_integrated(double x, int nf) {

    return - Pgg1sing(0., nf) * log(1. - x) ;

}

// Result

double C2_g1_x_Pgg1(double x, double m2Q2, int nf) {

    double regular, singular, local, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, nf};

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    gsl_function F;
    F.function = &C2_g1_x_Pgg1_reg_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &regular, &error);

    F.function = &C2_g1_x_Pgg1_sing_integrand;
    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &singular, &error);

    gsl_set_error_handler (old_handler);

    local = C2_g1(x, m2Q2) * (Pgg1loc(nf) - Pgg1sing_integrated(x, nf)) ;

    gsl_integration_workspace_free (w);

    return regular + singular + local ;

}

//==========================================================================================//
//  Convolution between first order massive gluon coefficient functions for F2 and splitting
//  function Pgg1
//------------------------------------------------------------------------------------------//

// Integrand of the regular part

double CL_g1_x_Pgg1_reg_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return CL_g1(z, m2Q2) * Pgg1reg(x / z, nf) / z;

}

// Integrand of the regular part

double CL_g1_x_Pgg1_sing_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return Pgg1sing(z, nf) * (CL_g1(x / z, m2Q2) / z - CL_g1(x, m2Q2)) ;

}

// Result

double CL_g1_x_Pgg1(double x, double m2Q2, int nf) {

    double regular, singular, local, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, nf};

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    gsl_function F;
    F.function = &CL_g1_x_Pgg1_reg_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &regular, &error);

    F.function = &CL_g1_x_Pgg1_sing_integrand;
    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &singular, &error);

    gsl_set_error_handler (old_handler);

    local = CL_g1(x, m2Q2) * (Pgg1loc(nf) - Pgg1sing_integrated(x, nf)) ;

    gsl_integration_workspace_free (w);

    return regular + singular + local ;

}

//==========================================================================================//
//  Convolution between the mu-independent term of the second order massive quark coefficient
//  functions for F2 and splitting function Pqg0
//------------------------------------------------------------------------------------------//

// Integrand

double C2_ps20_x_Pqg0_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return C2_ps20(z, m2Q2) * Pqg0(x / z, nf) / z ;

}

// Result

double C2_ps20_x_Pqg0(double x, double m2Q2, int nf) {

    double result, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, nf};

    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

    gsl_function F;
    F.function = &C2_ps20_x_Pqg0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &result, &error);

    gsl_set_error_handler (old_handler);

    gsl_integration_workspace_free (w);

    return result ;

}

//==========================================================================================//
//  Convolution between the mu-independent term of the second order massive quark coefficient
//  functions for FL and splitting function Pqg0
//------------------------------------------------------------------------------------------//

// Integrand

double CL_ps20_x_Pqg0_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return CL_ps20(z, m2Q2) * Pqg0(x / z, nf) / z ;

}

// Result

double CL_ps20_x_Pqg0(double x, double m2Q2, int nf) {

    double result, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, nf};

    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

    gsl_function F;
    F.function = &CL_ps20_x_Pqg0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &result, &error);

    gsl_set_error_handler (old_handler);

    gsl_integration_workspace_free (w);

    return result ;

}

//==========================================================================================//
//  Convolution between the mu-independent term of the second order massive gluon coefficient
//  functions for F2 and splitting function Pgg0
//------------------------------------------------------------------------------------------//

// Integrand of the regular part

double C2_g20_x_Pgg0_reg_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    return C2_g20(z, m2Q2) * Pgg0reg(x / z) / z;

}

// Integrand of the singular part

double C2_g20_x_Pgg0_sing_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    return Pgg0sing(z) * (C2_g20(x / z, m2Q2) / z - C2_g20(x , m2Q2)) ;

}

// Result

double C2_g20_x_Pgg0(double x, double m2Q2, int nf) {

    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

    double regular, singular, local, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, nf};

    gsl_function F;
    F.function = &C2_g20_x_Pgg0_reg_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &regular, &error);

    F.function = &C2_g20_x_Pgg0_sing_integrand;
    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &singular, &error);

    gsl_set_error_handler (old_handler);

    local = C2_g20(x, m2Q2) * (Pgg0loc(nf) - Pgg0sing_integrated(x)) ;

    gsl_integration_workspace_free (w);

    return regular + singular + local ;

}

//==========================================================================================//
//  Convolution between the mu-independent term of the second order massive gluon coefficient
//  functions for FL and splitting function Pgg0
//------------------------------------------------------------------------------------------//

// Integrand of the regular part

double CL_g20_x_Pgg0_reg_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    return CL_g20(z, m2Q2) * Pgg0reg(x / z) / z;

}

// Integrand of the singular part

double CL_g20_x_Pgg0_sing_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    return Pgg0sing(z) * (CL_g20(x / z, m2Q2) / z - CL_g20(x , m2Q2)) ;

}

// Result

double CL_g20_x_Pgg0(double x, double m2Q2, int nf) {

    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

    double regular, singular, local, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, nf};

    gsl_function F;
    F.function = &CL_g20_x_Pgg0_reg_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &regular, &error);

    F.function = &CL_g20_x_Pgg0_sing_integrand;
    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &singular, &error);

    gsl_set_error_handler (old_handler);

    local = CL_g20(x, m2Q2) * (Pgg0loc(nf) - Pgg0sing_integrated(x)) ;

    gsl_integration_workspace_free (w);

    return regular + singular + local ;

}

//==========================================================================================//
//  Analytical onvolution between the splitting functions Pqg0 and Pgq0
//------------------------------------------------------------------------------------------//

double Pqg0_x_Pgq0(double x, int nf) {

    return 4. * CF * nf * (
        1. + 4. / 3 / x - x - 4. * x * x / 3 + 2. * (1 + x) * log(x)
    ) ;

}

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions for F2 and the convolution
//  between the splitting functions Pqg0 and Pgq0
//------------------------------------------------------------------------------------------//

// Integrand

double C2_g1_x_Pqg0_x_Pgq0_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return C2_g1(z, m2Q2) * Pqg0_x_Pgq0(x / z , nf) / z ;

}

// Result

double C2_g1_x_Pqg0_x_Pgq0(double x, double m2Q2, int nf) {

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double result, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, nf};

    gsl_function F;
    F.function = &C2_g1_x_Pqg0_x_Pgq0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &result, &error);

    gsl_set_error_handler (old_handler);

    gsl_integration_workspace_free (w);

    return result ;

}

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions for FL and the convolution
//  between the splitting functions Pqg0 and Pgq0
//------------------------------------------------------------------------------------------//

// Integrand

double CL_g1_x_Pqg0_x_Pgq0_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return CL_g1(z, m2Q2) * Pqg0_x_Pgq0(x / z , nf) / z ;

}

// Result

double CL_g1_x_Pqg0_x_Pgq0(double x, double m2Q2, int nf) {

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double result, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, nf};

    gsl_function F;
    F.function = &CL_g1_x_Pqg0_x_Pgq0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &result, &error);

    gsl_set_error_handler (old_handler);

    gsl_integration_workspace_free (w);

    return result ;

}

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions for F2 and the convolution
//  between the splitting functions Pgg0 and Pgg0
//------------------------------------------------------------------------------------------//

// Integrand of the regular part

double C2_g1_x_Pgg0_x_Pgg0_reg_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return C2_g1_x_Pgg0(z, m2Q2, nf) * Pgg0reg(x / z) / z;

}

// Integrand of the singular part

double C2_g1_x_Pgg0_x_Pgg0_sing_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return Pgg0sing(z) * (C2_g1_x_Pgg0(x / z, m2Q2, nf) / z - C2_g1_x_Pgg0(x, m2Q2, nf)) ;

}

// Result

double C2_g1_x_Pgg0_x_Pgg0(double x, double m2Q2, int nf) {

    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

    double regular, singular, local, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, nf};

    double C2_g1xPgg0 = C2_g1_x_Pgg0(x, m2Q2, nf) ;

    gsl_function F;
    F.function = &C2_g1_x_Pgg0_x_Pgg0_reg_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &regular, &error);

    F.function = &C2_g1_x_Pgg0_x_Pgg0_sing_integrand;
    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &singular, &error);

    gsl_set_error_handler (old_handler);

    local = C2_g1xPgg0 * (Pgg0loc(nf) - Pgg0sing_integrated(x)) ;

    gsl_integration_workspace_free (w);

    return regular + singular + local ;

}

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions for FL and the convolution
//  between the splitting functions Pgg0 and Pgg0
//------------------------------------------------------------------------------------------//

// Integrand of the regular part

double CL_g1_x_Pgg0_x_Pgg0_reg_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return CL_g1_x_Pgg0(z, m2Q2, nf) * Pgg0reg(x / z) / z;

}

// Integrand of the singular part

double CL_g1_x_Pgg0_x_Pgg0_sing_integrand(double z, void * p) {

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return Pgg0sing(z) * (CL_g1_x_Pgg0(x / z, m2Q2, nf) / z - CL_g1_x_Pgg0(x , m2Q2, nf)) ;

}

// Result

double CL_g1_x_Pgg0_x_Pgg0(double x, double m2Q2, int nf) {

    gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

    double regular, singular, local, error, abserr = 0.001, relerr = 0.001;
    struct function_params params = {x, m2Q2, nf};

    double CL_g1xPgg0 = CL_g1_x_Pgg0(x, m2Q2, nf) ;

    gsl_function F;
    F.function = &CL_g1_x_Pgg0_x_Pgg0_reg_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &regular, &error);

    F.function = &CL_g1_x_Pgg0_x_Pgg0_sing_integrand;
    gsl_integration_qag(&F, x, 1., abserr, relerr, 1000, 4, w, &singular, &error);

    gsl_set_error_handler (old_handler);

    local = CL_g1xPgg0 * (Pgg0loc(nf) - Pgg0sing_integrated(x)) ;

    gsl_integration_workspace_free (w);

    return regular + singular + local ;

}

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions for F2 and the convolution
//  between the splitting functions Pgg0 and Pgg0 using monte carlo mathods
//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgg0_reg1_integrand(double z[], size_t dim, void * p) {

    if (dim != 2) {
        std::cout << "error: dim != 2" << std::endl ;
        exit(-1);
    }

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    double z1 = z[0], z2 = z[1] ;

    double int_bound = theta(z1 - x) * theta(z2 - z1) ;

    return 1. / (z1 * z2) * Pgg0reg(x / z1) * Pgg0reg(z1 / z2) * C2_g1(z2, m2Q2) * int_bound ;

}

//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgg0_reg2_integrand(double z[], size_t dim, void * p) {

    if (dim != 2) {
        std::cout << "error: dim != 2" << std::endl ;
        exit(-1);
    }

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    double z1 = z[0], z2 = z[1] ;

    double int_bound = theta(z1 - x) * theta(z2 - z1) ;

    return 1. / z1 * Pgg0reg(x / z1) * Pgg0sing(z2) * (C2_g1(z1 / z2, m2Q2) / z2 -  C2_g1(z1, m2Q2)) * int_bound ;

}

//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgg0_reg3_integrand(double z[], size_t dim, void * p) {

    if (dim != 2) {
        std::cout << "error: dim != 2" << std::endl ;
        exit(-1);
    }

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    double z1 = z[0], z2 = z[1] ;

    double int_bound = theta(z1 - x) * theta(z1 - z2) ;

    return - 1. / z1 * Pgg0reg(x / z1) * C2_g1(z1, m2Q2) * Pgg0sing(z2) * int_bound ;

}

//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgg0_reg(double x, double m2Q2, int nf, size_t calls) {

    struct function_params params = {x, m2Q2, nf} ;
    double xl[2] = {x, x};
    double xu[2] = {1, 1};

    double err, regular1, regular2, regular3, regular4 ;

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);

    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    gsl_monte_function F;

    F.f = &C2_g1_x_Pgg0_x_Pgg0_reg1_integrand;
    F.dim = 2;
    F.params = &params;

    gsl_monte_vegas_init (s);
    gsl_monte_vegas_integrate(&F, xl, xu, 2, calls, r, s, &regular1, &err);

    F.f = &C2_g1_x_Pgg0_x_Pgg0_reg2_integrand;
    gsl_monte_vegas_init (s);
    gsl_monte_vegas_integrate(&F, xl, xu, 2, calls, r, s, &regular2, &err);

    double xl_new[] = {x, 0} ;

    F.f = &C2_g1_x_Pgg0_x_Pgg0_reg3_integrand;
    gsl_monte_vegas_init (s);
    gsl_monte_vegas_integrate(&F, xl_new, xu, 2, calls, r, s, &regular3, &err);

    gsl_monte_vegas_free (s);

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000) ;

    double abserr = 0.001, relerr = 0.001;

    gsl_function f;
    f.function = &C2_g1_x_Pgg0_reg_integrand;
    f.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&f, x, 1, abserr, relerr, 1000, 4, w, &regular4, &err);

    gsl_set_error_handler (old_handler);

    regular4 *= Pgg0loc(nf) ;

    gsl_integration_workspace_free (w);

    return regular1 + regular2 + regular3 + regular4;

}

//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgg0_sing1_integrand(double z[], size_t dim, void * p) {

    if (dim != 2) {
        std::cout << "error: dim != 2" << std::endl ;
        exit(-1);
    }

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    double z1 = z[0], z2 = z[1] ;

    return 1. / z2 * theta(z1 - x) * Pgg0sing(z1) * (
        theta(z2 - x / z1) / z1 * Pgg0reg(x / (z1 * z2))
        - theta(z2 - x) * Pgg0reg(x / z2)
   ) * C2_g1(z2, m2Q2) ;

}

//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgg0_sing2_integrand(double z[], size_t dim, void * p) {

    if (dim != 2) {
        std::cout << "error: dim != 2" << std::endl ;
        exit(-1);
    }

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    double z1 = z[0], z2 = z[1] ;

    return theta(z1 - x) * Pgg0sing(z1) * (
        Pgg0sing(z2) / z1 * (C2_g1(x / (z1 * z2), m2Q2) / z2 - C2_g1(x / z1, m2Q2)) * theta(z2 - x / z1)
        - Pgg0sing(z2) * (C2_g1(x / z2, m2Q2) / z2 - C2_g1(x, m2Q2)) * theta(z2 - x)
    ) ;

}

//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgg0_sing3_integrand(double z[], size_t dim, void * p) {

    if (dim != 2) {
        std::cout << "error: dim != 2" << std::endl ;
        exit(-1);
    }

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    double z1 = z[0], z2 = z[1] ;

    return - (theta(z1 - x) * Pgg0sing(z1) * Pgg0sing(z2) * (C2_g1(x / z1 , m2Q2) / z1 * theta(x / z1 - z2) - C2_g1(x , m2Q2) * theta(x - z2)));

}

//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgg0_sing(double x, double m2Q2, int nf, size_t calls) {

    struct function_params params = {x, m2Q2, nf} ;
    double xl[2] = {x, 0};
    //double xl[2] = {x, x};
    double xu[2] = {1, 1};

    double err, singular1, singular2, singular3, singular4 ;

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);

    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    gsl_monte_function F;

    F.f = &C2_g1_x_Pgg0_x_Pgg0_sing1_integrand;
    F.dim = 2;
    F.params = &params;

    gsl_monte_vegas_integrate(&F, xl, xu, 2, calls, r, s, &singular1, &err);

    F.f = &C2_g1_x_Pgg0_x_Pgg0_sing2_integrand;
    gsl_monte_vegas_integrate(&F, xl, xu, 2, calls, r, s, &singular2, &err);

    //double xl_new[] = {x, 0} ;

    F.f = &C2_g1_x_Pgg0_x_Pgg0_sing3_integrand;
    gsl_monte_vegas_integrate(&F, xl, xu, 2, calls, r, s, &singular3, &err);

    gsl_monte_vegas_free (s);

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double abserr = 0.001, relerr = 0.001;

    gsl_function f;
    f.function = &C2_g1_x_Pgg0_sing_integrand;
    f.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&f, x, 1, abserr, relerr, 1000, 4, w, &singular4, &err);

    gsl_set_error_handler (old_handler);

    singular4 *= Pgg0loc(nf) ;

    gsl_integration_workspace_free (w);

    return singular1 + singular2 + singular3 + singular4;

}

//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgg0_MC(double x, double m2Q2, int nf, size_t calls) {

    return (
        C2_g1_x_Pgg0_x_Pgg0_reg(x, m2Q2, nf, calls)
        + C2_g1_x_Pgg0_x_Pgg0_sing(x, m2Q2, nf, calls)
        + C2_g1_x_Pgg0(x, m2Q2, nf) * (Pgg0loc(nf) - Pgg0sing_integrated(x))
    );

}

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions for FL and the convolution
//  between the splitting functions Pgg0 and Pgg0 using monte carlo mathods
//------------------------------------------------------------------------------------------//

double CL_g1_x_Pgg0_x_Pgg0_reg1_integrand(double z[], size_t dim, void * p) {

    if (dim != 2) {
        std::cout << "error: dim != 2" << std::endl ;
        exit(-1);
    }

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    double z1 = z[0], z2 = z[1] ;

    double int_bound = theta(z1 - x) * theta(z2 - z1) ;

    return 1. / (z1 * z2) * Pgg0reg(x / z1) * Pgg0reg(z1 / z2) * CL_g1(z2, m2Q2) * int_bound ;

}

//------------------------------------------------------------------------------------------//

double CL_g1_x_Pgg0_x_Pgg0_reg2_integrand(double z[], size_t dim, void * p) {

    if (dim != 2) {
        std::cout << "error: dim != 2" << std::endl ;
        exit(-1);
    }

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    double z1 = z[0], z2 = z[1] ;

    double int_bound = theta(z1 - x) * theta(z2 - z1) ;

    return 1. / z1 * Pgg0reg(x / z1) * Pgg0sing(z2) * (CL_g1(z1 / z2, m2Q2) / z2 -  CL_g1(z1, m2Q2)) * int_bound ;

}

//------------------------------------------------------------------------------------------//

double CL_g1_x_Pgg0_x_Pgg0_reg3_integrand(double z[], size_t dim, void * p) {

    if (dim != 2) {
        std::cout << "error: dim != 2" << std::endl ;
        exit(-1);
    }

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    double z1 = z[0], z2 = z[1] ;

    double int_bound = theta(z1 - x) * theta(z1 - z2) ;

    return - 1. / z1 * Pgg0reg(x / z1) * CL_g1(z1, m2Q2) * Pgg0sing(z2) * int_bound ;

}

//------------------------------------------------------------------------------------------//

double CL_g1_x_Pgg0_x_Pgg0_reg(double x, double m2Q2, int nf, size_t calls) {

    struct function_params params = {x, m2Q2, nf} ;
    double xl[2] = {x, x};
    double xu[2] = {1, 1};

    double err, regular1, regular2, regular3, regular4 ;

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);

    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    gsl_monte_function F;

    F.f = &CL_g1_x_Pgg0_x_Pgg0_reg1_integrand;
    F.dim = 2;
    F.params = &params;

    gsl_monte_vegas_init (s);
    gsl_monte_vegas_integrate(&F, xl, xu, 2, calls, r, s, &regular1, &err);

    F.f = &CL_g1_x_Pgg0_x_Pgg0_reg2_integrand;
    gsl_monte_vegas_init (s);
    gsl_monte_vegas_integrate(&F, xl, xu, 2, calls, r, s, &regular2, &err);

    double xl_new[] = {x, 0} ;

    F.f = &CL_g1_x_Pgg0_x_Pgg0_reg3_integrand;
    gsl_monte_vegas_init (s);
    gsl_monte_vegas_integrate(&F, xl_new, xu, 2, calls, r, s, &regular3, &err);

    gsl_monte_vegas_free (s);

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double abserr = 0.001, relerr = 0.001;

    gsl_function f;
    f.function = &CL_g1_x_Pgg0_reg_integrand;
    f.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&f, x, 1, abserr, relerr, 1000, 4, w, &regular4, &err);

    gsl_set_error_handler (old_handler);

    regular4 *= Pgg0loc(nf) ;

    gsl_integration_workspace_free (w);

    return regular1 + regular2 + regular3 + regular4;

}

//------------------------------------------------------------------------------------------//

double CL_g1_x_Pgg0_x_Pgg0_sing1_integrand(double z[], size_t dim, void * p) {

    if (dim != 2) {
        std::cout << "error: dim != 2" << std::endl ;
        exit(-1);
    }

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    double z1 = z[0], z2 = z[1] ;

    return 1. / z2 * theta(z1 - x) * Pgg0sing(z1) * (
            theta(z2 - x / z1) / z1 * Pgg0reg(x / (z1 * z2)) - theta(z2 - x) * Pgg0reg(x / z2)
        ) * CL_g1(z2, m2Q2) ;

}

//------------------------------------------------------------------------------------------//

double CL_g1_x_Pgg0_x_Pgg0_sing2_integrand(double z[], size_t dim, void * p) {

    if (dim != 2) {
        std::cout << "error: dim != 2" << std::endl ;
        exit(-1);
    }

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    double z1 = z[0], z2 = z[1] ;

    return theta(z1 - x) * Pgg0sing(z1) * (
        Pgg0sing(z2) / z1 * (CL_g1(x / (z1 * z2), m2Q2) / z2 - CL_g1(x / z1, m2Q2)) * theta(z2 - x / z1)
        - Pgg0sing(z2) * (CL_g1(x / z2, m2Q2) / z2 - CL_g1(x, m2Q2)) * theta(z2 - x)
    ) ;

}

//------------------------------------------------------------------------------------------//

double CL_g1_x_Pgg0_x_Pgg0_sing3_integrand(double z[], size_t dim, void * p) {

    if (dim != 2) {
        std::cout << "error: dim != 2" << std::endl ;
        exit(-1);
    }

    struct function_params * params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    //int nf = (params->nf);

    double z1 = z[0], z2 = z[1] ;

    return - (theta(z1 - x) * Pgg0sing(z1) * Pgg0sing(z2) * (CL_g1(x / z1 , m2Q2) / z1 * theta(x / z1 - z2) - CL_g1(x , m2Q2) * theta(x - z2)));

}

//------------------------------------------------------------------------------------------//

double CL_g1_x_Pgg0_x_Pgg0_sing(double x, double m2Q2, int nf, size_t calls) {

    struct function_params params = {x, m2Q2, nf} ;
    double xl[2] = {x, 0};
    //double xl[2] = {x, x};
    double xu[2] = {1, 1};

    double err, singular1, singular2, singular3, singular4 ;

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);

    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup ();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    gsl_monte_function F;

    F.f = &CL_g1_x_Pgg0_x_Pgg0_sing1_integrand;
    F.dim = 2;
    F.params = &params;

    gsl_monte_vegas_integrate(&F, xl, xu, 2, calls, r, s, &singular1, &err);

    F.f = &CL_g1_x_Pgg0_x_Pgg0_sing2_integrand;
    gsl_monte_vegas_integrate(&F, xl, xu, 2, calls, r, s, &singular2, &err);

    //xl_new[] = {x, 0} ;

    F.f = &CL_g1_x_Pgg0_x_Pgg0_sing3_integrand;
    gsl_monte_vegas_integrate(&F, xl, xu, 2, calls, r, s, &singular3, &err);

    gsl_monte_vegas_free (s);

    gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

    double abserr = 0.001, relerr = 0.001;

    gsl_function f;
    f.function = &CL_g1_x_Pgg0_sing_integrand;
    f.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(&f, x, 1, abserr, relerr, 1000, 4, w, &singular4, &err);

    gsl_set_error_handler (old_handler);

    singular4 *= Pgg0loc(nf) ;

    gsl_integration_workspace_free (w);

    return singular1 + singular2 + singular3 + singular4;

}

//------------------------------------------------------------------------------------------//

double CL_g1_x_Pgg0_x_Pgg0_MC(double x, double m2Q2, int nf, size_t calls) {

    return (
        CL_g1_x_Pgg0_x_Pgg0_reg(x, m2Q2, nf, calls)
        + CL_g1_x_Pgg0_x_Pgg0_sing(x, m2Q2, nf, calls)
        + CL_g1_x_Pgg0(x, m2Q2, nf) * (Pgg0loc(nf) - Pgg0sing_integrated(x))
    );

}
