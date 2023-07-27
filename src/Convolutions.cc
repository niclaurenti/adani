#include "adani/Convolutions.h"
#include "adani/Constants.h"
#include "adani/ExactCoefficientFunctions.h"
#include "adani/MasslessCoefficientFunctions.h"
#include "adani/MatchingConditions.h"
#include "adani/SpecialFunctions.h"
#include "adani/SplittingFunctions.h"

#include <cmath>
#include <iostream>

#define REL 0.001
#define ABS 0.001
#define CALLS 25000
#define DIM 1000

// TODO : in all the numerical convolutions the gsl default error is first
// switched off and then swithced on again: see if it can be done once for all

//==========================================================================================//
//  Convolution between first order massless quark coefficient functions for F2
//  and first order matching KQg1
//------------------------------------------------------------------------------------------//

double C2_Q1_x_K_Qg1(double x, double m2mu2) {

    double x2 = x * x;
    double x1 = 1 - x;
    double L = log(x);
    double L1 = log(x1);

    double res =
        -5. / 2. + 2. * (3. - 4. * x) * x + zeta2 * (-1. + 2. * x - 4. * x2)
        + L1 * L1 * (1. - 2. * x1 * x)
        - 0.5 * L1 * (7. + 4. * x * (-4. + 3. * x) + (4. - 8. * x1 * x) * L)
        + 0.5 * L * (-1. + 4. * x * (3. * x - 2.) + (1. - 2. * x + 4. * x2) * L)
        + (2. * x - 1.) * Li2(x);

    return 4. * CF * TR * res * log(1. / m2mu2);
}

//==========================================================================================//
//  Convolution between first order massless quark coefficient functions for FL
//  and first order matching KQg1
//------------------------------------------------------------------------------------------//

double CL_Q1_x_K_Qg1(double x, double m2mu2) {

    double x2 = x * x;
    double L = log(x);

    return 8. * CF * TR * (1. + x - 2. * x2 + 2. * x * L) * log(1. / m2mu2);
}

//==========================================================================================//
//  Convolution between first order massive gluon coefficient functions for F2
//  and splitting function Pgq0
//------------------------------------------------------------------------------------------//

//  Integrand

double C2_g1_x_Pgq0_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    return C2_g1(z, m2Q2) * Pgq0(x / z) / z;
}

// Result

double C2_g1_x_Pgq0(double x, double m2Q2) {

    double x_max = 1. / (1. + 4 * m2Q2);
    double result, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, static_cast<int>(nan("")) };
    // It is not dependent on nf so nf is put to nan

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    gsl_function F;
    F.function = &C2_g1_x_Pgq0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &result, &error
    );

    gsl_set_error_handler(old_handler);

    gsl_integration_workspace_free(w);

    return result;
}

//==========================================================================================//
//  Convolution between first order massive gluon coefficient functions for FL
//  and splitting function Pgq0
//------------------------------------------------------------------------------------------//

// Integrand

double CL_g1_x_Pgq0_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    return CL_g1(z, m2Q2) * Pgq0(x / z) / z;
}

// Result

double CL_g1_x_Pgq0(double x, double m2Q2) {

    double x_max = 1. / (1. + 4 * m2Q2);
    double result, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, static_cast<int>(nan("")) };
    // It is not dependent on nf so nf is put to nan

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    gsl_function F;
    F.function = &CL_g1_x_Pgq0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &result, &error
    );

    gsl_set_error_handler(old_handler);

    gsl_integration_workspace_free(w);

    return result;
}

//==========================================================================================//
//  Convolution between first order massive gluon coefficient functions for F2
//  and splitting function Pgg0
//------------------------------------------------------------------------------------------//

// Integrand of the regular part

double C2_g1_x_Pgg0_reg_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    return C2_g1(z, m2Q2) * Pgg0reg(x / z) / z;
}

// Integrand of the singular part

double C2_g1_x_Pgg0_sing_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    return Pgg0sing(z) * (C2_g1(x / z, m2Q2) / z - C2_g1(x, m2Q2));
}

// Needed to integrate singular part

double Pgg0sing_integrated(double x) { return -Pgg0sing(0.) * log(1. - x); }

// Result

double C2_g1_x_Pgg0(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double regular, singular, local, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, nf };

    gsl_function F;
    F.function = &C2_g1_x_Pgg0_reg_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &regular, &error
    );

    F.function = &C2_g1_x_Pgg0_sing_integrand;
    gsl_integration_qag(
        &F, x / x_max, 1., abserr, relerr, DIM, 4, w, &singular, &error
    );

    gsl_set_error_handler(old_handler);

    local = C2_g1(x, m2Q2) * (Pgg0loc(nf) - Pgg0sing_integrated(x / x_max));

    gsl_integration_workspace_free(w);

    return regular + singular + local;
}

//==========================================================================================//
//  Convolution between first order massive gluon coefficient functions for F2
//  and splitting function Pgg0
//------------------------------------------------------------------------------------------//

// Integrand of the regular part

double CL_g1_x_Pgg0_reg_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    return CL_g1(z, m2Q2) * Pgg0reg(x / z) / z;
}

// Integrand of the singular part

double CL_g1_x_Pgg0_sing_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    return Pgg0sing(z) * (CL_g1(x / z, m2Q2) / z - CL_g1(x, m2Q2));
}

// Result

double CL_g1_x_Pgg0(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double regular, singular, local, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, nf };

    gsl_function F;
    F.function = &CL_g1_x_Pgg0_reg_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &regular, &error
    );

    F.function = &CL_g1_x_Pgg0_sing_integrand;
    gsl_integration_qag(
        &F, x / x_max, 1., abserr, relerr, DIM, 4, w, &singular, &error
    );

    gsl_set_error_handler(old_handler);

    local = CL_g1(x, m2Q2) * (Pgg0loc(nf) - Pgg0sing_integrated(x / x_max));

    gsl_integration_workspace_free(w);

    return regular + singular + local;
}

//==========================================================================================//
//  Convolution between first order massive gluon coefficient functions for F2
//  and splitting function Pgq1
//------------------------------------------------------------------------------------------//

// Integrand

double C2_g1_x_Pgq1_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return C2_g1(z, m2Q2) * Pgq1(x / z, nf) / z;
}

// Result

double C2_g1_x_Pgq1(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double result, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, nf };

    gsl_function F;
    F.function = &C2_g1_x_Pgq1_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &result, &error
    );

    gsl_set_error_handler(old_handler);

    gsl_integration_workspace_free(w);

    return result;
}

//==========================================================================================//
//  Convolution between first order massive gluon coefficient functions for FL
//  and splitting function Pgq1
//------------------------------------------------------------------------------------------//

// Integrand

double CL_g1_x_Pgq1_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return CL_g1(z, m2Q2) * Pgq1(x / z, nf) / z;
}

// Result

double CL_g1_x_Pgq1(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double result, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, nf };

    gsl_function F;
    F.function = &CL_g1_x_Pgq1_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &result, &error
    );

    gsl_set_error_handler(old_handler);

    gsl_integration_workspace_free(w);

    return result;
}

//==========================================================================================//
//  Convolution between the mu-independent term of the second order massive
//  gluon coefficient functions for F2 and splitting function Pgq0
//------------------------------------------------------------------------------------------//

// Integrand

double C2_g20_x_Pgq0_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    return C2_g20(z, m2Q2) * Pgq0(x / z) / z;
}

// Result

double C2_g20_x_Pgq0(double x, double m2Q2) {

    double x_max = 1. / (1. + 4 * m2Q2);
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double result, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, static_cast<int>(nan("")) };

    gsl_function F;
    F.function = &C2_g20_x_Pgq0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &result, &error
    );

    gsl_set_error_handler(old_handler);

    gsl_integration_workspace_free(w);

    return result;
}

//==========================================================================================//
//  Convolution between the mu-independent term of the second order massive
//  gluon coefficient functions for FL and splitting function Pgq0
//------------------------------------------------------------------------------------------//

// Integrand

double CL_g20_x_Pgq0_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    return CL_g20(z, m2Q2) * Pgq0(x / z) / z;
}

// Result

double CL_g20_x_Pgq0(double x, double m2Q2) {

    double x_max = 1. / (1. + 4 * m2Q2);
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double result, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, static_cast<int>(nan("")) };

    gsl_function F;
    F.function = &CL_g20_x_Pgq0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &result, &error
    );

    gsl_set_error_handler(old_handler);

    gsl_integration_workspace_free(w);

    return result;
}

//==========================================================================================//
//  Convolution between the mu-independent term of the second order massive
//  quark coefficient functions for F2 and splitting function Pqq0
//------------------------------------------------------------------------------------------//

// Integrand of the regular part

double C2_ps20_x_Pqq0_reg_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    return C2_ps20(z, m2Q2) * Pqq0reg(x / z) / z;
}

// Integrand of the singular part

double C2_ps20_x_Pqq0_sing_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    return Pqq0sing(z) * (C2_ps20(x / z, m2Q2) / z - C2_ps20(x, m2Q2));
}

// Needed to integrate singular part

double Pqq0sing_integrated(double x) { return -Pqq0sing(0.) * log(1. - x); }

// Result

double C2_ps20_x_Pqq0(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double regular, singular, local, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, nf };

    gsl_function F;
    F.function = &C2_ps20_x_Pqq0_reg_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &regular, &error
    );

    F.function = &C2_ps20_x_Pqq0_sing_integrand;
    gsl_integration_qag(
        &F, x / x_max, 1., abserr, relerr, DIM, 4, w, &singular, &error
    );

    gsl_set_error_handler(old_handler);

    local = C2_ps20(x, m2Q2) * (Pqq0loc() - Pqq0sing_integrated(x / x_max));

    gsl_integration_workspace_free(w);

    return regular + singular + local;
}

//==========================================================================================//
//  Convolution between the mu-independent term of the second order massive
//  quark coefficient functions for FL and splitting function Pqq0
//------------------------------------------------------------------------------------------//

// Integrand of the regular part

double CL_ps20_x_Pqq0_reg_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    return CL_ps20(z, m2Q2) * Pqq0reg(x / z) / z;
}

// Integrand of the singular part

double CL_ps20_x_Pqq0_sing_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    return Pqq0sing(z) * (CL_ps20(x / z, m2Q2) / z - CL_ps20(x, m2Q2));
}

// Result

double CL_ps20_x_Pqq0(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double regular, singular, local, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, nf };

    gsl_function F;
    F.function = &CL_ps20_x_Pqq0_reg_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &regular, &error
    );

    F.function = &CL_ps20_x_Pqq0_sing_integrand;
    gsl_integration_qag(
        &F, x / x_max, 1., abserr, relerr, DIM, 4, w, &singular, &error
    );

    gsl_set_error_handler(old_handler);

    local = CL_ps20(x, m2Q2) * (Pqq0loc() - Pqq0sing_integrated(x / x_max));

    gsl_integration_workspace_free(w);

    return regular + singular + local;
}

//==========================================================================================//
//  Analytical onvolution between the splitting functions Pgg0 and Pgq0
//------------------------------------------------------------------------------------------//

double Pgg0_x_Pgq0(double x, int nf) {

    return (-4. * CF * nf * (2. + (-2. + x) * x)
            + 2. * CA * CF
                  * (-40. + x * (26. + x * (17. + 8. * x))
                     + 12. * (2. + (-2. + x) * x) * log(1. - x)
                     - 24. * (1. + x + x * x) * log(x)))
           / 3. / x;
}

//==========================================================================================//
//  Analytical onvolution between the splitting functions Pqq0 and Pgq0
//------------------------------------------------------------------------------------------//

double Pqq0_x_Pgq0(double x) {

    return (2. * CF * CF
            * (4. * (2. + (-2. + x) * x) * log(1. - x)
               - x * (-4. + x + 2. * (-2. + x) * log(x))))
           / x;
}

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions for
//  F2 and the convolution between the splitting functions Pgg0 and Pgq0
//------------------------------------------------------------------------------------------//

// Integrand

double C2_g1_x_Pgg0_x_Pgq0_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return C2_g1(z, m2Q2) * Pgg0_x_Pgq0(x / z, nf) / z;
}

// Result

double C2_g1_x_Pgg0_x_Pgq0(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double result, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, nf };

    gsl_function F;
    F.function = &C2_g1_x_Pgg0_x_Pgq0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &result, &error
    );

    gsl_set_error_handler(old_handler);

    gsl_integration_workspace_free(w);

    return result;
}

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions for
//  FL and the convolution between the splitting functions Pgg0 and Pgq0
//------------------------------------------------------------------------------------------//

// Integrand

double CL_g1_x_Pgg0_x_Pgq0_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return CL_g1(z, m2Q2) * Pgg0_x_Pgq0(x / z, nf) / z;
}

// Result

double CL_g1_x_Pgg0_x_Pgq0(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double result, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, nf };

    gsl_function F;
    F.function = &CL_g1_x_Pgg0_x_Pgq0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &result, &error
    );

    gsl_set_error_handler(old_handler);

    gsl_integration_workspace_free(w);

    return result;
}

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions for
//  F2 and the convolution between the splitting functions Pqq0 and Pgq0
//------------------------------------------------------------------------------------------//

// Integrand

double C2_g1_x_Pqq0_x_Pgq0_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    return C2_g1(z, m2Q2) * Pqq0_x_Pgq0(x / z) / z;
}

// Result

double C2_g1_x_Pqq0_x_Pgq0(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double result, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, nf };

    gsl_function F;
    F.function = &C2_g1_x_Pqq0_x_Pgq0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &result, &error
    );

    gsl_set_error_handler(old_handler);

    gsl_integration_workspace_free(w);

    return result;
}

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions for
//  FL and the convolution between the splitting functions Pqq0 and Pgq0
//------------------------------------------------------------------------------------------//

// Integrand

double CL_g1_x_Pqq0_x_Pgq0_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    return CL_g1(z, m2Q2) * Pqq0_x_Pgq0(x / z) / z;
}

// Result

double CL_g1_x_Pqq0_x_Pgq0(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double result, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, nf };

    gsl_function F;
    F.function = &CL_g1_x_Pqq0_x_Pgq0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &result, &error
    );

    gsl_set_error_handler(old_handler);

    gsl_integration_workspace_free(w);

    return result;
}

//==========================================================================================//
//  Convolution between first order massive gluon coefficient functions for F2
//  and splitting function Pgg1
//------------------------------------------------------------------------------------------//

// Integrand of the regular part

double C2_g1_x_Pgg1_reg_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return C2_g1(z, m2Q2) * Pgg1reg(x / z, nf) / z;
}

// Integrand of the singular part

double C2_g1_x_Pgg1_sing_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return Pgg1sing(z, nf) * (C2_g1(x / z, m2Q2) / z - C2_g1(x, m2Q2));
}

// Needed to integrate singular part

double Pgg1sing_integrated(double x, int nf) {

    return -Pgg1sing(0., nf) * log(1. - x);
}

// Result

double C2_g1_x_Pgg1(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    double regular, singular, local, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, nf };

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    gsl_function F;
    F.function = &C2_g1_x_Pgg1_reg_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &regular, &error
    );

    F.function = &C2_g1_x_Pgg1_sing_integrand;
    gsl_integration_qag(
        &F, x / x_max, 1., abserr, relerr, DIM, 4, w, &singular, &error
    );

    gsl_set_error_handler(old_handler);

    local = C2_g1(x, m2Q2) * (Pgg1loc(nf) - Pgg1sing_integrated(x / x_max, nf));

    gsl_integration_workspace_free(w);

    return regular + singular + local;
}

//==========================================================================================//
//  Convolution between first order massive gluon coefficient functions for F2
//  and splitting function Pgg1
//------------------------------------------------------------------------------------------//

// Integrand of the regular part

double CL_g1_x_Pgg1_reg_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return CL_g1(z, m2Q2) * Pgg1reg(x / z, nf) / z;
}

// Integrand of the regular part

double CL_g1_x_Pgg1_sing_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return Pgg1sing(z, nf) * (CL_g1(x / z, m2Q2) / z - CL_g1(x, m2Q2));
}

// Result

double CL_g1_x_Pgg1(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    double regular, singular, local, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, nf };

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    gsl_function F;
    F.function = &CL_g1_x_Pgg1_reg_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &regular, &error
    );

    F.function = &CL_g1_x_Pgg1_sing_integrand;
    gsl_integration_qag(
        &F, x / x_max, 1., abserr, relerr, DIM, 4, w, &singular, &error
    );

    gsl_set_error_handler(old_handler);

    local = CL_g1(x, m2Q2) * (Pgg1loc(nf) - Pgg1sing_integrated(x / x_max, nf));

    gsl_integration_workspace_free(w);

    return regular + singular + local;
}

//==========================================================================================//
//  Convolution between the mu-independent term of the second order massive
//  quark coefficient functions for F2 and splitting function Pqg0
//------------------------------------------------------------------------------------------//

// Integrand

double C2_ps20_x_Pqg0_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return C2_ps20(z, m2Q2) * Pqg0(x / z, nf) / z;
}

// Result

double C2_ps20_x_Pqg0(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    double result, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, nf };

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    gsl_function F;
    F.function = &C2_ps20_x_Pqg0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &result, &error
    );

    gsl_set_error_handler(old_handler);

    gsl_integration_workspace_free(w);

    return result;
}

//==========================================================================================//
//  Convolution between the mu-independent term of the second order massive
//  quark coefficient functions for FL and splitting function Pqg0
//------------------------------------------------------------------------------------------//

// Integrand

double CL_ps20_x_Pqg0_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return CL_ps20(z, m2Q2) * Pqg0(x / z, nf) / z;
}

// Result

double CL_ps20_x_Pqg0(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    double result, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, nf };

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    gsl_function F;
    F.function = &CL_ps20_x_Pqg0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &result, &error
    );

    gsl_set_error_handler(old_handler);

    gsl_integration_workspace_free(w);

    return result;
}

//==========================================================================================//
//  Convolution between the mu-independent term of the second order massive
//  gluon coefficient functions for F2 and splitting function Pgg0
//------------------------------------------------------------------------------------------//

// Integrand of the regular part

double C2_g20_x_Pgg0_reg_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    return C2_g20(z, m2Q2) * Pgg0reg(x / z) / z;
}

// Integrand of the singular part

double C2_g20_x_Pgg0_sing_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    return Pgg0sing(z) * (C2_g20(x / z, m2Q2) / z - C2_g20(x, m2Q2));
}

// Result

double C2_g20_x_Pgg0(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double regular, singular, local, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, nf };

    gsl_function F;
    F.function = &C2_g20_x_Pgg0_reg_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &regular, &error
    );

    F.function = &C2_g20_x_Pgg0_sing_integrand;
    gsl_integration_qag(
        &F, x / x_max, 1., abserr, relerr, DIM, 4, w, &singular, &error
    );

    gsl_set_error_handler(old_handler);

    local = C2_g20(x, m2Q2) * (Pgg0loc(nf) - Pgg0sing_integrated(x / x_max));

    gsl_integration_workspace_free(w);

    return regular + singular + local;
}

//==========================================================================================//
//  Convolution between the mu-independent term of the second order massive
//  gluon coefficient functions for FL and splitting function Pgg0
//------------------------------------------------------------------------------------------//

// Integrand of the regular part

double CL_g20_x_Pgg0_reg_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    return CL_g20(z, m2Q2) * Pgg0reg(x / z) / z;
}

// Integrand of the singular part

double CL_g20_x_Pgg0_sing_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    return Pgg0sing(z) * (CL_g20(x / z, m2Q2) / z - CL_g20(x, m2Q2));
}

// Result

double CL_g20_x_Pgg0(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double regular, singular, local, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, nf };

    gsl_function F;
    F.function = &CL_g20_x_Pgg0_reg_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &regular, &error
    );

    F.function = &CL_g20_x_Pgg0_sing_integrand;
    gsl_integration_qag(
        &F, x / x_max, 1., abserr, relerr, DIM, 4, w, &singular, &error
    );

    gsl_set_error_handler(old_handler);

    local = CL_g20(x, m2Q2) * (Pgg0loc(nf) - Pgg0sing_integrated(x / x_max));

    gsl_integration_workspace_free(w);

    return regular + singular + local;
}

//==========================================================================================//
//  Analytical onvolution between the splitting functions Pqg0 and Pgq0
//------------------------------------------------------------------------------------------//

double Pqg0_x_Pgq0(double x, int nf) {

    return 4. * CF * nf
           * (1. + 4. / 3 / x - x - 4. * x * x / 3 + 2. * (1 + x) * log(x));
}

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions for
//  F2 and the convolution between the splitting functions Pqg0 and Pgq0
//------------------------------------------------------------------------------------------//

// Integrand

double C2_g1_x_Pqg0_x_Pgq0_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return C2_g1(z, m2Q2) * Pqg0_x_Pgq0(x / z, nf) / z;
}

// Result

double C2_g1_x_Pqg0_x_Pgq0(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double result, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, nf };

    gsl_function F;
    F.function = &C2_g1_x_Pqg0_x_Pgq0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &result, &error
    );

    gsl_set_error_handler(old_handler);

    gsl_integration_workspace_free(w);

    return result;
}

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions for
//  FL and the convolution between the splitting functions Pqg0 and Pgq0
//------------------------------------------------------------------------------------------//

// Integrand

double CL_g1_x_Pqg0_x_Pgq0_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return CL_g1(z, m2Q2) * Pqg0_x_Pgq0(x / z, nf) / z;
}

// Result

double CL_g1_x_Pqg0_x_Pgq0(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double result, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, nf };

    gsl_function F;
    F.function = &CL_g1_x_Pqg0_x_Pgq0_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &result, &error
    );

    gsl_set_error_handler(old_handler);

    gsl_integration_workspace_free(w);

    return result;
}

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions for
//  F2 and the convolution between the splitting functions Pgg0 and Pgg0
//------------------------------------------------------------------------------------------//

// Integrand of the regular part

double C2_g1_x_Pgg0_x_Pgg0_reg_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return C2_g1_x_Pgg0(z, m2Q2, nf) * Pgg0reg(x / z) / z;
}

// Integrand of the singular part

double C2_g1_x_Pgg0_x_Pgg0_sing_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return Pgg0sing(z)
           * (C2_g1_x_Pgg0(x / z, m2Q2, nf) / z - C2_g1_x_Pgg0(x, m2Q2, nf));
}

// Result

double C2_g1_x_Pgg0_x_Pgg0(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double regular, singular, local, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, nf };

    double C2_g1xPgg0 = C2_g1_x_Pgg0(x, m2Q2, nf);

    gsl_function F;
    F.function = &C2_g1_x_Pgg0_x_Pgg0_reg_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &regular, &error
    );

    F.function = &C2_g1_x_Pgg0_x_Pgg0_sing_integrand;
    gsl_integration_qag(
        &F, x / x_max, 1., abserr, relerr, DIM, 4, w, &singular, &error
    );

    gsl_set_error_handler(old_handler);

    local = C2_g1xPgg0 * (Pgg0loc(nf) - Pgg0sing_integrated(x / x_max));

    gsl_integration_workspace_free(w);

    return regular + singular + local;
}

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions for
//  FL and the convolution between the splitting functions Pgg0 and Pgg0
//------------------------------------------------------------------------------------------//

// Integrand of the regular part

double CL_g1_x_Pgg0_x_Pgg0_reg_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return CL_g1_x_Pgg0(z, m2Q2, nf) * Pgg0reg(x / z) / z;
}

// Integrand of the singular part

double CL_g1_x_Pgg0_x_Pgg0_sing_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return Pgg0sing(z)
           * (CL_g1_x_Pgg0(x / z, m2Q2, nf) / z - CL_g1_x_Pgg0(x, m2Q2, nf));
}

// Result

double CL_g1_x_Pgg0_x_Pgg0(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double regular, singular, local, error, abserr = ABS, relerr = REL;
    struct function_params params = { x, m2Q2, nf };

    double CL_g1xPgg0 = CL_g1_x_Pgg0(x, m2Q2, nf);

    gsl_function F;
    F.function = &CL_g1_x_Pgg0_x_Pgg0_reg_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, DIM, 4, w, &regular, &error
    );

    F.function = &CL_g1_x_Pgg0_x_Pgg0_sing_integrand;
    gsl_integration_qag(
        &F, x / x_max, 1., abserr, relerr, DIM, 4, w, &singular, &error
    );

    gsl_set_error_handler(old_handler);

    local = CL_g1xPgg0 * (Pgg0loc(nf) - Pgg0sing_integrated(x / x_max));

    gsl_integration_workspace_free(w);

    return regular + singular + local;
}

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions for
//  F2 and the convolution between the splitting functions Pgg0 and Pgg0 using
//  monte carlo methods
//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgg0_reg1_integrand(double z[], size_t dim, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    double z1 = z[0], z2 = z[1];

    if (z2 > z1) {
        return 1. / (z1 * z2) * Pgg0reg(x / z1) * Pgg0reg(z1 / z2)
               * C2_g1(z2, m2Q2);
    } else {
        return 0.;
    }
}

//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgg0_reg2_integrand(double z[], size_t dim, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    double z1 = z[0], z2 = z[1];

    if (z2 > z1) {
        return 1. / (z1 * z2) * Pgg0reg(x / z1) * Pgg0sing(z1 / z2)
               * (C2_g1(z2, m2Q2) - z1 / z2 * C2_g1(z1, m2Q2));
    } else {
        return 0.;
    }
}

//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgg0_reg3_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);
    double x_max = 1. / (1. + 4 * m2Q2);

    return -1. / z * Pgg0reg(x / z) * C2_g1(z, m2Q2)
           * Pgg0sing_integrated(z / x_max);
}

//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgg0_reg(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    struct function_params params = { x, m2Q2, nf };
    double xl[2] = { x, x };
    double xu[2] = { x_max, x_max };

    double err, regular1, regular2, regular3, regular4;
    size_t calls = CALLS;

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);

    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_function F;

    F.f = &C2_g1_x_Pgg0_x_Pgg0_reg1_integrand;
    F.dim = 2;
    F.params = &params;

    gsl_monte_vegas_init(s);
    gsl_monte_vegas_integrate(&F, xl, xu, 2, calls, r, s, &regular1, &err);

    F.f = &C2_g1_x_Pgg0_x_Pgg0_reg2_integrand;
    gsl_monte_vegas_init(s);
    gsl_monte_vegas_integrate(&F, xl, xu, 2, calls, r, s, &regular2, &err);

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double abserr = ABS, relerr = REL;

    gsl_function f;
    f.function = &C2_g1_x_Pgg0_x_Pgg0_reg3_integrand;
    f.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &f, x, x_max, abserr, relerr, DIM, 4, w, &regular3, &err
    );

    f.function = &C2_g1_x_Pgg0_reg_integrand;

    gsl_integration_qag(
        &f, x, x_max, abserr, relerr, DIM, 4, w, &regular4, &err
    );

    gsl_set_error_handler(old_handler);

    regular4 *= Pgg0loc(nf);

    gsl_integration_workspace_free(w);

    return regular1 + regular2 + regular3 + regular4;
}

//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgg0_sing1_integrand(double z[], size_t dim, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    double z1 = z[0], z2 = z[1];

    double tmp;
    if (z2 > x / z1) {
        tmp = Pgg0reg(x / (z1 * z2)) / z1;
    } else {
        tmp = 0.;
    }

    return Pgg0sing(z1) * (tmp - Pgg0reg(x / z2)) * C2_g1(z2, m2Q2) / z2;
}

//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgg0_sing2_integrand(double z[], size_t dim, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    double x_max = 1. / (1. + 4 * m2Q2);

    double z1 = z[0], z2 = z[1];

    double tmp;
    if (z2 > x / (x_max * z1)) {
        tmp = (C2_g1(x / (z1 * z2), m2Q2) / z2 - C2_g1(x / z1, m2Q2)) / z1;
    } else {
        tmp = 0.;
    }

    return Pgg0sing(z1) * Pgg0sing(z2)
           * (tmp - (C2_g1(x / z2, m2Q2) / z2 - C2_g1(x, m2Q2)));
}

//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgg0_sing3_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    double x_max = 1. / (1. + 4 * m2Q2);

    return -(
        Pgg0sing(z)
        * (C2_g1(x / z, m2Q2) * Pgg0sing_integrated(x / (x_max * z)) / z
           - C2_g1(x, m2Q2) * Pgg0sing_integrated(x / x_max))
    );
}

//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgg0_sing(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    struct function_params params = { x, m2Q2, nf };

    double xl[2] = { x / x_max, x };
    double xu[2] = { 1, x_max };

    double err, singular1, singular2, singular3, singular4;
    size_t calls = CALLS;

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);

    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_function F;

    F.f = &C2_g1_x_Pgg0_x_Pgg0_sing1_integrand;
    F.dim = 2;
    F.params = &params;

    gsl_monte_vegas_init(s);
    gsl_monte_vegas_integrate(&F, xl, xu, 2, calls, r, s, &singular1, &err);

    xl[1] = x / x_max;
    xu[1] = 1;

    F.f = &C2_g1_x_Pgg0_x_Pgg0_sing2_integrand;
    gsl_monte_vegas_init(s);
    gsl_monte_vegas_integrate(&F, xl, xu, 2, calls, r, s, &singular2, &err);

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double abserr = ABS, relerr = REL;

    gsl_function f;
    f.function = &C2_g1_x_Pgg0_x_Pgg0_sing3_integrand;
    f.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &f, x / x_max, 1, abserr, relerr, DIM, 4, w, &singular3, &err
    );

    f.function = &C2_g1_x_Pgg0_sing_integrand;

    gsl_integration_qag(
        &f, x / x_max, 1, abserr, relerr, DIM, 4, w, &singular4, &err
    );

    gsl_set_error_handler(old_handler);

    singular4 *= Pgg0loc(nf);

    gsl_integration_workspace_free(w);

    return singular1 + singular2 + singular3 + singular4;
}

//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgg0_MC(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);

    return C2_g1_x_Pgg0_x_Pgg0_reg(x, m2Q2, nf)
           + C2_g1_x_Pgg0_x_Pgg0_sing(x, m2Q2, nf)
           + C2_g1_x_Pgg0(x, m2Q2, nf)
                 * (Pgg0loc(nf) - Pgg0sing_integrated(x / x_max));
}

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions for
//  FL and the convolution between the splitting functions Pgg0 and Pgg0 using
//  monte carlo mathods
//------------------------------------------------------------------------------------------//

double CL_g1_x_Pgg0_x_Pgg0_reg1_integrand(double z[], size_t dim, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    double z1 = z[0], z2 = z[1];

    if (z2 > z1) {
        return 1. / (z1 * z2) * Pgg0reg(x / z1) * Pgg0reg(z1 / z2)
               * CL_g1(z2, m2Q2);
    } else {
        return 0.;
    }
}

//------------------------------------------------------------------------------------------//

double CL_g1_x_Pgg0_x_Pgg0_reg2_integrand(double z[], size_t dim, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    double z1 = z[0], z2 = z[1];

    if (z2 > z1) {
        return 1. / (z1 * z2) * Pgg0reg(x / z1) * Pgg0sing(z1 / z2)
               * (CL_g1(z2, m2Q2) - z1 / z2 * CL_g1(z1, m2Q2));
    } else {
        return 0.;
    }
}

//------------------------------------------------------------------------------------------//

double CL_g1_x_Pgg0_x_Pgg0_reg3_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);
    double x_max = 1. / (1. + 4 * m2Q2);

    return -1. / z * Pgg0reg(x / z) * CL_g1(z, m2Q2)
           * Pgg0sing_integrated(z / x_max);
}

//------------------------------------------------------------------------------------------//

double CL_g1_x_Pgg0_x_Pgg0_reg(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    struct function_params params = { x, m2Q2, nf };
    double xl[2] = { x, x };
    double xu[2] = { x_max, x_max };

    double err, regular1, regular2, regular3, regular4;
    size_t calls = CALLS;

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);

    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_function F;

    F.f = &CL_g1_x_Pgg0_x_Pgg0_reg1_integrand;
    F.dim = 2;
    F.params = &params;

    gsl_monte_vegas_init(s);
    gsl_monte_vegas_integrate(&F, xl, xu, 2, calls, r, s, &regular1, &err);

    F.f = &CL_g1_x_Pgg0_x_Pgg0_reg2_integrand;
    gsl_monte_vegas_init(s);
    gsl_monte_vegas_integrate(&F, xl, xu, 2, calls, r, s, &regular2, &err);

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double abserr = ABS, relerr = REL;

    gsl_function f;
    f.function = &CL_g1_x_Pgg0_x_Pgg0_reg3_integrand;
    f.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &f, x, x_max, abserr, relerr, DIM, 4, w, &regular3, &err
    );

    f.function = &CL_g1_x_Pgg0_reg_integrand;

    gsl_integration_qag(
        &f, x, x_max, abserr, relerr, DIM, 4, w, &regular4, &err
    );

    gsl_set_error_handler(old_handler);

    regular4 *= Pgg0loc(nf);

    gsl_integration_workspace_free(w);

    return regular1 + regular2 + regular3 + regular4;
}

//------------------------------------------------------------------------------------------//

double CL_g1_x_Pgg0_x_Pgg0_sing1_integrand(double z[], size_t dim, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    double z1 = z[0], z2 = z[1];

    double tmp;
    if (z2 > x / z1) {
        tmp = Pgg0reg(x / (z1 * z2)) / z1;
    } else {
        tmp = 0.;
    }

    return Pgg0sing(z1) * (tmp - Pgg0reg(x / z2)) * CL_g1(z2, m2Q2) / z2;
}

//------------------------------------------------------------------------------------------//

double CL_g1_x_Pgg0_x_Pgg0_sing2_integrand(double z[], size_t dim, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    double x_max = 1. / (1. + 4 * m2Q2);

    double z1 = z[0], z2 = z[1];

    double tmp;
    if (z2 > x / (x_max * z1)) {
        tmp = (CL_g1(x / (z1 * z2), m2Q2) / z2 - CL_g1(x / z1, m2Q2)) / z1;
    } else {
        tmp = 0.;
    }

    return Pgg0sing(z1) * Pgg0sing(z2)
           * (tmp - (CL_g1(x / z2, m2Q2) / z2 - CL_g1(x, m2Q2)));
}

//------------------------------------------------------------------------------------------//

double CL_g1_x_Pgg0_x_Pgg0_sing3_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    // int nf = (params->nf);

    double x_max = 1. / (1. + 4 * m2Q2);

    return -(
        Pgg0sing(z)
        * (CL_g1(x / z, m2Q2) * Pgg0sing_integrated(x / (x_max * z)) / z
           - CL_g1(x, m2Q2) * Pgg0sing_integrated(x / x_max))
    );
}

//------------------------------------------------------------------------------------------//

double CL_g1_x_Pgg0_x_Pgg0_sing(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);
    struct function_params params = { x, m2Q2, nf };

    double xl[2] = { x / x_max, x };
    double xu[2] = { 1, x_max };

    double err, singular1, singular2, singular3, singular4;
    size_t calls = CALLS;

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);

    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_function F;

    F.f = &CL_g1_x_Pgg0_x_Pgg0_sing1_integrand;
    F.dim = 2;
    F.params = &params;

    gsl_monte_vegas_init(s);
    gsl_monte_vegas_integrate(&F, xl, xu, 2, calls, r, s, &singular1, &err);

    xl[1] = x / x_max;
    xu[1] = 1;

    F.f = &CL_g1_x_Pgg0_x_Pgg0_sing2_integrand;
    gsl_monte_vegas_init(s);
    gsl_monte_vegas_integrate(&F, xl, xu, 2, calls, r, s, &singular2, &err);

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(DIM);

    double abserr = ABS, relerr = REL;

    gsl_function f;
    f.function = &CL_g1_x_Pgg0_x_Pgg0_sing3_integrand;
    f.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &f, x / x_max, 1, abserr, relerr, DIM, 4, w, &singular3, &err
    );

    f.function = &CL_g1_x_Pgg0_sing_integrand;

    gsl_integration_qag(
        &f, x / x_max, 1, abserr, relerr, DIM, 4, w, &singular4, &err
    );

    gsl_set_error_handler(old_handler);

    singular4 *= Pgg0loc(nf);

    gsl_integration_workspace_free(w);

    return singular1 + singular2 + singular3 + singular4;
}

//------------------------------------------------------------------------------------------//

double CL_g1_x_Pgg0_x_Pgg0_MC(double x, double m2Q2, int nf) {

    double x_max = 1. / (1. + 4 * m2Q2);

    return CL_g1_x_Pgg0_x_Pgg0_reg(x, m2Q2, nf)
           + CL_g1_x_Pgg0_x_Pgg0_sing(x, m2Q2, nf)
           + CL_g1_x_Pgg0(x, m2Q2, nf)
                 * (Pgg0loc(nf) - Pgg0sing_integrated(x / x_max));
}
