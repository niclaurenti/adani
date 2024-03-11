#include "adani/Convolutions.h"
#include "adani/Constants.h"
#include "adani/ExactCoefficientFunctions.h"
#include "adani/MasslessCoefficientFunctions.h"
#include "adani/MatchingConditions.h"
#include "adani/SpecialFunctions.h"
#include "adani/SplittingFunctions.h"

#include <cmath>
#include <iostream>

// #define REL 0.001
// #define ABS 0.001
#define CALLS 25000
// #define DIM 1000

// TODO : in all the numerical convolutions the gsl default error is first
// switched off and then swithced on again: see if it can be done once for all

Convolution::Convolution(const CoefficientFunction& coefffunc, const SplittingFunction& splitfunc, const double& abserr, const double& relerr, const int& dim) {
    abserr_ = abserr;
    relerr_ = relerr;
    dim_ = dim;
    
    coefffunc_ = &coefffunc;
    splitfunc_ = &splitfunc;
}

Convolution::Convolution(const CoefficientFunction* coefffunc, const SplittingFunction* splitfunc, const double& abserr, const double& relerr, const int& dim) {
    abserr_ = abserr;
    relerr_ = relerr;
    dim_ = dim;
    
    coefffunc_ = coefffunc;
    splitfunc_ = splitfunc;
}

Convolution::~Convolution() {
    
    delete coefffunc_;
    delete splitfunc_;
    
}

void Convolution::SetAbserr(const double& abserr) {
    // check abserr
    if (abserr <= 0) {
        cout << "Error: abserr must be positive. Got " << abserr << endl;
        exit(-1);
    }
    abserr_ = abserr;
}

void Convolution::SetRelerr(const double& relerr) {
    // check relerr
    if (relerr <= 0) {
        cout << "Error: relerr must be positive. Got " << relerr << endl;
        exit(-1);
    }
    relerr_ = relerr;
}

void Convolution::SetDim(const int& dim) {
    // check dim
    if (dim <= 0) {
        cout << "Error: MCcalls must be positive. Got " << dim << endl;
        exit(-1);
    }
    dim_ = dim;
}

double Convolution::regular_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return coefffunc_ -> MuIndependentTerms(z, m2Q2, nf) * splitfunc_ -> Regular(x / z, nf) / z ;
}

double Convolution::singular_integrand(double z, void *p) {
    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return splitfunc_ -> Regular(z, nf) * (coefffunc_-> MuIndependentTerms(x / z, m2Q2, nf) - coefffunc_ -> MuIndependentTerms(x, m2Q2, nf) ) ;

}

//==========================================================================================//
//  WARNING: from now on all the convolutions are numerical.
//  In all the integrals we have switched off the ERROR
//  raised by gsl when the numerical precision is not
//  achieved. This is not ideal but as a wise man once said:
//  "If one has no better options, then he has sex even with
//  his own wife" (the identity of this enlightened person
//  will not be revealed)
//------------------------------------------------------------------------------------------//

double Convolution::convolute(double x, double m2Q2, int nf) const {
    
    double x_max = 1. / (1. + 4 * m2Q2);
    
    // TODO: consider moving w in the data members
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(dim_);
    
    double regular, singular, local, error;
    struct function_params params = { x, m2Q2, nf };

    gsl_function F;
    F.function = &regular_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr_, relerr_, dim_, 4, w, &regular, &error
    );

    F.function = &singular_integrand;
    gsl_integration_qag(
        &F, x / x_max, 1., abserr_, relerr_, dim_, 4, w, &singular, &error
    );

    gsl_set_error_handler(old_handler);

    local = coefffunc_->MuIndependentTerms(x, m2Q2, m2mu2, nf) * (splitfunc_->Local(nf) - splitfunc_->SingularIntegrated(x / x_max, nf));

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
