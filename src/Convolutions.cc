#include "adani/Convolutions.h"
#include "adani/Constants.h"
#include "adani/ExactCoefficientFunctions.h"
#include "adani/MasslessCoefficientFunctions.h"
#include "adani/MatchingConditions.h"
#include "adani/SpecialFunctions.h"
#include "adani/SplittingFunctions.h"

#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

// TODO : in all the numerical convolutions the gsl default error is first
// switched off and then swithced on again: see if it can be done once for all

// AbstractConvolution::AbstractConvolution(const CoefficientFunction& coefffunc, const SplittingFunction& splitfunc, const double& abserr, const double& relerr, const int& dim) {
//     SetAbserr(abserr);
//     SetRelerr(relerr);
//     SetDim(dim);

//     coefffunc_ = &coefffunc;
//     splitfunc_ = &splitfunc;
// }

AbstractConvolution::AbstractConvolution(CoefficientFunction* coefffunc, SplittingFunction* splitfunc, const double& abserr, const double& relerr, const int& dim) {
    abserr_ = abserr;
    relerr_ = relerr;
    dim_ = dim;

    coefffunc_ = coefffunc;
    splitfunc_ = splitfunc;
}

void AbstractConvolution::SetAbserr(const double& abserr) {
    // check abserr
    if (abserr <= 0) {
        cout << "Error: abserr must be positive. Got " << abserr << endl;
        exit(-1);
    }
    abserr_ = abserr;
}

void AbstractConvolution::SetRelerr(const double& relerr) {
    // check relerr
    if (relerr <= 0) {
        cout << "Error: relerr must be positive. Got " << relerr << endl;
        exit(-1);
    }
    relerr_ = relerr;
}

void AbstractConvolution::SetDim(const int& dim) {
    // check dim
    if (dim <= 0) {
        cout << "Error: dim must be positive. Got " << dim << endl;
        exit(-1);
    }
    dim_ = dim;
}

double AbstractConvolution::Convolute(double x, double m2Q2, int nf) const {
    return RegularPart(x, m2Q2, nf) + SingularPart(x, m2Q2, nf) + LocalPart(x, m2Q2, nf) ;
}


// double regular_integrand(double z, void *p) {

//     struct function_params *params = (struct function_params *)p;

//     double m2Q2 = (params->m2Q2);
//     double x = (params->x);
//     int nf = (params->nf);

//     return coefffunc_ -> MuIndependentTerms(z, m2Q2, nf) * splitfunc_ -> Regular(x / z, nf) / z ;
// }

// double singular_integrand(double z, void *p) {

//     struct function_params *params = (struct function_params *)p;

//     double m2Q2 = (params->m2Q2);
//     double x = (params->x);
//     int nf = (params->nf);

//     return splitfunc_ -> Regular(z, nf) * (coefffunc_-> MuIndependentTerms(x / z, m2Q2, nf) - coefffunc_ -> MuIndependentTerms(x, m2Q2, nf) ) ;

// }

double Convolution::regular_integrand(double z, void *p) const {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return coefffunc_ -> MuIndependentTerms(z, m2Q2, nf) * splitfunc_ -> Regular(x / z, nf) / z ;
}

double Convolution::singular_integrand(double z, void *p) const {
    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return splitfunc_ -> Singular(z, nf) * (coefffunc_-> MuIndependentTerms(x / z, m2Q2, nf) - coefffunc_ -> MuIndependentTerms(x, m2Q2, nf) ) ;

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

double Convolution::RegularPart(double x, double m2Q2, int nf) const {
    double x_max = 1. / (1. + 4 * m2Q2);

    double abserr = GetAbserr();
    double relerr = GetRelerr();
    int dim = GetDim();

    // TODO: consider moving w in the data members
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(dim);

    double regular, error;

    struct function_params params = { x, m2Q2, nf };

    gsl_function F;
    F.function = &regular_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, dim, 4, w, &regular, &error
    );

    gsl_integration_workspace_free(w);

    return regular;
}

double Convolution::SingularPart(double x, double m2Q2, int nf) const {

    double x_max = 1. / (1. + 4 * m2Q2);

    double abserr = GetAbserr();
    double relerr = GetRelerr();
    int dim = GetDim();

    // TODO: consider moving w in the data members
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(dim);

    double singular, error;

    struct function_params params = { x, m2Q2, nf };

    gsl_function F;
    F.function = &singular_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x / x_max, 1., abserr, relerr, dim, 4, w, &singular, &error
    );

    gsl_set_error_handler(old_handler);

    gsl_integration_workspace_free(w);

    return singular;
}

double Convolution::LocalPart(double x, double m2Q2, int nf) const {

    double x_max = 1. / (1. + 4 * m2Q2);

    double local, error;

    local = coefffunc_ -> MuIndependentTerms(x, m2Q2, nf) * (splitfunc_ -> Local(nf) - splitfunc_ -> SingularIntegrated(x / x_max, nf));

    return local;

}

MonteCarloDoubleConvolution::MonteCarloDoubleConvolution(CoefficientFunction* coefffunc, SplittingFunction* splitfunc, const double& abserr, const double& relerr, const int& dim, const int& MCcalls) : AbstractConvolution(coefffunc, splitfunc, abserr, relerr, dim) {

    SetMCcalls(MCcalls);
    convolution_ = new Convolution(coefffunc, splitfunc, abserr, relerr, dim);
}

MonteCarloDoubleConvolution::~MonteCarloDoubleConvolution() {

    delete convolution_;
}

void MonteCarloDoubleConvolution::SetMCcalls(const int& MCcalls) {
    // check dim
    if (MCcalls <= 0) {
        cout << "Error: MCcalls must be positive. Got " << MCcalls << endl;
        exit(-1);
    }
    MCcalls_ = MCcalls;
}

double MonteCarloDoubleConvolution::regular1_integrand(double z[], size_t dim, void *p) const {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    double z1 = z[0], z2 = z[1];

    if (z2 > z1) {
        return 1. / (z1 * z2) * splitfunc_ -> Regular(x / z1, nf) * splitfunc_ -> Regular(z1 / z2, nf)
               * coefffunc_-> MuIndependentTerms(z2, m2Q2, nf);
    } else {
        return 0.;
    }
}

double MonteCarloDoubleConvolution::regular2_integrand(double z[], size_t dim, void *p) const {
    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    double z1 = z[0], z2 = z[1];

    if (z2 > z1) {
        return 1. / (z1 * z2) * splitfunc_ -> Regular(x / z1, nf) * splitfunc_ -> Singular(z1 / z2, nf)
               * (coefffunc_-> MuIndependentTerms(z2, m2Q2, nf) - z1 / z2 * coefffunc_-> MuIndependentTerms(z1, m2Q2, nf));
    } else {
        return 0.;
    }
}

double MonteCarloDoubleConvolution::regular3_integrand(double z, void *p) const {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    double x_max = 1. / (1. + 4 * m2Q2);

    return -1. / z * splitfunc_ -> Regular(x / z, nf) * coefffunc_-> MuIndependentTerms(z, m2Q2, nf)
           * splitfunc_ -> SingularIntegrated(z / x_max, nf);
}

double MonteCarloDoubleConvolution::RegularPart(double x, double m2Q2, int nf) const {

    double x_max = 1. / (1. + 4 * m2Q2);
    struct function_params params = { x, m2Q2, nf };

    double xl[2] = { x, x };
    double xu[2] = { x_max, x_max };

    double err, regular1, regular2, regular3, regular4;

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);

    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_function F;

    F.f = &regular1_integrand;
    F.dim = 2;
    F.params = &params;

    gsl_monte_vegas_init(s);
    gsl_monte_vegas_integrate(&F, xl, xu, 2, MCcalls_, r, s, &regular1, &err);

    F.f = &regular2_integrand;
    gsl_monte_vegas_init(s);
    gsl_monte_vegas_integrate(&F, xl, xu, 2, MCcalls_, r, s, &regular2, &err);

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);

    double abserr = GetAbserr();
    double relerr = GetRelerr();
    int dim = GetDim();

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(dim);

    gsl_function f;
    f.function = &regular3_integrand;
    f.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &f, x, x_max, abserr, relerr, dim, 4, w, &regular3, &err
    );

    gsl_set_error_handler(old_handler);

    gsl_integration_workspace_free(w);

    regular4 = convolution_ -> RegularPart(x, m2Q2, nf) * splitfunc_ -> Local(nf);

    return regular1 + regular2 + regular3 + regular4;
}

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions for
//  F2 and the convolution between the splitting functions Pgg0 and Pgg0 using
//  monte carlo methods
//------------------------------------------------------------------------------------------//

double MonteCarloDoubleConvolution::singular1_integrand(double z[], size_t dim, void *p) const {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    double z1 = z[0], z2 = z[1];

    double tmp;
    if (z2 > x / z1) {
        tmp = splitfunc_ -> Regular(x / (z1 * z2), nf) / z1;
    } else {
        tmp = 0.;
    }

    return splitfunc_ -> Singular(z1, nf) * (tmp - splitfunc_ -> Regular(x / z2, nf)) * coefffunc_ -> MuIndependentTerms(z2, m2Q2, nf) / z2;
}

double MonteCarloDoubleConvolution::singular2_integrand(double z[], size_t dim, void *p) const {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    double x_max = 1. / (1. + 4 * m2Q2);

    double z1 = z[0], z2 = z[1];

    double tmp;
    if (z2 > x / (x_max * z1)) {
        tmp = (coefffunc_ -> MuIndependentTerms(x / (z1 * z2), m2Q2, nf) / z2 - coefffunc_ -> MuIndependentTerms(x / z1, m2Q2, nf)) / z1;
    } else {
        tmp = 0.;
    }

    return splitfunc_ -> Singular(z1, nf) * splitfunc_ -> Singular(z2, nf)
           * (tmp - (coefffunc_ -> MuIndependentTerms(x / z2, m2Q2, nf) / z2 - coefffunc_ -> MuIndependentTerms(x, m2Q2, nf)));
}

double MonteCarloDoubleConvolution::singular3_integrand(double z, void *p) const {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    double x_max = 1. / (1. + 4 * m2Q2);

    return -(
        splitfunc_ -> Singular(z, nf)
        * (coefffunc_ -> MuIndependentTerms(x / z, m2Q2, nf) * splitfunc_ -> SingularIntegrated(x / (x_max * z), nf) / z
           - coefffunc_ -> MuIndependentTerms(x, m2Q2, nf) * splitfunc_ -> SingularIntegrated(x / x_max, nf))
    );
}

double MonteCarloDoubleConvolution::SingularPart(double x, double m2Q2, int nf) const {

    double x_max = 1. / (1. + 4 * m2Q2);
    struct function_params params = { x, m2Q2, nf };

    double xl[2] = { x / x_max, x };
    double xu[2] = { 1, x_max };

    double err, singular1, singular2, singular3, singular4;

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);

    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_function F;

    F.f = &singular1_integrand;
    F.dim = 2;
    F.params = &params;

    gsl_monte_vegas_init(s);
    gsl_monte_vegas_integrate(&F, xl, xu, 2, MCcalls_, r, s, &singular1, &err);

    xl[1] = x / x_max;
    xu[1] = 1;

    F.f = &singular2_integrand;
    gsl_monte_vegas_init(s);
    gsl_monte_vegas_integrate(&F, xl, xu, 2, MCcalls_, r, s, &singular2, &err);

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);

    double abserr = GetAbserr();
    double relerr = GetRelerr();
    int dim = GetDim();

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(dim);

    gsl_function f;
    f.function = &singular3_integrand;
    f.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &f, x / x_max, 1, abserr, relerr, dim, 4, w, &singular3, &err
    );

    gsl_set_error_handler(old_handler);

    singular4 = convolution_ -> SingularPart(x, m2Q2, nf) * splitfunc_ -> Local(nf);

    gsl_integration_workspace_free(w);

    return singular1 + singular2 + singular3 + singular4;
}

double MonteCarloDoubleConvolution::LocalPart(double x, double m2Q2, int nf) const {

    double x_max = 1. / (1. + 4 * m2Q2);

    return convolution_ -> Convolute(x, m2Q2, nf) * (splitfunc_ -> Local(nf) - splitfunc_ -> SingularIntegrated(x / x_max, nf));
}

//------------------------------------------------------------------------------------------//
