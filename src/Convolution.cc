#include "adani/Convolution.h"

#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

// TODO : in all the numerical convolutions the gsl default error is first
// switched off and then swithced on again: see if it can be done once for all

//==========================================================================================//
//  AbstractConvolution: constructor
//------------------------------------------------------------------------------------------//

AbstractConvolution::AbstractConvolution(
    CoefficientFunction *coefffunc, AbstractSplittingFunction *splitfunc,
    const double &abserr, const double &relerr, const int &dim
) {
    abserr_ = abserr;
    relerr_ = relerr;
    dim_ = dim;

    w_ = gsl_integration_workspace_alloc(dim);

    coefffunc_ = coefffunc;
    splitfunc_ = splitfunc;
}

//==========================================================================================//
//  AbstractConvolution: destructor
//------------------------------------------------------------------------------------------//

AbstractConvolution::~AbstractConvolution() {

    gsl_integration_workspace_free(w_);
};

//==========================================================================================//
//  AbstractConvolution: set method for abserr
//------------------------------------------------------------------------------------------//

void AbstractConvolution::SetAbserr(const double &abserr) {
    // check abserr
    if (abserr <= 0) {
        cout << "Error: abserr must be positive. Got " << abserr << endl;
        exit(-1);
    }
    abserr_ = abserr;
}

//==========================================================================================//
//  AbstractConvolution: set method for relerr
//------------------------------------------------------------------------------------------//

void AbstractConvolution::SetRelerr(const double &relerr) {
    // check relerr
    if (relerr <= 0) {
        cout << "Error: relerr must be positive. Got " << relerr << endl;
        exit(-1);
    }
    relerr_ = relerr;
}

//==========================================================================================//
//  AbstractConvolution: set method for dim
//------------------------------------------------------------------------------------------//

void AbstractConvolution::SetDim(const int &dim) {
    // check dim
    if (dim <= 0) {
        cout << "Error: dim must be positive. Got " << dim << endl;
        exit(-1);
    }
    dim_ = dim;
}

//==========================================================================================//
//  AbstractConvolution: convolute splitting function with coefficient function
//------------------------------------------------------------------------------------------//

double AbstractConvolution::Convolute(double x, double m2Q2, int nf) const {
    return RegularPart(x, m2Q2, nf) + SingularPart(x, m2Q2, nf)
           + LocalPart(x, m2Q2, nf);
}

//==========================================================================================//
//  AbstractConvolution: integrand of the regular part
//------------------------------------------------------------------------------------------//

double Convolution::regular_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    CoefficientFunction *cf = (params->conv)->GetCoeffFunc();
    AbstractSplittingFunction *sf = (params->conv)->GetSplitFunc();

    return cf->MuIndependentTerms(z, m2Q2, nf) * sf->Regular(x / z, nf) / z;
}

//==========================================================================================//
//  Convolution: integrand of the singular part
//------------------------------------------------------------------------------------------//

double Convolution::singular_integrand(double z, void *p) {
    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    CoefficientFunction *cf = (params->conv)->GetCoeffFunc();
    AbstractSplittingFunction *sf = (params->conv)->GetSplitFunc();

    return sf->Singular(z, nf)
           * (cf->MuIndependentTerms(x / z, m2Q2, nf) / z
              - cf->MuIndependentTerms(x, m2Q2, nf));
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

//==========================================================================================//
//  Convolution: regular part
//------------------------------------------------------------------------------------------//

double Convolution::RegularPart(double x, double m2Q2, int nf) const {
    double x_max = 1. / (1. + 4 * m2Q2);

    double abserr = GetAbserr();
    double relerr = GetRelerr();
    int dim = GetDim();
    gsl_integration_workspace *w = GetWorkspace();

    double regular, error;

    struct function_params params = { x, m2Q2, nf, this };

    gsl_function F;
    F.function = &regular_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x, x_max, abserr, relerr, dim, 4, w, &regular, &error
    );

    gsl_set_error_handler(old_handler);

    return regular;
}

//==========================================================================================//
//  Convolution: singular part
//------------------------------------------------------------------------------------------//

double Convolution::SingularPart(double x, double m2Q2, int nf) const {

    double x_max = 1. / (1. + 4 * m2Q2);

    double abserr = GetAbserr();
    double relerr = GetRelerr();
    int dim = GetDim();
    gsl_integration_workspace *w = GetWorkspace();

    double singular, error;

    struct function_params params = { x, m2Q2, nf, this };

    gsl_function F;
    F.function = &singular_integrand;
    F.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &F, x / x_max, 1., abserr, relerr, dim, 4, w, &singular, &error
    );

    gsl_set_error_handler(old_handler);

    return singular;
}

//==========================================================================================//
//  Convolution: local part
//------------------------------------------------------------------------------------------//

double Convolution::LocalPart(double x, double m2Q2, int nf) const {

    double x_max = 1. / (1. + 4 * m2Q2);

    return coefffunc_->MuIndependentTerms(x, m2Q2, nf)
           * (splitfunc_->Local(nf)
              - splitfunc_->SingularIntegrated(x / x_max, nf));
}

//==========================================================================================//
//  ConvolutedCoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

ConvolutedCoefficientFunction::ConvolutedCoefficientFunction(
    CoefficientFunction *coefffunc, AbstractSplittingFunction *splitfunc,
    const double &abserr, const double &relerr, const int &dim
)
    : CoefficientFunction(coefffunc) {

    conv_ = new Convolution(coefffunc, splitfunc, abserr, relerr, dim);
}

//==========================================================================================//
//  ConvolutedCoefficientFunction: destructor
//------------------------------------------------------------------------------------------//

ConvolutedCoefficientFunction::~ConvolutedCoefficientFunction() {

    delete conv_;
}

//==========================================================================================//
//  ConvolutedCoefficientFunction: mu independent terms
//------------------------------------------------------------------------------------------//

double ConvolutedCoefficientFunction::MuIndependentTerms(
    double x, double m2Q2, int nf
) const {

    return conv_->Convolute(x, m2Q2, nf);
}

//==========================================================================================//
//  ConvolutedCoefficientFunction: central value of fx
//------------------------------------------------------------------------------------------//

double ConvolutedCoefficientFunction::fx(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return MuIndependentTerms(x, m2Q2, nf)
           + MuDependentTerms(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  ConvolutedCoefficientFunction: band of fx
//------------------------------------------------------------------------------------------//

Value ConvolutedCoefficientFunction::fxBand(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return Value(fx(x, m2Q2, m2mu2, nf));
}

//==========================================================================================//
//  DoubleConvolution: constructor
//------------------------------------------------------------------------------------------//

DoubleConvolution::DoubleConvolution(
    CoefficientFunction *coefffunc, AbstractSplittingFunction *splitfunc,
    const double &abserr, const double &relerr, const int &dim,
    const bool &MCintegral, const int &MCcalls
)
    : AbstractConvolution(coefffunc, splitfunc, abserr, relerr, dim) {

    SetMCintegral(MCintegral);
    SetMCcalls(MCcalls);

    if (MCintegral) {
        convolution_ =
            new Convolution(coefffunc, splitfunc, abserr, relerr, dim);
        conv_coeff_ = nullptr;

        s_ = gsl_monte_vegas_alloc(2);
        gsl_rng_env_setup();
        r_ = gsl_rng_alloc(gsl_rng_default);

    } else {
        conv_coeff_ = new ConvolutedCoefficientFunction(
            coefffunc, splitfunc, abserr, relerr, dim
        );
        convolution_ =
            new Convolution(conv_coeff_, splitfunc, abserr, relerr, dim);
    }
}

//==========================================================================================//
//  DoubleConvolution: destructor
//------------------------------------------------------------------------------------------//

DoubleConvolution::~DoubleConvolution() {

    if (MCintegral_) {
        gsl_monte_vegas_free(s_);
        gsl_rng_free(r_);
    }

    delete convolution_;
    delete conv_coeff_;
}

//==========================================================================================//
//  DoubleConvolution: set method for method_flag
//------------------------------------------------------------------------------------------//

void DoubleConvolution::SetMCintegral(const bool &MCintegral) {

    MCintegral_ = MCintegral;
}

//==========================================================================================//
//  DoubleConvolution: set method for MCcalls
//------------------------------------------------------------------------------------------//

void DoubleConvolution::SetMCcalls(const int &MCcalls) {
    // check dim
    if (MCcalls <= 0) {
        cout << "Error: MCcalls must be positive. Got " << MCcalls << endl;
        exit(-1);
    }
    MCcalls_ = MCcalls;
}

//==========================================================================================//
//  DoubleConvolution: integrand of the first regular part
//------------------------------------------------------------------------------------------//

double
DoubleConvolution::regular1_integrand(double z[], size_t /*dim*/, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    CoefficientFunction *coefffunc = (params->conv)->GetCoeffFunc();
    AbstractSplittingFunction *splitfunc = (params->conv)->GetSplitFunc();

    double z1 = z[0], z2 = z[1];

    if (z2 > z1) {
        return 1. / (z1 * z2) * splitfunc->Regular(x / z1, nf)
               * splitfunc->Regular(z1 / z2, nf)
               * coefffunc->MuIndependentTerms(z2, m2Q2, nf);
    } else {
        return 0.;
    }
}

//==========================================================================================//
//  DoubleConvolution: integrand of the second regular part
//------------------------------------------------------------------------------------------//

double
DoubleConvolution::regular2_integrand(double z[], size_t /*dim*/, void *p) {
    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    CoefficientFunction *coefffunc = (params->conv)->GetCoeffFunc();
    AbstractSplittingFunction *splitfunc = (params->conv)->GetSplitFunc();

    double z1 = z[0], z2 = z[1];

    if (z2 > z1) {
        return 1. / (z1 * z2) * splitfunc->Regular(x / z1, nf)
               * splitfunc->Singular(z1 / z2, nf)
               * (coefffunc->MuIndependentTerms(z2, m2Q2, nf)
                  - z1 / z2 * coefffunc->MuIndependentTerms(z1, m2Q2, nf));
    } else {
        return 0.;
    }
}

//==========================================================================================//
//  DoubleConvolution: integrand of the third regular part
//------------------------------------------------------------------------------------------//

double DoubleConvolution::regular3_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    CoefficientFunction *coefffunc = (params->conv)->GetCoeffFunc();
    AbstractSplittingFunction *splitfunc = (params->conv)->GetSplitFunc();

    double x_max = 1. / (1. + 4 * m2Q2);

    return -1. / z * splitfunc->Regular(x / z, nf)
           * coefffunc->MuIndependentTerms(z, m2Q2, nf)
           * splitfunc->SingularIntegrated(z / x_max, nf);
}

//==========================================================================================//
//  DoubleConvolution: regular part
//------------------------------------------------------------------------------------------//

double DoubleConvolution::RegularPart(double x, double m2Q2, int nf) const {

    if (!MCintegral_)
        return convolution_->RegularPart(x, m2Q2, nf);

    double x_max = 1. / (1. + 4 * m2Q2);
    struct function_params params = { x, m2Q2, nf, this };

    double xl[2] = { x, x };
    double xu[2] = { x_max, x_max };

    double err, regular1, regular2, regular3, regular4;

    gsl_monte_function F;

    F.f = &regular1_integrand;
    F.dim = 2;
    F.params = &params;

    gsl_monte_vegas_init(s_);
    gsl_monte_vegas_integrate(&F, xl, xu, 2, MCcalls_, r_, s_, &regular1, &err);

    F.f = &regular2_integrand;
    gsl_monte_vegas_init(s_);
    gsl_monte_vegas_integrate(&F, xl, xu, 2, MCcalls_, r_, s_, &regular2, &err);

    double abserr = GetAbserr();
    double relerr = GetRelerr();
    int dim = GetDim();
    gsl_integration_workspace *w = GetWorkspace();

    gsl_function f;
    f.function = &regular3_integrand;
    f.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &f, x, x_max, abserr, relerr, dim, 4, w, &regular3, &err
    );

    gsl_set_error_handler(old_handler);

    regular4 = convolution_->RegularPart(x, m2Q2, nf) * splitfunc_->Local(nf);

    return regular1 + regular2 + regular3 + regular4;
}

//==========================================================================================//
//  DoubleConvolution: integrand of the first singular part
//------------------------------------------------------------------------------------------//

double
DoubleConvolution::singular1_integrand(double z[], size_t /*dim*/, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    CoefficientFunction *coefffunc = (params->conv)->GetCoeffFunc();
    AbstractSplittingFunction *splitfunc = (params->conv)->GetSplitFunc();

    double z1 = z[0], z2 = z[1];

    double tmp;
    if (z2 > x / z1) {
        tmp = splitfunc->Regular(x / (z1 * z2), nf) / z1;
    } else {
        tmp = 0.;
    }

    return splitfunc->Singular(z1, nf) * (tmp - splitfunc->Regular(x / z2, nf))
           * coefffunc->MuIndependentTerms(z2, m2Q2, nf) / z2;
}

//==========================================================================================//
//  DoubleConvolution: integrand of the second singular part
//------------------------------------------------------------------------------------------//

double
DoubleConvolution::singular2_integrand(double z[], size_t /*dim*/, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    CoefficientFunction *coefffunc = (params->conv)->GetCoeffFunc();
    AbstractSplittingFunction *splitfunc = (params->conv)->GetSplitFunc();

    double x_max = 1. / (1. + 4 * m2Q2);

    double z1 = z[0], z2 = z[1];

    double tmp;
    if (z2 > x / (x_max * z1)) {
        tmp = (coefffunc->MuIndependentTerms(x / (z1 * z2), m2Q2, nf) / z2
               - coefffunc->MuIndependentTerms(x / z1, m2Q2, nf))
              / z1;
    } else {
        tmp = 0.;
    }

    return splitfunc->Singular(z1, nf) * splitfunc->Singular(z2, nf)
           * (tmp
              - (coefffunc->MuIndependentTerms(x / z2, m2Q2, nf) / z2
                 - coefffunc->MuIndependentTerms(x, m2Q2, nf)));
}

//==========================================================================================//
//  DoubleConvolution: integrand of the third singular part
//------------------------------------------------------------------------------------------//

double DoubleConvolution::singular3_integrand(double z, void *p) {

    struct function_params *params = (struct function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    CoefficientFunction *coefffunc = (params->conv)->GetCoeffFunc();
    AbstractSplittingFunction *splitfunc = (params->conv)->GetSplitFunc();

    double x_max = 1. / (1. + 4 * m2Q2);

    return -(
        splitfunc->Singular(z, nf)
        * (coefffunc->MuIndependentTerms(x / z, m2Q2, nf)
               * splitfunc->SingularIntegrated(x / (x_max * z), nf) / z
           - coefffunc->MuIndependentTerms(x, m2Q2, nf)
                 * splitfunc->SingularIntegrated(x / x_max, nf))
    );
}

//==========================================================================================//
//  DoubleConvolution: singular part
//------------------------------------------------------------------------------------------//

double DoubleConvolution::SingularPart(double x, double m2Q2, int nf) const {

    if (!MCintegral_)
        return convolution_->SingularPart(x, m2Q2, nf);

    double x_max = 1. / (1. + 4 * m2Q2);
    struct function_params params = { x, m2Q2, nf, this };

    double xl[2] = { x / x_max, x };
    double xu[2] = { 1, x_max };

    double err, singular1, singular2, singular3, singular4;

    gsl_monte_function F;

    F.f = &singular1_integrand;
    F.dim = 2;
    F.params = &params;

    gsl_monte_vegas_init(s_);
    gsl_monte_vegas_integrate(
        &F, xl, xu, 2, MCcalls_, r_, s_, &singular1, &err
    );

    xl[1] = x / x_max;
    xu[1] = 1;

    F.f = &singular2_integrand;
    gsl_monte_vegas_init(s_);
    gsl_monte_vegas_integrate(
        &F, xl, xu, 2, MCcalls_, r_, s_, &singular2, &err
    );

    double abserr = GetAbserr();
    double relerr = GetRelerr();
    int dim = GetDim();
    gsl_integration_workspace *w = GetWorkspace();

    gsl_function f;
    f.function = &singular3_integrand;
    f.params = &params;

    gsl_error_handler_t *old_handler = gsl_set_error_handler(NULL);
    gsl_set_error_handler_off();

    gsl_integration_qag(
        &f, x / x_max, 1, abserr, relerr, dim, 4, w, &singular3, &err
    );

    gsl_set_error_handler(old_handler);

    singular4 = convolution_->SingularPart(x, m2Q2, nf) * splitfunc_->Local(nf);

    return singular1 + singular2 + singular3 + singular4;
}

//==========================================================================================//
//  DoubleConvolution: local part
//------------------------------------------------------------------------------------------//

double DoubleConvolution::LocalPart(double x, double m2Q2, int nf) const {

    if (!MCintegral_)
        return convolution_->LocalPart(x, m2Q2, nf);

    double x_max = 1. / (1. + 4 * m2Q2);

    return convolution_->Convolute(x, m2Q2, nf)
           * (splitfunc_->Local(nf)
              - splitfunc_->SingularIntegrated(x / x_max, nf));
}
