#include "adani/Convolution.h"

#include <cmath>
#include <future>

//==========================================================================================//
//  AbstractConvolution: constructor
//------------------------------------------------------------------------------------------//

AbstractConvolution::AbstractConvolution(
    std::shared_ptr<const CoefficientFunction> coefffunc,
    std::shared_ptr<const AbstractSplittingFunction> splitfunc,
    const double &abserr, const double &relerr, const int &dim
)
    : dim_(dim) {

    try {
        SetAbserr(abserr);
        SetRelerr(relerr);
        SetDim(dim);
    } catch (NotValidException &e) {
        e.runtime_error();
    }

    coefffunc_ = coefffunc;
    splitfunc_ = splitfunc;
}

//==========================================================================================//
//  AbstractConvolution: destructor
//------------------------------------------------------------------------------------------//

AbstractConvolution::~AbstractConvolution() = default;

//==========================================================================================//
//  AbstractConvolution: copy operator
//------------------------------------------------------------------------------------------//

AbstractConvolution &
    AbstractConvolution::operator=(const AbstractConvolution &obj) {
    if (this != &obj) {
        abserr_ = obj.abserr_;
        relerr_ = obj.relerr_;
        dim_ = obj.dim_;

        coefffunc_ = obj.coefffunc_;
        splitfunc_ = obj.splitfunc_;
    }
    return *this;
}

//==========================================================================================//
//  AbstractConvolution: set method for abserr
//------------------------------------------------------------------------------------------//

void AbstractConvolution::SetAbserr(const double &abserr) {
    // check abserr
    if (abserr <= 0) {
        throw NotValidException(
            "abserr must be positive. Got abserr=" + to_string(abserr),
            __PRETTY_FUNCTION__, __LINE__
        );
    }
    abserr_ = abserr;
}

//==========================================================================================//
//  AbstractConvolution: set method for relerr
//------------------------------------------------------------------------------------------//

void AbstractConvolution::SetRelerr(const double &relerr) {
    // check relerr
    if (relerr <= 0) {
        throw NotValidException(
            "relerr must be positive. Got relerr=" + to_string(relerr),
            __PRETTY_FUNCTION__, __LINE__
        );
    }
    relerr_ = relerr;
}

//==========================================================================================//
//  AbstractConvolution: method for the allocation of workspace
//------------------------------------------------------------------------------------------//

void AbstractConvolution::SetDim(const int &dim) {
    // check dim
    if (dim <= 0) {
        throw NotValidException(
            "dim must be positive. Got dim=" + to_string(dim),
            __PRETTY_FUNCTION__, __LINE__
        );
    }
    dim_ = dim;
}

//==========================================================================================//
//  AbstractConvolution: convolute splitting function with coefficient function
//------------------------------------------------------------------------------------------//

double AbstractConvolution::Convolute(double x, double m2Q2, int nf) const {

    std::future<double> future_f1 = std::async(
        std::launch::async, &AbstractConvolution::RegularPart, this, x, m2Q2, nf
    );
    std::future<double> future_f2 = std::async(
        std::launch::async, &AbstractConvolution::SingularPart, this, x, m2Q2,
        nf
    );
    std::future<double> future_f3 = std::async(
        std::launch::async, &AbstractConvolution::LocalPart, this, x, m2Q2, nf
    );

    return future_f1.get() + future_f2.get() + future_f3.get();
}

//==========================================================================================//
//  struct function_params to be passed to gsl
//------------------------------------------------------------------------------------------//

struct function_params {
        double x;
        double m2Q2;
        int nf;
        const AbstractConvolution *conv;
};

//==========================================================================================//
//  AbstractConvolution: integrand of the regular part
//------------------------------------------------------------------------------------------//

double Convolution::regular_integrand(double z, void *p) {

    function_params *params = (function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return (params->conv)->GetCoeffFunc()->MuIndependentTerms(z, m2Q2, nf)
           * (params->conv)->GetSplitFunc()->Regular(x / z, nf) / z;
}

//==========================================================================================//
//  Convolution: initialize static data members
//------------------------------------------------------------------------------------------//

int Convolution::number_of_instances = 0;
gsl_error_handler_t *Convolution::old_handler_ = nullptr;

//==========================================================================================//
//  Convolution: constructor
//------------------------------------------------------------------------------------------//

Convolution::Convolution(
    std::shared_ptr<const CoefficientFunction> coefffunc,
    std::shared_ptr<const AbstractSplittingFunction> splitfunc,
    const double &abserr, const double &relerr, const int &dim
)
    : AbstractConvolution(coefffunc, splitfunc, abserr, relerr, dim) {
    number_of_instances++;
    if (number_of_instances == 1) {
        old_handler_ = gsl_set_error_handler(NULL);
        gsl_set_error_handler_off();
    }
};

//==========================================================================================//
//  Convolution: copy constructor
//------------------------------------------------------------------------------------------//

Convolution::Convolution(const Convolution &obj)
    : Convolution(
          obj.GetCoeffFunc(), obj.GetSplitFunc(), obj.GetAbsErr(),
          obj.GetRelErr(), obj.GetDim()
      ) {}

//==========================================================================================//
//  Convolution: destructor
//------------------------------------------------------------------------------------------//

Convolution::~Convolution() {
    number_of_instances--;
    if (number_of_instances == 0) {
        gsl_set_error_handler(old_handler_);
    }
}

//==========================================================================================//
//  Convolution: copy operator
//------------------------------------------------------------------------------------------//

Convolution &Convolution::operator=(const Convolution &obj) {
    if (this != &obj) {
        AbstractConvolution::operator=(obj);
    }

    return *this;
}

//==========================================================================================//
//  Convolution: integrand of the singular part
//------------------------------------------------------------------------------------------//

double Convolution::singular_integrand(double z, void *p) {
    function_params *params = (function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    return (params->conv)->GetSplitFunc()->Singular(z, nf)
           * ((params->conv)
                      ->GetCoeffFunc()
                      ->MuIndependentTerms(x / z, m2Q2, nf)
                  / z
              - (params->conv)->GetCoeffFunc()->MuIndependentTerms(x, m2Q2, nf)
           );
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
    double x_max = CoefficientFunction::xMax(m2Q2);

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(GetDim());

    double regular, error;

    function_params params = { x, m2Q2, nf, this };

    gsl_function F;
    F.function = &regular_integrand;
    F.params = &params;

    gsl_integration_qag(
        &F, x, x_max, GetAbsErr(), GetRelErr(), GetDim(), 4, w, &regular, &error
    );

    gsl_integration_workspace_free(w);

    return regular;
}

//==========================================================================================//
//  Convolution: singular part
//------------------------------------------------------------------------------------------//

double Convolution::SingularPart(double x, double m2Q2, int nf) const {

    double x_max = CoefficientFunction::xMax(m2Q2);

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(GetDim());

    double singular, error;

    function_params params = { x, m2Q2, nf, this };

    gsl_function F;
    F.function = &singular_integrand;
    F.params = &params;

    gsl_integration_qag(
        &F, x / x_max, 1., GetAbsErr(), GetRelErr(), GetDim(), 4, w, &singular,
        &error
    );

    gsl_integration_workspace_free(w);

    return singular;
}

//==========================================================================================//
//  Convolution: local part
//------------------------------------------------------------------------------------------//

double Convolution::LocalPart(double x, double m2Q2, int nf) const {

    double x_max = CoefficientFunction::xMax(m2Q2);

    return GetCoeffFunc()->MuIndependentTerms(x, m2Q2, nf)
           * (GetSplitFunc()->Local(nf)
              - GetSplitFunc()->SingularIntegrated(x / x_max, nf));
}

//==========================================================================================//
//  ConvolutedCoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

ConvolutedCoefficientFunction::ConvolutedCoefficientFunction(
    std::shared_ptr<const CoefficientFunction> coefffunc,
    std::shared_ptr<const AbstractSplittingFunction> splitfunc,
    const double &abserr, const double &relerr, const int &dim
)
    : CoefficientFunction(
          coefffunc->GetOrder(), coefffunc->GetKind(), coefffunc->GetChannel()
      ) {

    conv_ = std::make_shared<Convolution>(
        coefffunc, splitfunc, abserr, relerr, dim
    );
}

//==========================================================================================//
//  ConvolutedCoefficientFunction: copy constructor
//------------------------------------------------------------------------------------------//

ConvolutedCoefficientFunction::ConvolutedCoefficientFunction(
    const ConvolutedCoefficientFunction &obj
)
    : ConvolutedCoefficientFunction(
          obj.GetConv()->GetCoeffFunc(), obj.GetConv()->GetSplitFunc(),
          obj.GetConv()->GetAbsErr(), obj.GetConv()->GetRelErr(),
          obj.GetConv()->GetDim()
      ) {}

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
    std::shared_ptr<const CoefficientFunction> coefffunc,
    std::shared_ptr<const AbstractSplittingFunction> splitfunc,
    const double &abserr, const double &relerr, const int &dim,
    const bool &MCintegral, const int &MCcalls
)
    : AbstractConvolution(coefffunc, splitfunc, abserr, relerr, dim),
      MCintegral_(MCintegral) {

    try {
        SetMCcalls(MCcalls);
    } catch (NotValidException &e) {
        e.runtime_error();
    }

    if (MCintegral) {
        convolution_ = std::make_unique<const Convolution>(
            coefffunc, splitfunc, abserr, relerr, dim
        );
        conv_coeff_ = nullptr;

    } else {
        conv_coeff_ = std::make_shared<const ConvolutedCoefficientFunction>(
            coefffunc, splitfunc, abserr, relerr, dim
        );
        convolution_ = std::make_unique<const Convolution>(
            conv_coeff_, splitfunc, abserr, relerr, dim
        );
    }
}

//==========================================================================================//
//  Convolution: copy constructor
//------------------------------------------------------------------------------------------//

DoubleConvolution::DoubleConvolution(const DoubleConvolution &obj)
    : DoubleConvolution(
          obj.GetCoeffFunc(), obj.GetSplitFunc(), obj.GetAbsErr(),
          obj.GetRelErr(), obj.GetDim(), obj.GetMCintegral(), obj.GetMCcalls()
      ) {}

//==========================================================================================//
//  DoubleConvolution: copy operator
//------------------------------------------------------------------------------------------//

DoubleConvolution &DoubleConvolution::operator=(const DoubleConvolution &obj) {
    if (this != &obj) {
        AbstractConvolution::operator=(obj);

        MCintegral_ = obj.MCintegral_;
        MCcalls_ = obj.MCcalls_;

        conv_coeff_ = obj.conv_coeff_;
        convolution_ = std::make_unique<const Convolution>(*obj.convolution_);
    }

    return *this;
}

//==========================================================================================//
//  DoubleConvolution: set method for MCcalls
//------------------------------------------------------------------------------------------//

void DoubleConvolution::SetMCcalls(const int &MCcalls) {
    // check MCcalls
    if (MCcalls <= 0) {
        throw NotValidException(
            "MCcalls must be positive. Got MCcalls=" + to_string(MCcalls),
            __PRETTY_FUNCTION__, __LINE__
        );
    }
    MCcalls_ = MCcalls;
}

//==========================================================================================//
//  DoubleConvolution: set method for MCcalls
//------------------------------------------------------------------------------------------//

void DoubleConvolution::SetMCintegral(const bool &MCintegral) {
    if (MCintegral != MCintegral_) {
        MCintegral_ = MCintegral;

        if (MCintegral) {
            convolution_ = std::make_unique<const Convolution>(
                GetCoeffFunc(), GetSplitFunc(), GetAbsErr(), GetRelErr(),
                GetDim()
            );
            conv_coeff_ = nullptr;
        } else {
            conv_coeff_ = std::make_shared<const ConvolutedCoefficientFunction>(
                GetCoeffFunc(), GetSplitFunc(), GetAbsErr(), GetRelErr(),
                GetDim()
            );
            convolution_ = std::make_unique<const Convolution>(
                conv_coeff_, GetSplitFunc(), GetAbsErr(), GetRelErr(), GetDim()
            );
        }
    }
}

//==========================================================================================//
//  DoubleConvolution: integrand of the first regular part
//------------------------------------------------------------------------------------------//

double
    DoubleConvolution::regular1_integrand(double z[], size_t /*dim*/, void *p) {

    function_params *params = (function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    double z1 = z[0], z2 = z[1];

    if (z2 > z1) {
        return 1. / (z1 * z2)
               * (params->conv)->GetSplitFunc()->Regular(x / z1, nf)
               * (params->conv)->GetSplitFunc()->Regular(z1 / z2, nf)
               * (params->conv)
                     ->GetCoeffFunc()
                     ->MuIndependentTerms(z2, m2Q2, nf);
    } else {
        return 0.;
    }
}

//==========================================================================================//
//  DoubleConvolution: integrand of the second regular part
//------------------------------------------------------------------------------------------//

double
    DoubleConvolution::regular2_integrand(double z[], size_t /*dim*/, void *p) {
    function_params *params = (function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    double z1 = z[0], z2 = z[1];

    if (z2 > z1) {
        return 1. / (z1 * z2)
               * (params->conv)->GetSplitFunc()->Regular(x / z1, nf)
               * (params->conv)->GetSplitFunc()->Singular(z1 / z2, nf)
               * ((params->conv)
                      ->GetCoeffFunc()
                      ->MuIndependentTerms(z2, m2Q2, nf)
                  - z1 / z2
                        * (params->conv)
                              ->GetCoeffFunc()
                              ->MuIndependentTerms(z1, m2Q2, nf));
    } else {
        return 0.;
    }
}

//==========================================================================================//
//  DoubleConvolution: integrand of the third regular part
//------------------------------------------------------------------------------------------//

double DoubleConvolution::regular3_integrand(double z, void *p) {

    function_params *params = (function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    double x_max = CoefficientFunction::xMax(m2Q2);

    return -1. / z * (params->conv)->GetSplitFunc()->Regular(x / z, nf)
           * (params->conv)->GetCoeffFunc()->MuIndependentTerms(z, m2Q2, nf)
           * (params->conv)->GetSplitFunc()->SingularIntegrated(z / x_max, nf);
}

//==========================================================================================//
//  DoubleConvolution: regular part
//------------------------------------------------------------------------------------------//

double DoubleConvolution::RegularPart(double x, double m2Q2, int nf) const {

    if (!MCintegral_)
        return convolution_->RegularPart(x, m2Q2, nf);

    double x_max = CoefficientFunction::xMax(m2Q2);
    function_params params = { x, m2Q2, nf, this };

    double xl[2] = { x, x };
    double xu[2] = { x_max, x_max };

    double err, regular1, regular2, regular3, regular4;

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);
    gsl_rng_env_setup();
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

    gsl_monte_function F;

    F.f = &regular1_integrand;
    F.dim = 2;
    F.params = &params;

    gsl_monte_vegas_init(s);
    gsl_monte_vegas_integrate(&F, xl, xu, 2, MCcalls_, r, s, &regular1, &err);

    F.f = &regular2_integrand;
    gsl_monte_vegas_init(s);
    gsl_monte_vegas_integrate(&F, xl, xu, 2, MCcalls_, r, s, &regular2, &err);

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(GetDim());

    gsl_function f;
    f.function = &regular3_integrand;
    f.params = &params;

    gsl_integration_qag(
        &f, x, x_max, GetAbsErr(), GetRelErr(), GetDim(), 4, w, &regular3, &err
    );

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);
    gsl_integration_workspace_free(w);

    regular4 =
        convolution_->RegularPart(x, m2Q2, nf) * GetSplitFunc()->Local(nf);

    return regular1 + regular2 + regular3 + regular4;
}

//==========================================================================================//
//  DoubleConvolution: integrand of the first singular part
//------------------------------------------------------------------------------------------//

double DoubleConvolution::singular1_integrand(
    double z[], size_t /*dim*/, void *p
) {

    function_params *params = (function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    double z1 = z[0], z2 = z[1];

    double tmp;
    if (z2 > x / z1) {
        tmp = (params->conv)->GetSplitFunc()->Regular(x / (z1 * z2), nf) / z1;
    } else {
        tmp = 0.;
    }

    return (params->conv)->GetSplitFunc()->Singular(z1, nf)
           * (tmp - (params->conv)->GetSplitFunc()->Regular(x / z2, nf))
           * (params->conv)->GetCoeffFunc()->MuIndependentTerms(z2, m2Q2, nf)
           / z2;
}

//==========================================================================================//
//  DoubleConvolution: integrand of the second singular part
//------------------------------------------------------------------------------------------//

double DoubleConvolution::singular2_integrand(
    double z[], size_t /*dim*/, void *p
) {

    function_params *params = (function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    double x_max = CoefficientFunction::xMax(m2Q2);

    double z1 = z[0], z2 = z[1];

    double tmp;
    if (z2 > x / (x_max * z1)) {
        tmp = ((params->conv)
                       ->GetCoeffFunc()
                       ->MuIndependentTerms(x / (z1 * z2), m2Q2, nf)
                   / z2
               - (params->conv)
                     ->GetCoeffFunc()
                     ->MuIndependentTerms(x / z1, m2Q2, nf))
              / z1;
    } else {
        tmp = 0.;
    }

    return (params->conv)->GetSplitFunc()->Singular(z1, nf)
           * (params->conv)->GetSplitFunc()->Singular(z2, nf)
           * (tmp
              - ((params->conv)
                         ->GetCoeffFunc()
                         ->MuIndependentTerms(x / z2, m2Q2, nf)
                     / z2
                 - (params->conv)
                       ->GetCoeffFunc()
                       ->MuIndependentTerms(x, m2Q2, nf)));
}

//==========================================================================================//
//  DoubleConvolution: integrand of the third singular part
//------------------------------------------------------------------------------------------//

double DoubleConvolution::singular3_integrand(double z, void *p) {

    function_params *params = (function_params *)p;

    double m2Q2 = (params->m2Q2);
    double x = (params->x);
    int nf = (params->nf);

    double x_max = CoefficientFunction::xMax(m2Q2);

    return -(
        (params->conv)->GetSplitFunc()->Singular(z, nf)
        * ((params->conv)->GetCoeffFunc()->MuIndependentTerms(x / z, m2Q2, nf)
               * (params->conv)
                     ->GetSplitFunc()
                     ->SingularIntegrated(x / (x_max * z), nf)
               / z
           - (params->conv)->GetCoeffFunc()->MuIndependentTerms(x, m2Q2, nf)
                 * (params->conv)
                       ->GetSplitFunc()
                       ->SingularIntegrated(x / x_max, nf))
    );
}

//==========================================================================================//
//  DoubleConvolution: singular part
//------------------------------------------------------------------------------------------//

double DoubleConvolution::SingularPart(double x, double m2Q2, int nf) const {

    if (!MCintegral_)
        return convolution_->SingularPart(x, m2Q2, nf);

    double x_max = CoefficientFunction::xMax(m2Q2);
    function_params params = { x, m2Q2, nf, this };

    double xl[2] = { x / x_max, x };
    double xu[2] = { 1, x_max };

    double err, singular1, singular2, singular3, singular4;

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);
    gsl_rng_env_setup();
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

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

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(GetDim());

    gsl_function f;
    f.function = &singular3_integrand;
    f.params = &params;

    gsl_integration_qag(
        &f, x / x_max, 1, GetAbsErr(), GetRelErr(), GetDim(), 4, w, &singular3,
        &err
    );

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);
    gsl_integration_workspace_free(w);

    singular4 =
        convolution_->SingularPart(x, m2Q2, nf) * GetSplitFunc()->Local(nf);

    return singular1 + singular2 + singular3 + singular4;
}

//==========================================================================================//
//  DoubleConvolution: local part
//------------------------------------------------------------------------------------------//

double DoubleConvolution::LocalPart(double x, double m2Q2, int nf) const {

    if (!MCintegral_)
        return convolution_->LocalPart(x, m2Q2, nf);

    double x_max = CoefficientFunction::xMax(m2Q2);

    return convolution_->Convolute(x, m2Q2, nf)
           * (GetSplitFunc()->Local(nf)
              - GetSplitFunc()->SingularIntegrated(x / x_max, nf));
}
