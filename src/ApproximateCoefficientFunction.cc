#include "adani/ApproximateCoefficientFunction.h"

#include "adani/Constants.h"
#include "adani/HighEnergyCoefficientFunction.h"
#include "adani/HighScaleCoefficientFunction.h"
#include "adani/SpecialFunctions.h"

#include <cmath>
#include <vector>

using std::vector;

//==========================================================================================//
//  AbstractApproximate: constructor
//------------------------------------------------------------------------------------------//

AbstractApproximate::AbstractApproximate(
    const int &order, const char &kind, const char &channel,
    const double &abserr, const double &relerr, const int &dim
)
    : CoefficientFunction(order, kind, channel) {

    muterms_ =
        new ExactCoefficientFunction(order, kind, channel, abserr, relerr, dim);
}

//==========================================================================================//
//  AbstractApproximate: destructor
//------------------------------------------------------------------------------------------//

AbstractApproximate::~AbstractApproximate() { delete muterms_; }

//==========================================================================================//
//  AbstractApproximate: function that sets the double integral method
//------------------------------------------------------------------------------------------//

void AbstractApproximate::SetDoubleIntegralMethod(
    const string &double_int_method, const double &abserr, const double &relerr,
    const int &dim, const int &MCcalls
) {
    muterms_->SetDoubleIntegralMethod(
        double_int_method, abserr, relerr, dim, MCcalls
    );
}

//==========================================================================================//
//  AbstractApproximate: central value of the mu-independent terms
//------------------------------------------------------------------------------------------//

double AbstractApproximate::MuIndependentTerms(
    double x, double m2Q2, int nf
) const {
    return MuIndependentTermsBand(x, m2Q2, nf).GetCentral();
}

//==========================================================================================//
//  AbstractApproximate: mu dependent terms
//------------------------------------------------------------------------------------------//

double AbstractApproximate::MuDependentTerms(
    double x, double m2Q2, double m2mu2, int nf
) const {
    return muterms_->MuDependentTerms(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  AbstractApproximate: band of the total result
//------------------------------------------------------------------------------------------//

Value AbstractApproximate::fxBand(
    double x, double m2Q2, double m2mu2, int nf
) const {
    double x_max = 1. / (1. + 4 * m2Q2);
    if (x >= x_max || x <= 0)
        return Value(0.);

    return MuIndependentTermsBand(x, m2Q2, nf)
           + MuDependentTerms(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  ApproximateCoefficientFunction: parameters of the legacy approximation
//------------------------------------------------------------------------------------------//

#define a 2.5
#define b 5.

struct approximation_parameters C2_g1_params = { 0.2, 2.5, 2.5, 1.2 };
struct approximation_parameters CL_g1_params = { 20., 11., 3., 2. };

struct approximation_parameters C2_g2_params = { 1.7, 2.5, 2.5, 1.2 };
struct approximation_parameters CL_g2_params = { 20., 11., 3., 2. };
struct approximation_parameters C2_ps2_params = { 1.7, 2.5, 2.5, 1.2 };
struct approximation_parameters CL_ps2_params = { 20., 11., 3., 2. };

struct approximation_parameters C2_g3_params = { 0.3, 2.5, 2.5, 1.2 };
struct approximation_parameters CL_g3_params = { 10., 11., 3., 2. };
struct approximation_parameters C2_ps3_params = { 0.3, 2.5, 2.5, 1.2 };
struct approximation_parameters CL_ps3_params = { 20., 11., 3., 2. };

struct variation_parameters C2_var = { 3., 0.3 };
struct variation_parameters CL_var = { 2., 0.2 };

//==========================================================================================//
//  ApproximateCoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

ApproximateCoefficientFunction::ApproximateCoefficientFunction(
    const int &order, const char &kind, const char &channel, const bool &NLL,
    const string &highscale_version, const double &abserr, const double &relerr,
    const int &dim
)
    : AbstractApproximate(order, kind, channel, abserr, relerr, dim) {

    threshold_ = new ThresholdCoefficientFunction(order, kind, channel);
    asymptotic_ = new AsymptoticCoefficientFunction(
        order, kind, channel, NLL, highscale_version
    );

    fx_ = &ApproximateCoefficientFunction::Approximation;

    legacy_appr_ = false;
    approximation_ = nullptr;
    variation_ = nullptr;
}

//==========================================================================================//
//  ApproximateCoefficientFunction: destructor
//------------------------------------------------------------------------------------------//

ApproximateCoefficientFunction::~ApproximateCoefficientFunction() {
    delete threshold_;
    delete asymptotic_;
    delete approximation_;
    delete variation_;
}

//==========================================================================================//
//  ApproximateCoefficientFunction: set parameters of legacy approximation
//------------------------------------------------------------------------------------------//

void ApproximateCoefficientFunction::SetLegacyParameters() {
    try {
        if (GetOrder() == 1) {
            if (GetKind() == '2') {
                if (GetChannel() == 'g')
                    approximation_ = new approximation_parameters(C2_g1_params);
                else {
                    throw UnexpectedException(
                        "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                    );
                }
                variation_ = new variation_parameters(C2_var);
            } else if (GetKind() == 'L') {
                if (GetChannel() == 'g')
                    approximation_ = new approximation_parameters(CL_g2_params);
                else {
                    throw UnexpectedException(
                        "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                    );
                }
                variation_ = new variation_parameters(CL_var);
            } else {
                throw UnexpectedException(
                    "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                );
            }
        } else if (GetOrder() == 2) {
            if (GetKind() == '2') {
                if (GetChannel() == 'g')
                    approximation_ = new approximation_parameters(C2_g2_params);
                else if (GetChannel() == 'q')
                    approximation_ = new approximation_parameters(C2_ps2_params);
                else {
                    throw UnexpectedException(
                        "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                    );
                }
                variation_ = new variation_parameters(C2_var);
            } else if (GetKind() == 'L') {
                if (GetChannel() == 'g')
                    approximation_ = new approximation_parameters(CL_g2_params);
                else if (GetChannel() == 'q')
                    approximation_ = new approximation_parameters(CL_ps2_params);
                else {
                    throw UnexpectedException(
                        "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                    );
                }
                variation_ = new variation_parameters(CL_var);
            }
        } else if (GetOrder() == 3) {
            if (GetKind() == '2') {
                if (GetChannel() == 'g')
                    approximation_ = new approximation_parameters(C2_g3_params);
                else if (GetChannel() == 'q')
                    approximation_ = new approximation_parameters(C2_ps3_params);
                else {
                    throw UnexpectedException(
                        "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                    );
                }
                variation_ = new variation_parameters(C2_var);
            } else if (GetKind() == 'L') {
                if (GetChannel() == 'g')
                    approximation_ = new approximation_parameters(CL_g3_params);
                else if (GetChannel() == 'q')
                    approximation_ = new approximation_parameters(CL_ps3_params);
                else {
                    throw UnexpectedException(
                        "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                    );
                }
                variation_ = new variation_parameters(CL_var);
            }
        } else {
            throw UnexpectedException(
                "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
            );
        }
    } catch (UnexpectedException &e) {
        e.runtime_error();
    }
}

//==========================================================================================//
//  ApproximateCoefficientFunction: set legacy threshold behavior
//------------------------------------------------------------------------------------------//

void ApproximateCoefficientFunction::SetLegacyThreshold(
    const bool &legacy_threshold
) {
    threshold_->SetLegacyThreshold(legacy_threshold);
}

//==========================================================================================//
//  ApproximateCoefficientFunction: restore legacy behavior for power terms
//------------------------------------------------------------------------------------------//

void ApproximateCoefficientFunction::SetLegacyPowerTerms(const bool &legacy_pt
) {
    asymptotic_->SetLegacyPowerTerms(legacy_pt);
}

//==========================================================================================//
//  ApproximateCoefficientFunction: restore legacy approximation
//------------------------------------------------------------------------------------------//

void ApproximateCoefficientFunction::SetLegacyApproximation(const bool &legacy_appr) {
    try {
        if (legacy_appr == legacy_appr_) {
            throw NotValidException(
                "Setting legacy approximation identical to its previous value!",
                __PRETTY_FUNCTION__, __LINE__
            );
        }

        if (legacy_appr) {
            SetLegacyThreshold(true);
            SetLegacyPowerTerms(true);
            SetLegacyParameters();
            fx_ = &ApproximateCoefficientFunction::ApproximationLegacy;

        } else {
            SetLegacyThreshold(false);
            SetLegacyPowerTerms(false);
            fx_ = &ApproximateCoefficientFunction::Approximation;
            delete approximation_;
            delete variation_;
        }
    } catch (NotValidException &e) {
        e.warning();
    }
}

//==========================================================================================//
//  ApproximateCoefficientFunction: band of the approximate mu independent terms
//------------------------------------------------------------------------------------------//

Value ApproximateCoefficientFunction::MuIndependentTermsBand(
    double x, double m2Q2, int nf
) const {
    return (this->*fx_)(x, m2Q2, nf);
}

//==========================================================================================//
//  ApproximateCoefficientFunction: band of the approximate mu independent terms (new)
//------------------------------------------------------------------------------------------//

Value ApproximateCoefficientFunction::Approximation(
    double x, double m2Q2, int nf
) const {

    double x_max = 1. / (1 + 4 * m2Q2);
    if (x <= 0 || x > x_max)
        return Value(0.);

    double rho = 1.3, eta0 = 2.;

    Value thresh = threshold_->MuIndependentTermsBand(x, m2Q2, nf);
    Value asy = asymptotic_->MuIndependentTermsBand(x, m2Q2, nf);

    double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

    double damp_thr = 1. / (1. + pow(eta / eta0, rho));
    double damp_asy = 1. - damp_thr;

    return asy * damp_asy + thresh * damp_thr;
}

//==========================================================================================//
//  ApproximateCoefficientFunction: band of the approximate mu independent terms (legacy form)
//------------------------------------------------------------------------------------------//

Value ApproximateCoefficientFunction::ApproximationLegacy(
    double x, double m2Q2, int nf
) const {

    double x_max = 1. / (1 + 4 * m2Q2);
    if (x <= 0 || x > x_max)
        return Value(0.);

    double A = approximation_->A, B = approximation_->B, C = approximation_->C,
           D = approximation_->D;
    double var2 = variation_->var2, var1 = variation_->var1;

    double Amax = var1 * A, Amin = A / var1, Bmax = B * var1, Bmin = B / var1;
    double Cmax, Cmin, Dmax, Dmin;

    Cmax = (1. + var2) * C;
    Cmin = (1. - var2) * C;
    Dmax = (1. + var2) * D;
    Dmin = (1. - var2) * D;

    double Avec[3] = { A, Amax, Amin };
    double Bvec[3] = { B, Bmax, Bmin };
    double Cvec[3] = { C, Cmax, Cmin };
    double Dvec[3] = { D, Dmax, Dmin };

    vector<double> asy =
        (asymptotic_->MuIndependentTermsBand(x, m2Q2, nf)).ToVect();
    vector<double> thresh =
        (threshold_->MuIndependentTermsBand(x, m2Q2, nf)).ToVect();

    double central = ApproximationLegacyForm(x, m2Q2, asy[0], thresh[0], A, B, C, D);
    double higher = central, lower = central, tmp;

    for (int i = 0; i < int(asy.size()); i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    for (int m = 0; m < 3; m++) {
                        for (int n = 0; n < 3; n++) {
                            tmp = ApproximationLegacyForm(
                                x, m2Q2, asy[i], thresh[j], Avec[k], Bvec[l],
                                Cvec[m], Dvec[n]
                            );
                            if (tmp > higher)
                                higher = tmp;
                            if (tmp < lower)
                                lower = tmp;
                        }
                    }
                }
            }
        }
    }

    return Value(central, higher, lower);
}

//==========================================================================================//
//  ApproximateCoefficientFunction: functional form of the legacy approximation
//------------------------------------------------------------------------------------------//

double ApproximateCoefficientFunction::ApproximationLegacyForm(
    double x, double m2Q2, double asy, double thresh, double A, double B,
    double C, double D
) const {

    double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;
    double xi = 1. / m2Q2;

    double h = A + (B - A) / (1. + exp(a * (log(xi) - b)));
    double k = C + (D - C) / (1. + exp(a * (log(xi) - b)));

    double damp_thr = 1. / (1. + pow(eta / h, k));
    double damp_asy = 1. - damp_thr;

    return asy * damp_asy + thresh * damp_thr;
}

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: parameters of the approximation
//------------------------------------------------------------------------------------------//

struct klmv_params klmv_C2g2A = { 1., 42.5, 0, 0, 0 };
struct klmv_params klmv_C2g2B = { 0.8, 19.4, 0, 0, 0 };
struct klmv_params klmv_C2q2A = { 1., 42.5, 0, 0, 0 };
struct klmv_params klmv_C2q2B = { 0.8, 19.4, 0, 0, 0 };

struct klmv_params klmv_C2g3A = { 1., 20., 0.007, 4, 0.28 };
struct klmv_params klmv_C2g3B = { 0.8, 10.7, 0.055, 2, 0.423 };
struct klmv_params klmv_C2q3A = { 1., 20., 0.004, 4, 0.125 };
struct klmv_params klmv_C2q3B = { 0.8, 10.7, 0.0245, 2, 0.17 };
struct klmv_params klmv_C2g3B_lowxi = { 0.8, 10.7, 0.055125, 2, 0.3825 };

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: constructor
//------------------------------------------------------------------------------------------//

ApproximateCoefficientFunctionKLMV::ApproximateCoefficientFunctionKLMV(
    const int &order, const char &kind, const char &channel,
    const string &highscale_version, const bool &lowxi, const double &abserr,
    const double &relerr, const int &dim
)
    : AbstractApproximate(order, kind, channel, abserr, relerr, dim), lowxi_(lowxi) {
    try {
        if (GetOrder() == 1) {
            throw NotImplementedException(
                "KLMV approximation is not implemented at order 1! Got order="
                    + to_string(order),
                __PRETTY_FUNCTION__, __LINE__
            );
        }
        if (GetKind() == 'L') {
            throw NotImplementedException(
                "KLMV approximation is not implemented for kind = 'L'! Got "
                "kind='"
                    + string(1, kind) + "'",
                __PRETTY_FUNCTION__, __LINE__
            );
        }

        SetFunctions();
    } catch (const NotImplementedException &e) {
        e.runtime_error();
    } catch (UnexpectedException &e) {
        e.runtime_error();
    }

    threshold_ = new ThresholdCoefficientFunction(order, kind, channel);
    threshold_->SetLegacyThreshold(true);

    highscale_ = new HighScaleCoefficientFunction(
        order, kind, channel, highscale_version
    );

    highenergy_ =
        new HighEnergyCoefficientFunction(order, kind, channel, false);
}

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: destructor
//------------------------------------------------------------------------------------------//

ApproximateCoefficientFunctionKLMV::~ApproximateCoefficientFunctionKLMV() {
    delete threshold_;
    delete highscale_;
    delete highenergy_;
    delete params_A_;
    delete params_B_;
}

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: set lowxi
//------------------------------------------------------------------------------------------//

void ApproximateCoefficientFunctionKLMV::SetLowXi(const bool &lowxi) {
    try {
        if (lowxi == lowxi_) {
            throw NotValidException(
                "Setting lowxi approximation identical to its previous value!",
                __PRETTY_FUNCTION__, __LINE__
            );
        }

        delete params_B_;
        if (lowxi_)
            params_B_ = new klmv_params(klmv_C2g3B_lowxi);
        else
            params_B_ = new klmv_params(klmv_C2g3B);

    } catch (NotValidException &e) {
        e.warning();
    }
}

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: set function to be used
//------------------------------------------------------------------------------------------//

void ApproximateCoefficientFunctionKLMV::SetFunctions() {
    switch (GetOrder()) {
        case 2:
            fx_A_=&ApproximateCoefficientFunctionKLMV::Approximation2A;
            fx_B_=&ApproximateCoefficientFunctionKLMV::Approximation2B;
            switch (GetChannel()) {
                case 'g':
                    params_A_ = new klmv_params(klmv_C2g2A);
                    params_B_ = new klmv_params(klmv_C2g2B);
                    break;
                case 'q':
                    params_A_ = new klmv_params(klmv_C2q2A);
                    params_B_ = new klmv_params(klmv_C2q2B);
                    break;
                default:
                    throw UnexpectedException(
                        "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                    );
            }
            break;
        case 3:
            fx_A_=&ApproximateCoefficientFunctionKLMV::Approximation3A;
            fx_B_=&ApproximateCoefficientFunctionKLMV::Approximation3B;
            switch (GetChannel()) {
                case 'g':
                    params_A_ = new klmv_params(klmv_C2g3A);
                    if (lowxi_)
                        params_B_ = new klmv_params(klmv_C2g3B_lowxi);
                    else
                        params_B_ = new klmv_params(klmv_C2g3B);
                    break;
                case 'q':
                    params_A_ = new klmv_params(klmv_C2q3A);
                    params_B_ = new klmv_params(klmv_C2q3B);
                    break;
                default:
                    throw UnexpectedException(
                        "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                    );
            }
            break;
        default:
            throw UnexpectedException(
                "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
            );
    }
}

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: band of the mu independent terms
//------------------------------------------------------------------------------------------//

Value ApproximateCoefficientFunctionKLMV::MuIndependentTermsBand(
    double x, double m2Q2, int nf
) const {

    double xmax = 1. / (1. + 4 * m2Q2);
    if (x <= 0 || x >= xmax)
        return Value(0.);

    double res_A = (this->*fx_A_)(x, m2Q2, nf);
    double res_B = (this->*fx_B_)(x, m2Q2, nf);

    if (res_A > res_B)
        return Value(res_A, res_B);
    else
        return Value(res_B, res_A);
}

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: functional form of the approximation A at O(2)
//------------------------------------------------------------------------------------------//

double ApproximateCoefficientFunctionKLMV::Approximation2A(
    double x, double m2Q2, int nf
) const {
    double beta = sqrt(1. - 4 * m2Q2 * x / (1. - x));
    double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

    double thr = threshold_->MuIndependentTerms(x, m2Q2, nf);
    double thr_const = threshold_->BetaIndependentTerms(x, m2Q2, 1.);

    double he_ll = highenergy_->MuIndependentTerms(x, m2Q2, nf);

    double hs = highscale_->MuIndependentTerms(x, m2Q2, nf);

    double xi = 1. / m2Q2;
    double f = 1. / (1. + exp(2. * (xi - 4.)));

    double eta_gamma = pow(eta, params_A_->gamma);

    double beta3 = beta * beta * beta;

    return thr - thr_const + (1. - f) * beta * hs
           + f * beta3
                 * he_ll * eta_gamma / (params_A_->C + eta_gamma);
}

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: functional form of the approximation B at O(2)
//------------------------------------------------------------------------------------------//

double ApproximateCoefficientFunctionKLMV::Approximation2B(
    double x, double m2Q2, int nf
) const {
    double beta = sqrt(1. - 4 * m2Q2 * x / (1. - x));
    double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

    double thr = threshold_->MuIndependentTerms(x, m2Q2, nf);

    double he_ll = highenergy_->MuIndependentTerms(x, m2Q2, nf);

    double hs = highscale_->MuIndependentTerms(x, m2Q2, nf);

    double xi = 1. / m2Q2;
    double f = 1. / (1. + exp(2. * (xi - 4.)));

    double eta_delta = pow(eta, params_B_->gamma);

    double beta3 = beta * beta * beta;

    return thr + (1. - f) * beta3 * hs
           + f * beta3
                 * he_ll * eta_delta / (params_B_->C + eta_delta);
}

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: functional form of the approximation A at O(3)
//------------------------------------------------------------------------------------------//

double ApproximateCoefficientFunctionKLMV::Approximation3A(
    double x, double m2Q2, int nf
) const {
    double beta = sqrt(1. - 4 * m2Q2 * x / (1. - x));
    double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

    double thr = threshold_->MuIndependentTerms(x, m2Q2, nf);
    double thr_const = threshold_->BetaIndependentTerms(x, m2Q2, 1.);

    double he_ll = highenergy_->MuIndependentTerms(x, m2Q2, nf);
    double he_nll = ApproximateNLL(x, m2Q2).GetHigher();

    double hs = highscale_->fxBand_NotOrdered(x, m2Q2, 1., nf)[1];

    double xi = 1. / m2Q2;
    double f = 1. / (1. + exp(2. * (xi - 4.)));

    double eta_gamma = pow(eta, params_A_->gamma);

    double beta3 = beta * beta * beta;

    return thr - thr_const + (1. - f) * beta * hs
           + f * beta3
                 * (-log(eta) / log(x) * he_ll
                    + he_nll * eta_gamma / (params_A_->C + eta_gamma));
}

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: functional form of the approximation B at O(3)
//------------------------------------------------------------------------------------------//

double ApproximateCoefficientFunctionKLMV::Approximation3B(
    double x, double m2Q2, int nf
) const {
    double beta = sqrt(1. - 4 * m2Q2 * x / (1. - x));
    double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

    double thr = threshold_->MuIndependentTerms(x, m2Q2, nf);
    double thr_const = threshold_->BetaIndependentTerms(x, m2Q2, 1.);

    double he_ll = highenergy_->MuIndependentTerms(x, m2Q2, nf);
    double he_nll = ApproximateNLL(x, m2Q2).GetLower();

    double hs = highscale_->fxBand_NotOrdered(x, m2Q2, 1., nf)[2];

    double xi = 1. / m2Q2;
    double f = 1. / (1. + exp(2. * (xi - 4.)));

    double eta_delta = pow(eta, params_B_->gamma);

    double beta3 = beta * beta * beta;

    return thr - thr_const + 2 * f * thr_const + (1. - f) * beta3 * hs
           + f * beta3
                 * (-log(eta) / log(x) * he_ll
                    + he_nll * eta_delta / (params_B_->C + eta_delta));
}

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: NLL coefficient of the small-x
//  approximation
//------------------------------------------------------------------------------------------//

Value ApproximateCoefficientFunctionKLMV::ApproximateNLL(
    double x, double m2Q2
) const {
    double pi3 = M_PI * M_PI * M_PI;
    double tmp_A =
        (64. * pi3)
        * (params_A_->log_coeff * pow(log(1. / m2Q2) / log(5), params_A_->log_pow)
           - params_A_->const_coeff)
        * 4. / m2Q2 / x;
    double tmp_B =
        (64. * pi3)
        * (params_B_->log_coeff * pow(log(1. / m2Q2) / log(5), params_B_->log_pow)
           - params_B_->const_coeff)
        * 4. / m2Q2 / x;
    return Value(tmp_A, tmp_B);
}
