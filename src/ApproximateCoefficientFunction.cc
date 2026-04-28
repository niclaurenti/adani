#include "adani/ApproximateCoefficientFunction.h"

#include "adani/Constants.h"
#include "adani/HighEnergyCoefficientFunction.h"
#include "adani/HighScaleCoefficientFunction.h"
#include "adani/SpecialFunctions.h"

#include <array>
#include <cmath>

using std::array;

//==========================================================================================//
//  AbstractApproximate: constructor
//------------------------------------------------------------------------------------------//

AbstractApproximate::AbstractApproximate(
    const int &order, const char &kind, const char &channel,
    const double &abserr, const double &relerr, const int &dim
)
    : CoefficientFunction(order, kind, channel) {

    muterms_ = std::make_unique<ExactCoefficientFunction>(
        order, kind, channel, abserr, relerr, dim
    );
}

//==========================================================================================//
//  AbstractApproximate: get method for abserr
//------------------------------------------------------------------------------------------//

double AbstractApproximate::GetAbsErr() const { return muterms_->GetAbsErr(); }

//==========================================================================================//
//  AbstractApproximate: get method for relerr
//------------------------------------------------------------------------------------------//

double AbstractApproximate::GetRelErr() const { return muterms_->GetRelErr(); }

//==========================================================================================//
//  ExactCoefficientFunction: get method for dim
//------------------------------------------------------------------------------------------//

int AbstractApproximate::GetDim() const { return muterms_->GetDim(); }

//==========================================================================================//
//  AbstractApproximate: function that sets the double integral method
//------------------------------------------------------------------------------------------//

void AbstractApproximate::SetDoubleIntegralMethod(
    const DoubleIntegralMethod &double_int_method, const double &abserr,
    const double &relerr, const int &dim, const int &MCcalls
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
    if (x >= xMax(m2Q2) || x <= 0)
        return Value(0.);

    return MuIndependentTermsBand(x, m2Q2, nf)
           + MuDependentTerms(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  ApproximateCoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

ApproximateCoefficientFunction::ApproximateCoefficientFunction(
    const int &order, const char &kind, const char &channel, const int &damp_power, const bool &NLL,
    const HighScaleVersion &highscale_version, const double &abserr,
    const double &relerr, const int &dim
)
    : AbstractApproximate(order, kind, channel, abserr, relerr, dim) {

    threshold_ =
        std::make_unique<ThresholdCoefficientFunction>(order, kind, channel);
    asymptotic_ = std::make_unique<AsymptoticCoefficientFunction>(
        order, kind, channel, NLL, highscale_version
    );

    approximation_ = nullptr;
    variation_ = nullptr;
}

//==========================================================================================//
//  ApproximateCoefficientFunction: copy constructor
//------------------------------------------------------------------------------------------//

ApproximateCoefficientFunction::ApproximateCoefficientFunction(
    const ApproximateCoefficientFunction &obj
)
    : ApproximateCoefficientFunction(
          obj.GetOrder(), obj.GetKind(), obj.GetChannel(), obj.GetDampPower(), obj.GetNLL(),
          obj.GetHighScaleVersion(), obj.GetAbsErr(), obj.GetRelErr(),
          obj.GetDim()
      ) {}

//==========================================================================================//
//  ApproximateCoefficientFunction: band of the approximate mu independent terms
//------------------------------------------------------------------------------------------//

Value ApproximateCoefficientFunction::MuIndependentTermsBand(
    double x, double m2Q2, int nf
) const {
    return Approximation(x, m2Q2, nf);
}

//==========================================================================================//
//  ApproximateCoefficientFunction: band of the approximate mu independent terms
//  (new)
//------------------------------------------------------------------------------------------//

Value ApproximateCoefficientFunction::Approximation(
    double x, double m2Q2, int nf
) const {

    if (x <= 0 || x > xMax(m2Q2))
        return Value(0.);

    double rho = 1.3, eta0 = 2.;

    Value thresh = threshold_->MuIndependentTermsBand(x, m2Q2, nf);
    Value asy = asymptotic_->MuIndependentTermsBand(x, m2Q2, nf);

    double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

    double damp_thr = 1. / (1. + pow(eta / eta0, rho));
    double damp_asy = 1. - damp_thr;

    return (asy * damp_asy + thresh * damp_thr) / std::pow(1. + m2Q2, damp_power_);
}

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: parameters of the approximation
//------------------------------------------------------------------------------------------//

klmv_params klmv_C2g2A = { 1., 42.5, 0, 0, 0 };
klmv_params klmv_C2g2B = { 0.8, 19.4, 0, 0, 0 };
klmv_params klmv_C2q2A = { 1., 42.5, 0, 0, 0 };
klmv_params klmv_C2q2B = { 0.8, 19.4, 0, 0, 0 };

klmv_params klmv_C2g3A = { 1., 20., 0.007, 4, 0.28 };
klmv_params klmv_C2g3B = { 0.8, 10.7, 0.055, 2, 0.423 };
klmv_params klmv_C2q3A = { 1., 20., 0.004, 4, 0.125 };
klmv_params klmv_C2q3B = { 0.8, 10.7, 0.0245, 2, 0.17 };
klmv_params klmv_C2g3B_lowxi = { 0.8, 10.7, 0.055125, 2, 0.3825 };

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: constructor
//------------------------------------------------------------------------------------------//

ApproximateCoefficientFunctionKLMV::ApproximateCoefficientFunctionKLMV(
    const int &order, const char &kind, const char &channel,
    const HighScaleVersion &highscale_version, const bool &lowxi,
    const double &abserr, const double &relerr, const int &dim
)
    : AbstractApproximate(order, kind, channel, abserr, relerr, dim),
      lowxi_(lowxi) {
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

    threshold_ =
        std::make_unique<ThresholdCoefficientFunction>(order, kind, channel);
    threshold_->SetPlainThreshold(true);

    highscale_ = std::make_unique<HighScaleCoefficientFunction>(
        order, kind, channel, highscale_version
    );

    highenergy_ = std::make_unique<HighEnergyCoefficientFunction>(
        order, kind, channel, false
    );
}

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: copy constructor
//------------------------------------------------------------------------------------------//

ApproximateCoefficientFunctionKLMV::ApproximateCoefficientFunctionKLMV(
    const ApproximateCoefficientFunctionKLMV &obj
)
    : ApproximateCoefficientFunctionKLMV(
          obj.GetOrder(), obj.GetKind(), obj.GetChannel(),
          obj.GetHighScaleVersion(), obj.GetAbsErr(), obj.GetRelErr(),
          obj.GetDim()
      ) {}

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

        if (lowxi_)
            params_B_ = std::make_unique<klmv_params>(klmv_C2g3B_lowxi);
        else
            params_B_ = std::make_unique<klmv_params>(klmv_C2g3B);

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
        fx_ = &ApproximateCoefficientFunctionKLMV::Order2;
        switch (GetChannel()) {
        case 'g':
            params_A_ = std::make_unique<klmv_params>(klmv_C2g2A);
            params_B_ = std::make_unique<klmv_params>(klmv_C2g2B);
            break;
        case 'q':
            params_A_ = std::make_unique<klmv_params>(klmv_C2q2A);
            params_B_ = std::make_unique<klmv_params>(klmv_C2q2B);
            break;
        default:
            throw UnexpectedException(
                "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
            );
        }
        break;
    case 3:
        fx_ = &ApproximateCoefficientFunctionKLMV::Order3;
        switch (GetChannel()) {
        case 'g':
            params_A_ = std::make_unique<klmv_params>(klmv_C2g3A);
            if (lowxi_)
                params_B_ = std::make_unique<klmv_params>(klmv_C2g3B_lowxi);
            else
                params_B_ = std::make_unique<klmv_params>(klmv_C2g3B);
            break;
        case 'q':
            params_A_ = std::make_unique<klmv_params>(klmv_C2q3A);
            params_B_ = std::make_unique<klmv_params>(klmv_C2q3B);
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

    return (this->*fx_)(x, m2Q2, nf);
}

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: functional form of the approximation A
//  at O(2)
//
//  Eq. (4.10,4.11) of Ref. [arXiv:1205.5727].
//------------------------------------------------------------------------------------------//

Value ApproximateCoefficientFunctionKLMV::Order2(
    double x, double m2Q2, int nf
) const {

    double xmax = 1. / (1. + 4 * m2Q2);
    if (x <= 0 || x >= xmax)
        return Value(0.);

    double beta = sqrt(1. - 4 * m2Q2 * x / (1. - x));
    double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

    double thr = threshold_->MuIndependentTerms(x, m2Q2, nf);
    double thr_const = threshold_->BetaIndependentTerms(x, m2Q2, 1.);

    double he_ll = highenergy_->MuIndependentTerms(x, m2Q2, nf);

    double hs = highscale_->MuIndependentTerms(x, m2Q2, nf);

    double xi = 1. / m2Q2;
    double f = 1. / (1. + exp(2. * (xi - 4.)));

    double eta_gamma = pow(eta, params_A_->eta_exponent);
    double eta_delta = pow(eta, params_B_->eta_exponent);

    double beta3 = beta * beta * beta;

    double res_A =
        thr - thr_const + (1. - f) * beta * hs
        + f * beta3 * he_ll * eta_gamma / (params_A_->shift + eta_gamma);
    double res_B =
        thr + (1. - f) * beta3 * hs
        + f * beta3 * he_ll * eta_delta / (params_B_->shift + eta_delta);

    if (res_A > res_B)
        return Value(res_A, res_B);
    else
        return Value(res_B, res_A);
}

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: functional form of the approximation A
//  at O(3)
//
//  Eq. (4.17,4.18,4.21,4.22) of Ref. [arXiv:1205.5727].
//------------------------------------------------------------------------------------------//

Value ApproximateCoefficientFunctionKLMV::Order3(
    double x, double m2Q2, int nf
) const {
    double xmax = 1. / (1. + 4 * m2Q2);
    if (x <= 0 || x >= xmax)
        return Value(0.);

    double beta = sqrt(1. - 4 * m2Q2 * x / (1. - x));
    double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

    double thr = threshold_->MuIndependentTerms(x, m2Q2, nf);
    double thr_const = threshold_->BetaIndependentTerms(x, m2Q2, 1.);

    double he_ll = highenergy_->MuIndependentTerms(x, m2Q2, nf);
    Value he_nll = ApproximateNLL(x, m2Q2);

    array<double, 3> hs = highscale_->fxBand_NotOrdered(x, m2Q2, 1., nf);

    double xi = 1. / m2Q2;
    double f = 1. / (1. + exp(2. * (xi - 4.)));

    double eta_gamma = pow(eta, params_A_->eta_exponent);
    double eta_delta = pow(eta, params_B_->eta_exponent);

    double beta3 = beta * beta * beta;

    double res_A = thr - thr_const + (1. - f) * beta * hs[1]
                   + f * beta3
                         * (-log(eta) / log(x) * he_ll
                            + he_nll.GetHigher() * eta_gamma
                                  / (params_A_->shift + eta_gamma));
    double res_B = thr - thr_const + 2 * f * thr_const
                   + (1. - f) * beta3 * hs[2]
                   + f * beta3
                         * (-log(eta) / log(x) * he_ll
                            + he_nll.GetLower() * eta_delta
                                  / (params_B_->shift + eta_delta));

    if (res_A > res_B)
        return Value(res_A, res_B);
    else
        return Value(res_B, res_A);
}

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: NLL coefficient of the small-x
//  approximation
//
//  Eq. (4.19,4.20,4.23,4.24) of Ref. [arXiv:1205.5727].
//------------------------------------------------------------------------------------------//

Value ApproximateCoefficientFunctionKLMV::ApproximateNLL(
    double x, double m2Q2
) const {
    double pi3 = M_PI * M_PI * M_PI;
    double tmp_A = (64. * pi3)
                   * (params_A_->log_coeff
                          * pow(log(1. / m2Q2) / log(5), params_A_->log_pow)
                      - params_A_->const_coeff)
                   * 4. / m2Q2 / x;
    double tmp_B = (64. * pi3)
                   * (params_B_->log_coeff
                          * pow(log(1. / m2Q2) / log(5), params_B_->log_pow)
                      - params_B_->const_coeff)
                   * 4. / m2Q2 / x;
    return Value(tmp_A, tmp_B);
}
