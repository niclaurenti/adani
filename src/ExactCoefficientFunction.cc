#include "adani/ExactCoefficientFunction.h"
#include "adani/Constants.h"
#include "adani/SpecialFunctions.h"
#include "adani/ThresholdCoefficientFunction.h"

#include <cmath>
#include <future>

//==========================================================================================//
//  ExactCoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

ExactCoefficientFunction::ExactCoefficientFunction(
    const int &order, const char &kind, const char &channel,
    const double &abserr, const double &relerr, const int &dim
)
    : CoefficientFunction(order, kind, channel) {

    gluon_as1_ = nullptr;
    gluon_as2_ = nullptr;
    quark_as2_ = nullptr;

    Pgq0_ = nullptr;
    Pgg0_ = nullptr;
    Pgq1_ = nullptr;
    Pqq0_ = nullptr;
    Pgg0Pgq0_ = nullptr;
    Pqq0Pgq0_ = nullptr;
    Pgg1_ = nullptr;
    Pqg0_ = nullptr;
    Pgq0Pqg0_ = nullptr;
    Pgg0Pgg0_ = nullptr;

    delta_ = nullptr;

    asy_ = nullptr;
    thr_ = nullptr;

    double_int_meth_ = DoubleIntegralMethod::Analytical;

    try {

        if (GetOrder() > 1) {
            // needed in both channels
            gluon_as1_ = std::make_shared<const ExactCoefficientFunction>(
                1, GetKind(), 'g'
            );
            delta_ = std::make_shared<const Delta>();

            switch (GetChannel()) {
            case 'q':
                Pgq0_ = std::make_shared<const SplittingFunction>(0, 'g', 'q');
                break;
            case 'g':
                Pgg0_ = std::make_shared<const SplittingFunction>(0, 'g', 'g');
                break;
            default:
                throw UnexpectedException(
                    "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                );
            }
        }
        if (GetOrder() > 2) {
            // needed in both channels
            gluon_as2_ = std::make_shared<const ExactCoefficientFunction>(
                2, GetKind(), 'g'
            );
            quark_as2_ = std::make_shared<const ExactCoefficientFunction>(
                2, GetKind(), 'q'
            );

            switch (GetChannel()) {
            case 'q':
                Pgq1_ = std::make_shared<const SplittingFunction>(1, 'g', 'q');
                Pqq0_ = std::make_shared<const SplittingFunction>(0, 'q', 'q');
                Pgg0Pgq0_ =
                    std::make_shared<const ConvolutedSplittingFunctions>(
                        0, 'g', 'g', 0, 'g', 'q'
                    );
                Pqq0Pgq0_ =
                    std::make_shared<const ConvolutedSplittingFunctions>(
                        0, 'q', 'q', 0, 'g', 'q'
                    );
                break;
            case 'g':
                Pgg1_ = std::make_shared<const SplittingFunction>(1, 'g', 'g');
                Pqg0_ = std::make_shared<const SplittingFunction>(0, 'q', 'g');
                Pgq0Pqg0_ =
                    std::make_shared<const ConvolutedSplittingFunctions>(
                        0, 'g', 'q', 0, 'q', 'g'
                    );
                Pgg0Pgg0_ =
                    std::make_shared<const ConvolutedSplittingFunctions>(
                        0, 'g', 'g', 0, 'g', 'g'
                    );
                break;
            default:
                throw UnexpectedException(
                    "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                );
            }
        }

        switch (GetOrder()) {
        case 1:
            break;
        case 2:
            asy_ = std::make_unique<AsymptoticCoefficientFunction>(
                order, kind, channel
            );
            thr_ = std::make_unique<ThresholdCoefficientFunction>(
                order, kind, channel
            );
            switch (GetChannel()) {
            case 'q':
                convolutions_lmu1_.push_back(std::make_unique<Convolution>(
                    gluon_as1_, Pgq0_, abserr, relerr, dim
                ));
                break;
            case 'g':
                convolutions_lmu1_.push_back(std::make_unique<Convolution>(
                    gluon_as1_, Pgg0_, abserr, relerr, dim
                ));
                convolutions_lmu1_.push_back(
                    std::make_unique<Convolution>(gluon_as1_, delta_)
                );
                break;
            default:
                throw UnexpectedException(
                    "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                );
            }
            break;
        case 3:
            switch (GetChannel()) {
            case 'q':
                convolutions_lmu1_.push_back(std::make_unique<Convolution>(
                    gluon_as1_, Pgq1_, abserr, relerr, dim
                ));
                convolutions_lmu1_.push_back(std::make_unique<Convolution>(
                    gluon_as2_, Pgq0_, abserr, relerr, dim
                ));
                convolutions_lmu1_.push_back(std::make_unique<Convolution>(
                    quark_as2_, Pqq0_, abserr, relerr, dim
                ));
                convolutions_lmu1_.push_back(
                    std::make_unique<Convolution>(quark_as2_, delta_)
                );

                convolutions_lmu2_.push_back(std::make_unique<Convolution>(
                    gluon_as1_, Pgg0Pgq0_, abserr, relerr, dim
                ));
                convolutions_lmu2_.push_back(std::make_unique<Convolution>(
                    gluon_as1_, Pqq0Pgq0_, abserr, relerr, dim
                ));
                convolutions_lmu2_.push_back(std::make_unique<Convolution>(
                    gluon_as1_, Pgq0_, abserr, relerr, dim
                ));
                break;
            case 'g':
                convolutions_lmu1_.push_back(std::make_unique<Convolution>(
                    gluon_as1_, Pgg1_, abserr, relerr, dim
                ));
                convolutions_lmu1_.push_back(
                    std::make_unique<Convolution>(gluon_as1_, delta_)
                );
                convolutions_lmu1_.push_back(std::make_unique<Convolution>(
                    quark_as2_, Pqg0_, abserr, relerr, dim
                ));
                convolutions_lmu1_.push_back(std::make_unique<Convolution>(
                    gluon_as2_, Pgg0_, abserr, relerr, dim
                ));
                convolutions_lmu1_.push_back(
                    std::make_unique<Convolution>(gluon_as2_, delta_)
                );

                // by default option I integrate with analytical double integral
                // method
                convolutions_lmu2_.push_back(std::make_unique<Convolution>(
                    gluon_as1_, Pgg0Pgg0_, abserr, relerr, dim
                ));
                convolutions_lmu2_.push_back(std::make_unique<Convolution>(
                    gluon_as1_, Pgq0Pqg0_, abserr, relerr, dim
                ));
                convolutions_lmu2_.push_back(std::make_unique<Convolution>(
                    gluon_as1_, Pgg0_, abserr, relerr, dim
                ));
                convolutions_lmu2_.push_back(
                    std::make_unique<Convolution>(gluon_as1_, delta_)
                );
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

        SetFunctions();
    } catch (UnexpectedException &e) {
        e.runtime_error();
    }
}

//==========================================================================================//
//  ExactCoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

ExactCoefficientFunction::ExactCoefficientFunction(ExactCoefficientFunction &obj
)
    : ExactCoefficientFunction(
          obj.GetOrder(), obj.GetKind(), obj.GetChannel(), obj.GetAbsErr(),
          obj.GetRelErr(), obj.GetDim()
      ) {
    SetDoubleIntegralMethod(
        obj.GetDoubleIntegralMethod(), obj.GetAbsErr(), obj.GetRelErr(),
        obj.GetDim(), obj.GetMCcalls()
    );
}

//==========================================================================================//
//  ExactCoefficientFunction: get method for abserr
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::GetAbsErr() const {
    if (GetOrder() > 1)
        return convolutions_lmu1_[0]->GetAbsErr();
    else
        return 1e-3;
}

//==========================================================================================//
//  ExactCoefficientFunction: get method for relerr
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::GetRelErr() const {
    if (GetOrder() > 1)
        return convolutions_lmu1_[0]->GetRelErr();
    else
        return 1e-3;
}

//==========================================================================================//
//  ExactCoefficientFunction: get method for dim
//------------------------------------------------------------------------------------------//

int ExactCoefficientFunction::GetDim() const {
    if (GetOrder() > 1)
        return convolutions_lmu1_[0]->GetDim();
    else
        return 1000;
}

//==========================================================================================//
//  ExactCoefficientFunction: get method for dim
//------------------------------------------------------------------------------------------//

int ExactCoefficientFunction::GetMCcalls() const {
    if (auto ptr =
            dynamic_cast<DoubleConvolution *>(convolutions_lmu2_[0].get())) {
        return ptr->GetMCcalls();
    } else {
        return 25000;
    }
}

//==========================================================================================//
//  ExactCoefficientFunction: function that sets the pointer for mu_indep_ and
//  mu_dep_
//------------------------------------------------------------------------------------------//

void ExactCoefficientFunction::SetFunctions() {
    switch (GetOrder()) {
    case 1:
        switch (GetKind()) {
        case '2':
            mu_indep_ = &ExactCoefficientFunction::C2_g1;
            break;
        case 'L':
            mu_indep_ = &ExactCoefficientFunction::CL_g1;
            break;
        default:
            throw UnexpectedException(
                "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
            );
        }
        mu_dep_ = &ExactCoefficientFunction::ZeroFunction;
        break;
    case 2:
        switch (GetChannel()) {
        case 'q':
            switch (GetKind()) {
            case '2':
                mu_indep_ = &ExactCoefficientFunction::C2_ps20;
                break;
            case 'L':
                mu_indep_ = &ExactCoefficientFunction::CL_ps20;
                break;
            default:
                throw UnexpectedException(
                    "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                );
            }
            mu_dep_ = &ExactCoefficientFunction::C_ps2_MuDep;
            break;
        case 'g':
            switch (GetKind()) {
            case '2':
                mu_indep_ = &ExactCoefficientFunction::C2_g20;
                break;
            case 'L':
                mu_indep_ = &ExactCoefficientFunction::CL_g20;
                break;
            default:
                throw UnexpectedException(
                    "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                );
            }
            mu_dep_ = &ExactCoefficientFunction::C_g2_MuDep;
            break;
        default:
            throw UnexpectedException(
                "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
            );
        }
        break;
    case 3:
        switch (GetChannel()) {
        case 'q':
            mu_dep_ = &ExactCoefficientFunction::C_ps3_MuDep;
            break;
        case 'g':
            mu_dep_ = &ExactCoefficientFunction::C_g3_MuDep;
            break;
        default:
            throw UnexpectedException(
                "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
            );
        }
        mu_indep_ = &ExactCoefficientFunction::WarningFunc;
        break;
    default:
        throw UnexpectedException(
            "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
        );
    }
}

//==========================================================================================//
//  ExactCoefficientFunction: function that sets the double integral method
//------------------------------------------------------------------------------------------//

void ExactCoefficientFunction::SetDoubleIntegralMethod(
    const DoubleIntegralMethod &double_int_method, const double &abserr,
    const double &relerr, const int &dim, const int &MCcalls
) {
    try {

        if (double_int_meth_ == double_int_method) {
            throw NotValidException(
                "Setting double integral method identical to its previous "
                "value!",
                __PRETTY_FUNCTION__, __LINE__
            );
        }

        if (GetOrder() < 3) {
            throw NotPresentException(
                "Double Integration is not needed at order="
                    + to_string(GetOrder()),
                __PRETTY_FUNCTION__, __LINE__
            );
        }

        if (GetChannel() == 'q') {
            throw NotPresentException(
                "Double Integration is not present for channel='q'",
                __PRETTY_FUNCTION__, __LINE__
            );
        }

        // at this point I must be in the g channel at order 3

        double_int_meth_ = double_int_method;

        switch (double_int_method) {
        case DoubleIntegralMethod::MonteCarlo:
            convolutions_lmu2_[0] = std::make_unique<DoubleConvolution>(
                gluon_as1_, Pgg0_, abserr, relerr, dim, true, MCcalls
            );
            break;
        case DoubleIntegralMethod::DoubleNumerical:
            convolutions_lmu2_[0] = std::make_unique<DoubleConvolution>(
                gluon_as1_, Pgg0_, abserr, relerr, dim, false, MCcalls
            );
            break;
        case DoubleIntegralMethod::Analytical:
            if (Pgg0Pgg0_ == nullptr) {
                // TODO: this if is probably useless
                Pgg0Pgg0_ =
                    std::make_shared<const ConvolutedSplittingFunctions>(
                        0, 'g', 'g', 0, 'g', 'g'
                    );
            }
            convolutions_lmu2_[0] = std::make_unique<Convolution>(
                gluon_as1_, Pgg0Pgg0_, abserr, relerr, dim
            );
        }

    } catch (const NotValidException &e) {
        e.runtime_error();
    } catch (const NotPresentException &e) {
        e.warning();
    }
}

//==========================================================================================//
//  ExactCoefficientFunction: central value of the exact coefficient function
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::fx(
    double x, double m2Q2, double m2mu2, int nf
) const {
    return MuIndependentTerms(x, m2Q2, nf)
           + MuDependentTerms(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  ExactCoefficientFunction: mu-independent terms
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::MuIndependentTerms(
    double x, double m2Q2, int nf
) const {
    double res;
    try {
        res = (this->*mu_indep_)(x, m2Q2, nf);
    } catch (NotValidException &e) {
        e.runtime_error();
    } catch (NotKnownException &e) {
        e.runtime_error();
    }
    return res;
}

//==========================================================================================//
//  ExactCoefficientFunction: mu dependent terms
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::MuDependentTerms(
    double x, double m2Q2, double m2mu2, int nf
) const {

    if (x <= 0 || x > xMax(m2Q2))
        return 0.;

    return (this->*mu_dep_)(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  ExactCoefficientFunction: implemented only because it is pure virtual in
//  base class. Returning three identical values
//------------------------------------------------------------------------------------------//

Value ExactCoefficientFunction::fxBand(
    double x, double m2Q2, double m2mu2, int nf
) const {
    return Value(fx(x, m2Q2, m2mu2, nf));
}

//==========================================================================================//
//  Exact massive gluon coefficient functions for F2 at O(as)
//
//  Eq. (50) from Ref. [arXiv:1001.2312]
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C2_g1(double x, double m2Q2, int /*nf*/)
    const { // m2Q2=m^2/Q^2

    double beta = sqrt(1. - 4. * m2Q2 * x / (1 - x));
    double x2 = x * x;
    double m4Q4 = m2Q2 * m2Q2;
    double L = log((1. + beta) / (1. - beta));

    return 4. * TR
           * (L
                  * (-8. * x2 * m4Q4 - 4. * x * m2Q2 * (3. * x - 1) + 2. * x2
                     - 2. * x + 1.)
              + beta * (4. * m2Q2 * x * (x - 1.) - (8. * x2 - 8. * x + 1.)));
}

//==========================================================================================//
//  Exact massive gluon coefficient functions for FL at O(as)
//
//  Eq. (2.9) from Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double
    ExactCoefficientFunction::CL_g1(double x, double m2Q2, int /*nf*/) const {

    double beta = sqrt(1. - 4. * m2Q2 * x / (1. - x));
    double x2 = x * x;
    double L = log((1. + beta) / (1. - beta));

    return 16. * TR * (x * (1. - x) * beta - 2. * x2 * m2Q2 * L);
}

//==========================================================================================//
//  OBSERVATION: in the O(as^2) exact coefficeint functions the
//  mu-independent part is an interpolation in a certain grid. When this
//  function is called for a (x,Q) value outside this grid, te returned value is
//  the appropriate limit.
//------------------------------------------------------------------------------------------//

/// @cond UNNECESSARY
/**
 * @name Fortran massive coefficient functions
 * Fortran functions for the O(as^2)
 * coefficient functions from 'src/hqcoef.f'.
 */
///@{
extern "C" {
    // double c2log_(double *wr,double *xi);
    double c2nlog_(double *wr, double *xi);
    double clnlog_(double *wr, double *xi);
    double c2nloq_(double *wr, double *xi);
    double clnloq_(double *wr, double *xi);
    // double c2nlobarg_(double *wr, double *xi);
    // double clnlobarg_(double *wr, double *xi);
    // double c2nlobarq_(double *wr, double *xi);
    // double clnlobarq_(double *wr, double *xi);
}
///@}
/// \endcond

//==========================================================================================//
//  mu independent part of the exact massive quark coefficient functions for F2
//  at O(as^2)
//
//  Exact (but numerical) result from [arXiv:hep-ph/9411431].
//  Taken from the Fortran code 'src/hqcoef.f'
//------------------------------------------------------------------------------------------//

double
    ExactCoefficientFunction::C2_ps20(double x, double m2Q2, int /*nf*/) const {

    double xi = 1. / m2Q2;
    double eta = 0.25 * xi * (1 - x) / x - 1.;

    // using nf = nan since at O(as2) the coefficient functions don't depend on
    // nf
    int nf = static_cast<int>(nan(""));

    if (eta < 1e-6)
        return thr_->MuIndependentTerms(x, m2Q2, nf);
    if (eta > 1e6 || xi > 1e5)
        return asy_->MuIndependentTerms(x, m2Q2, nf);
    if (xi < 1e-3) {
        throw NotValidException(
            "max value of m2Q2 is 1e3. Got m2Q2=" + to_string(m2Q2),
            __PRETTY_FUNCTION__, __LINE__
        );
    }

    return 16 * M_PI * xi * c2nloq_(&eta, &xi) / x;
}

//==========================================================================================//
//  Exact massive quark coefficient functions at O(as^2):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.1) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C_ps21(double x, double m2Q2) const {

    int nf = static_cast<int>(nan(""));

    return -convolutions_lmu1_[0]->Convolute(x, m2Q2, nf);
    // The minus sign comes from the fact that in [arXiv:1205.5727]
    // the expansion is performed in terms of log(m^2/mu^2)
    // (even if it says the opposite) but we are
    // expanding in terms of log(mu^2/m^2)
}

//==========================================================================================//
//  mu independent part of the exact massive quark coefficient functions for FL
//  at O(as)
//
//  Exact (but numerical) result from [arXiv:hep-ph/9411431].
//  Taken from the Fortran code 'src/hqcoef.f'
//------------------------------------------------------------------------------------------//

double
    ExactCoefficientFunction::CL_ps20(double x, double m2Q2, int /*nf*/) const {

    double xi = 1. / m2Q2;
    double eta = 0.25 * xi * (1 - x) / x - 1.;

    // using nf = nan since at O(as2) the coefficient functions don't depend on
    // nf
    int nf = static_cast<int>(nan(""));

    if (eta < 1e-6)
        return thr_->MuIndependentTerms(x, m2Q2, nf);
    if (eta > 1e6 || xi > 1e5)
        return asy_->MuIndependentTerms(x, m2Q2, nf);
    if (xi < 1e-3) {
        throw NotValidException(
            "max value of m2Q2 is 1e3. Got m2Q2=" + to_string(m2Q2),
            __PRETTY_FUNCTION__, __LINE__
        );
    }

    return 16 * M_PI * xi * clnloq_(&eta, &xi) / x;
}

//==========================================================================================//
//  mu independent part of the exact massive gluon coefficient functions for F2
//  at O(as^2)
//
//  Exact (but numerical) result from [arXiv:hep-ph/9411431].
//  Taken from the Fortran code 'src/hqcoef.f'
//------------------------------------------------------------------------------------------//

double
    ExactCoefficientFunction::C2_g20(double x, double m2Q2, int /*nf*/) const {

    double xi = 1. / m2Q2;
    double eta = 0.25 * xi * (1 - x) / x - 1.;

    // using nf = nan since at O(as2) the coefficient functions don't depend on
    // nf
    int nf = static_cast<int>(nan(""));

    if (eta < 1e-6)
        return thr_->MuIndependentTerms(x, m2Q2, nf);
    if (eta > 1e6 || xi > 1e5)
        return asy_->MuIndependentTerms(x, m2Q2, nf);
    if (xi < 1e-3) {
        throw NotValidException(
            "max value of m2Q2 is 1e3. Got m2Q2=" + to_string(m2Q2),
            __PRETTY_FUNCTION__, __LINE__
        );
    }

    return 16 * M_PI * xi * c2nlog_(&eta, &xi) / x;
}

//==========================================================================================//
//  Exact massive gluon coefficient functions at O(as^2):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.4) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C_g21(double x, double m2Q2) const {

    int nf_one = 1;
    // Put nf to 1 since the nf contribution cancels for any value of nf

    std::future<double> future_f0 = std::async(
        std::launch::async, &AbstractConvolution::Convolute,
        std::ref(*convolutions_lmu1_[0]), x, m2Q2, nf_one
    );
    std::future<double> future_f1 = std::async(
        std::launch::async, &AbstractConvolution::Convolute,
        std::ref(*convolutions_lmu1_[1]), x, m2Q2, nf_one
    );

    return -(future_f0.get() - future_f1.get() * beta0(nf_one));
}

//==========================================================================================//
//  mu independent part of the exact massive gluon coefficient functions for FL
//  at O(as^2)
//
//  Exact (but numerical) result from [arXiv:hep-ph/9411431].
//  Taken from the Fortran code 'src/hqcoef.f'
//------------------------------------------------------------------------------------------//

double
    ExactCoefficientFunction::CL_g20(double x, double m2Q2, int /*nf*/) const {

    double xi = 1. / m2Q2;
    double eta = 0.25 * xi * (1 - x) / x - 1.;

    // using nf = nan since at O(as2) the coefficient functions don't depend on
    // nf
    int nf = static_cast<int>(nan(""));

    if (eta < 1e-6)
        return thr_->MuIndependentTerms(x, m2Q2, nf);
    if (eta > 1e6 || xi > 1e5)
        return asy_->MuIndependentTerms(x, m2Q2, nf);
    if (xi < 1e-3) {
        throw NotValidException(
            "max value of m2Q2 is 1e3. Got m2Q2=" + to_string(m2Q2),
            __PRETTY_FUNCTION__, __LINE__
        );
    }

    return 16 * M_PI * xi * clnlog_(&eta, &xi) / x;
}

//==========================================================================================//
//  ExactCoefficientFunction: mu dependent part of the exact massive gluon
//  coefficient function at O(as^2):
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::
    C_g2_MuDep(double x, double m2Q2, double m2mu2, int /*nf*/) const {

    return C_g21(x, m2Q2) * log(1. / m2mu2);
}

//==========================================================================================//
//  ExactCoefficientFunction: mu dependent part of the exact massive quark
//  coefficient function at O(as^2):
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::
    C_ps2_MuDep(double x, double m2Q2, double m2mu2, int /*nf*/) const {

    return C_ps21(x, m2Q2) * log(1. / m2mu2);
}

//==========================================================================================//
//  Exact massive quark coefficient functions at O(as^3):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.2) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C_ps31(double x, double m2Q2, int nf) const {

    std::future<double> future_f0 = std::async(
        std::launch::async, &AbstractConvolution::Convolute,
        std::ref(*convolutions_lmu1_[0]), x, m2Q2, nf
    );
    std::future<double> future_f1 = std::async(
        std::launch::async, &AbstractConvolution::Convolute,
        std::ref(*convolutions_lmu1_[1]), x, m2Q2, nf
    );
    std::future<double> future_f2 = std::async(
        std::launch::async, &AbstractConvolution::Convolute,
        std::ref(*convolutions_lmu1_[2]), x, m2Q2, nf
    );
    std::future<double> future_f3 = std::async(
        std::launch::async, &AbstractConvolution::Convolute,
        std::ref(*convolutions_lmu1_[3]), x, m2Q2, nf
    );

    return -(
        future_f0.get() + future_f1.get() + future_f2.get()
        - 2. * beta0(nf) * future_f3.get()
    );
}

//==========================================================================================//
//  Exact massive quark coefficient functions at O(as^3):
//  Term proportional to log(mu^2/m^2)^2
//
//  Eq. (4.3) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C_ps32(double x, double m2Q2, int nf) const {

    std::future<double> future_f0 = std::async(
        std::launch::async, &AbstractConvolution::Convolute,
        std::ref(*convolutions_lmu2_[0]), x, m2Q2, nf
    );
    std::future<double> future_f1 = std::async(
        std::launch::async, &AbstractConvolution::Convolute,
        std::ref(*convolutions_lmu2_[1]), x, m2Q2, nf
    );
    std::future<double> future_f2 = std::async(
        std::launch::async, &AbstractConvolution::Convolute,
        std::ref(*convolutions_lmu2_[2]), x, m2Q2, nf
    );

    return 0.5 * (future_f0.get() + future_f1.get())
           - 3. / 2 * beta0(nf) * future_f2.get();
}

//==========================================================================================//
//  Exact massive gluon coefficient functions for F2 at O(as^3):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.5) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C_g31(double x, double m2Q2, int nf) const {

    std::future<double> future_f0 = std::async(
        std::launch::async, &AbstractConvolution::Convolute,
        std::ref(*convolutions_lmu1_[0]), x, m2Q2, nf
    );
    std::future<double> future_f1 = std::async(
        std::launch::async, &AbstractConvolution::Convolute,
        std::ref(*convolutions_lmu1_[1]), x, m2Q2, nf
    );
    std::future<double> future_f2 = std::async(
        std::launch::async, &AbstractConvolution::Convolute,
        std::ref(*convolutions_lmu1_[2]), x, m2Q2, nf
    );
    std::future<double> future_f3 = std::async(
        std::launch::async, &AbstractConvolution::Convolute,
        std::ref(*convolutions_lmu1_[3]), x, m2Q2, nf
    );
    std::future<double> future_f4 = std::async(
        std::launch::async, &AbstractConvolution::Convolute,
        std::ref(*convolutions_lmu1_[4]), x, m2Q2, nf
    );

    return -(
        future_f0.get() - beta1(nf) * future_f1.get() + future_f2.get()
        + future_f3.get() - 2. * beta0(nf) * future_f4.get()
    );
}

//==========================================================================================//
//  Exact massive gluon coefficient functions at O(as^3):
//  Term proportional to log(mu^2/m^2)^2
//
//  Eq. (4.6) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C_g32(double x, double m2Q2, int nf) const {

    double beta_0 = beta0(nf);

    std::future<double> future_f0 = std::async(
        std::launch::async, &AbstractConvolution::Convolute,
        std::ref(*convolutions_lmu2_[0]), x, m2Q2, nf
    );
    std::future<double> future_f1 = std::async(
        std::launch::async, &AbstractConvolution::Convolute,
        std::ref(*convolutions_lmu2_[1]), x, m2Q2, nf
    );
    std::future<double> future_f2 = std::async(
        std::launch::async, &AbstractConvolution::Convolute,
        std::ref(*convolutions_lmu2_[2]), x, m2Q2, nf
    );
    std::future<double> future_f3 = std::async(
        std::launch::async, &AbstractConvolution::Convolute,
        std::ref(*convolutions_lmu2_[3]), x, m2Q2, nf
    );
    return 0.5 * (future_f0.get() + future_f1.get())
           - 3. / 2 * beta_0 * future_f2.get()
           + beta_0 * beta_0 * future_f3.get();
}

//==========================================================================================//
//  ExactCoefficientFunction: mu dependent part of the exact massive gluon
//  coefficient function at O(as^3):
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C_g3_MuDep(
    double x, double m2Q2, double m2mu2, int nf
) const {

    double lmu = log(1. / m2mu2);

    std::future<double> future_f1 = std::async(
        std::launch::async, &ExactCoefficientFunction::C_g31, this, x, m2Q2, nf
    );
    std::future<double> future_f2 = std::async(
        std::launch::async, &ExactCoefficientFunction::C_g32, this, x, m2Q2, nf
    );

    return future_f1.get() * lmu + future_f2.get() * lmu * lmu;
}

//==========================================================================================//
//  ExactCoefficientFunction: mu dependent part of the exact massive quark
//  coefficient function at O(as^3):
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C_ps3_MuDep(
    double x, double m2Q2, double m2mu2, int nf
) const {

    double lmu = log(1. / m2mu2);

    std::future<double> future_f1 = std::async(
        std::launch::async, &ExactCoefficientFunction::C_ps31, this, x, m2Q2, nf
    );
    std::future<double> future_f2 = std::async(
        std::launch::async, &ExactCoefficientFunction::C_ps32, this, x, m2Q2, nf
    );

    return future_f1.get() * lmu + future_f2.get() * lmu * lmu;
}

//==========================================================================================//
//  ExactCoefficientFunction: print warning for mu independent part of O(as^3)
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::
    WarningFunc(double /*x*/, double /*m2Q2*/, int /*nf*/) const {
    throw NotKnownException(
        "mu-independent terms of the exact coefficient function at order=3 are "
        "not known!",
        __PRETTY_FUNCTION__, __LINE__
    );
    return 0.;
}
