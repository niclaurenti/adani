#include "adani/HighEnergyCoefficientFunction.h"
#include "adani/Constants.h"
#include "adani/SpecialFunctions.h"
#include <cmath>

//==========================================================================================//
//  AbstractHighEnergyCoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

AbstractHighEnergyCoefficientFunction::AbstractHighEnergyCoefficientFunction(
    const int &order, const char &kind, const char &channel, const bool &NLL
)
    : CoefficientFunction(order, kind, channel), NLL_(NLL) {
    switch (order) {
    case 1:
        fx_ = &AbstractHighEnergyCoefficientFunction::ZeroFunctionBand;
        break;
    case 2:
        fx_ = &AbstractHighEnergyCoefficientFunction::Order2;
        break;
    case 3:
        if (NLL) {
            fx_ = &AbstractHighEnergyCoefficientFunction::Order3;
        } else {
            fx_ = &AbstractHighEnergyCoefficientFunction::Order3LL;
        }
        break;
    default:
        throw UnexpectedException(
            "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
        );
    }
};

//==========================================================================================//
//  AbstractHighEnergyCoefficientFunction: band of the high energy coefficient
//  function
//------------------------------------------------------------------------------------------//

Value AbstractHighEnergyCoefficientFunction::Order2(
    double x, double m2Q2, double m2mu2, int /*nf*/
) const {
    return Value(LL(m2Q2, m2mu2) / x);
}

//==========================================================================================//
//  AbstractHighEnergyCoefficientFunction: band of the high energy coefficient
//  function
//------------------------------------------------------------------------------------------//

Value AbstractHighEnergyCoefficientFunction::Order3(
    double x, double m2Q2, double m2mu2, int nf
) const {
    return (LL(m2Q2, m2mu2) * log(x) + NLL(m2Q2, m2mu2, nf)) / x;
}

//==========================================================================================//
//  AbstractHighEnergyCoefficientFunction: band of the high energy coefficient
//  function
//------------------------------------------------------------------------------------------//

Value AbstractHighEnergyCoefficientFunction::Order3LL(
    double x, double m2Q2, double m2mu2, int /*nf*/
) const {
    return Value(LL(m2Q2, m2mu2) * log(x) / x);
}

//==========================================================================================//
//  AbstractHighEnergyCoefficientFunction: band of the high energy coefficient
//  function
//------------------------------------------------------------------------------------------//

Value AbstractHighEnergyCoefficientFunction::fxBand(
    double x, double m2Q2, double m2mu2, int nf
) const {
    return (this->*fx_)(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  HighEnergyCoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

HighEnergyCoefficientFunction::HighEnergyCoefficientFunction(
    const int &order, const char &kind, const char &channel, const bool &NLL
)
    : AbstractHighEnergyCoefficientFunction(order, kind, channel, NLL) {
    try {
        SetFunctions();
    } catch (UnexpectedException &e) {
        e.runtime_error();
    }
}

//==========================================================================================//
//  HighEnergyCoefficientFunction: copy constructor
//------------------------------------------------------------------------------------------//

HighEnergyCoefficientFunction::HighEnergyCoefficientFunction(
    const HighEnergyCoefficientFunction &obj
)
    : HighEnergyCoefficientFunction(
          obj.GetOrder(), obj.GetKind(), obj.GetChannel(), obj.GetNLL()
      ) {}

//==========================================================================================//
//  HighEnergyCoefficientFunction: function that sets the pointer for LL_ and
//  NLL_
//------------------------------------------------------------------------------------------//

void HighEnergyCoefficientFunction::SetFunctions() {

    switch (GetOrder()) {
    case 1:
        LL_ = nullptr;
        NLL_ = nullptr;
        break;
    case 2:
        NLL_ = nullptr;
        switch (GetKind()) {
        case '2':
            switch (GetChannel()) {
            case 'g':
                LL_ = &HighEnergyCoefficientFunction::C2_g2_highenergyLL;
                break;
            case 'q':
                LL_ = &HighEnergyCoefficientFunction::C2_ps2_highenergyLL;
                break;
            default:
                throw UnexpectedException(
                    "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                );
            }
            break;
        case 'L':
            switch (GetChannel()) {
            case 'g':
                LL_ = &HighEnergyCoefficientFunction::CL_g2_highenergyLL;
                break;
            case 'q':
                LL_ = &HighEnergyCoefficientFunction::CL_ps2_highenergyLL;
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
        break;
    case 3:
        switch (GetKind()) {
        case '2':
            switch (GetChannel()) {
            case 'g':
                LL_ = &HighEnergyCoefficientFunction::C2_g3_highenergyLL;
                NLL_ = &HighEnergyCoefficientFunction::C2_g3_highenergyNLL;
                break;
            case 'q':
                LL_ = &HighEnergyCoefficientFunction::C2_ps3_highenergyLL;
                NLL_ = &HighEnergyCoefficientFunction::C2_ps3_highenergyNLL;
                break;
            default:
                throw UnexpectedException(
                    "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                );
            }
            break;
        case 'L':
            switch (GetChannel()) {
            case 'g':
                LL_ = &HighEnergyCoefficientFunction::CL_g3_highenergyLL;
                NLL_ = &HighEnergyCoefficientFunction::CL_g3_highenergyNLL;
                break;
            case 'q':
                LL_ = &HighEnergyCoefficientFunction::CL_ps3_highenergyLL;
                NLL_ = &HighEnergyCoefficientFunction::CL_ps3_highenergyNLL;
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
        break;
    default:
        throw UnexpectedException(
            "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
        );
    }

    if (!GetNLL())
        NLL_ = &HighEnergyCoefficientFunction::ThrowException;
}

//==========================================================================================//
//  HighEnergyCoefficientFunction: band of the high energy coefficient function
//------------------------------------------------------------------------------------------//

double HighEnergyCoefficientFunction::LL(double m2Q2, double m2mu2) const {
    return (this->*LL_)(m2Q2, m2mu2);
}

//==========================================================================================//
//  HighEnergyCoefficientFunction: band of the high energy coefficient function
//------------------------------------------------------------------------------------------//

Value HighEnergyCoefficientFunction::NLL(
    double m2Q2, double m2mu2, int nf
) const {
    return (this->*NLL_)(m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  HighEnergyCoefficientFunction: function to be used when NLL_=false
//------------------------------------------------------------------------------------------//

Value HighEnergyCoefficientFunction::
    ThrowException(double /*m2Q2*/, double /*m2mu2*/, int /*nf*/) const {
    throw NotValidException(
        "Called HighEnergyCoefficientFunction::NLL with NLL_=false!",
        __PRETTY_FUNCTION__, __LINE__
    );
}

//==========================================================================================//
//  HighEnergyHighScaleCoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

HighEnergyHighScaleCoefficientFunction::HighEnergyHighScaleCoefficientFunction(
    const int &order, const char &kind, const char &channel, const bool &NLL
)
    : AbstractHighEnergyCoefficientFunction(order, kind, channel, NLL) {
    try {
        SetFunctions();
    } catch (UnexpectedException &e) {
        e.runtime_error();
    }
}

//==========================================================================================//
//  HighEnergyHighScaleCoefficientFunction: copy constructor
//------------------------------------------------------------------------------------------//

HighEnergyHighScaleCoefficientFunction::HighEnergyHighScaleCoefficientFunction(
    const HighEnergyHighScaleCoefficientFunction &obj
)
    : HighEnergyHighScaleCoefficientFunction(
          obj.GetOrder(), obj.GetKind(), obj.GetChannel(), obj.GetNLL()
      ) {}

//==========================================================================================//
//  HighEnergyHighScaleCoefficientFunction: function that sets the pointer for
//  LL_ and NLL_
//------------------------------------------------------------------------------------------//

void HighEnergyHighScaleCoefficientFunction::SetFunctions() {

    switch (GetOrder()) {
    case 1:
        LL_ = nullptr;
        NLL_ = nullptr;
        break;
    case 2:
        NLL_ = nullptr;
        switch (GetKind()) {
        case '2':
            switch (GetChannel()) {
            case 'g':
                LL_ = &HighEnergyHighScaleCoefficientFunction::
                          C2_g2_highenergy_highscaleLL;
                break;
            case 'q':
                LL_ = &HighEnergyHighScaleCoefficientFunction::
                          C2_ps2_highenergy_highscaleLL;
                break;
            default:
                throw UnexpectedException(
                    "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                );
            }
            break;
        case 'L':
            switch (GetChannel()) {
            case 'g':
                LL_ = &HighEnergyHighScaleCoefficientFunction::
                          CL_g2_highenergy_highscaleLL;
                break;
            case 'q':
                LL_ = &HighEnergyHighScaleCoefficientFunction::
                          CL_ps2_highenergy_highscaleLL;
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
        break;
    case 3:
        switch (GetKind()) {
        case '2':
            switch (GetChannel()) {
            case 'g':
                LL_ = &HighEnergyHighScaleCoefficientFunction::
                          C2_g3_highenergy_highscaleLL;
                NLL_ = &HighEnergyHighScaleCoefficientFunction::
                           C2_g3_highenergy_highscaleNLL;
                break;
            case 'q':
                LL_ = &HighEnergyHighScaleCoefficientFunction::
                          C2_ps3_highenergy_highscaleLL;
                NLL_ = &HighEnergyHighScaleCoefficientFunction::
                           C2_ps3_highenergy_highscaleNLL;
                break;
            default:
                throw UnexpectedException(
                    "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
                );
            }
            break;
        case 'L':
            switch (GetChannel()) {
            case 'g':
                LL_ = &HighEnergyHighScaleCoefficientFunction::
                          CL_g3_highenergy_highscaleLL;
                NLL_ = &HighEnergyHighScaleCoefficientFunction::
                           CL_g3_highenergy_highscaleNLL;
                break;
            case 'q':
                LL_ = &HighEnergyHighScaleCoefficientFunction::
                          CL_ps3_highenergy_highscaleLL;
                NLL_ = &HighEnergyHighScaleCoefficientFunction::
                           CL_ps3_highenergy_highscaleNLL;
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
        break;
    default:
        throw UnexpectedException(
            "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
        );
    }

    if (!GetNLL())
        NLL_ = &HighEnergyHighScaleCoefficientFunction::ThrowException;
}

//==========================================================================================//
//  HighEnergyCoefficientFunction: band of the high energy coefficient function
//------------------------------------------------------------------------------------------//

double HighEnergyHighScaleCoefficientFunction::LL(
    double m2Q2, double m2mu2
) const {
    return (this->*LL_)(m2Q2, m2mu2);
    ;
}

//==========================================================================================//
//  HighEnergyCoefficientFunction: band of the high energy coefficient function
//------------------------------------------------------------------------------------------//

Value HighEnergyHighScaleCoefficientFunction::NLL(
    double m2Q2, double m2mu2, int nf
) const {
    return (this->*NLL_)(m2Q2, m2mu2, nf);
    ;
}

//==========================================================================================//
//  HighEnergyHighScaleCoefficientFunction: function to be used when NLL_=false
//------------------------------------------------------------------------------------------//

Value HighEnergyHighScaleCoefficientFunction::
    ThrowException(double /*m2Q2*/, double /*m2mu2*/, int /*nf*/) const {
    throw NotValidException(
        "Called HighEnergyHighScaleCoefficientFunction::NLL with NLL_=false!",
        __PRETTY_FUNCTION__, __LINE__
    );
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for F2 at O(as^2).
//
//  Eq. (3.38) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double HighEnergyCoefficientFunction::C2_g2_highenergyLL(
    double m2Q2, double m2mu2
) const {
    return a_11() * Coff2_1(m2Q2,m2mu2);
}

//==========================================================================================//
//  High energy limit of the quark coefficient function for F2 at O(as^2).
//------------------------------------------------------------------------------------------//

double HighEnergyCoefficientFunction::C2_ps2_highenergyLL(
    double m2Q2, double m2mu2
) const {

    return CF / CA * C2_g2_highenergyLL(m2Q2, m2mu2);
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for FL at O(as^2).
//
//  Eq. (16, 21) of Ref. [arXiv:hep-ph/9411431]
//------------------------------------------------------------------------------------------//

double HighEnergyCoefficientFunction::CL_g2_highenergyLL(
    double m2Q2, double m2mu2
) const {
    return a_11() * CoffL_1(m2Q2,m2mu2);
}

//==========================================================================================//
//  High energy limit of the quark coefficient function for FL at O(as^2).
//------------------------------------------------------------------------------------------//

double HighEnergyCoefficientFunction::CL_ps2_highenergyLL(
    double m2Q2, double m2mu2
) const {

    return CF / CA * CL_g2_highenergyLL(m2Q2, m2mu2);
}

//==========================================================================================//
//  Q^2>>m^2 limit of the high energy limit of the gluon coefficient function
//  for F2 at O(as^2).
//
//  Eq. (3.40) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double HighEnergyHighScaleCoefficientFunction::C2_g2_highenergy_highscaleLL(
    double m2Q2, double m2mu2
) const {
    return a_11() * Coff2_1(m2Q2,m2mu2);
}

//==========================================================================================//
//  Q^2>>m^2 limit of the high energy limit of the quark coefficient function
//  for F2 at O(as^2).
//------------------------------------------------------------------------------------------//

double HighEnergyHighScaleCoefficientFunction::C2_ps2_highenergy_highscaleLL(
    double m2Q2, double m2mu2
) const {

    return CF / CA * C2_g2_highenergy_highscaleLL(m2Q2, m2mu2);
}

//==========================================================================================//
//  Q^2>>m^2 limit of the high energy limit of the gluon coefficient function
//  for FL at O(as^2).
//
//  Q^2>>m^2 limit of Eq. (16, 21) of Ref. [arXiv:hep-ph/9411431]
//------------------------------------------------------------------------------------------//

double HighEnergyHighScaleCoefficientFunction::CL_g2_highenergy_highscaleLL(
    double m2Q2, double m2mu2
) const {
    return a_11() * CoffL_1(m2Q2,m2mu2);
}

//==========================================================================================//
//  Q^2>>m^2 limit of the high energy limit of the quark coefficient function
//  for FL at O(as^2).
//------------------------------------------------------------------------------------------//

double HighEnergyHighScaleCoefficientFunction::CL_ps2_highenergy_highscaleLL(
    double m2Q2, double m2mu2
) const {

    return CF / CA * CL_g2_highenergy_highscaleLL(m2Q2, m2mu2);
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for F2 at O(as^3)
//  at leading log.
//
//  Eq. (3.39) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double HighEnergyCoefficientFunction::C2_g3_highenergyLL(
    double m2Q2, double m2mu2
) const {
    return -a_11()*a_11() * Coff2_2(m2Q2,m2mu2);
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for F2 at O(as^3)
//  at next to leading log.
//
//  Eq. (3.39) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

Value HighEnergyCoefficientFunction::C2_g3_highenergyNLL(
    double m2Q2, double m2mu2, int nf
) const {

    double a11 = a_11();
    double a21 = a_21(nf);
    double a21_new = a_21_new(nf);
    double a10 = a_10(nf);
    double beta_0 = beta0(nf);

    double central_value =
        C2_g3_highenergyNLL(m2Q2, m2mu2, a11, a10, a21, beta_0);
    double error = C2_g3_highenergyNLL(m2Q2, m2mu2, a11, a10, a21_new, beta_0);

    double delta = std::abs(central_value - error);

    double higher = central_value + delta;
    double lower = central_value - delta;

    return Value(central_value, higher, lower);
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for F2 at O(as^3)
//  at next to leading log.
//
//  Eq. (3.39) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double HighEnergyCoefficientFunction::C2_g3_highenergyNLL(
    double m2Q2, double m2mu2, double a11, double a10, double a21, double b0
) const {
    //return 0;
    return (2*a11*a10 - a11*b0) * Coff2_2(m2Q2,m2mu2) + a21 * Coff2_1(m2Q2,m2mu2);
}

//==========================================================================================//
//  High energy limit of the quark coefficient function for F2 at O(as^3)
//  at leading log.
//------------------------------------------------------------------------------------------//

double HighEnergyCoefficientFunction::C2_ps3_highenergyLL(
    double m2Q2, double m2mu2
) const {

    return CF / CA * C2_g3_highenergyLL(m2Q2, m2mu2);
}

//==========================================================================================//
//  High energy limit of the quark coefficient function for F2 at O(as^3)
//  at next-to-leading log.
//------------------------------------------------------------------------------------------//

Value HighEnergyCoefficientFunction::C2_ps3_highenergyNLL(
    double m2Q2, double m2mu2, int nf
) const {

    return CF / CA * C2_g3_highenergyNLL(m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for FL at O(as^3)
//  at leading log.
//------------------------------------------------------------------------------------------//

double HighEnergyCoefficientFunction::CL_g3_highenergyLL(
    double m2Q2, double m2mu2
) const {
    return -a_11()*a_11() * CoffL_2(m2Q2,m2mu2);
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for FL at O(as^3)
//  at next to leading log.
//------------------------------------------------------------------------------------------//

Value HighEnergyCoefficientFunction::CL_g3_highenergyNLL(
    double m2Q2, double m2mu2, int nf
) const {

    double a11 = a_11();
    double a21 = a_21(nf);
    double a21_new = a_21_new(nf);
    double a10 = a_10(nf);
    double beta_0 = beta0(nf);

    double central_value =
        CL_g3_highenergyNLL(m2Q2, m2mu2, a11, a10, a21, beta_0);
    double error = CL_g3_highenergyNLL(m2Q2, m2mu2, a11, a10, a21_new, beta_0);

    double delta = std::abs(central_value - error);

    double higher = central_value + delta;
    double lower = central_value - delta;

    return Value(central_value, higher, lower);
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for FL at O(as^3)
//  at next to leading log.
//------------------------------------------------------------------------------------------//

double HighEnergyCoefficientFunction::CL_g3_highenergyNLL(
    double m2Q2, double m2mu2, double a11, double a10, double a21, double b0
) const {
    //return 0;
    return (2*a11*a10 - a11*b0) * CoffL_2(m2Q2,m2mu2) + a21 * CoffL_1(m2Q2,m2mu2);
}

//==========================================================================================//
//  High energy limit of the quark coefficient function for FL at O(as^3)
//  at leading log.
//------------------------------------------------------------------------------------------//

double HighEnergyCoefficientFunction::CL_ps3_highenergyLL(
    double m2Q2, double m2mu2
) const {

    return CF / CA * CL_g3_highenergyLL(m2Q2, m2mu2);
}

//==========================================================================================//
//  High energy limit of the quark coefficient function for FL at O(as^3)
//  at next-to-leading log.
//------------------------------------------------------------------------------------------//

Value HighEnergyCoefficientFunction::CL_ps3_highenergyNLL(
    double m2Q2, double m2mu2, int nf
) const {

    return CF / CA * CL_g3_highenergyNLL(m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the gluon coefficient function
//  for F2 at O(as^3) at leading log.
//
//  Eq. (3.41) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double HighEnergyHighScaleCoefficientFunction::C2_g3_highenergy_highscaleLL(
    double m2Q2, double m2mu2
) const {
    return -a_11()*a_11() * Coff2_2(m2Q2,m2mu2);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the gluon coefficient function
//  for F2 at O(as^3) at next to leading log.
//------------------------------------------------------------------------------------------//

Value HighEnergyHighScaleCoefficientFunction::C2_g3_highenergy_highscaleNLL(
    double m2Q2, double m2mu2, int nf
) const {

    double a11 = a_11();
    double a21 = a_21(nf);
    double a21_new = a_21_new(nf);
    double a10 = a_10(nf);
    double beta_0 = beta0(nf);

    double central_value =
        C2_g3_highenergy_highscaleNLL(m2Q2, m2mu2, a11, a10, a21, beta_0);
    double error =
        C2_g3_highenergy_highscaleNLL(m2Q2, m2mu2, a11, a10, a21_new, beta_0);

    double delta = std::abs(central_value - error);

    double higher = central_value + delta;
    double lower = central_value - delta;

    return Value(central_value, higher, lower);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the gluon coefficient function
//  for F2 at O(as^3) at next to leading log.
//------------------------------------------------------------------------------------------//

double HighEnergyHighScaleCoefficientFunction::C2_g3_highenergy_highscaleNLL(
    double m2Q2, double m2mu2, double a11, double a10, double a21, double b0
) const {
    //return 0;
    return (2*a11*a10 - a11*b0) * Coff2_2(m2Q2,m2mu2) + a21 * Coff2_1(m2Q2,m2mu2);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the quark coefficient function
//  for F2 at O(as^3) at leading log.
//------------------------------------------------------------------------------------------//

double HighEnergyHighScaleCoefficientFunction::C2_ps3_highenergy_highscaleLL(
    double m2Q2, double m2mu2
) const {

    return CF / CA * C2_g3_highenergy_highscaleLL(m2Q2, m2mu2);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the quark coefficient function
//  for F2 at O(as^3) at next-to-leading log.
//------------------------------------------------------------------------------------------//

Value HighEnergyHighScaleCoefficientFunction::C2_ps3_highenergy_highscaleNLL(
    double m2Q2, double m2mu2, int nf
) const {

    return CF / CA * C2_g3_highenergy_highscaleNLL(m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the gluon coefficient function
//  for FL at O(as^3) at leading log.
//------------------------------------------------------------------------------------------//

double HighEnergyHighScaleCoefficientFunction::CL_g3_highenergy_highscaleLL(
    double m2Q2, double m2mu2
) const {
    return -a_11()*a_11() * CoffL_2(m2Q2,m2mu2);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the gluon coefficient function
//  for FL at O(as^3) at next to leading log.
//------------------------------------------------------------------------------------------//

Value HighEnergyHighScaleCoefficientFunction::CL_g3_highenergy_highscaleNLL(
    double m2Q2, double m2mu2, int nf
) const {

    double a11 = a_11();
    double a21 = a_21(nf);
    double a21_new = a_21_new(nf);
    double a10 = a_10(nf);
    double beta_0 = beta0(nf);

    double central_value =
        CL_g3_highenergy_highscaleNLL(m2Q2, m2mu2, a11, a10, a21, beta_0);
    double error =
        CL_g3_highenergy_highscaleNLL(m2Q2, m2mu2, a11, a10, a21_new, beta_0);

    double delta = std::abs(central_value - error);

    double higher = central_value + delta;
    double lower = central_value - delta;

    return Value(central_value, higher, lower);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the gluon coefficient function
//  for FL at O(as^3) at next to leading log.
//------------------------------------------------------------------------------------------//

double HighEnergyHighScaleCoefficientFunction::CL_g3_highenergy_highscaleNLL(
    double m2Q2, double m2mu2, double a11, double a10, double a21, double b0
) const {
    //return 0;
    return (2*a11*a10 - a11*b0) * CoffL_2(m2Q2,m2mu2) + a21 * CoffL_1(m2Q2,m2mu2);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the quark coefficient function
//  for F2 at O(as^3) at leading log.
//------------------------------------------------------------------------------------------//

double HighEnergyHighScaleCoefficientFunction::CL_ps3_highenergy_highscaleLL(
    double m2Q2, double m2mu2
) const {

    return CF / CA * CL_g3_highenergy_highscaleLL(m2Q2, m2mu2);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the quark coefficient function
//  for F2 at O(as^3) at next-to-leading log.
//------------------------------------------------------------------------------------------//

Value HighEnergyHighScaleCoefficientFunction::CL_ps3_highenergy_highscaleNLL(
    double m2Q2, double m2mu2, int nf
) const {

    return CF / CA * CL_g3_highenergy_highscaleNLL(m2Q2, m2mu2, nf);
}

//==========================================================================================//
//                  Color factors O(as^3)
//------------------------------------------------------------------------------------------//

// eq D.11-D.13 of our paper, normalized to as/4/pi
double AbstractHighEnergyCoefficientFunction::a_11() const {
    return 4*CA;
}

double AbstractHighEnergyCoefficientFunction::a_10(int nf) const {
    return -(11. * CA + 2. * nf * (1. - 2. * CF / CA)) / 3.;
}

double AbstractHighEnergyCoefficientFunction::a_21(int nf) const {
    return nf * (26. * CF - 23. * CA) / 9. * 4.;
}

double AbstractHighEnergyCoefficientFunction::a_21_new(int nf) const {
    return beta0(nf) * a_11() * (21. / 8 * zeta3 - 4 * ln2);
}



// eq D.20-D.25 of our paper, normalized to as/4/pi
double HighEnergyCoefficientFunction::Coff2_0(double m2Q2) const {
  double z = sqrt(1. / (1. + 4. * m2Q2));
  double L = log((1. + z) / (1. - z));
  double J = 4 * z * L;
  return (2 + (1-m2Q2)*J) /3.;
}
double HighEnergyCoefficientFunction::Coff2_1(double m2Q2) const {
  double z = sqrt(1. / (1. + 4. * m2Q2));
  double L = log((1. + z) / (1. - z));
  double Hmp = H_11(z) + H_1m1(z) - H_m11(z) - H_m1m1(z);
  double II = 4 * z * Hmp;
  double J  = 4 * z * L;
  return 2./3. * (5./3. + (1-m2Q2)/2.*II + (13-10*m2Q2)/12.*J);
}
double HighEnergyCoefficientFunction::Coff2_2(double m2Q2) const {
  double z = sqrt(1. / (1. + 4. * m2Q2));
  double L = log((1. + z) / (1. - z));
  double Hmp = H_11(z) + H_1m1(z) - H_m11(z) - H_m1m1(z);
  double Hmpm = H_111(z) - H_11m1(z) + H_1m11(z) - H_1m1m1(z) - H_m111(z) + H_m11m1(z) - H_m1m11(z) + H_m1m1m1(z);
  double II = 4 * z * Hmp;
  double J  = 4 * z * L;
  double K  = 4 * z * Hmpm;
  return 2./3. * (46./9. - (1-m2Q2)/4.*K + ((13-10*m2Q2)/12.-(1-m2Q2)/4.*log(4.*m2Q2/(1+4*m2Q2)))*II + (71-92*m2Q2)/36.*J);
}
double HighEnergyCoefficientFunction::CoffL_0(double m2Q2) const {
  double z = sqrt(1. / (1. + 4. * m2Q2));
  double L = log((1. + z) / (1. - z));
  double J = 4 * z * L;
  return 4./3. * (1+6*m2Q2 - (1+3*m2Q2)*m2Q2*J) / (1+4*m2Q2);
}
double HighEnergyCoefficientFunction::CoffL_1(double m2Q2) const {
  double z = sqrt(1. / (1. + 4. * m2Q2));
  double L = log((1. + z) / (1. - z));
  double Hmp = H_11(z) + H_1m1(z) - H_m11(z) - H_m1m1(z);
  double II = 4 * z * Hmp;
  double J  = 4 * z * L;
  return 2./3. * (-2./3.+8*m2Q2 - 2*(1+3*m2Q2)*m2Q2*II + (1./2.+2./3.*m2Q2-4*m2Q2*m2Q2)*J) / (1+4*m2Q2);
}
double HighEnergyCoefficientFunction::CoffL_2(double m2Q2) const {
  double z = sqrt(1. / (1. + 4. * m2Q2));
  double L = log((1. + z) / (1. - z));
  double Hmp = H_11(z) + H_1m1(z) - H_m11(z) - H_m1m1(z);
  double Hmpm = H_111(z) - H_11m1(z) + H_1m11(z) - H_1m1m1(z) - H_m111(z) + H_m11m1(z) - H_m1m11(z) + H_m1m1m1(z);
  double II = 4 * z * Hmp;
  double J  = 4 * z * L;
  double K  = 4 * z * Hmpm;
  return 2./3. * (68./9.+160./3.*m2Q2 + (1+3*m2Q2)*m2Q2*K + (1./2.+2./3.*m2Q2-4*m2Q2*m2Q2+(1+3*m2Q2)*m2Q2*log(4.*m2Q2/(1+4*m2Q2)))*II - (1./6.+68./9.*m2Q2+80./3.*m2Q2*m2Q2)*J) / (1+4*m2Q2);
}
double HighEnergyCoefficientFunction::Coff2_1(double m2Q2, double m2mu2) const {
  double Lmu = log(m2mu2);
  return Coff2_1(m2Q2) + Coff2_0(m2Q2)*Lmu;
}
double HighEnergyCoefficientFunction::Coff2_2(double m2Q2, double m2mu2) const {
  double Lmu = log(m2mu2);
  return Coff2_2(m2Q2) + Coff2_1(m2Q2)*Lmu + Coff2_0(m2Q2)*Lmu*Lmu/2.;
}
double HighEnergyCoefficientFunction::CoffL_1(double m2Q2, double m2mu2) const {
  double Lmu = log(m2mu2);
  return CoffL_1(m2Q2) + CoffL_0(m2Q2)*Lmu;
}
double HighEnergyCoefficientFunction::CoffL_2(double m2Q2, double m2mu2) const {
  double Lmu = log(m2mu2);
  return CoffL_2(m2Q2) + CoffL_1(m2Q2)*Lmu + CoffL_0(m2Q2)*Lmu*Lmu/2.;
}


// eq D.29-D.34 of our paper, normalized to as/4/pi
double HighEnergyHighScaleCoefficientFunction::Coff2_0(double m2Q2) const {
  return 2./3. * (1 - 2*log(m2Q2));
}
double HighEnergyHighScaleCoefficientFunction::Coff2_1(double m2Q2) const {
  double L = log(m2Q2);
  return 2./3. * (5./3. - 2*zeta2 - 13./3.*L + L*L);
}
double HighEnergyHighScaleCoefficientFunction::Coff2_2(double m2Q2) const {
  double L = log(m2Q2);
  return 2./3. * (46./9. - 13./3.*zeta2 + 4*zeta3 + (-71./9.+2*zeta2)*L + 13./6.*L*L - L*L*L/3.);
}
double HighEnergyHighScaleCoefficientFunction::CoffL_0() const {
  return 4./3.;
}
double HighEnergyHighScaleCoefficientFunction::CoffL_1(double m2Q2) const {
  return 2./3. * (-2./3. - 2*log(m2Q2));
}
double HighEnergyHighScaleCoefficientFunction::CoffL_2(double m2Q2) const {
  double L = log(m2Q2);
  return 2./3. * (68./9. - 2*zeta2 + 2./3.*L + L*L);
}
double HighEnergyHighScaleCoefficientFunction::Coff2_1(double m2Q2, double m2mu2) const {
  double Lmu = log(m2mu2);
  return Coff2_1(m2Q2) + Coff2_0(m2Q2)*Lmu;
}
double HighEnergyHighScaleCoefficientFunction::Coff2_2(double m2Q2, double m2mu2) const {
  double Lmu = log(m2mu2);
  return Coff2_2(m2Q2) + Coff2_1(m2Q2)*Lmu + Coff2_0(m2Q2)*Lmu*Lmu/2.;
}
double HighEnergyHighScaleCoefficientFunction::CoffL_1(double m2Q2, double m2mu2) const {
  double Lmu = log(m2mu2);
  return CoffL_1(m2Q2) + CoffL_0()*Lmu;
}
double HighEnergyHighScaleCoefficientFunction::CoffL_2(double m2Q2, double m2mu2) const {
  double Lmu = log(m2mu2);
  return CoffL_2(m2Q2) + CoffL_1(m2Q2)*Lmu + CoffL_0()*Lmu*Lmu/2.;
}
