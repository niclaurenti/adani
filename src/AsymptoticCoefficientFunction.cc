#include "adani/AsymptoticCoefficientFunction.h"

#include <cmath>

//==========================================================================================//
//  AsymptoticCoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

AsymptoticCoefficientFunction::AsymptoticCoefficientFunction(
    const int &order, const char &kind, const char &channel, const bool &NLL,
    const HighScaleVersion &highscale_version
)
    : CoefficientFunction(order, kind, channel), legacy_pt_(false) {

    try {
        highscale_ = std::make_unique<HighScaleCoefficientFunction>(
            GetOrder(), GetKind(), GetChannel(), highscale_version
        );
        highenergy_ = std::make_unique<HighEnergyCoefficientFunction>(
            GetOrder(), GetKind(), GetChannel(), NLL
        );
        highenergyhighscale_ =
            std::make_unique<HighEnergyHighScaleCoefficientFunction>(
                GetOrder(), GetKind(), GetChannel(), NLL
            );

        SetFunctions();
    } catch (const UnexpectedException &e) {
        e.runtime_error();
    }
}

//==========================================================================================//
//  AsymptoticCoefficientFunction: copy constructor
//------------------------------------------------------------------------------------------//

AsymptoticCoefficientFunction::AsymptoticCoefficientFunction(
    const AsymptoticCoefficientFunction &obj
)
    : AsymptoticCoefficientFunction(
          obj.GetOrder(), obj.GetKind(), obj.GetChannel(), obj.GetNLL(),
          obj.GetHighScaleVersion()
      ) {}

//==========================================================================================//
//  AsymptoticCoefficientFunction: SetFunctions
//------------------------------------------------------------------------------------------//

void AsymptoticCoefficientFunction::SetFunctions() {
    switch (GetOrder()) {
    case 1:
        fx_ = &AsymptoticCoefficientFunction::PlainAdditiveMatching;
        break;
    case 2:
        switch (GetKind()) {
        case '2':
            fx_ = &AsymptoticCoefficientFunction::C2_2_asymptotic;
            break;
        case 'L':
            fx_ = &AsymptoticCoefficientFunction::CL_2_asymptotic;
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
            fx_ = &AsymptoticCoefficientFunction::C2_3_asymptotic;
            break;
        case 'L':
            fx_ = &AsymptoticCoefficientFunction::CL_3_asymptotic;
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
//  AsymptoticCoefficientFunction: restore legacy behavior for power terms
//------------------------------------------------------------------------------------------//

void AsymptoticCoefficientFunction::SetLegacyPowerTerms(const bool &legacy_pt) {
    try {
        if (legacy_pt == legacy_pt_) {
            throw NotValidException(
                "Setting legacy power terms identical to its previous value!",
                __PRETTY_FUNCTION__, __LINE__
            );
        }
        legacy_pt_ = legacy_pt;
        switch (GetOrder()) {
        case 1:
            throw NotPresentException(
                "For order='1' legacy power terms are identical to the "
                "current ones!",
                __PRETTY_FUNCTION__, __LINE__
            );
            break;
        default:
            if (legacy_pt) {
                fx_ = &AsymptoticCoefficientFunction::PlainAdditiveMatching;
            } else {
                SetFunctions();
            }
        }
    } catch (const NotPresentException &e) {
        e.warning();
    } catch (const NotValidException &e) {
        e.warning();
    } catch (UnexpectedException &e) {
        e.runtime_error();
    }
}

//==========================================================================================//
//  AsymptoticCoefficientFunction: band of fx
//------------------------------------------------------------------------------------------//

Value AsymptoticCoefficientFunction::fxBand(
    double x, double m2Q2, double m2mu2, int nf
) const {
    return (this->*fx_)(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  AsymptoticCoefficientFunction: band of fx with additive matching
//------------------------------------------------------------------------------------------//

Value AsymptoticCoefficientFunction::PlainAdditiveMatching(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return highscale_->fxBand(x, m2Q2, m2mu2, nf)
           + (highenergy_->fxBand(x, m2Q2, m2mu2, nf)
              - highenergyhighscale_->fxBand(x, m2Q2, m2mu2, nf));
}

//==========================================================================================//
//  AsymptoticCoefficientFunction: band of fx with pure LL multiplicative
//  matching
//------------------------------------------------------------------------------------------//

Value AsymptoticCoefficientFunction::PlainMultiplicativeMatching(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return highscale_->fxBand(x, m2Q2, m2mu2, nf) * highenergy_->LL(m2Q2, m2mu2)
           / highenergyhighscale_->LL(m2Q2, m2mu2);
}

//==========================================================================================//
//  AsymptoticCoefficientFunction: band of fx with multiplicative matching at
//  NLL (version 1)
//------------------------------------------------------------------------------------------//

Value AsymptoticCoefficientFunction::ModifiedMultiplicativeMatching1(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return (highscale_->fxBand(x, m2Q2, m2mu2, nf)
            + (highenergy_->NLL(m2Q2, m2mu2, nf)
               - highenergyhighscale_->NLL(m2Q2, m2mu2, nf))
                  / x)
           * highenergy_->LL(m2Q2, m2mu2)
           / highenergyhighscale_->LL(m2Q2, m2mu2);
}

//==========================================================================================//
//  AsymptoticCoefficientFunction: band of fx with multiplicative matching at
//  NLL (version 2)
//------------------------------------------------------------------------------------------//

Value AsymptoticCoefficientFunction::ModifiedMultiplicativeMatching2(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return highscale_->fxBand(x, m2Q2, m2mu2, nf) * highenergy_->LL(m2Q2, m2mu2)
               / highenergyhighscale_->LL(m2Q2, m2mu2)
           + (highenergy_->NLL(m2Q2, m2mu2, nf)
              - highenergyhighscale_->NLL(m2Q2, m2mu2, nf))
                 / x;
}

//==========================================================================================//
//  AsymptoticCoefficientFunction: asymptotic coefficient function for C2 at
//  O(as^2)
//------------------------------------------------------------------------------------------//

Value AsymptoticCoefficientFunction::C2_2_asymptotic(
    double x, double m2Q2, double m2mu2, int nf
) const {

    Value central = PlainAdditiveMatching(x, m2Q2, m2mu2, nf);
    Value variation = PlainMultiplicativeMatching(x, m2Q2, m2mu2, nf);

    return Delta2(central, variation, m2Q2, m2mu2);
}

//==========================================================================================//
//  AsymptoticCoefficientFunction: asymptotic coefficient function for CL at
//  O(as^2)
//------------------------------------------------------------------------------------------//

Value AsymptoticCoefficientFunction::CL_2_asymptotic(
    double x, double m2Q2, double m2mu2, int nf
) const {

    Value central = PlainMultiplicativeMatching(x, m2Q2, m2mu2, nf);
    Value variation = PlainAdditiveMatching(x, m2Q2, m2mu2, nf);

    return Delta2(central, variation, m2Q2, m2mu2);
}

//==========================================================================================//
//  AsymptoticCoefficientFunction: asymptotic coefficient function for C2 at
//  O(as^3)
//------------------------------------------------------------------------------------------//

Value AsymptoticCoefficientFunction::C2_3_asymptotic(
    double x, double m2Q2, double m2mu2, int nf
) const {

    Value central = PlainAdditiveMatching(x, m2Q2, m2mu2, nf);
    Value variation1 = ModifiedMultiplicativeMatching1(x, m2Q2, m2mu2, nf);
    Value variation2 = ModifiedMultiplicativeMatching2(x, m2Q2, m2mu2, nf);

    return Delta3(central, variation1, variation2, m2Q2, m2mu2);
}

//==========================================================================================//
//  AsymptoticCoefficientFunction: asymptotic coefficient function for CL at
//  O(as^3)
//------------------------------------------------------------------------------------------//

Value AsymptoticCoefficientFunction::CL_3_asymptotic(
    double x, double m2Q2, double m2mu2, int nf
) const {

    Value central = ModifiedMultiplicativeMatching1(x, m2Q2, m2mu2, nf);
    Value variation1 = PlainAdditiveMatching(x, m2Q2, m2mu2, nf);
    Value variation2 = ModifiedMultiplicativeMatching2(x, m2Q2, m2mu2, nf);

    return Delta3(central, variation1, variation2, m2Q2, m2mu2);
}

//==========================================================================================//
//  AsymptoticCoefficientFunction: compute error of asymptotic cefficient
//  function starting from the central value and the two variations
//------------------------------------------------------------------------------------------//

Value AsymptoticCoefficientFunction::Delta2(
    Value central, Value variation, double m2Q2, double m2mu2
) const {
    double central_c = central.GetHigher();
    double var_c = variation.GetHigher();

    // double central_delta = central.GetAvgDelta();
    // double var_delta = variation.GetAvgDelta();
    double damp = 1.
                  / (1.
                     + 2
                           * std::abs(
                               log(highenergy_->LL(m2Q2, m2mu2)
                                   / highenergyhighscale_->LL(m2Q2, m2mu2))
                           ));
    double delta = damp * std::abs(central_c - var_c);
    // double delta = damp * sqrt(tmp * tmp + central_delta * central_delta +
    // var_delta * var_delta);

    return Value(central_c, central_c + delta, central_c - delta);
}

//==========================================================================================//
//  AsymptoticCoefficientFunction: compute error of asymptotic cefficient
//  function starting from the central value and the two variations
//------------------------------------------------------------------------------------------//

Value AsymptoticCoefficientFunction::Delta3(
    Value central, Value variation1, Value variation2, double m2Q2, double m2mu2
) const {
    double central_c = central.GetHigher();
    double var1_c = variation1.GetHigher();
    double var2_c = variation2.GetHigher();

    // double central_delta = central.GetAvgDelta();
    // double var1_delta = variation1.GetAvgDelta();
    // double var2_delta = variation2.GetAvgDelta();

    double damp = 1.
                  / (1.
                     + 2
                           * std::abs(
                               log(highenergy_->LL(m2Q2, m2mu2)
                                   / highenergyhighscale_->LL(m2Q2, m2mu2))
                           ));

    double tmp1 = central_c - var1_c;
    double tmp2 = central_c - var2_c;

    double delta = damp * sqrt(tmp1 * tmp1 + tmp2 * tmp2);
    // double delta = damp * sqrt(tmp1 * tmp1 + tmp2 * tmp2 + central_delta *
    // central_delta + var1_delta * var1_delta + var2_delta * var2_delta);

    return Value(central_c, central_c + delta, central_c - delta);
}
