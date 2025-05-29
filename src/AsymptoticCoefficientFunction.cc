#include "adani/AsymptoticCoefficientFunction.h"

#include <cmath>

//==========================================================================================//
//  AsymptoticCoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

AsymptoticCoefficientFunction::AsymptoticCoefficientFunction(
    const int &order, const char &kind, const char &channel, const bool &NLL,
    const string &highscale_version
)
    : AbstractHighEnergyCoefficientFunction(order, kind, channel, NLL) {

    highscale_ = new HighScaleCoefficientFunction(
        GetOrder(), GetKind(), GetChannel(), highscale_version
    );
    if (GetKind() == '2') {
        powerterms_ = new PowerTermsCoefficientFunction(
            GetOrder(), GetKind(), GetChannel(), GetNLL()
        );
        fx_ = &AsymptoticCoefficientFunction::AdditiveMatching;
    } else {
        powerterms_ = new MultiplicativeAsymptotic(
            GetOrder(), GetKind(), GetChannel(), GetNLL()
        );
        fx_ = &AsymptoticCoefficientFunction::MultiplicativeMatching;
    }
}

//==========================================================================================//
//  AsymptoticCoefficientFunction: destructor
//------------------------------------------------------------------------------------------//

AsymptoticCoefficientFunction::~AsymptoticCoefficientFunction() {
    delete highscale_;
    delete powerterms_;
}

//==========================================================================================//
//  AsymptoticCoefficientFunction: restore legacy behavior for power terms
//------------------------------------------------------------------------------------------//

void AsymptoticCoefficientFunction::SetLegacyPowerTerms(const bool &legacy_pt) {
    try {
        if (legacy_pt) {
            if (GetKind() == '2') {
                throw NotPresentException(
                    "For kind='2' legacy power terms are identical to the current ones!",
                    __PRETTY_FUNCTION__, __LINE__
                );
            } else {
                if (legacy_pt) {
                    delete powerterms_;
                    powerterms_ = new PowerTermsCoefficientFunction(
                        GetOrder(), GetKind(), GetChannel(), GetNLL()
                    );
                    fx_ = &AsymptoticCoefficientFunction::AdditiveMatching;
                } else {
                    delete powerterms_;
                    powerterms_ = new MultiplicativeAsymptotic(
                        GetOrder(), GetKind(), GetChannel(), GetNLL()
                    );
                    fx_ = &AsymptoticCoefficientFunction::MultiplicativeMatching;
                }
            }
        }
    } catch (const NotPresentException &e) {
        e.warning();
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

Value AsymptoticCoefficientFunction::AdditiveMatching(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return highscale_->fxBand(x, m2Q2, m2mu2, nf)
           + powerterms_->fxBand(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  AsymptoticCoefficientFunction: band of fx with multiplicative matching
//------------------------------------------------------------------------------------------//

Value AsymptoticCoefficientFunction::MultiplicativeMatching(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return highscale_->fxBand(x, m2Q2, m2mu2, nf)
           * powerterms_->fxBand(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  AsymptoticCoefficientFunction: all possible variation (3x3=9) of the
//  combination between high scale and power terms
//------------------------------------------------------------------------------------------//

vector<double> AsymptoticCoefficientFunction::AllVariations(
    double x, double m2Q2, double m2mu2, int nf
) const {

    vector<double> hs_vec = (highscale_->fxBand(x, m2Q2, m2mu2, nf)).ToVect();
    vector<double> pt_vec = (powerterms_->fxBand(x, m2Q2, m2mu2, nf)).ToVect();

    vector<double> res;

    for (double hs : hs_vec) {
        for (double pt : pt_vec) {
            res.push_back(hs + pt);
        }
    }

    return res;
}
