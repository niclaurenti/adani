#include "adani/AsymptoticCoefficientFunction.h"

#include <cmath>

//==========================================================================================//
//  ExactCoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

AsymptoticCoefficientFunction::AsymptoticCoefficientFunction(
    const int &order, const char &kind, const char &channel, const bool &NLL,
    const bool &exact_highscale, const bool &revised_approx_highscale
)
    : AbstractHighEnergyCoefficientFunction(order, kind, channel, NLL) {

    highscale_ = new HighScaleCoefficientFunction(
        GetOrder(), GetKind(), GetChannel(), exact_highscale,
        revised_approx_highscale
    );
    powerterms_ = new PowerTermsCoefficientFunction(
        GetOrder(), GetKind(), GetChannel(), GetNLL()
    );
}

//==========================================================================================//
//  ExactCoefficientFunction: destructor
//------------------------------------------------------------------------------------------//

AsymptoticCoefficientFunction::~AsymptoticCoefficientFunction() {
    delete highscale_;
    delete powerterms_;
}

//==========================================================================================//
//  ExactCoefficientFunction: central value of fx
//------------------------------------------------------------------------------------------//

double AsymptoticCoefficientFunction::fx(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return highscale_->fx(x, m2Q2, m2mu2, nf)
           + powerterms_->fx(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  ExactCoefficientFunction: band of fx
//------------------------------------------------------------------------------------------//

Value AsymptoticCoefficientFunction::fxBand(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return highscale_->fxBand(x, m2Q2, m2mu2, nf)
           + powerterms_->fxBand(x, m2Q2, m2mu2, nf);
}
