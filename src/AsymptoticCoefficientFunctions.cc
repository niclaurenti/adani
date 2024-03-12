#include "adani/AsymptoticCoefficientFunctions.h"
#include "adani/HighEnergyCoefficientFunctions.h"
#include "adani/HighScaleCoefficientFunctions.h"
#include <cmath>

AsymptoticCoefficientFunction::AsymptoticCoefficientFunction(const int& order, const char& kind, const char& channel, const bool& NLL) : AbstractHighEnergyCoefficientFunction(order, kind, channel, NLL) {

    highscale_ = new HighScaleCoefficientFunction(GetOrder(), GetKind(), GetChannel());
    powerterms_ = new PowerTermsCoefficientFunction(GetOrder(), GetKind(), GetChannel(), GetNLL());
}

AsymptoticCoefficientFunction::~AsymptoticCoefficientFunction() {
    delete highscale_;
    delete powerterms_;
}

double AsymptoticCoefficientFunction::fx(double x, double m2Q2, double m2mu2, int nf) const {

    return highscale_->fx(x, m2Q2, m2mu2, nf) + powerterms_->fx(x, m2Q2, m2mu2, nf) ;
}
