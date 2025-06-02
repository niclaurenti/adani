#include "adani/PowerTermsCoefficientFunction.h"
#include <iostream>

//==========================================================================================//
//  AbstractPowerTerms: constructor
//------------------------------------------------------------------------------------------//

AbstractPowerTerms::AbstractPowerTerms(
    const int &order, const char &kind, const char &channel, const bool &NLL
)
    : CoefficientFunction(order, kind, channel) {
    highenergy_ = new HighEnergyCoefficientFunction(
        GetOrder(), GetKind(), GetChannel(), NLL
    );
    highenergyhighscale_ = new HighEnergyHighScaleCoefficientFunction(
        GetOrder(), GetKind(), GetChannel(), NLL
    );
}

//==========================================================================================//
//  AbstractPowerTerms: destructor
//------------------------------------------------------------------------------------------//

AbstractPowerTerms::~AbstractPowerTerms() {
    delete highenergy_;
    delete highenergyhighscale_;
}

//==========================================================================================//
//  PowerTermsCoefficientFunction: band of the power terms
//------------------------------------------------------------------------------------------//

Value PowerTermsCoefficientFunction::fxBand(
    double x, double m2Q2, double m2mu2, int nf
) const {
    return GetHighEnergy()->fxBand(x, m2Q2, m2mu2, nf)
           - GetHighEnergyHighScale()->fxBand(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  MultiplicativeAsymptotic: constructor
//------------------------------------------------------------------------------------------//

MultiplicativeAsymptotic::MultiplicativeAsymptotic(
    const int &order, const char &kind, const char &channel,
    const bool &NLL
) : AbstractPowerTerms(order, kind, channel, NLL) {
    if (order == 1) {
        fx_ = &MultiplicativeAsymptotic::OneFunction;
    } else if (order == 2) {
        fx_ = &MultiplicativeAsymptotic::PlainRatio;
    } else if (order == 3) {
        fx_ = &MultiplicativeAsymptotic::RegoularizedRatio;
    }
};

//==========================================================================================//
//  MultiplicativeAsymptotic: band of the power terms
//------------------------------------------------------------------------------------------//

Value MultiplicativeAsymptotic::fxBand(
    double x, double m2Q2, double m2mu2, int nf
) const {
    return (this->*fx_)(x, m2Q2, m2mu2, nf);;
}

//==========================================================================================//
//  MultiplicativeAsymptotic: band of the power terms
//------------------------------------------------------------------------------------------//

Value MultiplicativeAsymptotic::PlainRatio(
    double x, double m2Q2, double m2mu2, int nf
) const {
    return GetHighEnergy()->fxBand(x, m2Q2, m2mu2, nf)
           / GetHighEnergyHighScale()->fxBand(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  MultiplicativeAsymptotic: band of the power terms
//------------------------------------------------------------------------------------------//

Value MultiplicativeAsymptotic::RegoularizedRatio(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return Value(
        GetHighEnergy()->LL(m2Q2, m2mu2)
        / GetHighEnergyHighScale()->LL(m2Q2, m2mu2)
    );

    // return GetHighEnergy()->fxBand(x, m2Q2, m2mu2, nf)
    //        / (GetHighEnergyHighScale()->LL(m2Q2, m2mu2) * log(x))
    //        * (1 - GetHighEnergyHighScale()->NLL(m2Q2, m2mu2, nf).GetCentral() / (GetHighEnergyHighScale()->LL(m2Q2, m2mu2) * log(x))) ;
    // return GetHighEnergy()->fxBand(x, m2Q2, m2mu2, nf)
    //        / (GetHighEnergyHighScale()->LL(m2Q2, m2mu2) * log(x)) ;

}
