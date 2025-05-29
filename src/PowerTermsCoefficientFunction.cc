#include "adani/PowerTermsCoefficientFunction.h"

//==========================================================================================//
//  AbstractPowerTerms: constructor
//------------------------------------------------------------------------------------------//

AbstractPowerTerms::AbstractPowerTerms(
    const int &order, const char &kind, const char &channel, const bool &NLL
)
    : AbstractHighEnergyCoefficientFunction(order, kind, channel, NLL) {
    highenergy_ = new HighEnergyCoefficientFunction(
        GetOrder(), GetKind(), GetChannel(), GetNLL()
    );
    highenergyhighscale_ = new HighEnergyHighScaleCoefficientFunction(
        GetOrder(), GetKind(), GetChannel(), GetNLL()
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

// In this way the error is enormous

// Value PowerTermsCoefficientFunction::fxBand(
//     double x, double m2Q2, double m2mu2, int nf
// ) const {

//     double central = (highenergy_->fx(x, m2Q2, m2mu2, nf))
//                      - (highenergyhighscale_->fx(x, m2Q2, m2mu2, nf));

//     Value tmp1 = highenergy_->fxBand(x, m2Q2, m2mu2, nf);
//     Value tmp2 = highenergyhighscale_->fxBand(x, m2Q2, m2mu2, nf);

//     double delta_he_up = tmp1.GetHigher() - tmp1.GetCentral();
//     double delta_he_down = tmp1.GetCentral() - tmp1.GetLower();

//     double delta_hehs_up = tmp2.GetHigher() - tmp2.GetCentral();
//     double delta_hehs_down = tmp2.GetCentral() - tmp2.GetLower();

//     double err_up = sqrt(delta_he_up*delta_he_up +
//     delta_hehs_up*delta_hehs_up); double err_down =
//     sqrt(delta_he_down*delta_he_down + delta_hehs_down*delta_hehs_down);

//     return Value(central, central + err_up, central - err_down);
// }

// In this way the error is enormous

// Value PowerTermsCoefficientFunction::fxBand(
//     double x, double m2Q2, double m2mu2, int nf
// ) const {

//     double central = (highenergy_->fx(x, m2Q2, m2mu2, nf))
//                      - (highenergyhighscale_->fx(x, m2Q2, m2mu2, nf));

//     vector<double> tmp1 = highenergy_->fxBand(x, m2Q2, m2mu2, nf).ToVect();
//     vector<double> tmp2 = highenergyhighscale_->fxBand(x, m2Q2, m2mu2,
//     nf).ToVect();

//     double tmp, higher = central, lower = central;
//     for(double he : tmp1) {
//         for (double hehs : tmp2) {
//             tmp = he - hehs;

//             if(tmp > higher) higher = tmp;
//             if(tmp < lower) lower = tmp;
//         }
//     }

//     return Value(central, higher, lower);
// }

//==========================================================================================//
//  MultiplicativeAsymptotic: band of the power terms
//------------------------------------------------------------------------------------------//

Value MultiplicativeAsymptotic::fxBand(
    double x, double m2Q2, double m2mu2, int nf
) const {
    return GetHighEnergy()->fxBand(x, m2Q2, m2mu2, nf)
           / GetHighEnergyHighScale()->fxBand(x, m2Q2, m2mu2, nf);
}
