#include "adani/HighEnergyCoefficientFunction.h"
#include "adani/Constants.h"
#include "adani/SpecialFunctions.h"
#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

//==========================================================================================//
//  AbstractHighEnergyCoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

AbstractHighEnergyCoefficientFunction::AbstractHighEnergyCoefficientFunction(
    const int &order, const char &kind, const char &channel, const bool &NLL
)
    : CoefficientFunction(order, kind, channel) {

    SetNLL(NLL);
}

//==========================================================================================//
//  HighEnergyCoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

HighEnergyCoefficientFunction::HighEnergyCoefficientFunction(
    const int &order, const char &kind, const char &channel, const bool &NLL
)
    : AbstractHighEnergyCoefficientFunction(order, kind, channel, NLL) {
    SetFunctions();
}

//==========================================================================================//
//  HighEnergyCoefficientFunction: band of the high energy coefficient function
//------------------------------------------------------------------------------------------//

Value HighEnergyCoefficientFunction::fxBand(
    double x, double m2Q2, double m2mu2, int nf
) const {
    return (this->*LL_)(x, m2Q2, m2mu2) + (this->*NLL_)(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  HighEnergyCoefficientFunction: function that sets the pointer for LL_ and
//  NLL_
//------------------------------------------------------------------------------------------//

void HighEnergyCoefficientFunction::SetFunctions() {

    if (GetOrder() == 1) {
        LL_ = &HighEnergyCoefficientFunction::ZeroFunction;
        NLL_ = &HighEnergyCoefficientFunction::ZeroFunctionBand;
    } else if (GetOrder() == 2) {

        NLL_ = &HighEnergyCoefficientFunction::ZeroFunctionBand;
        if (GetKind() == '2' && GetChannel() == 'g')
            LL_ = &HighEnergyCoefficientFunction::C2_g2_highenergy;
        else if (GetKind() == '2' && GetChannel() == 'q')
            LL_ = &HighEnergyCoefficientFunction::C2_ps2_highenergy;
        else if (GetKind() == 'L' && GetChannel() == 'g')
            LL_ = &HighEnergyCoefficientFunction::CL_g2_highenergy;
        else if (GetKind() == 'L' && GetChannel() == 'q')
            LL_ = &HighEnergyCoefficientFunction::CL_ps2_highenergy;

    } else if (GetOrder() == 3) {
        if (GetKind() == '2' && GetChannel() == 'g') {
            LL_ = &HighEnergyCoefficientFunction::C2_g3_highenergyLL;
            NLL_ = &HighEnergyCoefficientFunction::C2_g3_highenergyNLL;
        } else if (GetKind() == '2' && GetChannel() == 'q') {
            LL_ = &HighEnergyCoefficientFunction::C2_ps3_highenergyLL;
            NLL_ = &HighEnergyCoefficientFunction::C2_ps3_highenergyNLL;
        } else if (GetKind() == 'L' && GetChannel() == 'g') {
            LL_ = &HighEnergyCoefficientFunction::CL_g3_highenergyLL;
            NLL_ = &HighEnergyCoefficientFunction::CL_g3_highenergyNLL;
        } else if (GetKind() == 'L' && GetChannel() == 'q') {
            LL_ = &HighEnergyCoefficientFunction::CL_ps3_highenergyLL;
            NLL_ = &HighEnergyCoefficientFunction::CL_ps3_highenergyNLL;
        } else {
            cout << "Error: something has gone wrong in "
                    "HighEnergyCoefficientFunction::SetFunctions!"
                 << endl;
            exit(-1);
        }
    } else {
        cout << "Error: something has bone wrong in "
                "HighEnergyCoefficientFunction::SetFunctions!"
             << endl;
        exit(-1);
    }
    if (!GetNLL())
        NLL_ = &HighEnergyCoefficientFunction::ZeroFunctionBand;
}

//==========================================================================================//
//  HighEnergyHighScaleCoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

HighEnergyHighScaleCoefficientFunction::HighEnergyHighScaleCoefficientFunction(
    const int &order, const char &kind, const char &channel, const bool &NLL
)
    : AbstractHighEnergyCoefficientFunction(order, kind, channel, NLL) {
    SetFunctions();
}

//==========================================================================================//
//  HighEnergyHighScaleCoefficientFunction: band of the high scale limit of the
//  high energy
//------------------------------------------------------------------------------------------//

Value HighEnergyHighScaleCoefficientFunction::fxBand(
    double x, double m2Q2, double m2mu2, int nf
) const {
    return (this->*LL_)(x, m2Q2, m2mu2) + (this->*NLL_)(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  HighEnergyHighScaleCoefficientFunction: function that sets the pointer for
//  LL_ and NLL_
//------------------------------------------------------------------------------------------//

void HighEnergyHighScaleCoefficientFunction::SetFunctions() {

    if (GetOrder() == 1) {
        LL_ = &HighEnergyHighScaleCoefficientFunction::ZeroFunction;
        NLL_ = &HighEnergyHighScaleCoefficientFunction::ZeroFunctionBand;
    } else if (GetOrder() == 2) {

        NLL_ = &HighEnergyHighScaleCoefficientFunction::ZeroFunctionBand;
        if (GetKind() == '2' && GetChannel() == 'g')
            LL_ = &HighEnergyHighScaleCoefficientFunction::
                      C2_g2_highenergy_highscale;
        else if (GetKind() == '2' && GetChannel() == 'q')
            LL_ = &HighEnergyHighScaleCoefficientFunction::
                      C2_ps2_highenergy_highscale;
        else if (GetKind() == 'L' && GetChannel() == 'g')
            LL_ = &HighEnergyHighScaleCoefficientFunction::
                      CL_g2_highenergy_highscale;
        else if (GetKind() == 'L' && GetChannel() == 'q')
            LL_ = &HighEnergyHighScaleCoefficientFunction::
                      CL_ps2_highenergy_highscale;

    } else if (GetOrder() == 3) {
        if (GetKind() == '2' && GetChannel() == 'g') {
            LL_ = &HighEnergyHighScaleCoefficientFunction::
                      C2_g3_highenergy_highscaleLL;
            NLL_ = &HighEnergyHighScaleCoefficientFunction::
                       C2_g3_highenergy_highscaleNLL;
        } else if (GetKind() == '2' && GetChannel() == 'q') {
            LL_ = &HighEnergyHighScaleCoefficientFunction::
                      C2_ps3_highenergy_highscaleLL;
            NLL_ = &HighEnergyHighScaleCoefficientFunction::
                       C2_ps3_highenergy_highscaleNLL;
        } else if (GetKind() == 'L' && GetChannel() == 'g') {
            LL_ = &HighEnergyHighScaleCoefficientFunction::
                      CL_g3_highenergy_highscaleLL;
            NLL_ = &HighEnergyHighScaleCoefficientFunction::
                       CL_g3_highenergy_highscaleNLL;
        } else if (GetKind() == 'L' && GetChannel() == 'q') {
            LL_ = &HighEnergyHighScaleCoefficientFunction::
                      CL_ps3_highenergy_highscaleLL;
            NLL_ = &HighEnergyHighScaleCoefficientFunction::
                       CL_ps3_highenergy_highscaleNLL;
        } else {
            cout << "Error: something has gone wrong in "
                    "HighEnergyHighScaleCoefficientFunction::SetFunctions!"
                 << endl;
            exit(-1);
        }
    } else {
        cout << "Error: something has bone wrong in "
                "HighEnergyHighScaleCoefficientFunction::SetFunctions!"
             << endl;
        exit(-1);
    }
    if (!GetNLL())
        NLL_ = &HighEnergyHighScaleCoefficientFunction::ZeroFunctionBand;
}

//==========================================================================================//
//  PowerTermsCoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

PowerTermsCoefficientFunction::PowerTermsCoefficientFunction(
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
//  PowerTermsCoefficientFunction: destructor
//------------------------------------------------------------------------------------------//

PowerTermsCoefficientFunction::~PowerTermsCoefficientFunction() {
    delete highenergy_;
    delete highenergyhighscale_;
}

//==========================================================================================//
//  PowerTermsCoefficientFunction: band of the power terms
//------------------------------------------------------------------------------------------//

Value PowerTermsCoefficientFunction::fxBand(
    double x, double m2Q2, double m2mu2, int nf
) const {
    // TODO: in this way the error is very small: should take all the
    // combinations?
    double central = (highenergy_->fx(x, m2Q2, m2mu2, nf))
                     - (highenergyhighscale_->fx(x, m2Q2, m2mu2, nf));
    double higher =
        (highenergy_->fxBand(x, m2Q2, m2mu2, nf)).GetHigher()
        - (highenergyhighscale_->fxBand(x, m2Q2, m2mu2, nf)).GetHigher();
    double lower =
        (highenergy_->fxBand(x, m2Q2, m2mu2, nf)).GetLower()
        - (highenergyhighscale_->fxBand(x, m2Q2, m2mu2, nf)).GetLower();

    if (higher > lower)
        return Value(central, higher, lower);
    else
        return Value(central, lower, higher);
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
//  High energy limit of the gluon coefficient function for F2 at O(as^2).
//
//  Eq. (3.38) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double HighEnergyCoefficientFunction::C2_g2_highenergy(
    double x, double m2Q2, double m2mu2
) const {

    double z = sqrt(1. / (1. + 4. * m2Q2));

    double L = log((1. + z) / (1. - z));

    double Hmp = H_11(z) + H_1m1(z) - H_m11(z) - H_m1m1(z);

    double I = 4. * z * Hmp;
    double J = 4. * z * L;

    double Lmu = log(m2mu2);

    double c_const =
        10. / 3. + (1. - m2Q2) * I + (13. / 6. - 5. / 3. * m2Q2) * J;

    double c_Lmu = 2. + (1. - m2Q2) * J;

    return 4. / 3. * CA * (c_const + c_Lmu * Lmu) / x;
}

//==========================================================================================//
//  High energy limit of the quark coefficient function for F2 at O(as^2).
//------------------------------------------------------------------------------------------//

double HighEnergyCoefficientFunction::C2_ps2_highenergy(
    double x, double m2Q2, double m2mu2
) const {

    return CF / CA * C2_g2_highenergy(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for FL at O(as^2).
//
//  Eq. (16, 21) of Ref. [arXiv:hep-ph/9411431]
//------------------------------------------------------------------------------------------//

double HighEnergyCoefficientFunction::CL_g2_highenergy(
    double x, double m2Q2, double m2mu2
) const {

    double z = sqrt(1. / (1. + 4. * m2Q2));

    double L = log((1. + z) / (1. - z));

    double Hmp = H_11(z) + H_1m1(z) - H_m11(z) - H_m1m1(z);

    double I = 4. * z * Hmp;
    double J = 4. * z * L;

    double Lmu = log(1. / m2mu2);

    double c_const = (4. * (-1. + 12. * m2Q2) / (3. + 12. * m2Q2)
                      + (5. - 12. * m2Q2 + 1. / (1. + 4. * m2Q2)) * J / 6.
                      - 4. * m2Q2 * (1. + 3. * m2Q2) / (1. + 4. * m2Q2) * I)
                     / 3.;

    double c_Lmu = (-4. * (1. + 6. * m2Q2) / (1. + 4. * m2Q2)
                    + 4. * m2Q2 * (1. + 3. * m2Q2) / (1. + 4. * m2Q2) * J)
                   / 3.;

    return 8 * CA * TR * (c_const + c_Lmu * Lmu) / x;
}

//==========================================================================================//
//  High energy limit of the quark coefficient function for FL at O(as^2).
//------------------------------------------------------------------------------------------//

double HighEnergyCoefficientFunction::CL_ps2_highenergy(
    double x, double m2Q2, double m2mu2
) const {

    return CF / CA * CL_g2_highenergy(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  Q^2>>m^2 limit of the high energy limit of the gluon coefficient function
//  for F2 at O(as^2).
//
//  Eq. (3.40) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double HighEnergyHighScaleCoefficientFunction::C2_g2_highenergy_highscale(
    double x, double m2Q2, double m2mu2
) const {

    double LQ = log(1. / m2Q2);
    double L2Q = LQ * LQ;

    double Lm = log(m2mu2);

    return CA
           * (8. / 3. * L2Q + 104. / 9. * LQ + 40. / 9. - 16. / 3. * zeta2
              + (16. / 3. * LQ + 8. / 3.) * Lm)
           / x;
}

//==========================================================================================//
//  Q^2>>m^2 limit of the high energy limit of the quark coefficient function
//  for F2 at O(as^2).
//------------------------------------------------------------------------------------------//

double HighEnergyHighScaleCoefficientFunction::C2_ps2_highenergy_highscale(
    double x, double m2Q2, double m2mu2
) const {

    return CF / CA * C2_g2_highenergy_highscale(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  Q^2>>m^2 limit of the high energy limit of the gluon coefficient function
//  for FL at O(as^2).
//
//  Q^2>>m^2 limit of Eq. (16, 21) of Ref. [arXiv:hep-ph/9411431]
//------------------------------------------------------------------------------------------//

double HighEnergyHighScaleCoefficientFunction::CL_g2_highenergy_highscale(
    double x, double m2Q2, double m2mu2
) const {

    double LQ = log(1. / m2Q2);

    double Lm = log(1. / m2mu2);

    double c_const = -2. / 9. * (1. - 3. * LQ);

    double c_log = -2. / 3.;

    return 16 * CA * TR / x * (c_const + c_log * Lm);
}

//==========================================================================================//
//  Q^2>>m^2 limit of the high energy limit of the quark coefficient function
//  for FL at O(as^2).
//------------------------------------------------------------------------------------------//

double HighEnergyHighScaleCoefficientFunction::CL_ps2_highenergy_highscale(
    double x, double m2Q2, double m2mu2
) const {

    return CF / CA * CL_g2_highenergy_highscale(x, m2Q2, m2mu2);
}

// //==========================================================================================//
// //  Power terms in the small x limit of the gluon coefficient function for F2
// at
// //  O(as^2).
// //------------------------------------------------------------------------------------------//

// double PowerTermsCoefficientFunction::C2_g2_power_terms(double x, double
// m2Q2, double m2mu2) {

//     highenergy = HighEnergyCoefficientFunction(order_, kind_, channel_, NLL_)

//     return C2_g2_highenergy(x, m2Q2, m2mu2)
//            - C2_g2_highenergy_highscale(x, m2Q2, m2mu2);
// }

// //==========================================================================================//
// //  Power terms in the small x limit of the quark coefficient function for F2
// at
// //  O(as^2).
// //------------------------------------------------------------------------------------------//

// double PowerTermsCoefficientFunction::C2_ps2_power_terms(double x, double
// m2Q2, double m2mu2) {

//     return CF / CA * C2_g2_power_terms(x, m2Q2, m2mu2);
// }

// //==========================================================================================//
// //  Power terms in the small x limit the gluon coefficient function for FL at
// //  O(as^2).
// //------------------------------------------------------------------------------------------//

// double PowerTermsCoefficientFunction::CL_g2_power_terms(double x, double
// m2Q2, double m2mu2) {

//     return CL_g2_highenergy(x, m2Q2, m2mu2)
//            - CL_g2_highenergy_highscale(x, m2Q2, m2mu2);
// }

// //==========================================================================================//
// //  Power terms in the small x limit the quark coefficient function for FL at
// //  O(as^2).
// //------------------------------------------------------------------------------------------//

// double PowerTermsCoefficientFunction::CL_ps2_power_terms(double x, double
// m2Q2, double m2mu2) {

//     return CF / CA * CL_g2_power_terms(x, m2Q2, m2mu2);
// }

//==========================================================================================//
//  High energy limit of the gluon coefficient function for F2 at O(as^3)
//  at leading log.
//
//  Eq. (3.39) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double HighEnergyCoefficientFunction::C2_g3_highenergyLL(
    double x, double m2Q2, double m2mu2
) const {

    double z = sqrt(1. / (1. + 4. * m2Q2));

    double L = log((1. + z) / (1. - z));

    double Hmp = H_11(z) + H_1m1(z) - H_m11(z) - H_m1m1(z);

    double Hmpm = H_111(z) - H_11m1(z) + H_1m11(z) - H_1m1m1(z) - H_m111(z)
                  + H_m11m1(z) - H_m1m11(z) + H_m1m1m1(z);

    double I = 4 * z * Hmp;
    double J = 4 * z * L;
    double K = 4 * z * Hmpm;

    double Logxi = log(1. + 1. / (4. * m2Q2));

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    return CA * CA
           * (-1472. / 27 - 8. / 3 * K * (-1. + m2Q2)
              + 8. / 27 * J * (-71. + 92. * m2Q2)
              + I
                    * (8. / 3 * Logxi * (-1. + m2Q2)
                       + 8. / 9 * (-13. + 10. * m2Q2))
              + (-160. / 9 + 16. / 3 * I * (-1. + m2Q2)
                 + 8. / 9 * J * (-13. + 10. * m2Q2))
                    * Lmu
              + (-16. / 3 + 8. / 3 * J * (-1. + m2Q2)) * Lmu2)
           * log(x) / x;
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for F2 at O(as^3)
//  at next to leading log.
//
//  Eq. (3.39) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

Value HighEnergyCoefficientFunction::C2_g3_highenergyNLL(
    double x, double m2Q2, double m2mu2, int nf
) const {

    double z = sqrt(1. / (1. + 4. * m2Q2));

    double L = log((1. + z) / (1. - z));

    double Hmp = H_11(z) + H_1m1(z) - H_m11(z) - H_m1m1(z);

    double Hmpm = H_111(z) - H_11m1(z) + H_1m11(z) - H_1m1m1(z) - H_m111(z)
                  + H_m11m1(z) - H_m1m11(z) + H_m1m1m1(z);

    double II = 4 * z * Hmp;
    double J = 4 * z * L;
    double K = 4 * z * Hmpm;

    double a11 = a_11();
    double a21 = a_21(nf);
    double a10 = a_10(nf);

    double beta0 = beta(0, nf);

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    double Logxi = log(1. + 1. / (4. * m2Q2));

    double central_value =
        (a21
             * (160. / 9 - 16. / 3 * II * (-1. + m2Q2)
                - 8. / 9 * J * (-13. + 10. * m2Q2))
         + a10 * a11
               * (2944. / 27 + 16. / 3 * K * (-1. + m2Q2)
                  - 16. / 27 * J * (-71. + 92. * m2Q2)
                  + II
                        * (-16. / 3 * Logxi * (-1. + m2Q2)
                           - 16. / 9 * (-13. + 10. * m2Q2)))
         + a11 * beta0
               * (-1472. / 27 - 8. / 3 * K * (-1. + m2Q2)
                  + 8. / 27 * J * (-71. + 92. * m2Q2)
                  + II
                        * (8. / 3 * Logxi * (-1. + m2Q2)
                           + 8. / 9 * (-13. + 10. * m2Q2)))
         + (a21 * (32. / 3 - 16. / 3 * J * (-1. + m2Q2))
            + a10 * a11
                  * (320. / 9 - 32. / 3 * II * (-1. + m2Q2)
                     - 16. / 9 * J * (-13. + 10. * m2Q2))
            + a11 * beta0
                  * (-160. / 9 + 16. / 3 * II * (-1. + m2Q2)
                     + 8. / 9 * J * (-13. + 10. * m2Q2)))
               * Lmu
         + (a10 * a11 * (32. / 3 - 16. / 3 * J * (-1. + m2Q2))
            + a11 * beta0 * (-16. / 3 + 8. / 3 * J * (-1. + m2Q2)))
               * Lmu2)
        / x;

    double error = (a10 * a11
                        * (2944. / 27 + 16. / 3 * K * (-1. + m2Q2)
                           - 16. / 27 * J * (-71 + 92 * m2Q2)
                           + II
                                 * (-16. / 3 * Logxi * (-1. + m2Q2)
                                    - 16. / 9 * (-13 + 10 * m2Q2)))
                    + (a10 * a11 * (32. / 3 - 16. / 3 * J * (-1. + m2Q2))
                       + a11 * beta0 * (-16. / 3 + 8. / 3 * J * (-1 + m2Q2)))
                          * Lmu2
                    + a11 * beta0
                          * (-1472. / 27 - 8. / 3 * K * (-1 + m2Q2)
                             - 640. * ln2 / 9 + 140. * zeta3 / 3
                             + II
                                   * (8. / 3 * Logxi * (-1. + m2Q2)
                                      + 8. / 9 * (-13 + 10 * m2Q2)
                                      + 64. / 3 * (-1 + m2Q2) * ln2
                                      - 14. * (-1 + m2Q2) * zeta3)
                             + J
                                   * (8. / 27 * (-71 + 92 * m2Q2)
                                      + 32. / 9 * (-13 + 10 * m2Q2) * ln2
                                      - 7. / 3 * (-13 + 10 * m2Q2) * zeta3))
                    + Lmu
                          * (a10 * a11
                                 * (320. / 9 - 32. / 3 * II * (-1 + m2Q2)
                                    - 16. / 9 * J * (-13 + 10 * m2Q2))
                             + a11 * beta0
                                   * (-160. / 9 + 16. / 3 * II * (-1 + m2Q2)
                                      - 128. * ln2 / 3 + 28 * zeta3
                                      + J
                                            * (8. / 9 * (-13 + 10 * m2Q2)
                                               + 64. / 3 * (-1 + m2Q2) * ln2
                                               - 14 * (-1 + m2Q2) * zeta3))))
                   / x;

    double delta = fabs(central_value - error);

    double higher = central_value + delta;
    double lower = central_value - delta;

    return Value(central_value, higher, lower);
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for F2 at O(as^3).
//
//  Eq. (3.39) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

Value HighEnergyCoefficientFunction::C2_g3_highenergy(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return C2_g3_highenergyLL(x, m2Q2, m2mu2)
           + C2_g3_highenergyNLL(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  High energy limit of the quark coefficient function for F2 at O(as^3)
//  at leading log.
//------------------------------------------------------------------------------------------//

double HighEnergyCoefficientFunction::C2_ps3_highenergyLL(
    double x, double m2Q2, double m2mu2
) const {

    return CF / CA * C2_g3_highenergyLL(x, m2Q2, m2mu2);
}

Value HighEnergyCoefficientFunction::C2_ps3_highenergyNLL(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return CF / CA * C2_g3_highenergyNLL(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  High energy limit of the quark coefficient function for F2 at O(as^3).
//------------------------------------------------------------------------------------------//

Value HighEnergyCoefficientFunction::C2_ps3_highenergy(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return CF / CA * C2_g3_highenergy(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for FL at O(as^3)
//  at leading log.
//------------------------------------------------------------------------------------------//

double HighEnergyCoefficientFunction::CL_g3_highenergyLL(
    double x, double m2Q2, double m2mu2
) const {

    double m4Q4 = m2Q2 * m2Q2;

    double z = sqrt(1. / (1. + 4. * m2Q2));

    double L = log((1. + z) / (1. - z));

    double Hmp = H_11(z) + H_1m1(z) - H_m11(z) - H_m1m1(z);

    double Hmpm = H_111(z) - H_11m1(z) + H_1m11(z) - H_1m1m1(z) - H_m111(z)
                  + H_m11m1(z) - H_m1m11(z) + H_m1m1m1(z);

    double II = 4 * z * Hmp;
    double J = 4 * z * L;
    double K = 4 * z * Hmpm;

    double a11 = CA;

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    double Logxi = log(1 + 1. / (4 * m2Q2));

    return a11 * a11
           * (-32. / 3 * K * m2Q2 * (1. + 3. * m2Q2)
              - 128. / 27 * (17. + 120. * m2Q2)
              + 16. / 27 * J * (3. + 136. * m2Q2 + 480. * m4Q4)
              + II
                    * (32. / 3 * Logxi * m2Q2 * (1. + 3. * m2Q2)
                       + 16. / 9 * (-3. - 4. * m2Q2 + 24. * m4Q4))
              + (64. / 3 * II * m2Q2 * (1. + 3. * m2Q2)
                 - 64. / 9 * (-1. + 12 * m2Q2)
                 + 16. / 9 * J * (-3. - 4. * m2Q2 + 24. * m4Q4))
                    * Lmu
              + (32. / 3 * J * m2Q2 * (1. + 3. * m2Q2)
                 - 32. / 3 * (1. + 6. * m2Q2))
                    * Lmu2)
           * log(x) / x / (1. + 4. * m2Q2);
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for FL at O(as^3)
//  at next to leading log.
//------------------------------------------------------------------------------------------//

Value HighEnergyCoefficientFunction::CL_g3_highenergyNLL(
    double x, double m2Q2, double m2mu2, int nf
) const {

    double m4Q4 = m2Q2 * m2Q2;

    double z = sqrt(1. / (1. + 4. * m2Q2));

    double L = log((1. + z) / (1. - z));

    double Hmp = H_11(z) + H_1m1(z) - H_m11(z) - H_m1m1(z);

    double Hmpm = H_111(z) - H_11m1(z) + H_1m11(z) - H_1m1m1(z) - H_m111(z)
                  + H_m11m1(z) - H_m1m11(z) + H_m1m1m1(z);

    double II = 4 * z * Hmp;
    double J = 4 * z * L;
    double K = 4 * z * Hmpm;

    double a11 = a_11();
    double a21 = a_21(nf);
    double a10 = a_10(nf);

    double beta0 = beta(0, nf);

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    double Logxi = log(1. + 1. / (4. * m2Q2));

    double central =
        (a21
             * (-64. / 3 * II * m2Q2 * (1. + 3. * m2Q2)
                + 64. / 9 * (-1. + 12. * m2Q2)
                - 16. / 9 * J * (-3. - 4. * m2Q2 + 24. * m4Q4))
         + a10 * a11
               * (64. / 3 * K * m2Q2 * (1. + 3. * m2Q2)
                  + 256. / 27 * (17. + 120. * m2Q2)
                  - 32. / 27 * J * (3. + 136. * m2Q2 + 480. * m4Q4)
                  + II
                        * (-64. / 3 * Logxi * m2Q2 * (1. + 3. * m2Q2)
                           - 32. / 9 * (-3. - 4. * m2Q2 + 24. * m4Q4)))
         + a11 * beta0
               * (-32. / 3 * K * m2Q2 * (1. + 3. * m2Q2)
                  - 128. / 27 * (17. + 120. * m2Q2)
                  + 16. / 27 * J * (3. + 136. * m2Q2 + 480. * m4Q4)
                  + II
                        * (32. / 3 * Logxi * m2Q2 * (1. + 3. * m2Q2)
                           + 16. / 9 * (-3. - 4. * m2Q2 + 24. * m4Q4)))
         + (a21
                * (-64. / 3 * J * m2Q2 * (1. + 3. * m2Q2)
                   + 64. / 3 * (1. + 6. * m2Q2))
            + a10 * a11
                  * (-128. / 3 * II * m2Q2 * (1. + 3. * m2Q2)
                     + 128. / 9 * (-1. + 12. * m2Q2)
                     - 32. / 9 * J * (-3. - 4. * m2Q2 + 24. * m4Q4))
            + a11 * beta0
                  * (64. / 3 * II * m2Q2 * (1. + 3. * m2Q2)
                     - 64. / 9 * (-1. + 12. * m2Q2)
                     + 16. / 9 * J * (-3. - 4. * m2Q2 + 24. * m4Q4)))
               * Lmu
         + (a11 * beta0
                * (32. / 3 * J * m2Q2 * (1. + 3. * m2Q2)
                   - 32. / 3 * (1. + 6. * m2Q2))
            + a10 * a11
                  * (-64. / 3 * J * m2Q2 * (1. + 3. * m2Q2)
                     + 64. / 3 * (1. + 6. * m2Q2)))
               * Lmu2)
        / x / (1. + 4. * m2Q2);

    double error =
        (a10 * a11
             * (64. / 3 * K * m2Q2 * (1 + 3 * m2Q2)
                + 256. / 27 * (17 + 120 * m2Q2)
                - 32. / 27 * J * (3 + 136 * m2Q2 + 480 * m4Q4)
                + II
                      * (-64. / 3 * Logxi * m2Q2 * (1 + 3 * m2Q2)
                         - 32. / 9 * (-3 - 4 * m2Q2 + 24 * m4Q4)))
         + (a11 * beta0
                * (32. / 3 * J * m2Q2 * (1 + 3 * m2Q2)
                   - 32. / 3 * (1 + 6 * m2Q2))
            + a10 * a11
                  * (-64. / 3 * J * m2Q2 * (1 + 3 * m2Q2)
                     + 64. / 3 * (1 + 6 * m2Q2)))
               * Lmu2
         + a11 * beta0
               * (-32. / 3 * K * m2Q2 * (1 + 3 * m2Q2)
                  - 128. / 27 * (17 + 120 * m2Q2)
                  - 256. / 9 * (-1 + 12 * m2Q2) * ln2
                  + 56. / 3 * (-1 + 12 * m2Q2) * zeta3
                  + II
                        * (32. / 3 * Logxi * m2Q2 * (1 + 3 * m2Q2)
                           + 16. / 9 * (-3 - 4 * m2Q2 + 24 * m4Q4)
                           + 256. / 3 * m2Q2 * (1 + 3 * m2Q2) * ln2
                           - 56 * m2Q2 * (1 + 3 * m2Q2) * zeta3)
                  + J
                        * (16. / 27 * (3 + 136 * m2Q2 + 480 * m4Q4)
                           + 64. / 9 * (-3 - 4 * m2Q2 + 24 * m4Q4) * ln2
                           - 14. / 3 * (-3 - 4 * m2Q2 + 24 * m4Q4) * zeta3))
         + Lmu
               * (a10 * a11
                      * (-128. / 3 * II * m2Q2 * (1 + 3 * m2Q2)
                         + 128. / 9 * (-1 + 12 * m2Q2)
                         - 32. / 9 * J * (-3 - 4 * m2Q2 + 24 * m4Q4))
                  + a11 * beta0
                        * (64. / 3 * II * m2Q2 * (1 + 3 * m2Q2)
                           - 64. / 9 * (-1 + 12 * m2Q2)
                           - 256. / 3 * (1 + 6 * m2Q2) * ln2
                           + 56 * (1 + 6 * m2Q2) * zeta3
                           + J
                                 * (16. / 9 * (-3 - 4 * m2Q2 + 24 * m4Q4)
                                    + 256. / 3 * m2Q2 * (1 + 3 * m2Q2) * ln2
                                    - 56 * m2Q2 * (1 + 3 * m2Q2) * zeta3))))
        / x / (1. + 4. * m2Q2);

    double delta = fabs(central - error);

    double higher = central + delta;
    double lower = central - delta;

    return Value(central, higher, lower);
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for FL at O(as^3).
//------------------------------------------------------------------------------------------//

Value HighEnergyCoefficientFunction::CL_g3_highenergy(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return CL_g3_highenergyLL(x, m2Q2, m2mu2)
           + CL_g3_highenergyNLL(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  High energy limit of the quark coefficient function for FL at O(as^3)
//  at leading log.
//------------------------------------------------------------------------------------------//

double HighEnergyCoefficientFunction::CL_ps3_highenergyLL(
    double x, double m2Q2, double m2mu2
) const {

    return CF / CA * CL_g3_highenergyLL(x, m2Q2, m2mu2);
}

Value HighEnergyCoefficientFunction::CL_ps3_highenergyNLL(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return CF / CA * CL_g3_highenergyNLL(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  High energy limit of the quark coefficient function for FL at O(as^3).
//------------------------------------------------------------------------------------------//

Value HighEnergyCoefficientFunction::CL_ps3_highenergy(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return CF / CA * CL_g3_highenergy(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the gluon coefficient function
//  for F2 at O(as^3) at leading log.
//
//  Eq. (3.41) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double HighEnergyHighScaleCoefficientFunction::C2_g3_highenergy_highscaleLL(
    double x, double m2Q2, double m2mu2
) const {

    double LQ = log(m2Q2);
    double LQ2 = LQ * LQ;
    double LQ3 = LQ * LQ2;

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    return CA * CA * log(x)
           * (-32. / 27 * (-71. + 18 * zeta2) * LQ - 208. / 9 * LQ2
              + 32. / 9 * LQ3 + Lmu2 * (-16. / 3 + 32. / 3 * LQ)
              + Lmu
                    * (32. / 9 * (-5. + 6 * zeta2) + 416. / 9 * LQ
                       - 32. / 3 * LQ2)
              + 16. / 27 * (-92. + 78. * zeta2 - 72. * zeta3))
           / x;
}

//==========================================================================================//
//  High scale limit of the high energy limit of the gluon coefficient function
//  for F2 at O(as^3) at next to leading log.
//------------------------------------------------------------------------------------------//

Value HighEnergyHighScaleCoefficientFunction::C2_g3_highenergy_highscaleNLL(
    double x, double m2Q2, double m2mu2, int nf
) const {

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    double pi2 = M_PI * M_PI;

    double LQ = log(m2Q2);
    double LQ2 = LQ * LQ;
    double LQ3 = LQ2 * LQ;

    double a11 = a_11();
    double a21 = a_21(nf);
    double a10 = a_10(nf);

    double beta0 = beta(0, nf);

    double central =
        (-32. / 9 * a21 * (-5. + pi2)
         + (-416. * a21 / 9 + 64. / 27 * a10 * a11 * (-71. + 3. * pi2)
            - 32. / 27 * a11 * beta0 * (-71. + 3. * pi2))
               * LQ
         + (416. * a10 * a11 / 9 + 32. * a21 / 3 - 208. * a11 * beta0 / 9) * LQ2
         + (-64. * a10 * a11 / 9 + 32. * a11 * beta0 / 9) * LQ3
         + Lmu2
               * (32. * a10 * a11 / 3 - 16. * a11 * beta0 / 3
                  + (-64 * a10 * a11 / 3 + 32. * a11 * beta0 / 3) * LQ)
         + Lmu
               * (32. * a21 / 3 - 64. / 9 * a10 * a11 * (-5. + pi2)
                  + 32. / 9 * a11 * beta0 * (-5. + pi2)
                  + (-832. * a10 * a11 / 9 - 64 * a21 / 3
                     + 416. * a11 * beta0 / 9)
                        * LQ
                  + (64. * a10 * a11 / 3 - 32. * a11 * beta0 / 3) * LQ2)
         - 32. / 27 * a10 * a11 * (-92. + 13. * pi2 - 72. * zeta3)
         + 16. / 27 * a11 * beta0 * (-92. + 13. * pi2 - 72. * zeta3))
        / x;

    double error =
        ((-64 * a10 * a11 / 9 + 32 * a11 * beta0 / 9) * LQ3
         + Lmu2
               * (32 * a10 * a11 / 3 - 16 * a11 * beta0 / 3
                  + (-64 * a10 * a11 / 3 + 32 * a11 * beta0 / 3) * LQ)
         + LQ2
               * (416 * a10 * a11 / 9
                  - 4. / 9 * a11 * beta0 * (52 + 96 * ln2 - 63 * zeta3))
         - 32. / 27 * a10 * a11 * (-92 + 13 * pi2 - 72 * zeta3)
         + 4. / 27 * a11 * beta0
               * (-368 + 52 * pi2 - 480 * ln2 + 96 * pi2 * ln2 + 27 * zeta3
                  - 63 * pi2 * zeta3)
         + Lmu
               * (-64. / 9 * a10 * a11 * (-5 + pi2)
                  + (64 * a10 * a11 / 3 - 32 * a11 * beta0 / 3) * LQ2
                  + LQ
                        * (-832 * a10 * a11 / 9
                           + 8. / 9 * a11 * beta0 * (52 + 96 * ln2 - 63 * zeta3)
                        )
                  + 4. / 9 * a11 * beta0
                        * (-40 + 8 * pi2 - 96 * ln2 + 63 * zeta3))
         + LQ
               * (64. / 27 * a10 * a11 * (-71 + 3 * pi2)
                  - 4. / 27 * a11 * beta0
                        * (-568 + 24 * pi2 - 1248 * ln2 + 819 * zeta3)))
        / x;

    double delta = fabs(central - error);

    double higher = central + delta;
    double lower = central - delta;

    return Value(central, higher, lower);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the gluon coefficient function
//  for F2 at O(as^3).
//------------------------------------------------------------------------------------------//

Value HighEnergyHighScaleCoefficientFunction::C2_g3_highenergy_highscale(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return C2_g3_highenergy_highscaleLL(x, m2Q2, m2mu2)
           + C2_g3_highenergy_highscaleNLL(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the quark coefficient function
//  for F2 at O(as^3).
//------------------------------------------------------------------------------------------//

Value HighEnergyHighScaleCoefficientFunction::C2_ps3_highenergy_highscale(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return CF / CA * C2_g3_highenergy_highscale(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the quark coefficient function
//  for F2 at O(as^3) at leading log.
//------------------------------------------------------------------------------------------//

double HighEnergyHighScaleCoefficientFunction::C2_ps3_highenergy_highscaleLL(
    double x, double m2Q2, double m2mu2
) const {

    return CF / CA * C2_g3_highenergy_highscaleLL(x, m2Q2, m2mu2);
}

Value HighEnergyHighScaleCoefficientFunction::C2_ps3_highenergy_highscaleNLL(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return CF / CA * C2_g3_highenergy_highscaleNLL(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the gluon coefficient function
//  for FL at O(as^3) at leading log.
//------------------------------------------------------------------------------------------//

double HighEnergyHighScaleCoefficientFunction::CL_g3_highenergy_highscaleLL(
    double x, double m2Q2, double m2mu2
) const {

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    double LQ = log(m2Q2);
    double LQ2 = LQ * LQ;

    return CA * CA
           * (32. / 27 * (-68. + 18. * zeta2) - 32. / 3 * Lmu2 - 64. / 9 * LQ
              - 32. / 3 * LQ2 + Lmu * (64. / 9 + 64. / 3 * LQ))
           * log(x) / x;
}

//==========================================================================================//
//  High scale limit of the high energy limit of the gluon coefficient function
//  for FL at O(as^3) at next to leading log.
//------------------------------------------------------------------------------------------//

Value HighEnergyHighScaleCoefficientFunction::CL_g3_highenergy_highscaleNLL(
    double x, double m2Q2, double m2mu2, int nf
) const {

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    double pi2 = M_PI * M_PI;

    double LQ = log(m2Q2);
    double LQ2 = LQ * LQ;

    double a11 = a_11();
    double a21 = a_21(nf);
    double a10 = a_10(nf);

    double beta0 = beta(0, nf);

    double central =
        (-64. * a21 / 9 - 64. / 27 * a10 * a11 * (-68. + 3 * pi2)
         + 32. / 27 * a11 * beta0 * (-68. + 3 * pi2)
         + (64. * a10 * a11 / 3 - 32. * a11 * beta0 / 3) * Lmu2
         + (128. * a10 * a11 / 9 - 64. * a21 / 3 - 64. * a11 * beta0 / 9) * LQ
         + (64. * a10 * a11 / 3 - 32. * a11 * beta0 / 3) * LQ2
         + Lmu
               * (-128. * a10 * a11 / 9 + 64. * a21 / 3 + 64. * a11 * beta0 / 9
                  + (-128. * a10 * a11 / 3 + 64. * a11 * beta0 / 3) * LQ))
        / x;

    double error =
        (-64. / 27 * a10 * a11 * (-68 + 3 * pi2)
         + (64 * a10 * a11 / 3 - 32 * a11 * beta0 / 3) * Lmu2
         + (64 * a10 * a11 / 3 - 32 * a11 * beta0 / 3) * LQ2
         + Lmu
               * (-128 * a10 * a11 / 9
                  + (-128 * a10 * a11 / 3 + 64 * a11 * beta0 / 3) * LQ
                  - 8. / 9 * a11 * beta0 * (-8 + 96 * ln2 - 63 * zeta3))
         + LQ
               * (128 * a10 * a11 / 9
                  + 8. / 9 * a11 * beta0 * (-8 + 96 * ln2 - 63 * zeta3))
         + 8. / 27 * a11 * beta0 * (-272 + 12 * pi2 + 96 * ln2 - 63 * zeta3))
        / x;

    double delta = fabs(central - error);

    double higher = central + delta;
    double lower = central - delta;

    return Value(central, higher, lower);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the gluon coefficient function
//  for FL at O(as^3).
//------------------------------------------------------------------------------------------//

Value HighEnergyHighScaleCoefficientFunction::CL_g3_highenergy_highscale(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return CL_g3_highenergy_highscaleLL(x, m2Q2, m2mu2)
           + CL_g3_highenergy_highscaleNLL(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the quark coefficient function
//  for F2 at O(as^3).
//------------------------------------------------------------------------------------------//

Value HighEnergyHighScaleCoefficientFunction::CL_ps3_highenergy_highscale(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return CF / CA * CL_g3_highenergy_highscale(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  High scale limit of the high energy limit of the quark coefficient function
//  for F2 at O(as^3) at leading log.
//------------------------------------------------------------------------------------------//

double HighEnergyHighScaleCoefficientFunction::CL_ps3_highenergy_highscaleLL(
    double x, double m2Q2, double m2mu2
) const {

    return CF / CA * CL_g3_highenergy_highscaleLL(x, m2Q2, m2mu2);
}

Value HighEnergyHighScaleCoefficientFunction::CL_ps3_highenergy_highscaleNLL(
    double x, double m2Q2, double m2mu2, int nf
) const {

    return CF / CA * CL_g3_highenergy_highscaleNLL(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//                  Color factors O(as^3)
//------------------------------------------------------------------------------------------//

double AbstractHighEnergyCoefficientFunction::a_10(int nf) const {

    return -(11. * CA + 2. * nf * (1. - 2. * CF / CA)) / 12.;
}

double AbstractHighEnergyCoefficientFunction::a_11() const { return CA; }

double AbstractHighEnergyCoefficientFunction::a_21(int nf) const {
    return nf * (26. * CF - 23. * CA) / 36.;
}
