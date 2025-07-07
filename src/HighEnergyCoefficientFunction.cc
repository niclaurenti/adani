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
    : CoefficientFunction(order, kind, channel), NLL_(NLL){
        if(order == 1) {
            fx_ = &AbstractHighEnergyCoefficientFunction::ZeroFunctionBand;
        } else if (order == 2) {
            fx_ = &AbstractHighEnergyCoefficientFunction::Order2;
        } else {
            if (NLL) {
                fx_ = &AbstractHighEnergyCoefficientFunction::Order3;
            } else {
                fx_ = &AbstractHighEnergyCoefficientFunction::Order3LL;
            }
        }
    };

//==========================================================================================//
//  AbstractHighEnergyCoefficientFunction: band of the high energy coefficient function
//------------------------------------------------------------------------------------------//

Value AbstractHighEnergyCoefficientFunction::Order2(
    double x, double m2Q2, double m2mu2, int /*nf*/
) const {
    return Value(LL(m2Q2, m2mu2) / x);
}

//==========================================================================================//
//  AbstractHighEnergyCoefficientFunction: band of the high energy coefficient function
//------------------------------------------------------------------------------------------//

Value AbstractHighEnergyCoefficientFunction::Order3(
    double x, double m2Q2, double m2mu2, int nf
) const {
    return (LL(m2Q2, m2mu2) * log(x) + NLL(m2Q2, m2mu2, nf)) / x;
}

//==========================================================================================//
//  AbstractHighEnergyCoefficientFunction: band of the high energy coefficient function
//------------------------------------------------------------------------------------------//

Value AbstractHighEnergyCoefficientFunction::Order3LL(
    double x, double m2Q2, double m2mu2, int /*nf*/
) const {
    return Value(LL(m2Q2, m2mu2) * log(x) / x);
}

//==========================================================================================//
//  AbstractHighEnergyCoefficientFunction: band of the high energy coefficient function
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
//  HighEnergyCoefficientFunction: function that sets the pointer for LL_ and
//  NLL_
//------------------------------------------------------------------------------------------//

void HighEnergyCoefficientFunction::SetFunctions() {

    if (GetOrder() == 1) {

        LL_ = nullptr;
        NLL_ = nullptr;

    } else if (GetOrder() == 2) {

        NLL_ = nullptr;
        if (GetKind() == '2' && GetChannel() == 'g')
            LL_ = &HighEnergyCoefficientFunction::C2_g2_highenergyLL;
        else if (GetKind() == '2' && GetChannel() == 'q')
            LL_ = &HighEnergyCoefficientFunction::C2_ps2_highenergyLL;
        else if (GetKind() == 'L' && GetChannel() == 'g')
            LL_ = &HighEnergyCoefficientFunction::CL_g2_highenergyLL;
        else if (GetKind() == 'L' && GetChannel() == 'q')
            LL_ = &HighEnergyCoefficientFunction::CL_ps2_highenergyLL;
        else {
            throw UnexpectedException(
                "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
            );
        }

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
            throw UnexpectedException(
                "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
            );
        }
    } else {
        throw UnexpectedException(
            "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
        );
    }

    if (!GetNLL())
        NLL_ = nullptr;
}

//==========================================================================================//
//  HighEnergyCoefficientFunction: band of the high energy coefficient function
//------------------------------------------------------------------------------------------//

double HighEnergyCoefficientFunction::LL(
    double m2Q2, double m2mu2
) const {
    return (this->*LL_)(m2Q2, m2mu2);;
}

//==========================================================================================//
//  HighEnergyCoefficientFunction: band of the high energy coefficient function
//------------------------------------------------------------------------------------------//

Value HighEnergyCoefficientFunction::NLL(
    double m2Q2, double m2mu2, int nf
) const {
    return (this->*NLL_)(m2Q2, m2mu2, nf);;
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
//  HighEnergyHighScaleCoefficientFunction: function that sets the pointer for
//  LL_ and NLL_
//------------------------------------------------------------------------------------------//

void HighEnergyHighScaleCoefficientFunction::SetFunctions() {

    if (GetOrder() == 1) {

        LL_ = nullptr;
        NLL_ = nullptr;

    } else if (GetOrder() == 2) {

        NLL_ = nullptr;
        if (GetKind() == '2' && GetChannel() == 'g')
            LL_ = &HighEnergyHighScaleCoefficientFunction::
                      C2_g2_highenergy_highscaleLL;
        else if (GetKind() == '2' && GetChannel() == 'q')
            LL_ = &HighEnergyHighScaleCoefficientFunction::
                      C2_ps2_highenergy_highscaleLL;
        else if (GetKind() == 'L' && GetChannel() == 'g')
            LL_ = &HighEnergyHighScaleCoefficientFunction::
                      CL_g2_highenergy_highscaleLL;
        else if (GetKind() == 'L' && GetChannel() == 'q')
            LL_ = &HighEnergyHighScaleCoefficientFunction::
                      CL_ps2_highenergy_highscaleLL;
        else {
            throw UnexpectedException(
                "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
            );
        }

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
            throw UnexpectedException(
                "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
            );
        }
    } else {
        throw UnexpectedException(
            "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
        );
    }
    if (!GetNLL())
        NLL_ = nullptr;
}

//==========================================================================================//
//  HighEnergyCoefficientFunction: band of the high energy coefficient function
//------------------------------------------------------------------------------------------//

double HighEnergyHighScaleCoefficientFunction::LL(
    double m2Q2, double m2mu2
) const {
    return (this->*LL_)(m2Q2, m2mu2);;
}

//==========================================================================================//
//  HighEnergyCoefficientFunction: band of the high energy coefficient function
//------------------------------------------------------------------------------------------//

Value HighEnergyHighScaleCoefficientFunction::NLL(
    double m2Q2, double m2mu2, int nf
) const {
    return (this->*NLL_)(m2Q2, m2mu2, nf);;
}

//==========================================================================================//
//  High energy limit of the gluon coefficient function for F2 at O(as^2).
//
//  Eq. (3.38) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double HighEnergyCoefficientFunction::C2_g2_highenergyLL(
    double m2Q2, double m2mu2
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

    return 4. / 3. * CA * (c_const + c_Lmu * Lmu);
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

    return 8 * CA * TR * (c_const + c_Lmu * Lmu);
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

    double LQ = log(1. / m2Q2);
    double L2Q = LQ * LQ;

    double Lm = log(m2mu2);

    return CA
           * (8. / 3. * L2Q + 104. / 9. * LQ + 40. / 9. - 16. / 3. * zeta2
              + (16. / 3. * LQ + 8. / 3.) * Lm);
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

    double LQ = log(1. / m2Q2);

    double Lm = log(1. / m2mu2);

    double c_const = -2. / 9. * (1. - 3. * LQ);

    double c_log = -2. / 3.;

    return 16 * CA * TR * (c_const + c_log * Lm);
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

    double a11 = a_11();

    return a11 * a11
           * (-1472. / 27 - 8. / 3 * K * (-1. + m2Q2)
              + 8. / 27 * J * (-71. + 92. * m2Q2)
              + I
                    * (8. / 3 * Logxi * (-1. + m2Q2)
                       + 8. / 9 * (-13. + 10. * m2Q2))
              + (-160. / 9 + 16. / 3 * I * (-1. + m2Q2)
                 + 8. / 9 * J * (-13. + 10. * m2Q2))
                    * Lmu
              + (-16. / 3 + 8. / 3 * J * (-1. + m2Q2)) * Lmu2);
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

    double central_value = C2_g3_highenergyNLL(m2Q2, m2mu2, a11, a10, a21, beta_0);
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
    double m2Q2, double m2mu2, double a11, double a10, double a21, double beta0
) const {

    double z = sqrt(1. / (1. + 4. * m2Q2));

    double L = log((1. + z) / (1. - z));

    double Hmp = H_11(z) + H_1m1(z) - H_m11(z) - H_m1m1(z);

    double Hmpm = H_111(z) - H_11m1(z) + H_1m11(z) - H_1m1m1(z) - H_m111(z)
                  + H_m11m1(z) - H_m1m11(z) + H_m1m1m1(z);

    double II = 4 * z * Hmp;
    double J = 4 * z * L;
    double K = 4 * z * Hmpm;

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    double Logxi = log(1. + 1. / (4. * m2Q2));

    double res =
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
               * Lmu2);

    return res;
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
           / (1. + 4. * m2Q2);
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

    double central_value = CL_g3_highenergyNLL(m2Q2, m2mu2, a11, a10, a21, beta_0);
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
    double m2Q2, double m2mu2, double a11, double a10, double a21, double beta0
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

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    double Logxi = log(1. + 1. / (4. * m2Q2));

    double res =
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
        / (1. + 4. * m2Q2);

    return res;
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

    double LQ = log(m2Q2);
    double LQ2 = LQ * LQ;
    double LQ3 = LQ * LQ2;

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    double a11 = a_11();

    return a11 * a11
           * (-32. / 27 * (-71. + 18 * zeta2) * LQ - 208. / 9 * LQ2
              + 32. / 9 * LQ3 + Lmu2 * (-16. / 3 + 32. / 3 * LQ)
              + Lmu
                    * (32. / 9 * (-5. + 6 * zeta2) + 416. / 9 * LQ
                       - 32. / 3 * LQ2)
              + 16. / 27 * (-92. + 78. * zeta2 - 72. * zeta3));
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

    double central_value = C2_g3_highenergy_highscaleNLL(m2Q2, m2mu2, a11, a10, a21, beta_0);
    double error = C2_g3_highenergy_highscaleNLL(m2Q2, m2mu2, a11, a10, a21_new, beta_0);

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
    double m2Q2, double m2mu2, double a11, double a10, double a21, double beta0
) const {

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    double pi2 = M_PI * M_PI;

    double LQ = log(m2Q2);
    double LQ2 = LQ * LQ;
    double LQ3 = LQ2 * LQ;

    double res =
        -32. / 9 * a21 * (-5. + pi2)
         + (-416. * a21 / 9 + 64. / 27 * a10 * a11 * (-71. + 3. * pi2)
            - 32. / 27 * a11 * beta0 * (-71. + 3. * pi2))
               * LQ
         + (416. * a10 * a11 / 9 + 32. * a21 / 3 - 208. * a11 * beta0 / 9)
               * LQ2
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
         + 16. / 27 * a11 * beta0 * (-92. + 13. * pi2 - 72. * zeta3);

    return res;
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

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    double LQ = log(m2Q2);
    double LQ2 = LQ * LQ;

    return CA * CA
           * (32. / 27 * (-68. + 18. * zeta2) - 32. / 3 * Lmu2 - 64. / 9 * LQ
              - 32. / 3 * LQ2 + Lmu * (64. / 9 + 64. / 3 * LQ));
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

    double central_value = CL_g3_highenergy_highscaleNLL(m2Q2, m2mu2, a11, a10, a21, beta_0);
    double error = CL_g3_highenergy_highscaleNLL(m2Q2, m2mu2, a11, a10, a21_new, beta_0);

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
    double m2Q2, double m2mu2, double a11, double a10, double a21, double beta0
) const {

    double Lmu = log(m2mu2);
    double Lmu2 = Lmu * Lmu;

    double pi2 = M_PI * M_PI;

    double LQ = log(m2Q2);
    double LQ2 = LQ * LQ;

    double res =
        (-64. * a21 / 9 - 64. / 27 * a10 * a11 * (-68. + 3 * pi2)
         + 32. / 27 * a11 * beta0 * (-68. + 3 * pi2)
         + (64. * a10 * a11 / 3 - 32. * a11 * beta0 / 3) * Lmu2
         + (128. * a10 * a11 / 9 - 64. * a21 / 3 - 64. * a11 * beta0 / 9) * LQ
         + (64. * a10 * a11 / 3 - 32. * a11 * beta0 / 3) * LQ2
         + Lmu
               * (-128. * a10 * a11 / 9 + 64. * a21 / 3 + 64. * a11 * beta0 / 9
                  + (-128. * a10 * a11 / 3 + 64. * a11 * beta0 / 3) * LQ));

    return res;
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

double AbstractHighEnergyCoefficientFunction::a_10(int nf) const {

    return -(11. * CA + 2. * nf * (1. - 2. * CF / CA)) / 12.;
}

double AbstractHighEnergyCoefficientFunction::a_11() const { return CA; }

double AbstractHighEnergyCoefficientFunction::a_21(int nf) const {
    return nf * (26. * CF - 23. * CA) / 36.;
}

double AbstractHighEnergyCoefficientFunction::a_21_new(int nf) const {
    return beta0(nf) * a_11() * (21. / 8 * zeta3 - 4 * ln2);
}
