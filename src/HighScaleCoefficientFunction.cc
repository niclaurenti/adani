#include "adani/HighScaleCoefficientFunction.h"
#include "adani/Constants.h"
#include "adani/SpecialFunctions.h"

#include <cmath>

//==========================================================================================//
//  HighScaleCoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

HighScaleCoefficientFunction::HighScaleCoefficientFunction(
    const int &order, const char &kind, const char &channel,
    const string &version
)
    : CoefficientFunction(order, kind, channel) {
    massless_as1_ = nullptr;
    massless_as2_ = nullptr;
    massless_as3_ = nullptr;
    a_muindep_ = nullptr;

    if (GetChannel() == 'g')
        massless_as1_ =
            new MasslessCoefficientFunction(1, GetKind(), GetChannel());

    if (order > 1)
        massless_as2_ =
            new MasslessCoefficientFunction(2, GetKind(), GetChannel());
    if (order > 2)
        massless_as3_ =
            new MasslessCoefficientFunction(3, GetKind(), GetChannel());

    if (GetOrder() == 3 && GetKind() == '2') {
        a_muindep_ = new MatchingCondition(3, 'Q', GetChannel(), version);
    }

    try {

        if (GetOrder() < 3 && version != "exact") {
            throw NotValidException(
                "HighScaleCoefficientFunction at orders 1 and 2 are only "
                "'exact'!",
                __PRETTY_FUNCTION__, __LINE__
            );
        }

        SetFunctions();
    } catch (NotPresentException &e) {
        e.runtime_error();
    } catch (UnexpectedException &e) {
        e.runtime_error();
    } catch (NotImplementedException &e) {
        e.runtime_error();
    }
}

//==========================================================================================//
//  HighScaleCoefficientFunction: destructor
//------------------------------------------------------------------------------------------//

HighScaleCoefficientFunction::~HighScaleCoefficientFunction() {
    delete massless_as1_;
    delete massless_as2_;
    delete massless_as3_;
    delete a_muindep_;
}

//==========================================================================================//
//  HighScaleCoefficientFunction: central value of the highscale coefficient
//  function
//------------------------------------------------------------------------------------------//

double HighScaleCoefficientFunction::fx(
    double x, double m2Q2, double m2mu2, int nf
) const {
    return fxBand(x, m2Q2, m2mu2, nf).GetCentral();
}

//==========================================================================================//
//  HighScaleCoefficientFunction: band of the highscale coefficient function
//------------------------------------------------------------------------------------------//

Value HighScaleCoefficientFunction::fxBand(
    double x, double m2Q2, double m2mu2, int nf
) const {
    return (this->*fx_)(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  HighScaleCoefficientFunction: band of the highscale coefficient function
//  without ordering the upper and lower bands
//------------------------------------------------------------------------------------------//

vector<double> HighScaleCoefficientFunction::fxBand_NotOrdered(
    double x, double m2Q2, double m2mu2, int nf
) const {

    if (GetOrder() == 2)
        return fxBand(x, m2Q2, m2mu2, nf).ToVect();

    double central, higher, lower;
    central = fxBand(x, m2Q2, m2mu2, nf).GetCentral();

    try {
        higher = fxBand(x, m2Q2, m2mu2, nf).GetHigher()
                 - a_muindep_->MuIndependentNfIndependentTerm(x).GetHigher()
                 + (a_muindep_->NotOrdered(x))[1];
        lower = fxBand(x, m2Q2, m2mu2, nf).GetLower()
                - a_muindep_->MuIndependentNfIndependentTerm(x).GetLower()
                + (a_muindep_->NotOrdered(x))[2];
    } catch (NotImplementedException &e) {
        e.runtime_error();
    } catch (UnexpectedException &e) {
        e.runtime_error();
    }

    return { central, higher, lower };
}

//==========================================================================================//
//  HighScaleCoefficientFunction: function that sets the pointer for fxBand
//------------------------------------------------------------------------------------------//

void HighScaleCoefficientFunction::SetFunctions() {

    if (GetOrder() == 1) {

        if (GetKind() == '2' && GetChannel() == 'g')
            fx_ = &HighScaleCoefficientFunction::C2_g1_highscale;
        else if (GetKind() == 'L' && GetChannel() == 'g')
            fx_ = &HighScaleCoefficientFunction::CL_g1_highscale;
        else {
            throw NotPresentException(
                "quark coefficient function is not present at order 1! Got "
                "order="
                    + to_string(GetOrder()),
                __PRETTY_FUNCTION__, __LINE__
            );
        }

    } else if (GetOrder() == 2) {

        if (GetKind() == '2' && GetChannel() == 'g')
            fx_ = &HighScaleCoefficientFunction::C2_g2_highscale;
        else if (GetKind() == '2' && GetChannel() == 'q')
            fx_ = &HighScaleCoefficientFunction::C2_ps2_highscale;
        else if (GetKind() == 'L' && GetChannel() == 'g')
            fx_ = &HighScaleCoefficientFunction::CL_g2_highscale;
        else if (GetKind() == 'L' && GetChannel() == 'q')
            fx_ = &HighScaleCoefficientFunction::CL_ps2_highscale;
        else {
            throw UnexpectedException(
                "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
            );
        }
    } else if (GetOrder() == 3) {

        if (GetKind() == '2' && GetChannel() == 'g')
            fx_ = &HighScaleCoefficientFunction::C2_g3_highscale;
        else if (GetKind() == '2' && GetChannel() == 'q')
            fx_ = &HighScaleCoefficientFunction::C2_ps3_highscale;
        else if (GetKind() == 'L' && GetChannel() == 'g')
            fx_ = &HighScaleCoefficientFunction::CL_g3_highscale;
        else if (GetKind() == 'L' && GetChannel() == 'q')
            fx_ = &HighScaleCoefficientFunction::CL_ps3_highscale;
        else {
            throw UnexpectedException(
                "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
            );
        }
    } else {
        throw UnexpectedException(
            "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
        );
    }
}

//==========================================================================================//
// OBSERVATION: the highscale coefficient functions depend on the massless
// coefficient function through C_massless(nf + 1)/(nf + 1). Since up to O(as^2)
// the massless coefficient functions are proportional to nf, the (nf + 1)
// dependence cancels and therefore it is equal to putting (nf + 1) = 1 .
//
//------------------------------------------------------------------------------------------//

//==========================================================================================//
//  High scale (Q^2 >> m^2) limit of the gluon coefficient functions for F2 at
//  O(as) expanded in terms of \alpha_s^{[nf]}
//
//  Eq. (B.4) of Ref. [arXiv:1205.5727].
//------------------------------------------------------------------------------------------//

Value HighScaleCoefficientFunction::
    C2_g1_highscale(double x, double m2Q2, double /*m2mu2*/, int /*nf*/) const {

    return Value(D2_g1_highscale(x, m2Q2));
}

//==========================================================================================//
//  High scale (Q^2 >> m^2) limit of the gluon coefficient functions for FL at
//  O(as) expanded in terms of \alpha_s^{[nf]}
//------------------------------------------------------------------------------------------//

Value HighScaleCoefficientFunction::
    CL_g1_highscale(double x, double /*m2Q2*/, double /*m2mu2*/, int /*nf*/)
        const {
    return Value(DL_g1_highscale(x));
}

//==========================================================================================//
//  High scale (Q^2 >> m^2) limit of the gluon coefficient functions for F2 at
//  O(as) expanded in terms of \alpha_s^{[nf+1]}
//
//  Eq. (B.4) of Ref. [arXiv:1205.5727].
//------------------------------------------------------------------------------------------//

double
    HighScaleCoefficientFunction::D2_g1_highscale(double x, double m2Q2) const {

    return 4 * TR * (x * x + (x - 1) * (x - 1)) * log(1. / m2Q2)
           + massless_as1_->MuIndependentTerms(x, 1);
}

//==========================================================================================//
//  High scale (Q^2 >> m^2) limit of the gluon coefficient functions for FL at
//  O(as) expanded in terms of \alpha_s^{[nf+1]}
//------------------------------------------------------------------------------------------//

double HighScaleCoefficientFunction::DL_g1_highscale(double x) const {
    return massless_as1_->MuIndependentTerms(x, 1);
}

//==========================================================================================//
//  High scale (Q^2 >> m^2) limit of the gluon coefficient functions for F2 at
//  O(as^2) expanded in terms of \alpha_s^{[nf]}
//
//  Eq. (B.6) of Ref. [arXiv:1205.5727].
//------------------------------------------------------------------------------------------//

Value HighScaleCoefficientFunction::
    C2_g2_highscale(double x, double m2Q2, double m2mu2, int /*nf*/) const {

    double Lmu = log(1. / m2mu2);

    double tmp = D2_g2_highscale(x, m2Q2, m2mu2)
                 + 2. / 3 * Lmu * D2_g1_highscale(x, m2Q2);
    return Value(tmp);
}

//==========================================================================================//
//  High scale (Q^2 >> m^2) limit of the quark coefficient functions for F2 at
//  O(as^2) expanded in terms of \alpha_s^{[nf]}
//
//  Eq. (B.9) of Ref. [arXiv:1205.5727].
//------------------------------------------------------------------------------------------//

double HighScaleCoefficientFunction::D2_ps2_highscale(
    double z, double m2Q2, double m2mu2
) const {

    double z2 = z * z;
    double z3 = z2 * z;

    double L_M = log(m2mu2);
    double L_M2 = L_M * L_M;
    double L_Q = log(1. / m2Q2) + L_M;
    double L_Q2 = L_Q * L_Q;

    double H0 = H_0(z);
    double H1 = H_1(z);
    double Hm1 = H_m1(z);
    double H01 = H_01(z);
    double H0m1 = H_0m1(z);
    double H001 = H_001(z);
    double H011 = H_011(z);

    return CF * TR
           * ((32. / 3 * H0m1 - 32. / 3 * Hm1 * H0) * (z + 1) * (z + 1)
                  * (z + 1) / z
              + (16. / 3 * H0 * H0 * H0 + 32 * H01 * H0 - 32 * zeta2 * H0
                 - 32 * H001 + 16 * H011 + 16 * zeta3)
                    * (z + 1)
              - 8 * z * (2 * z - 5) * H0 * H0
              + (4 * (z - 1) * (4 * z2 + 7 * z + 4) / (3 * z) - 8 * (z + 1) * H0
                ) * L_M2
              + 16 * (z - 1) * (52 * z2 - 24 * z - 5) / (9 * z)
              + 32. / 3 * (3 * z3 - 3 * z2 - 1) * zeta2 / z
              - 8. / 9 * (88 * z2 + 99 * z - 105) * H0
              + L_Q2
                    * (8 * (z + 1) * H0
                       - 4 * (z - 1) * (4 * z2 + 7 * z + 4) / (3 * z))
              + 16 * (z - 1) * (4 * z2 - 26 * z + 13) * H1 / (9 * z)
              + (z - 1) * (4 * z2 + 7 * z + 4) / z
                    * (-4. / 3 * H1 * H1 - 16. / 3 * H0 * H1)
              - 16 * (2 * z3 - 3 * z2 + 3 * z + 4) * H01 / (3 * z)
              + (8 * (z + 1) * H0 * H0 - 8. / 3 * (8 * z2 + 15 * z + 3) * H0
                 + 16 * (z - 1) * (28 * z2 + z + 10) / (9 * z))
                    * L_M
              + L_Q
                    * (32 * H0 * z2
                       + (z + 1) * (-16 * H0 * H0 - 16 * H01 + 16 * zeta2)
                       - 16 * (z - 1) * (4 * z2 - 26 * z + 13) / (9 * z)
                       + 8 * (z - 1) * (4 * z2 + 7 * z + 4) * H1 / (3 * z)));
}

//==========================================================================================//
//  High scale (Q^2 >> m^2) limit of the gluon coefficient functions for FL at
//  O(as^2) expanded in terms of \alpha_s^{[nf]}
//------------------------------------------------------------------------------------------//

Value HighScaleCoefficientFunction::
    CL_g2_highscale(double x, double m2Q2, double m2mu2, int /*nf*/) const {

    double Lmu = log(1. / m2mu2);

    double tmp =
        DL_g2_highscale(x, m2Q2, m2mu2) + 2. / 3 * Lmu * DL_g1_highscale(x);

    return Value(tmp);
}

//==========================================================================================//
//  High scale (Q^2 >> m^2) limit of the gluon coefficient functions for FL at
//  O(as^2) expanded in terms of \alpha_s^{[nf]}
//------------------------------------------------------------------------------------------//

Value HighScaleCoefficientFunction::
    CL_ps2_highscale(double x, double m2Q2, double m2mu2, int /*nf*/) const {

    return Value(DL_ps2_highscale(x, m2Q2, m2mu2));
}

//==========================================================================================//
//  High scale (Q^2 >> m^2) limit of the gluon coefficient functions for F2 at
//  O(as^2) expanded in terms of \alpha_s^{[nf+1]}
//
//  Eq. (B.5) of Ref. [arXiv:1205.5727].
//------------------------------------------------------------------------------------------//

double HighScaleCoefficientFunction::D2_g2_highscale(
    double x, double m2Q2, double m2mu2
) const {

    double x2 = x * x;

    double H0 = H_0(x);
    double H1 = H_1(x);
    double Hm1 = H_m1(x);
    double H00 = H_00(x);
    double H01 = H_01(x);
    double H11 = H_11(x);
    double H10 = H_10(x);
    double H011 = H_011(x);
    double H111 = H_111(x);
    double H100 = H_100(x);
    double H000 = H_000(x);
    double H010 = H_010(x);
    double H101 = H_101(x);
    double Hm1m10 = H_m1m10(x);
    double Hm10 = H_m10(x);
    double Hm100 = H_m100(x);
    double H110 = H_110(x);

    double LQm = log(1. / m2Q2);
    double LQm2 = LQm * LQm;

    double Lmmu = log(m2mu2);

    return +Lmmu
               * (-4. / 3 + 32. / 3 * x - 32. / 3 * x2 - 4. / 3 * H0
                  + 8. / 3 * H0 * x - 8. / 3 * H0 * x2 - 4. / 3 * H1
                  + 8. / 3 * H1 * x - 8. / 3 * H1 * x2)
           + LQm * Lmmu * (4. / 3 - 8. / 3 * x + 8. / 3 * x2)
           + CF
                 * (13 - 41 * x + 40 * x2 - 4 * zeta3 + 8 * zeta3 * x
                    + 8 * zeta3 * x2 - 4 * zeta2 - 24 * zeta2 * x
                    + 12 * zeta2 * x2 - 8 * H0 - 9 * H0 * x - 24 * H0 * x2 - H00
                    - 12 * H00 * x + 20 * H00 * x2 + 2 * H000 - 4 * H000 * x
                    + 8 * H000 * x2 + 26 * H1 * x - 24 * H1 * x2 + 2 * H10
                    - 24 * H10 * x + 20 * H10 * x2 + 4 * H100 - 8 * H100 * x
                    + 8 * H100 * x2 + 4 * H11 + 8 * H11 * x - 12 * H11 * x2
                    - 4 * H111 + 8 * H111 * x - 8 * H111 * x2 + 4 * H01
                    + 24 * H01 * x - 12 * H01 * x2 - 4 * H010 + 8 * H010 * x
                    - 4 * H011 + 8 * H011 * x - 8 * H011 * x2)
           + CF * LQm
                 * (18 - 34 * x + 8 * x2 - 12 * zeta2 + 24 * zeta2 * x
                    - 32 * zeta2 * x2 + 4 * H0 - 24 * H0 * x + 40 * H0 * x2
                    + 8 * H00 - 16 * H00 * x + 32 * H00 * x2 + 14 * H1
                    - 48 * H1 * x + 40 * H1 * x2 + 16 * H10 - 32 * H10 * x
                    + 32 * H10 * x2 + 16 * H11 - 32 * H11 * x + 32 * H11 * x2
                    + 12 * H01 - 24 * H01 * x + 32 * H01 * x2)
           + CF * LQm2
                 * (-1 + 4 * x - 2 * H0 + 4 * H0 * x - 8 * H0 * x2 - 4 * H1
                    + 8 * H1 * x - 8 * H1 * x2)
           + CA
                 * (-2. / 3 - 224. / 27 / x - 314. / 3 * x + 3176. / 27 * x2
                    + 16 * zeta3 + 56 * zeta3 * x + 8 * zeta2 * x
                    - 2 * zeta2 * x2 - 4 * Hm1 * zeta2 - 8 * Hm1 * zeta2 * x
                    - 8 * Hm1 * zeta2 * x2 - 8 * Hm1m10 - 16 * Hm1m10 * x
                    - 16 * Hm1m10 * x2 + 8 * Hm10 * x + 8 * Hm10 * x2
                    + 4 * Hm100 + 8 * Hm100 * x + 8 * Hm100 * x2 - 28. / 3 * H0
                    - 86. / 3 * H0 * x - 800. / 9 * H0 * x2 + 2 * H00
                    + 8 * H00 * x + 46. / 3 * H00 * x2 - 4 * H000 - 8 * H000 * x
                    - 2 * H1 - 8 * H1 * x + 8 * H1 * x2 + 6 * H10
                    + 16. / 3 * H10 / x + 32 * H10 * x - 130. / 3 * H10 * x2
                    - 2 * H11 - 8 * H11 * x + 10 * H11 * x2 - 4 * H110
                    + 8 * H110 * x - 8 * H110 * x2 + 4 * H111 - 8 * H111 * x
                    + 8 * H111 * x2 - 4 * H101 + 8 * H101 * x - 8 * H101 * x2
                    + 2 * H01 * x2 + 8 * H010 + 32 * H010 * x)
           + CA * Lmmu
                 * (-86. / 3 + 8. / 3 / x - 484. / 3 * x + 562. / 3 * x2
                    + 48 * zeta2 * x - 16 * zeta2 * x2 - 4 * H0 - 128 * H0 * x
                    + 124. / 3 * H0 * x2 - 8 * H00 - 32 * H00 * x + 4 * H1
                    - 16. / 3 * H1 / x - 96 * H1 * x + 316. / 3 * H1 * x2
                    + 8 * H10 - 16 * H10 * x + 16 * H10 * x2 + 16 * H11
                    - 32 * H11 * x + 32 * H11 * x2 - 48 * H01 * x
                    + 16 * H01 * x2)
           + CA * LQm
                 * (-110. / 3 + 104. / 9 / x - 184. / 3 * x + 814. / 9 * x2
                    + 32 * zeta2 * x - 16 * zeta2 * x2 - 8 * Hm10
                    - 16 * Hm10 * x - 16 * Hm10 * x2 - 96 * H0 * x
                    + 100 * H0 * x2 - 16 * H00 - 48 * H00 * x + 4 * H1
                    - 16. / 3 * H1 / x - 80 * H1 * x + 268. / 3 * H1 * x2
                    + 8 * H10 - 16 * H10 * x + 16 * H10 * x2 + 8 * H11
                    - 16 * H11 * x + 16 * H11 * x2 - 48 * H01 * x
                    + 16 * H01 * x2)
           + CA * LQm * Lmmu
                 * (4 + 16. / 3 / x + 32 * x - 124. / 3 * x2 + 8 * H0
                    + 32 * H0 * x - 8 * H1 + 16 * H1 * x - 16 * H1 * x2)
           + CA * LQm2
                 * (2 + 8. / 3 / x + 16 * x - 62. / 3 * x2 + 4 * H0
                    + 16 * H0 * x - 4 * H1 + 8 * H1 * x - 8 * H1 * x2)
           + massless_as2_->MuIndependentTerms(x, 1);
}

//==========================================================================================//
//  High scale (Q^2 >> m^2) limit of the quark coefficient functions for F2 at
//  O(as^2) expanded in terms of \alpha_s^{[nf+1]}
//
//  Eq. (B.6) of Ref. [arXiv:1205.5727].
//------------------------------------------------------------------------------------------//

Value HighScaleCoefficientFunction::
    C2_ps2_highscale(double x, double m2Q2, double m2mu2, int /*nf*/) const {

    return Value(D2_ps2_highscale(x, m2Q2, m2mu2));
}

//==========================================================================================//
//  High scale (Q^2 >> m^2) limit of the quark coefficient functions for FL at
//  O(as^2) expanded in terms of \alpha_s^{[nf+1]}
//------------------------------------------------------------------------------------------//

double HighScaleCoefficientFunction::DL_g2_highscale(
    double z, double m2Q2, double m2mu2
) const {

    double z2 = z * z;
    double z3 = z2 * z;
    double z4 = z3 * z;

    double L_M = log(m2mu2);
    double L_Q = log(1. / m2Q2) + L_M;

    double H0 = H_0(z);
    double H1 = H_1(z);
    double Hm1 = H_m1(z);
    double H01 = H_01(z);
    double H0m1 = H_0m1(z);

    return TR * TR * (-64. / 3 * (z - 1.) * z * L_M)
           + CA * TR
                 * (-32. * (z - 1.) * (53. * z2 + 2. * z - 1.) / (9. * z)
                    - 32. * (13. * z2 - 8. * z - 1.) * H0
                    - 32. / 3 * (z - 1.) * (29. * z2 + 2. * z - 1.) * H1 / z
                    + L_Q
                          * (32. / 3 * (z - 1.) * (17. * z2 + 2. * z - 1.) / z
                             - 128. * z * H0 + 64. * (z - 1.) * z * H1)
                    + (z - 1.) * z * (-32. * H1 * H1 - 64. * H0 * H1)
                    + z * (z + 1.) * (64. * Hm1 * H0 - 64. * H0m1)
                    + z * (96. * H0 * H0 + 128. * H01)
                    + 64. * (z - 2.) * z * zeta2)
           + CF * TR
                 * (-64. / 15 * z * (3. * z2 + 5.) * H0 * H0
                    + 16. * (36. * z3 - 78. * z2 - 13. * z - 4.) * H0 / 15. / z
                    + 32. * (z - 1.) * (63. * z2 + 6. * z - 2.) / 15. / z
                    + L_Q * (32. * z * H0 - 16. * (z - 1.) * (2. * z + 1.))
                    + 16. * (z - 1.) * (4. * z + 1.) * H1
                    + (z + 1.) * (6. * z4 - 6. * z3 + z2 - z + 1.) / z2
                          * (64. / 15 * Hm1 * H0 - 64. / 15 * H0m1)
                    - 32. * z * H01
                    + (16. * (z - 1.) * (2. * z + 1.) - 32. * z * H0) * L_M
                    + 32. / 15 * z * (12. * z2 + 5.) * zeta2);
}

//==========================================================================================//
//  High scale (Q^2 >> m^2) limit of the quark coefficient functions for FL at
//  O(as^2) expanded in terms of \alpha_s^{[nf+1]}
//------------------------------------------------------------------------------------------//

double HighScaleCoefficientFunction::DL_ps2_highscale(
    double z, double m2Q2, double m2mu2
) const {

    double z2 = z * z;

    double L_M = log(m2mu2);
    double L_Q = log(1. / m2Q2) + L_M;

    double H0 = H_0(z);
    double H1 = H_1(z);
    double H01 = H_01(z);

    return CF * TR
           * (32. / 9. / z * (z - 1.) * (10. * z2 - 2. * z + 1.)
              - 32. * (z + 1.) * (2. * z - 1.) * H0
              - 32. / 3. / z * (z - 1.) * (2. * z2 + 2. * z - 1.) * H1
              + L_Q
                    * (32. * (z - 1.) * (2. * z2 + 2. * z - 1.) / 3. / z
                       - 32. * z * H0)
              + z * (32. * H0 * H0 + 32. * H01 - 32. * zeta2));
}

//==========================================================================================//
//  High scale (Q^2 >> m^2) limit of the gluon coefficient functions for F2 at
//  O(as^3) expanded in terms of \alpha_s^{[nf]}
//
//  Eq. (B.8) of Ref. [arXiv:1205.5727].
//------------------------------------------------------------------------------------------//

Value HighScaleCoefficientFunction::C2_g3_highscale(
    double x, double m2Q2, double m2mu2, int nf
) const {

    double Lmu = log(m2mu2);
    double L2mu = Lmu * Lmu;

    return D2_g3_highscale(x, m2Q2, m2mu2, nf)
           - 4. / 3 * Lmu * D2_g2_highscale(x, m2Q2, m2mu2)
           - ((16. / 9 * CA - 15. / 2 * CF) + (10. / 3 * CA + 2 * CF) * Lmu
              - 4. / 9 * L2mu)
                 * D2_g1_highscale(x, m2Q2);
}

//==========================================================================================//
//  High scale (Q^2 >> m^2) limit of the quark coefficient functions for F2 at
//  O(as^3) expanded in terms of \alpha_s^{[nf]}
//
//  Eq. (B.11) of Ref. [arXiv:1205.5727].
//------------------------------------------------------------------------------------------//

Value HighScaleCoefficientFunction::C2_ps3_highscale(
    double x, double m2Q2, double m2mu2, int nf
) const {

    double Lmu = log(m2mu2);

    return D2_ps3_highscale(x, m2Q2, m2mu2, nf)
           - 4. / 3 * Lmu * D2_ps2_highscale(x, m2Q2, m2mu2);
}

//==========================================================================================//
//  High scale (Q^2 >> m^2) limit of the gluon coefficient functions for FL at
//  O(as^3) expanded in terms of \alpha_s^{[nf+1]}
//------------------------------------------------------------------------------------------//

double HighScaleCoefficientFunction::DL_g3_highscale(
    double z, double m2Q2, double m2mu2, int nf
) const {

    double z2 = z * z;
    double z3 = z2 * z;
    double z4 = z3 * z;
    double z5 = z4 * z;

    double L_M = log(m2mu2);
    double L_M2 = L_M * L_M;
    double L_Q = log(1. / m2Q2) + L_M;
    double L_Q2 = L_Q * L_Q;

    // Allocate pointers for the harmonic polylogs
    double wx = z;
    int nw = 4;
    int n1 = -1;
    int n2 = 1;
    int sz = n2 - n1 + 1;
    double *Hr1 = new double[sz];
    double *Hr2 = new double[sz * sz];
    double *Hr3 = new double[sz * sz * sz];
    double *Hr4 = new double[sz * sz * sz * sz];
    double *Hr5 = new double[sz * sz * sz * sz * sz];

    // Call polylogs
    apf_hplog_(&wx, &nw, Hr1, Hr2, Hr3, Hr4, Hr5, &n1, &n2);

    // weight 1
    const double Hm1 = Hr1[0];
    const double H0 = Hr1[1];
    const double H1 = Hr1[2];

    // weight 2
    const double H0m1 = Hr2[1];
    const double H01 = Hr2[7];

    // weight 3
    const double H0m1m1 = Hr3[1];
    const double H00m1 = Hr3[4];
    const double H01m1 = Hr3[7];
    const double H0m11 = Hr3[19];
    const double H001 = Hr3[22];
    const double H011 = Hr3[25];

    // weight 4
    const double H000m1 = Hr4[13];
    const double H0001 = Hr4[67];
    const double H0011 = Hr4[76];
    const double H0111 = Hr4[79];

    delete[] Hr1;
    delete[] Hr2;
    delete[] Hr3;
    delete[] Hr4;
    delete[] Hr5;

    return -TR * TR * TR * 256. / 9 * (z - 1) * z * L_M2
           + CA * TR * TR
                 * ((64. * (z - 1) * (17. * z2 + 2. * z - 1.) / (9. * z)
                     - 256. / 3 * z * H0 + 128. / 3 * (z - 1.) * z * H1)
                        * L_Q2
                    + (-64. * (z - 1.) * (461. * z2 + 11. * z - 25.) / (27. * z)
                       - 128. / 9 * z * (26. * z - 59.) * H0
                       - 128. * (z - 1.) * (39. * z2 + 2. * z - 1.) * H1
                             / (9. * z)
                       + (z - 1.) * z
                             * (-128. / 3 * H1 * H1 - 256. / 3 * H0 * H1)
                       + z * (z + 1.) * (256. / 3 * Hm1 * H0 - 256. / 3 * H0m1)
                       + z * (512. / 3 * H0 * H0 + 512. / 3 * H01)
                       + (128. * (z - 1.) * (17. * z2 + 2. * z - 1.) / (9. * z)
                          - 512. / 3 * z * H0 + 256. / 3 * (z - 1.) * z * H1)
                             * L_M
                       + 256. / 3 * (z - 2.) * z * zeta2)
                          * L_Q
                    - 32. / 9 * (28. * z - 3.) * H0 * H0
                    + (64. * (z - 1.) * (17. * z2 + 2. * z - 1.) / (9. * z)
                       - 256. / 3 * z * H0 + 128. / 3 * (z - 1.) * z * H1)
                          * L_M2
                    + 32. * (z - 1.) * (2714. * z2 - 106. * z - 139.)
                          / (81. * z)
                    - 64. / 27 * (110. * z2 + 277. * z - 33.) * H0
                    + 4160. / 27 * (z - 1.) * z * H1
                    + z
                          * (-64. / 9 * H0 * H0 * H0 - 64. / 3 * H01
                             + 64. * zeta2 / 3.)
                    + L_M
                          * (64. * (z - 1.) * (68. * z2 + z - 7.) / (9. * z)
                             - 128. / 9 * (4. * z - 1.) * (13. * z + 6.) * H0
                             - 128. * (z - 1.) * (19. * z2 + 2. * z - 1.) * H1
                                   / (9. * z)
                             + (z - 1.) * z
                                   * (-128. / 3 * H1 * H1 - 256. / 3 * H0 * H1)
                             + z * (z + 1.)
                                   * (256. / 3 * Hm1 * H0 - 256. / 3 * H0m1)
                             + z * (256. / 3 * H0 * H0 + 512. / 3 * H01)
                             + 256. / 3 * (z - 2.) * z * zeta2))
           + CA * TR * TR * nf
                 * ((64. * (z - 1.) * (17. * z2 + 2. * z - 1.) / (9. * z)
                     - 256. / 3 * z * H0 + 128. / 3 * (z - 1.) * z * H1)
                        * L_Q2
                    + (-64. * (z - 1.) * (461. * z2 + 11. * z - 25.) / (27. * z)
                       - 128. / 9 * z * (26. * z - 59.) * H0
                       - 128. * (z - 1.) * (39. * z2 + 2. * z - 1.) * H1
                             / (9. * z)
                       + (z - 1.) * z
                             * (-128. / 3 * H1 * H1 - 256. / 3 * H0 * H1)
                       + z * (z + 1.) * (256. / 3 * Hm1 * H0 - 256. / 3 * H0m1)
                       + z * (512. / 3 * H0 * H0 + 512. / 3 * H01)
                       + 256. / 3 * (z - 2.) * z * zeta2)
                          * L_Q)
           + CA * CA * TR
                 * ((-16. * (z - 1.) * (1033. * z2 - 26. * z - 65.) / (9. * z)
                     - 64. * (6. * z3 - 47. * z2 + 3. * z + 1.) * H0 / (3. * z)
                     - 32. * (z - 1.) * (79. * z2 + 8. * z - 4.) * H1 / (3. * z)
                     + (z - 1.) * z * (-128. * H1 * H1 - 128. * H0 * H1)
                     + 128. * z * (z + 3.) * H01
                     + z * (256. * H0 * H0 - 512. * zeta2))
                        * L_Q2
                    + (+(64. / 3 * (18. * z2 - 91. * z + 6.) * H0 * H0
                         + 32. * (2713. * z3 - 1405. * z2 - 60. * z + 4.) * H0
                               / (9. * z)
                         + 64. * (z - 1.) * (137. * z2 + 12. * z - 6.) * H1 * H0
                               / (3. * z)
                         + 128. * z * (3. * z - 5.) * H0m1 * H0
                         - 128. * z * (z + 11.) * H01 * H0
                         + 32. * (z - 1.) * (161. * z2 + 12. * z - 6.) * H1 * H1
                               / (3. * z)
                         + 32. * (z - 1.) * (680. * z2 - 60. * z - 13.)
                               / (9. * z)
                         + 32. * (z - 1.) * (1919. * z2 + 30. * z - 93.) * H1
                               / (9. * z)
                         + (z + 1.) * (79. * z2 - 8. * z - 4.)
                               * (64. / 3 * H0m1 - 64. / 3 * Hm1 * H0) / z
                         - 128. * (z3 + 53. * z2 - 6. * z + 1.) * H01 / (3. * z)
                         - 128. * z * (3. * z - 13.) * H00m1
                         - 128. * (z - 3.) * z * H001
                         - 256. * z * (z + 5.) * H011
                         - 64. / 3 * (135. * z2 - 160. * z - 6.) * zeta2
                         + z * (z + 1.)
                               * (256. * H0 * Hm1 * Hm1
                                  + (-192. * H0 * H0 - 512. * H0m1 - 256. * H01)
                                        * Hm1
                                  + 512. * zeta2 * Hm1 + 512. * H0m1m1
                                  + 256. * H0m11 + 256. * H01m1)
                         + (z - 1.) * z
                               * (128. * H1 * H1 * H1 + 384. * H0 * H1 * H1
                                  + 192. * H0 * H0 * H1 - 512. * zeta2 * H1)
                         + z
                               * (-1024. / 3 * H0 * H0 * H0 + 1792. * zeta2 * H0
                                  - 128. * zeta3))
                       * L_Q))
           + CF * CF * TR
                 * (-8. / 3 * (4 * z2 - 4. * z - 1.) * H0 * H0 * H0
                    - 4. * (20. * z2 - 11. * z - 1.) * H0 * H0
                    + 8. * (24. * z2 + 37. * z - 7.) * H0
                    - 16. * (z - 1.) * (10. * z - 1.) * H1 * H0
                    + 32. * (2. * z2 + 5. * z - 2.) * H01 * H0
                    + 48. * (z - 1.) * z * H1 * H1
                    - 16. * (z - 1.) * (20. * z + 3.)
                    + 32. * (z - 1.) * (6. * z - 1.) * H1
                    + 16. * (16. * z2 - 24. * z + 3.) * H01
                    - 32. * (2. * z2 + 17. * z - 3.) * H001
                    + (z - 1.) * (2. * z + 1.)
                          * (16. / 3 * H1 * H1 * H1 - 16. * H0 * H0 * H1
                             + 32. * H011)
                    - 16. * (z - 2.) * (6. * z - 1) * zeta2
                    + L_Q2
                          * (32. * (2. * z + 1.) * H1 * (z - 1.)
                             + 24. * (z - 1.)
                             + 16. * (2. * z - 1.) * (2. * z + 1.) * H0
                             + z * (-16. * H0 * H0 - 64. * H01 + 64. * zeta2))
                    + L_M2
                          * (32. * (2. * z + 1.) * H1 * (z - 1.)
                             + 24. * (z - 1.)
                             + 16. * (2. * z - 1.) * (2. * z + 1.) * H0
                             + z * (-16. * H0 * H0 - 64. * H01 + 64. * zeta2))
                    - 32. * (2. * z2 - 19. * z + 1.) * zeta3
                    + z
                          * (4. / 3 * H0 * H0 * H0 * H0 + 32. * H01 * H0 * H0
                             - 192. * H001 * H0 + 192. * zeta2 * H0
                             - 64. * zeta3 * H0 - 608. * zeta2 * zeta2 / 5.
                             + 384. * H0001 - 64. * H0011 - 64. * H0111)
                    + L_M
                          * (32. / 15 * (24. * z3 + 90. * z2 - 95. * z - 15.)
                                 * H0 * H0
                             + 32. * (78. * z3 + 141. * z2 - 34. * z + 8.) * H0
                                   / (15. * z)
                             + 128. * (2. * z2 - 3. * z - 1.) * H01 * H0
                             - 8. * (z - 1.) * (6. * z + 1.) * (153. * z - 32.)
                                   / (15. * z)
                             + 16. * (z - 1.) * (4. * z - 3.) * H1
                             + ((z + 1.)
                                * (12. * z4 + 3. * z3 - 73. * z2 - 2. * z + 2.)
                                * (128. / 15 * H0m1 - 128. / 15 * Hm1 * H0))
                                   / z2
                             + 32. * (6. * z + 1.) * H01
                             - 64. * (4. * z2 - 5. * z - 2.) * H001
                             - 32. / 15
                                   * (48. * z3 + 120. * z2 - 250. * z - 45.)
                                   * zeta2
                             + (z + 1.) * (2. * z - 1.)
                                   * (128. * H0 * Hm1 * Hm1
                                      + (-64. * H0 * H0 - 256. * H0m1) * Hm1
                                      + 128. * zeta2 * Hm1 - 128. * H0 * H0m1
                                      + 256. * H0m1m1 + 384. * H00m1)
                             + (z - 1.) * (2. * z + 1.)
                                   * (-64. * H1 * H0 * H0 + 128. * H1 * H0
                                      + 64. * H1 * H1 + 128. * H1 * zeta2)
                             + z
                                   * (32. / 3 * H0 * H0 * H0
                                      + (128. * H01 - 128. * H0m1) * H0 * H0
                                      + (512. * H0m1m1 + 512. * H00m1
                                         - 512. * H001)
                                            * H0
                                      - 256. * H0m1 * H0m1
                                      + 768. * zeta2 * zeta2 / 5. - 256. * H011
                                      - 768. * H000m1 + 768. * H0001
                                      + (64. * H0 + 256. * H0m1 - 256. * H01)
                                            * zeta2
                                      + (-512. * H0 - 576.) * zeta3))
                    + L_Q
                          * (-32. / 15 * (24. * z3 + 90. * z2 - 95. * z - 15.)
                                 * H0 * H0
                             - 32. * (78. * z3 + 141. * z2 - 34. * z + 8.) * H0
                                   / (15. * z)
                             - 128. * (2. * z2 - 3. * z - 1.) * H01 * H0
                             + 8. * (z - 1.) * (6. * z + 1.) * (153. * z - 32.)
                                   / (15. * z)
                             - 16. * (z - 1.) * (4. * z - 3.) * H1
                             + (z + 1.
                               ) * (12. * z4 + 3 * z3 - 73. * z2 - 2. * z + 2.)
                                   * (128. / 15 * Hm1 * H0 - 128. / 15 * H0m1)
                                   / z2
                             - 32. * (6. * z + 1) * H01
                             + 64. * (4. * z2 - 5. * z - 2.) * H001
                             + L_M
                                   * (-64. * (2. * z + 1.) * H1 * (z - 1.)
                                      - 48. * (z - 1.)
                                      - 32. * (2. * z - 1.) * (2. * z + 1.) * H0
                                      + z
                                            * (32. * H0 * H0 + 128. * H01
                                               - 128. * zeta2))
                             + 32. / 15
                                   * (48. * z3 + 120. * z2 - 250. * z - 45.)
                                   * zeta2
                             + (z + 1.) * (2. * z - 1.)
                                   * (-128. * H0 * Hm1 * Hm1
                                      + (64. * H0 * H0 + 256. * H0m1) * Hm1
                                      - 128. * zeta2 * Hm1 + 128. * H0 * H0m1
                                      - 256. * H0m1m1 - 384. * H00m1)
                             + (z - 1.) * (2. * z + 1.)
                                   * (64. * H1 * H0 * H0 - 128. * H1 * H0
                                      - 64. * H1 * H1 - 128. * H1 * zeta2)
                             + z
                                   * (-32. / 3 * H0 * H0 * H0
                                      + (128. * H0m1 - 128. * H01) * H0 * H0
                                      + (-512. * H0m1m1 - 512. * H00m1
                                         + 512. * H001)
                                            * H0
                                      + 256. * H0m1 * H0m1
                                      - (768. * zeta2 * zeta2) / 5.
                                      + 256. * H011 + 768. * H000m1
                                      - 768. * H0001
                                      + (-64. * H0 - 256. * H0m1 + 256. * H01)
                                            * zeta2
                                      + (512. * H0 + 576.) * zeta3)))
           + CF * TR * TR
                 * (-16. / 3 * z * H0 * H0 * H0 * H0
                    - 32. / 3 * (7. * z - 1.) * H0 * H0 * H0
                    - 32. * (19. * z - 3.) * H0 * H0
                    - 64. * (6. * z2 + 7. * z - 8.) * H0
                    + (-64. * z * H0 * H0
                       - 64. / 3 * (4. * z2 + 5. * z - 3.) * H0
                       + 64. * (z - 1.) * (40. * z2 - 17. * z - 2.) / (9. * z))
                          * L_M2
                    + 16. * (z - 1.) * (343. * z2 - 242. * z + 4.) / (3. * z)
                    + L_Q2
                          * (-64. * z * H0 * H0
                             - 64. / 3 * (z + 1.) * (4. * z - 3.) * H0
                             + 64. * (z - 1.) * (28. * z2 - 23. * z - 2.)
                                   / (9. * z))
                    + L_Q
                          * (-128. * (25. * z + 2.) * H1 * (z - 1.) * (z - 1.)
                                 / (9. * z)
                             - 32. * (2474. * z2 - 4897. * z + 44.) * (z - 1.)
                                   / (135. * z)
                             - 64. / 45
                                   * (12. * z3 - 180. * z2 - 265. * z + 90.)
                                   * H0 * H0
                             - 64. * (354. * z3 - 397. * z2 + 388. * z + 4.)
                                   * H0 / (45. * z)
                             + (z + 1.) * (6. * z4 - 6. * z3 + z2 - z + 1.) / z2
                                   * (256. / 45 * Hm1 * H0 - 256. / 45 * H0m1)
                             + 128. / 3 * (z + 1.) * (4. * z - 3.) * H01
                             + (128. * z * H0 * H0
                                + 128. / 3 * (4. * z2 + 3. * z - 3.) * H0
                                - 256. * (z - 1.) * (17. * z2 - 10. * z - 1.)
                                      / (9. * z))
                                   * L_M
                             + 128. / 45 * (12. * z3 - 60. * z2 - 25. * z + 45.)
                                   * zeta2
                             + z
                                   * (128. * H0 * H0 * H0 - 256. * zeta2 * H0
                                      + 256. * H001 - 256. * zeta3))
                    + L_M
                          * (-64. / 45 * (12. * z3 + 180. * z2 + 305. * z - 90.)
                                 * H0 * H0
                             + 64. * (426. * z3 - 553. * z2 + 362. * z - 4.)
                                   * H0 / (45. * z)
                             + 32. * (z - 1.) * (3716. * z2 - 4753. * z - 4.)
                                   / (135. * z)
                             + 128. * (z - 1.) * (37. * z2 - 20. * z - 2.) * H1
                                   / (9. * z)
                             + (z + 1.) * (6. * z4 - 6. * z3 + z2 - z + 1.)
                                   * (256. / 45 * Hm1 * H0 - 256. / 45 * H0m1)
                                   / z2
                             - 128. / 3 * (4. * z2 + 3. * z - 3.) * H01
                             + 128. / 45 * (12. * z3 + 60. * z2 + 35. * z - 45.)
                                   * zeta2
                             + z
                                   * (-128. * H0 * H0 * H0 + 256. * zeta2 * H0
                                      - 256. * H001 + 256. * zeta3)))
           + CF * nf * TR * TR
                 * ((-64. * z * H0 * H0
                     - 64. / 3 * (z + 1.) * (4. * z - 3.) * H0
                     + 64. * (z - 1.) * (28. * z2 - 23. * z - 2.) / (9. * z))
                        * L_Q2
                    + (-128. * (25. * z + 2.) * H1 * (z - 1.) * (z - 1.)
                           / (9. * z)
                       - 32. * (2474. * z2 - 4897. * z + 44.) * (z - 1.)
                             / (135. * z)
                       - 64. / 45 * (12. * z3 - 180. * z2 - 265. * z + 90.) * H0
                             * H0
                       - 64. * (354. * z3 - 397. * z2 + 388. * z + 4.) * H0
                             / (45. * z)
                       + (z + 1.) * (6. * z4 - 6. * z3 + z2 - z + 1.) / z2
                             * (256. / 45 * Hm1 * H0 - 256. / 45 * H0m1)
                       + 128. / 3 * (z + 1.) * (4. * z - 3.) * H01
                       + (64. / 3 * (z - 1.) * (2. * z + 1.) - 128. / 3 * z * H0
                         ) * L_M
                       + 128. / 45 * (12. * z3 - 60. * z2 - 25. * z + 45.)
                             * zeta2
                       + z
                             * (128. * H0 * H0 * H0 - 256. * zeta2 * H0
                                + 256. * H001 - 256. * zeta3))
                          * L_Q
                    + L_M
                          * (-32. / 9 * (z - 1.) * (68. * z + 25.)
                             - 64. / 9 * (6. * z2 - 31. * z - 3.) * H0
                             - 64. / 3 * (z - 1.) * (2. * z + 1.) * H1
                             + z
                                   * (128. / 3 * H0 * H0 + 128. / 3 * H01
                                      - 128. * zeta2 / 3.)))
           + CA * CF * TR
                 * (-16. / 3 * (3. * z + 1.) * H0 * H0 * H0
                    - 8. / 3 * (11. * z2 - 18. * z + 3.) * H0 * H0
                    + 16. / 9 * (772. * z2 + 480. * z - 39.) * H0
                    + 16. * (z - 1.) * (77. * z2 - 25. * z - 4.) * H1 * H0
                          / (3. * z)
                    - 32. * (7. * z - 2.) * H01 * H0
                    - 8. * (z - 1.) * (9. * z - 1.) * H1 * H1
                    - 32. * (z - 1.) * (2168. * z2 - 91. * z - 28.) / (27. * z)
                    - 16. * (z - 1.) * (16. * z - 1.) * H1
                    + (z - 1.) * (2. * z + 1.)
                          * (16. * H0 * H1 * H1 - 16. / 3 * H1 * H1 * H1)
                    + z * (z + 1.) * (192. * H0m1 - 192. * Hm1 * H0)
                    - 16. * (68. * z3 - 117. * z2 + 21. * z + 4.) * H01
                          / (3. * z)
                    - 32. * (2. * z2 - 2. * z - 1.) * H011
                    + L_M2
                          * (-16. * (z - 1.) * (43. * z2 - 11. * z - 2.)
                                 / (3. * z)
                             + 32. * (3. * z - 1.) * H0
                             - 32. * (z - 1.) * (2. * z + 1.) * H1
                             + z * (64. * H0 * H0 + 64. * H01 - 64. * zeta2))
                    - 16. * z * (3. * z + 17.) * zeta2
                    + (z + 1.) * (2. * z - 1.)
                          * (32. * H0 * Hm1 * Hm1
                             + (-16. * H0 * H0 - 64. * H0m1) * Hm1
                             + 32. * zeta2 * Hm1 + 32. * H0 * H0m1
                             + 64. * H0m1m1 - 32. * H00m1)
                    + L_Q2
                          * (16. * (z - 1.) * (65. * z2 - 2.) / (3. * z)
                             - 32. / 3 * (20. * z - 3.) * H0
                             + 32. * (z - 1.) * (2. * z + 1.) * H1
                             + z * (-64. * H0 * H0 - 64. * H01 + 64. * zeta2))
                    + (15. * z - 4.) * (32. * H001 - 32. * zeta3)
                    + z
                          * (8. / 3 * H0 * H0 * H0 * H0 - 32. * H0m1 * H0 * H0
                             + (128. * H0m1m1 + 128. * H00m1 - 256. * H001
                                - 64. * H011)
                                   * H0
                             - 448. * zeta3 * H0 - 64. * H0m1 * H0m1
                             - 1472. * zeta2 * zeta2 / 5. - 192. * H000m1
                             + 768. * H0001 + 128. * H0011 + 64. * H0111
                             + (64. * H0m1 - 32. * H0) * zeta2)
                    + L_Q
                          * (128. / 45 * z * (6. * z2 + 35.) * H0 * H0 * H0
                             + 64.
                                   * (84. * z4 - 9. * z3 + 272. * z2 - 48. * z
                                      + 6.)
                                   * H0 * H0 / (45. * z)
                             - 32. * (z + 1.)
                                   * (24. * z4 + 6. * z3 - 11. * z2 - 4. * z
                                      + 4.)
                                   * Hm1 * H0 * H0 / (15. * z2)
                             - 32. * (z - 1.) * (2. * z + 1.) * H1 * H0 * H0
                             - 16. * (4668. * z3 - 5233. * z2 - 130. * z - 64.)
                                   * H0 / (45. * z)
                             - 64. * (z - 1.) * (4. * z + 1.) * H1 * H0
                             - 64. * (30. * z4 - 35. * z3 - 15. * z2 - 4.)
                                   * H0m1 * H0 / (15. * z2)
                             + 64. * (z + 1.) * (2. * z - 1.) * H01 * H0
                             - 256. / 5 * z * (4. * z2 + 5.) * zeta2 * H0
                             - 32. * (z - 1.) * (6. * z + 1.) * H1 * H1
                             - 8. * (z - 1.) * (11062. * z2 + 1335. * z - 168.)
                                   / (45. * z)
                             - 64. / 45
                                   * (168. * z4 - 108. * z3 + 343. * z2 - 6. * z
                                      + 24.)
                                   * zeta2 / z
                             + 64. / 5 * (z + 1.)
                                   * (12. * z4 - 2. * z3 - 3. * z2 - 2. * z + 2.
                                   )
                                   * zeta2 / z2 * Hm1
                             - 64. * (z - 1.) * (371. * z2 + 37. * z - 14.) * H1
                                   / (15. * z)
                             - 64. / 15 * (z - 1.)
                                   * (12. * z4 - 18. * z3 - 13. * z2 + 2. * z
                                      + 2.)
                                   * zeta2 / z2 * H1
                             + (z + 1.)
                                   * (42. * z4 - 69. * z3 - 35. * z2 - 4. * z
                                      + 7.)
                                   * (256. / 45 * H0m1 - 256. / 45 * Hm1 * H0)
                                   / z2
                             + 64. * (24. * z3 + 208. * z2 - 17. * z + 4.) * H01
                                   / (15. * z)
                             + (z + 1.)
                                   * (12. * z4 + 18. * z3 - 13. * z2 - 2. * z
                                      + 2.)
                                   / z2
                                   * (64. / 15 * H0 * Hm1 * Hm1
                                      - 128. / 15 * H0m1 * Hm1
                                      + 128. / 15 * H0m1m1)
                             + 64.
                                   * (24. * z5 + 90. * z4 - 75. * z3 - 45. * z2
                                      - 4.)
                                   * H00m1 / (15. * z2)
                             + 64. / 15 * (24. * z3 - 30. * z2 + 55. * z + 15.)
                                   * H001
                             + (z + 1.) * (6. * z4 - 6. * z3 + z2 - z + 1.) / z2
                                   * (-256. / 15 * Hm1 * H01 + 256. / 15 * H0m11
                                      + 256. / 15 * H01m1)
                             + (352. / 3 * z * H0
                                - 176. / 3 * (z - 1.) * (2. * z + 1.))
                                   * L_M
                             - 256. / 3 * z * (3. * z2 + 2.) * zeta3
                             + z
                                   * ((64. * H01 - 64. * H0m1) * H0 * H0
                                      + (256. * H0m1m1 + 256. * H00m1
                                         - 256. * H001)
                                            * H0
                                      - 256. * zeta3 * H0 - 128. * H0m1 * H0m1
                                      + 384. * zeta2 * zeta2 / 5. + 128. * H011
                                      - 384. * H000m1 + 384. * H0001
                                      + (128. * H0m1 - 128. * H01) * zeta2))
                    + L_M
                          * (-32. / 5 * (4. * z3 + 5. * z2 - 10. * z - 5.) * H0
                                 * H0
                             + 16. * (2226. * z3 - 43. * z2 - 63. * z - 24.)
                                   * H0 / (45. * z)
                             - 8. * (z - 1.) * (3758. * z2 - 1299. * z - 152.)
                                   / (45. * z)
                             - 16. / 3 * (z - 1.) * (2. * z - 23.) * H1
                             + ((z + 1.)
                                * (12. * z4 - 27. * z3 - 58. * z2 - 2. * z + 2.)
                               ) / z2
                                   * (64. / 15 * Hm1 * H0 - 64. / 15 * H0m1)
                             - 64. * (6. * z2 - z - 3.) * H00m1
                             + 32. / 15 * z * (24. * z2 - 85.) * zeta2
                             + (z + 1.) * (2. * z - 1.)
                                   * (-64. * H0 * Hm1 * Hm1
                                      + (32. * H0 * H0 + 128. * H0m1) * Hm1
                                      - 64. * zeta2 * Hm1 - 128. * H0m1m1)
                             + (z - 1.) * (2. * z + 1.)
                                   * (32. * H1 * H0 * H0
                                      + (64. * H0m1 - 64. * H01) * H0
                                      - 32. * H1 * H1 + 64. * H001
                                      - 64. * H1 * zeta2)
                             + z
                                   * (-128. / 3 * H0 * H0 * H0
                                      + (64. * H0m1 - 64. * H01) * H0 * H0
                                      + (-256. * H0m1m1 - 256. * H00m1
                                         + 256. * H001)
                                            * H0
                                      + 256. * zeta3 * H0 + 128. * H0m1 * H0m1
                                      - 384. * zeta2 * zeta2 / 5.
                                      - 544. / 3 * H01 + 128. * H011
                                      + 384. * H000m1 - 384. * H0001
                                      + (128. * H01 - 128. * H0m1) * zeta2)))
           + massless_as3_->MuIndependentTerms(z, nf + 1) / (nf + 1.);
}

//==========================================================================================//
//  High scale (Q^2 >> m^2) limit of the gluon coefficient functions for FL at
//  O(as^3) expanded in terms of \alpha_s^{[nf]}
//------------------------------------------------------------------------------------------//

Value HighScaleCoefficientFunction::CL_g3_highscale(
    double x, double m2Q2, double m2mu2, int nf
) const {

    double Lmu = log(m2mu2);
    double L2mu = Lmu * Lmu;

    double tmp = DL_g3_highscale(x, m2Q2, m2mu2, nf)
                 - 4. / 3 * Lmu * DL_g2_highscale(x, m2Q2, m2mu2)
                 - ((16. / 9 * CA - 15. / 2 * CF)
                    + (10. / 3 * CA + 2 * CF) * Lmu - 4. / 9 * L2mu)
                       * DL_g1_highscale(x);
    return Value(tmp);
}

//==========================================================================================//
//  High scale (Q^2 >> m^2) limit of the quark coefficient functions for FL at
//  O(as^3) expanded in terms of \alpha_s^{[nf+1]}
//------------------------------------------------------------------------------------------//

double HighScaleCoefficientFunction::DL_ps3_highscale(
    double z, double m2Q2, double m2mu2, int nf
) const {

    double z2 = z * z;
    double z3 = z2 * z;
    double z4 = z3 * z;

    double L_M = log(m2mu2);
    double L_M2 = L_M * L_M;
    double L_Q = log(1. / m2Q2) + L_M;
    double L_Q2 = L_Q * L_Q;

    // Allocate pointers for the harmonic polylogs
    double wx = z;
    int nw = 4;
    int n1 = -1;
    int n2 = 1;
    int sz = n2 - n1 + 1;
    double *Hr1 = new double[sz];
    double *Hr2 = new double[sz * sz];
    double *Hr3 = new double[sz * sz * sz];
    double *Hr4 = new double[sz * sz * sz * sz];
    double *Hr5 = new double[sz * sz * sz * sz * sz];

    // Call polylogs
    apf_hplog_(&wx, &nw, Hr1, Hr2, Hr3, Hr4, Hr5, &n1, &n2);

    // weight 1
    const double Hm1 = Hr1[0];
    const double H0 = Hr1[1];
    const double H1 = Hr1[2];

    // weight 2
    const double H0m1 = Hr2[1];
    const double H01 = Hr2[7];

    // weight 3
    const double H00m1 = Hr3[4];
    const double H001 = Hr3[22];
    const double H011 = Hr3[25];

    // weight 4
    const double H0001 = Hr4[67];

    delete[] Hr1;
    delete[] Hr2;
    delete[] Hr3;
    delete[] Hr4;
    delete[] Hr5;

    return CF * CF * TR
               * (-8. / 3. * (5. * z + 2.) * H0 * H0 * H0
                  - 8. / 3. * (8. * z2 + 3.) * H0 * H0
                  + 16. / 9. * (160. * z2 + 93. * z - 39.) * H0
                  + (16. * z * H0 * H0 - 16. * (z + 2.) * H0
                     - 16. * (z - 1.) * (4. * z2 - 11. * z - 2.) / (3. * z))
                        * L_M2
                  - 32. * (z - 1) * (440. * z2 - 91. * z - 28.) / (27. * z)
                  + (z - 1.) * (4. * z2 - 11. * z - 2.)
                        * (32. / 3 * H0 * H1 - 32. / 3. * H01) / z
                  + (-32. / 3 * z * H0 * H0 * H0 + 16. * (5. * z + 2.) * H0 * H0
                     + 32. / 3 * (8. * z2 + 18. * z + 3.) * H0
                     - 32. * (z - 1.) * (80. * z2 + 17. * z - 10.) / (9. * z))
                        * L_M
                  + L_Q2
                        * (-64. / 3 * (2. * z2 - 3.) * H0
                           + 32. * (z - 1) * (2. * z2 - 9. * z - 1.) / (3. * z)
                           - 64. * (z - 1.) * (2. * z2 + 2. * z - 1.) * H1
                                 / (3. * z)
                           + z * (64. * H01 - 64. * zeta2))
                  + (z + 2.) * (32. * H0 * H01 - 64. * H001 + 64. * zeta3)
                  + z
                        * (4. / 3 * H0 * H0 * H0 * H0 - 64. * H001 * H0
                           - 128. * zeta3 * H0 - 384. * zeta2 * zeta2 / 5.
                           + 192. * H0001)
                  + L_Q
                        * (-32. * (z - 1.) * (86. * z2 - 33. * z + 6.)
                               / (45. * z)
                           - 32. * (56. * z3 - 813. * z2 + 142. * z + 16.) * H0
                                 / (45. * z)
                           + 64. * (z - 1.) * (4. * z2 + 55. * z - 5.) * H1
                                 / (9. * z)
                           + (z - 1.) * (2. * z2 + 2. * z - 1.)
                                 * (64. / 3 * H1 * H1 + 128. / 3. * H0 * H1) / z
                           + ((z + 1.)
                              * (6. * z4 - 6. * z3 + 11. * z2 + 4. * z - 4.)
                              * (128. / 45 * H0m1 - 128. / 45 * Hm1 * H0))
                                 / z2
                           + (128. * (4. * z3 - 6. * z2 - 3. * z - 1.) * H01)
                                 / (3. * z)
                           + (6. * z3 + 90. * z2 - 85. * z - 90.)
                                 * (64. / 45 * H0 * H0 - (128. * zeta2) / 45.)
                           + z
                                 * (-64. / 9 * H0 * H0 * H0
                                    + (128. / 3 * H0m1 - 128. * H01) * H0
                                    + 896. / 3 * zeta2 * H0 - 256. / 3 * H00m1
                                    - 128. * H011 + 192. * zeta3)))
           + CF * TR * TR * nf
                 * ((128. * (z - 1.) * (2. * z2 + 2. * z - 1.) / (9. * z)
                     - 128. / 3 * z * H0)
                        * L_Q2
                    + (256. / 3 * z * H0 * H0
                       - 256. / 9 * (4. * z2 - 8. * z - 3.) * H0
                       - 256. * (z - 1.) * (3. * z2 + 6. * z - 2.) / (9. * z))
                          * L_Q)
           + CF * TR * TR
                 * (((128. * (z - 1.) * (2. * z2 + 2. * z - 1.)) / (9. * z)
                     - 128. / 3 * z * H0)
                        * L_Q2
                    + (256. / 3 * z * H0 * H0
                       - 256. / 9 * (4. * z2 - 8. * z - 3) * H0
                       - 256 * (z - 1.) * (3. * z2 + 6. * z - 2.) / (9. * z))
                          * L_Q
                    + 64. * (z - 1) * (2. * z2 + 2. * z - 1.) * H1 * H1
                          / (9. * z)
                    + (128. * (z - 1.) * (2. * z2 + 2. * z - 1.) / (9. * z)
                       - 128. / 3 * z * H0)
                          * L_M2
                    + (256. * (z - 1.) * (55. * z2 + 43. * z - 14.)) / (81. * z)
                    - 128. / 27 * z * (19. * z + 67.) * H0
                    - 128. * (z - 1.) * (19. * z2 + 16. * z - 5.) * H1
                          / (27. * z)
                    + L_M
                          * ((256. * (z - 1.) * (19. * z2 + 16. * z - 5.))
                                 / (27. * z)
                             - 256. / 9 * z * (2. * z + 11.) * H0
                             - (256. * (z - 1.) * (2. * z2 + 2. * z - 1.) * H1)
                                   / (9. * z)
                             + z * (256. / 3 * H01 - 256. * zeta2 / 3.))
                    + z * (2. * z + 11.) * (128. / 9 * H01 - 128. * zeta2 / 9.)
                    + z * (128. * zeta3 / 3. - 128. / 3 * H011))
           + CF * CA * TR
                 * ((-16. * (z - 1.) * (46. * z2 - z - 21.) / (3. * z)
                     + 64. * (7. * z2 - 3. * z - 1.) * H0 / (3. * z)
                     - 64. * (z - 1.) * (2. * z2 + 2. * z - 1.) * H1 / (3. * z)
                     + z * (64. * H0 * H0 + 64. * H01 - 64. * zeta2))
                        * L_Q2
                    + (-32. / 3 * (19. * z - 12.) * H0 * H0
                       + 32. * (422. * z3 - 137. * z2 - 114. * z + 4.) * H0
                             / (9. * z)
                       - 32. * (z - 1.) * (670. * z2 - 245. * z + 46.)
                             / (27. * z)
                       + 32. * (z - 1.) * (106. * z2 - 23. * z - 65.) * H1
                             / (9. * z)
                       + (z - 1.) * (2. * z2 + 2. * z - 1.)
                             * (128. / 3 * H1 * H1 + 256. / 3 * H0 * H1) / z
                       + (z + 1.) * (2. * z2 - 2. * z - 1.)
                             * (256. / 3 * H0m1 - 256. / 3 * Hm1 * H0) / z
                       - 64. * (z - 4.) * H01
                       - 64. / 3 * z * (8. * z - 3.) * zeta2
                       + z
                             * (-256. / 3 * H0 * H0 * H0
                                + (-256. * H0m1 - 256. * H01) * H0
                                + 256. * zeta2 * H0 + 512. * H00m1 - 256. * H011
                                - 128. * zeta3))
                          * L_Q)
           + massless_as3_->MuIndependentTerms(z, nf + 1) / (nf + 1.);
}

//==========================================================================================//
//  High scale (Q^2 >> m^2) limit of the quark coefficient functions for FL at
//  O(as^3) expanded in terms of \alpha_s^{[nf]}
//------------------------------------------------------------------------------------------//

Value HighScaleCoefficientFunction::CL_ps3_highscale(
    double x, double m2Q2, double m2mu2, int nf
) const {

    double Lmu = log(m2mu2);

    double tmp = DL_ps3_highscale(x, m2Q2, m2mu2, nf)
                 - 4. / 3 * Lmu * DL_ps2_highscale(x, m2Q2, m2mu2);

    return Value(tmp);
}

//==========================================================================================//
//  High scale (Q^2 >> m^2) limit of the gluon coefficient functions for F2 at
//  O(as^3) expanded in terms of \alpha_s^{[nf+1]}. This function uses the
//  approximation for the currrently unknown term a_Qg_30 given in Eq. (3.49) of
//  Ref. [arXiv:1205.5727] and in Eq. (16) Ref. of [arXiv:1701.05838].
//
//  Eq. (B.11) of Ref. [arXiv:1205.5727].
//------------------------------------------------------------------------------------------//

Value HighScaleCoefficientFunction::D2_g3_highscale(
    double x, double m2Q2, double m2mu2, int nf
) const {

    double x2 = x * x;
    double x3 = x2 * x;

    double LQm = log(1. / m2Q2);
    double LQm2 = LQm * LQm;
    double LQm3 = LQm2 * LQm;

    double Lmmu = log(m2mu2);
    double Lmmu2 = Lmmu * Lmmu;
    // double Lmmu3=Lmmu2*Lmmu;

    // Allocate pointers for the harmonic polylogs
    double wx = x;
    int nw = 5;
    int n1 = -1;
    int n2 = 1;
    int sz = n2 - n1 + 1;
    double *Hr1 = new double[sz];
    double *Hr2 = new double[sz * sz];
    double *Hr3 = new double[sz * sz * sz];
    double *Hr4 = new double[sz * sz * sz * sz];
    double *Hr5 = new double[sz * sz * sz * sz * sz];

    // Call polylogs
    apf_hplog_(&wx, &nw, Hr1, Hr2, Hr3, Hr4, Hr5, &n1, &n2);

    // weight 1
    const double Hm1 = Hr1[0];
    const double H0 = Hr1[1];
    const double H1 = Hr1[2];

    // weight 2
    const double Hm1m1 = Hr2[0];
    const double H0m1 = Hr2[1];
    const double Hm10 = Hr2[3];
    const double H00 = Hr2[4];
    const double H10 = Hr2[5];
    const double H01 = Hr2[7];
    const double H11 = Hr2[8];

    // weight 3
    const double Hm1m1m1 = Hr3[0];
    const double H0m1m1 = Hr3[1];
    const double Hm10m1 = Hr3[3];
    const double H00m1 = Hr3[4];
    const double Hm1m10 = Hr3[9];
    const double H0m10 = Hr3[10];
    const double Hm100 = Hr3[12];
    const double H000 = Hr3[13];
    const double H100 = Hr3[14];
    const double H010 = Hr3[16];
    const double H110 = Hr3[17];
    const double Hm101 = Hr3[21];
    const double H001 = Hr3[22];
    const double H101 = Hr3[23];
    const double H011 = Hr3[25];
    const double H111 = Hr3[26];

    // weight 4
    const double Hm1m1m10 = Hr4[27];
    const double H0m1m10 = Hr4[28];
    const double Hm10m10 = Hr4[30];
    const double H00m10 = Hr4[31];
    const double H10m10 = Hr4[32];
    const double Hm1m100 = Hr4[36];
    const double H0m100 = Hr4[37];
    const double Hm1000 = Hr4[39];
    const double H0000 = Hr4[40];
    const double H1000 = Hr4[41];
    const double H0100 = Hr4[43];
    const double H1100 = Hr4[44];
    const double Hm1010 = Hr4[48];
    const double H0010 = Hr4[49];
    const double H1010 = Hr4[50];
    const double H0110 = Hr4[52];
    const double H1110 = Hr4[53];
    const double Hm1m101 = Hr4[63];
    const double H0m101 = Hr4[64];
    const double Hm1001 = Hr4[66];
    const double H0001 = Hr4[67];
    const double H1001 = Hr4[68];
    const double H0101 = Hr4[70];
    const double H1101 = Hr4[71];
    const double Hm1011 = Hr4[75];
    const double H0011 = Hr4[76];
    const double H1011 = Hr4[77];
    const double H0111 = Hr4[79];
    const double H1111 = Hr4[80];

    //  weight 5
    const double Hm1m1m1m10 = Hr5[81];
    const double H0m1m1m10 = Hr5[82];
    const double Hm10m1m10 = Hr5[84];
    const double H00m1m10 = Hr5[85];
    const double Hm1m10m10 = Hr5[90];
    const double H0m10m10 = Hr5[91];
    const double Hm100m10 = Hr5[93];
    const double Hm1m1m100 = Hr5[108];
    const double H0m1m100 = Hr5[109];
    const double Hm10m100 = Hr5[111];
    const double H00m100 = Hr5[112];
    const double Hm1m1000 = Hr5[117];
    const double H0m1000 = Hr5[118];
    const double Hm10000 = Hr5[120];
    const double H00000 = Hr5[121];
    const double H10000 = Hr5[122];
    const double H01000 = Hr5[124];
    const double H11000 = Hr5[125];
    const double Hm10100 = Hr5[129];
    const double H00100 = Hr5[130];
    const double H10100 = Hr5[131];
    const double H01100 = Hr5[133];
    const double H11100 = Hr5[134];
    const double Hm1m1010 = Hr5[144];
    const double H0m1010 = Hr5[145];
    const double Hm10010 = Hr5[147];
    const double H00010 = Hr5[148];
    const double H10010 = Hr5[149];
    const double H01010 = Hr5[151];
    const double H11010 = Hr5[152];
    const double Hm10110 = Hr5[156];
    const double H00110 = Hr5[157];
    const double H10110 = Hr5[158];
    const double H01110 = Hr5[160];
    const double H11110 = Hr5[161];
    const double Hm1m1m101 = Hr5[189];
    const double H0m1m101 = Hr5[190];
    const double Hm10m101 = Hr5[192];
    const double Hm1m1001 = Hr5[198];
    const double H0m1001 = Hr5[199];
    const double Hm10001 = Hr5[201];
    const double H00001 = Hr5[202];
    const double H10001 = Hr5[203];
    const double H01001 = Hr5[205];
    const double H11001 = Hr5[206];
    const double Hm10101 = Hr5[210];
    const double H00101 = Hr5[211];
    const double H10101 = Hr5[212];
    const double H01101 = Hr5[214];
    const double H11101 = Hr5[215];
    const double Hm1m1011 = Hr5[225];
    const double Hm10011 = Hr5[228];
    const double H00011 = Hr5[229];
    const double H10011 = Hr5[230];
    const double H01011 = Hr5[232];
    const double H11011 = Hr5[233];
    const double H00111 = Hr5[238];
    const double H10111 = Hr5[239];
    const double H01111 = Hr5[241];
    const double H11111 = Hr5[242];

    delete[] Hr1;
    delete[] Hr2;
    delete[] Hr3;
    delete[] Hr4;
    delete[] Hr5;

    return +Lmmu2
               * (-8. / 9 + 64. / 9 * x - 64. / 9 * x2 - 8. / 9 * H0
                  + 16. / 9 * H0 * x - 16. / 9 * H0 * x2 - 8. / 9 * H1
                  + 16. / 9 * H1 * x - 16. / 9 * H1 * x2)
           + LQm * Lmmu2 * (8. / 9 - 16. / 9 * x + 16. / 9 * x2)
           + CF
                 * (313. / 6 - 104. / 9. / x - 1435. / 3 * x + 3497. / 9 * x2
                    - 8 * zeta5 + 16 * zeta5 * x + 32. / 27 * zeta3 / x
                    - 700. / 9 * zeta3 + 1504. / 9 * zeta3 * x
                    - 680. / 27 * zeta3 * x2 - 128. / 27 * zeta2 / x
                    - 4081. / 18 * zeta2 + 1286. / 9 * zeta2 * x
                    + 992. / 27 * zeta2 * x2 - 364. / 15 * zeta2 * zeta2
                    + 776. / 15 * zeta2 * zeta2 * x
                    - 872. / 15 * zeta2 * zeta2 * x2 + 7 * H0 - 526 * H0 * x
                    + 538. / 3 * H0 * x2 - 296. / 9 * H0 * zeta3
                    + 592. / 9 * H0 * zeta3 * x - 32 * H0 * zeta3 * x2
                    - 1063. / 9 * H0 * zeta2 + 1124. / 9 * H0 * zeta2 * x
                    + 8. / 3 * H0 * zeta2 * x2 - 16. / 5 * H0 * zeta2 * zeta2
                    + 32. / 5 * H0 * zeta2 * zeta2 * x + 88. / 3 * H00
                    - 100 * H00 * x - 272. / 3 * H00 * x2
                    - 32. / 3 * H00 * zeta3 + 64. / 3 * H00 * zeta3 * x
                    - 118. / 3 * H00 * zeta2 + 224. / 3 * H00 * zeta2 * x
                    + 64. / 3 * H00 * zeta2 * x2 + 32. / 3 * H000
                    - 168 * H000 * x + 80. / 3 * H000 * x2 - 20 * H000 * zeta2
                    + 40 * H000 * zeta2 * x + 8. / 3 * H0000
                    - 160. / 3 * H0000 * x + 32. / 3 * H0000 * x2
                    + 793. / 3 * H1 + 16. / 3 * H1 / x - 1558. / 3 * H1 * x
                    + 270 * H1 * x2 - 160. / 9 * H1 * zeta3
                    + 320. / 9 * H1 * zeta3 * x - 320. / 9 * H1 * zeta3 * x2
                    - 20. / 9 * H1 * zeta2 - 80. / 9 * H1 * zeta2 * x
                    + 80. / 9 * H1 * zeta2 * x2 - 80. / 3 * H10 + 72 * H10 * x
                    - 128. / 3 * H10 * x2 + 16. / 3 * H10 * zeta2
                    - 32. / 3 * H10 * zeta2 * x + 32. / 3 * H10 * zeta2 * x2
                    + 8. / 3 * H100 - 32 * H100 * x + 80. / 3 * H100 * x2
                    + 16. / 3 * H1000 - 32. / 3 * H1000 * x
                    + 32. / 3 * H1000 * x2 + 32. / 3 * H11 - 136. / 3 * H11 * x
                    + 32 * H11 * x2 + 52. / 3 * H11 * zeta2
                    - 104. / 3 * H11 * zeta2 * x + 104. / 3 * H11 * zeta2 * x2
                    - 16. / 3 * H111 - 32. / 3 * H111 * x + 16 * H111 * x2
                    + 16. / 3 * H1111 - 32. / 3 * H1111 * x
                    + 32. / 3 * H1111 * x2 - 32. / 3 * H1101
                    + 64. / 3 * H1101 * x - 64. / 3 * H1101 * x2
                    - 32. / 3 * H1010 + 64. / 3 * H1010 * x
                    - 64. / 3 * H1010 * x2 + 160 * H01 - 236. / 3 * H01 * x
                    - 16 * H01 * x2 + 52. / 3 * H01 * zeta2
                    - 104. / 3 * H01 * zeta2 * x + 104. / 3 * H01 * zeta2 * x2
                    - 8 * H010 - 32. / 3 * H010 * x - 16. / 3 * H0100
                    + 32. / 3 * H0100 * x - 16. / 3 * H011 - 32 * H011 * x
                    + 16 * H011 * x2 + 16. / 3 * H0111 - 32. / 3 * H0111 * x
                    + 32. / 3 * H0111 * x2 - 32. / 3 * H0101
                    + 64. / 3 * H0101 * x - 64. / 3 * H0101 * x2 + 56 * H001
                    - 140 * H001 * x - 16. / 3 * H0010 + 32. / 3 * H0010 * x
                    - 64. / 3 * H0010 * x2 + 20 * H0001 - 48 * H0001 * x
                    + 8 * H00001 - 16 * H00001 * x)
           + CF * Lmmu
                 * (-1988. / 45 + 32. / 45. / x + 616. / 15 * x + 176. / 15 * x2
                    + 112. / 3 * zeta3 + 32. / 3 * zeta3 * x
                    + 320. / 3 * zeta3 * x2 + 16 * zeta2 - 704. / 9 * zeta2 * x
                    + 112 * zeta2 * x2 + 128. / 5 * zeta2 * x3
                    + 128. / 3 * H0m10 + 128. / 3 * H0m10 * x2
                    - 64. / 3 * Hm1 * zeta2 - 128. / 3 * Hm1 * zeta2 * x
                    - 64. / 3 * Hm1 * zeta2 * x2 - 128. / 3 * Hm1m10
                    - 256. / 3 * Hm1m10 * x - 128. / 3 * Hm1m10 * x2 + 64 * Hm10
                    + 32. / 45 * Hm10 / x2 + 256. / 9 * Hm10 * x
                    + 128. / 5 * Hm10 * x3 + 64. / 3 * Hm100
                    + 128. / 3 * Hm100 * x + 64. / 3 * Hm100 * x2
                    - 1604. / 45 * H0 - 32. / 45 * H0 / x + 392. / 15 * H0 * x
                    - 488. / 5 * H0 * x2 + 64. / 3 * H0 * zeta2
                    - 128. / 3 * H0 * zeta2 * x + 64 * H0 * zeta2 * x2
                    - 16. / 3 * H00 + 32. / 9 * H00 * x - 208. / 3 * H00 * x2
                    - 128. / 5 * H00 * x3 - 32. / 3 * H000 + 64. / 3 * H000 * x
                    - 128. / 3 * H000 * x2 - 68. / 3 * H1 + 96 * H1 * x
                    - 72 * H1 * x2 + 32. / 3 * H1 * zeta2
                    - 64. / 3 * H1 * zeta2 * x + 128. / 3 * H1 * zeta2 * x2
                    - 32 * H10 + 224. / 3 * H10 * x - 208. / 3 * H10 * x2
                    - 64. / 3 * H100 * x2 - 88. / 3 * H11 + 352. / 3 * H11 * x
                    - 112 * H11 * x2 - 64. / 3 * H110 + 128. / 3 * H110 * x
                    - 128. / 3 * H110 * x2 - 32 * H111 + 64 * H111 * x
                    - 64 * H111 * x2 - 32 * H101 + 64 * H101 * x
                    - 64 * H101 * x2 - 16 * H01 + 320. / 3 * H01 * x
                    - 112 * H01 * x2 - 64. / 3 * H010 + 128. / 3 * H010 * x
                    - 128. / 3 * H010 * x2 - 80. / 3 * H011
                    + 160. / 3 * H011 * x - 64 * H011 * x2 - 64. / 3 * H001
                    + 128. / 3 * H001 * x - 64 * H001 * x2)
           + CF * LQm
                 * (-189287. / 135 + 11744. / 405. / x + 276332. / 135 * x
                    - 269954. / 405 * x2 + 352. / 3 * zeta3 - 224 * zeta3 * x
                    + 128. / 3 * zeta3 * x2 - 64. / 9 * zeta2 / x
                    + 812. / 3 * zeta2 - 1120. / 9 * zeta2 * x
                    + 1280. / 9 * zeta2 * x2 + 64. / 5 * zeta2 * x3
                    + 8. / 5 * zeta2 * zeta2 - 16. / 5 * zeta2 * zeta2 * x
                    + 160. / 3 * H0m10 - 64 * H0m10 * x + 128. / 3 * H0m10 * x2
                    - 32. / 3 * Hm1 * zeta2 - 64. / 3 * Hm1 * zeta2 * x
                    - 32. / 3 * Hm1 * zeta2 * x2 - 64. / 3 * Hm1m10
                    - 128. / 3 * Hm1m10 * x - 64. / 3 * Hm1m10 * x2 + 80 * Hm10
                    + 16. / 45 * Hm10 / x2 - 64. / 9 * Hm10 / x
                    + 128. / 9 * Hm10 * x - 496. / 9 * Hm10 * x2
                    + 64. / 5 * Hm10 * x3 + 32. / 3 * Hm100
                    + 64. / 3 * Hm100 * x + 32. / 3 * Hm100 * x2
                    - 12284. / 15 * H0 - 16. / 45 * H0 / x
                    + 22528. / 45 * H0 * x - 5776. / 45 * H0 * x2
                    + 32 * H0 * zeta3 - 64 * H0 * zeta3 * x + 88 * H0 * zeta2
                    - 224 * H0 * zeta2 * x - 32. / 3 * H0 * zeta2 * x2
                    - 4760. / 9 * H00 + 4256. / 9 * H00 * x
                    - 896. / 9 * H00 * x2 - 64. / 5 * H00 * x3
                    + 48 * H00 * zeta2 - 96 * H00 * zeta2 * x - 160 * H000
                    + 480 * H000 * x + 64 * H000 * x2 - 96 * H0000
                    + 192 * H0000 * x - 1166. / 3 * H1 - 512. / 27 * H1 / x
                    + 5632. / 9 * H1 * x - 6160. / 27 * H1 * x2
                    + 16 * H1 * zeta2 - 32 * H1 * zeta2 * x
                    + 128. / 3 * H1 * zeta2 * x2 - 1120. / 9 * H10
                    + 64. / 9 * H10 / x + 2264. / 9 * H10 * x - 176 * H10 * x2
                    - 64. / 3 * H100 + 128. / 3 * H100 * x
                    - 160. / 3 * H100 * x2 - 1048. / 9 * H11 + 64. / 9 * H11 / x
                    + 2024. / 9 * H11 * x - 448. / 3 * H11 * x2 - 16 * H110
                    + 32 * H110 * x - 32 * H110 * x2 - 16 * H111 + 32 * H111 * x
                    - 32 * H111 * x2 - 80. / 3 * H101 + 160. / 3 * H101 * x
                    - 160. / 3 * H101 * x2 - 812. / 3 * H01 + 416. / 3 * H01 * x
                    - 1280. / 9 * H01 * x2 - 104. / 3 * H010
                    + 64. / 3 * H010 * x - 32. / 3 * H010 * x2 - 104. / 3 * H011
                    + 64. / 3 * H011 * x - 32. / 3 * H011 * x2 - 88 * H001
                    + 160 * H001 * x + 32. / 3 * H001 * x2 - 16 * H0010
                    + 32 * H0010 * x - 16 * H0011 + 32 * H0011 * x - 48 * H0001
                    + 96 * H0001 * x)
           + CF * LQm * Lmmu
                 * (28 - 160. / 3 * x + 56. / 3 * x2 - 16 * zeta2
                    + 32 * zeta2 * x - 128. / 3 * zeta2 * x2 + 16. / 3 * H0
                    - 32 * H0 * x + 160. / 3 * H0 * x2 + 32. / 3 * H00
                    - 64. / 3 * H00 * x + 128. / 3 * H00 * x2 + 56. / 3 * H1
                    - 64 * H1 * x + 160. / 3 * H1 * x2 + 64. / 3 * H10
                    - 128. / 3 * H10 * x + 128. / 3 * H10 * x2 + 64. / 3 * H11
                    - 128. / 3 * H11 * x + 128. / 3 * H11 * x2 + 16 * H01
                    - 32 * H01 * x + 128. / 3 * H01 * x2)
           + CF * LQm2
                 * (553. / 3 + 256. / 27. / x - 848. / 3 * x + 2156. / 27 * x2
                    - 8 * zeta3 + 16 * zeta3 * x - 52. / 3 * zeta2
                    + 32. / 3 * zeta2 * x - 16. / 3 * zeta2 * x2
                    + 1142. / 9 * H0 - 472. / 9 * H0 * x + 488. / 9 * H0 * x2
                    - 8 * H0 * zeta2 + 16 * H0 * zeta2 * x + 124. / 3 * H00
                    - 224. / 3 * H00 * x - 32. / 3 * H00 * x2 + 24 * H000
                    - 48 * H000 * x + 448. / 9 * H1 - 32. / 9 * H1 / x
                    - 860. / 9 * H1 * x + 520. / 9 * H1 * x2 + 32. / 3 * H10
                    - 64. / 3 * H10 * x + 64. / 3 * H10 * x2 + 8 * H11
                    - 16 * H11 * x + 16 * H11 * x2 + 52. / 3 * H01
                    - 32. / 3 * H01 * x + 16. / 3 * H01 * x2 + 8 * H001
                    - 16 * H001 * x)
           + CF * LQm2 * Lmmu
                 * (-4. / 3 + 16. / 3 * x - 8. / 3 * H0 + 16. / 3 * H0 * x
                    - 32. / 3 * H0 * x2 - 16. / 3 * H1 + 32. / 3 * H1 * x
                    - 32. / 3 * H1 * x2)
           + CF * LQm3
                 * (-112. / 9 + 32. / 27. / x + 196. / 9 * x - 248. / 27 * x2
                    - 44. / 9 * H0 + 16. / 9 * H0 * x - 8. / 3 * H00
                    + 16. / 3 * H00 * x - 16. / 9 * H1 + 32. / 9 * H1 * x
                    - 32. / 9 * H1 * x2)
           + CF * nf
                 * (-521477. / 972 - 30608. / 729. / x + 168241. / 486 * x
                    + 109745. / 729 * x2 - 64 * zeta5 + 128 * zeta5 * x
                    + 256. / 27 * zeta3 / x + 4952. / 27 * zeta3
                    + 6752. / 27 * zeta3 * x + 256. / 3 * zeta3 * x2
                    + 2104. / 81 * zeta2 - 3752. / 81 * zeta2 * x
                    + 4076. / 81 * zeta2 * x2 - 8. / 9 * zeta2 * zeta2
                    - 1648. / 45 * zeta2 * zeta2 * x
                    + 496. / 45 * zeta2 * zeta2 * x2 - 93808. / 243 * H0
                    - 21970. / 243 * H0 * x - 10432. / 27 * H0 * x2
                    + 608. / 9 * H0 * zeta3 - 640. / 9 * H0 * zeta3 * x
                    - 128. / 3 * H0 * zeta3 * x2 - 13990. / 81 * H00
                    - 12154. / 81 * H00 * x + 23164. / 81 * H00 * x2
                    + 128. / 3 * H00 * zeta3 - 256. / 3 * H00 * zeta3 * x
                    - 3394. / 27 * H000 - 676. / 27 * H000 * x
                    - 1328. / 9 * H000 * x2 - 376. / 9 * H0000
                    + 32. / 9 * H0000 * x + 352. / 9 * H0000 * x2 - 32 * H00000
                    + 64 * H00000 * x + 1420. / 243 * H1 + 1108. / 243 * H1 * x
                    - 2800. / 243 * H1 * x2 + 64. / 9 * H1 * zeta3
                    - 128. / 9 * H1 * zeta3 * x + 128. / 9 * H1 * zeta3 * x2
                    + 12826. / 81 * H10 + 416. / 27 * H10 / x
                    - 7376. / 81 * H10 * x - 7012. / 81 * H10 * x2
                    + 3616. / 27 * H100 - 128. / 9 * H100 / x
                    - 5264. / 27 * H100 * x + 1904. / 27 * H100 * x2
                    - 80. / 9 * H1000 + 160. / 9 * H1000 * x
                    - 160. / 9 * H1000 * x2 - 916. / 81 * H11
                    + 2780. / 81 * H11 * x - 2132. / 81 * H11 * x2
                    - 16. / 3 * H11 * zeta2 + 32. / 3 * H11 * zeta2 * x
                    - 32. / 3 * H11 * zeta2 * x2 + 304. / 27 * H111
                    - 80. / 27 * H111 * x - 64. / 27 * H111 * x2
                    - 40. / 9 * H1111 + 80. / 9 * H1111 * x
                    - 80. / 9 * H1111 * x2 + 16. / 3 * H1101
                    - 32. / 3 * H1101 * x + 32. / 3 * H1101 * x2
                    - 1132. / 81 * H01 + 1808. / 81 * H01 * x
                    - 2132. / 81 * H01 * x2 - 16. / 3 * H01 * zeta2
                    + 32. / 3 * H01 * zeta2 * x - 32. / 3 * H01 * zeta2 * x2
                    + 1300. / 9 * H010 + 520. / 9 * H010 * x
                    + 704. / 9 * H010 * x2 + 160. / 3 * H0100
                    - 32. / 3 * H0100 * x - 128. / 3 * H0100 * x2
                    + 304. / 27 * H011 + 496. / 27 * H011 * x
                    - 64. / 27 * H011 * x2 - 40. / 9 * H0111
                    + 80. / 9 * H0111 * x - 80. / 9 * H0111 * x2
                    + 16. / 3 * H0101 - 32. / 3 * H0101 * x
                    + 32. / 3 * H0101 * x2 + 48 * H0010 - 32 * H0010 * x
                    - 128. / 3 * H0010 * x2 + 32 * H00100 - 64 * H00100 * x
                    + 32 * H00010 - 64 * H00010 * x)
           + CF * nf * Lmmu
                 * (-96932. / 135 + 1088. / 135. / x + 146732. / 135 * x
                    - 49988. / 135 * x2 + 200. / 3 * zeta3
                    - 656. / 3 * zeta3 * x + 256. / 3 * zeta3 * x2
                    - 64. / 9 * zeta2 / x + 248 * zeta2 - 1000. / 9 * zeta2 * x
                    + 104 * zeta2 * x2 + 64. / 5 * zeta2 * x3
                    - 88. / 5 * zeta2 * zeta2 + 176. / 5 * zeta2 * zeta2 * x
                    + 160. / 3 * H0m10 - 64 * H0m10 * x + 128. / 3 * H0m10 * x2
                    - 32. / 3 * Hm1 * zeta2 - 64. / 3 * Hm1 * zeta2 * x
                    - 32. / 3 * Hm1 * zeta2 * x2 - 64. / 3 * Hm1m10
                    - 128. / 3 * Hm1m10 * x - 64. / 3 * Hm1m10 * x2 + 80 * Hm10
                    + 16. / 45 * Hm10 / x2 - 64. / 9 * Hm10 / x
                    + 128. / 9 * Hm10 * x - 496. / 9 * Hm10 * x2
                    + 64. / 5 * Hm10 * x3 + 32. / 3 * Hm100
                    + 64. / 3 * Hm100 * x + 32. / 3 * Hm100 * x2
                    - 6824. / 15 * H0 - 16. / 45 * H0 / x + 3192. / 5 * H0 * x
                    - 41408. / 135 * H0 * x2 + 248. / 3 * H0 * zeta2
                    - 640. / 3 * H0 * zeta2 * x - 32 * H0 * zeta2 * x2
                    - 872. / 3 * H00 + 3544. / 9 * H00 * x + 88 * H00 * x2
                    - 64. / 5 * H00 * x3 + 48 * H00 * zeta2
                    - 96 * H00 * zeta2 * x - 232. / 3 * H000
                    + 1184. / 3 * H000 * x + 128. / 3 * H000 * x2 - 48 * H0000
                    + 96 * H0000 * x - 3172. / 9 * H1 - 512. / 27 * H1 / x
                    + 5116. / 9 * H1 * x - 5248. / 27 * H1 * x2
                    + 16. / 3 * H1 * zeta2 - 32. / 3 * H1 * zeta2 * x
                    + 64. / 3 * H1 * zeta2 * x2 - 160 * H10 + 128. / 9 * H10 / x
                    + 832. / 3 * H10 * x - 1304. / 9 * H10 * x2
                    - 32. / 3 * H100 * x2 - 260. / 3 * H11 + 64. / 9 * H11 / x
                    + 536. / 3 * H11 * x - 1000. / 9 * H11 * x2 - 32. / 3 * H110
                    + 64. / 3 * H110 * x - 64. / 3 * H110 * x2 - 16 * H111
                    + 32 * H111 * x - 32 * H111 * x2 - 16 * H101 + 32 * H101 * x
                    - 32 * H101 * x2 - 248 * H01 + 376. / 3 * H01 * x
                    - 104 * H01 * x2 - 176. / 3 * H010 + 64. / 3 * H010 * x
                    + 64. / 3 * H010 * x2 - 112. / 3 * H011 + 80. / 3 * H011 * x
                    - 32. / 3 * H011 * x2 - 248. / 3 * H001
                    + 448. / 3 * H001 * x + 32 * H001 * x2 - 32 * H0010
                    + 64 * H0010 * x - 16 * H0011 + 32 * H0011 * x - 48 * H0001
                    + 96 * H0001 * x)
           + CF * nf * Lmmu2
                 * (788. / 9 + 16. / 9. / x - 1928. / 9 * x + 1124. / 9 * x2
                    - 8 * zeta3 + 16 * zeta3 * x - 12 * zeta2
                    + 32. / 3 * zeta2 * x2 + 48 * H0 - 68 * H0 * x
                    - 136. / 9 * H0 * x2 - 8 * H0 * zeta2 + 16 * H0 * zeta2 * x
                    + 12 * H00 - 48 * H00 * x - 32. / 3 * H00 * x2 + 8 * H000
                    - 16 * H000 * x + 36 * H1 - 32. / 9 * H1 / x - 60 * H1 * x
                    + 248. / 9 * H1 * x2 + 12 * H01 - 32. / 3 * H01 * x2
                    + 8 * H001 - 16 * H001 * x)
           + CF * nf * LQm
                 * (-153242. / 135 + 13904. / 405. / x + 212342. / 135 * x
                    - 182204. / 405 * x2 + 352. / 3 * zeta3 - 224 * zeta3 * x
                    + 128. / 3 * zeta3 * x2 - 64. / 9 * zeta2 / x
                    + 812. / 3 * zeta2 - 1120. / 9 * zeta2 * x
                    + 1280. / 9 * zeta2 * x2 + 64. / 5 * zeta2 * x3
                    + 8. / 5 * zeta2 * zeta2 - 16. / 5 * zeta2 * zeta2 * x
                    + 160. / 3 * H0m10 - 64 * H0m10 * x + 128. / 3 * H0m10 * x2
                    - 32. / 3 * Hm1 * zeta2 - 64. / 3 * Hm1 * zeta2 * x
                    - 32. / 3 * Hm1 * zeta2 * x2 - 64. / 3 * Hm1m10
                    - 128. / 3 * Hm1m10 * x - 64. / 3 * Hm1m10 * x2 + 80 * Hm10
                    + 16. / 45 * Hm10 / x2 - 64. / 9 * Hm10 / x
                    + 128. / 9 * Hm10 * x - 496. / 9 * Hm10 * x2
                    + 64. / 5 * Hm10 * x3 + 32. / 3 * Hm100
                    + 64. / 3 * Hm100 * x + 32. / 3 * Hm100 * x2
                    - 10124. / 15 * H0 - 16. / 45 * H0 / x
                    + 21628. / 45 * H0 * x - 7936. / 45 * H0 * x2
                    + 32 * H0 * zeta3 - 64 * H0 * zeta3 * x + 88 * H0 * zeta2
                    - 224 * H0 * zeta2 * x - 32. / 3 * H0 * zeta2 * x2
                    - 4256. / 9 * H00 + 2996. / 9 * H00 * x
                    - 896. / 9 * H00 * x2 - 64. / 5 * H00 * x3
                    + 48 * H00 * zeta2 - 96 * H00 * zeta2 * x - 140 * H000
                    + 432 * H000 * x + 64 * H000 * x2 - 88 * H0000
                    + 176 * H0000 * x - 1166. / 3 * H1 - 512. / 27 * H1 / x
                    + 5632. / 9 * H1 * x - 6160. / 27 * H1 * x2
                    + 16 * H1 * zeta2 - 32 * H1 * zeta2 * x
                    + 128. / 3 * H1 * zeta2 * x2 - 1120. / 9 * H10
                    + 64. / 9 * H10 / x + 2264. / 9 * H10 * x - 176 * H10 * x2
                    - 64. / 3 * H100 + 128. / 3 * H100 * x
                    - 160. / 3 * H100 * x2 - 1048. / 9 * H11 + 64. / 9 * H11 / x
                    + 2024. / 9 * H11 * x - 448. / 3 * H11 * x2 - 16 * H110
                    + 32 * H110 * x - 32 * H110 * x2 - 16 * H111 + 32 * H111 * x
                    - 32 * H111 * x2 - 80. / 3 * H101 + 160. / 3 * H101 * x
                    - 160. / 3 * H101 * x2 - 812. / 3 * H01 + 416. / 3 * H01 * x
                    - 1280. / 9 * H01 * x2 - 104. / 3 * H010
                    + 64. / 3 * H010 * x - 32. / 3 * H010 * x2 - 104. / 3 * H011
                    + 64. / 3 * H011 * x - 32. / 3 * H011 * x2 - 88 * H001
                    + 160 * H001 * x + 32. / 3 * H001 * x2 - 16 * H0010
                    + 32 * H0010 * x - 16 * H0011 + 32 * H0011 * x - 48 * H0001
                    + 96 * H0001 * x)
           + CF * nf * LQm * Lmmu
                 * (3196. / 9 + 512. / 27. / x - 4924. / 9 * x + 4528. / 27 * x2
                    - 16 * zeta3 + 32 * zeta3 * x - 32 * zeta2 + 16 * zeta2 * x
                    + 728. / 3 * H0 - 88 * H0 * x + 224. / 3 * H0 * x2
                    - 16 * H0 * zeta2 + 32 * H0 * zeta2 * x + 232. / 3 * H00
                    - 416. / 3 * H00 * x - 128. / 3 * H00 * x2 + 48 * H000
                    - 96 * H000 * x + 244. / 3 * H1 - 64. / 9 * H1 / x
                    - 152 * H1 * x + 736. / 9 * H1 * x2 + 32. / 3 * H10
                    - 64. / 3 * H10 * x + 64. / 3 * H10 * x2 + 32. / 3 * H11
                    - 64. / 3 * H11 * x + 64. / 3 * H11 * x2 + 32 * H01
                    - 16 * H01 * x + 16 * H001 - 32 * H001 * x)
           + CF * nf * LQm * Lmmu2
                 * (-36 + 32. / 9. / x + 60 * x - 248. / 9 * x2 - 12 * H0
                    + 32. / 3 * H0 * x2 - 8 * H00 + 16 * H00 * x)
           + CF * nf * LQm2
                 * (553. / 3 + 256. / 27. / x - 848. / 3 * x + 2156. / 27 * x2
                    - 8 * zeta3 + 16 * zeta3 * x - 52. / 3 * zeta2
                    + 32. / 3 * zeta2 * x - 16. / 3 * zeta2 * x2
                    + 1142. / 9 * H0 - 472. / 9 * H0 * x + 488. / 9 * H0 * x2
                    - 8 * H0 * zeta2 + 16 * H0 * zeta2 * x + 124. / 3 * H00
                    - 224. / 3 * H00 * x - 32. / 3 * H00 * x2 + 24 * H000
                    - 48 * H000 * x + 448. / 9 * H1 - 32. / 9 * H1 / x
                    - 860. / 9 * H1 * x + 520. / 9 * H1 * x2 + 32. / 3 * H10
                    - 64. / 3 * H10 * x + 64. / 3 * H10 * x2 + 8 * H11
                    - 16 * H11 * x + 16 * H11 * x2 + 52. / 3 * H01
                    - 32. / 3 * H01 * x + 16. / 3 * H01 * x2 + 8 * H001
                    - 16 * H001 * x)
           + CF * nf * LQm2 * Lmmu
                 * (-110. / 3 + 32. / 9. / x + 188. / 3 * x - 248. / 9 * x2
                    - 40. / 3 * H0 + 8. / 3 * H0 * x + 16. / 3 * H0 * x2
                    - 8 * H00 + 16 * H00 * x - 8. / 3 * H1 + 16. / 3 * H1 * x
                    - 16. / 3 * H1 * x2)
           + CF * nf * LQm3
                 * (-112. / 9 + 32. / 27. / x + 196. / 9 * x - 248. / 27 * x2
                    - 44. / 9 * H0 + 16. / 9 * H0 * x - 8. / 3 * H00
                    + 16. / 3 * H00 * x - 16. / 9 * H1 + 32. / 9 * H1 * x
                    - 32. / 9 * H1 * x2)
           + CF * CF
                 * (-77. / 2 - 42 * x - 39 * x2 - 4 * zeta5 + 8 * zeta5 * x
                    - 224 * zeta5 * x2 + 26. / 3 * zeta3 - 68. / 3 * zeta3 * x
                    + 168 * zeta3 * x2 + 375. / 4 * zeta2 + 255. / 2 * zeta2 * x
                    - 58 * zeta2 * x2 - 192 * zeta2 * ln2
                    + 384 * zeta2 * ln2 * x - 384 * zeta2 * ln2 * x2
                    + 8. / 3 * zeta2 * zeta3 - 112. / 3 * zeta2 * zeta3 * x
                    - 160. / 3 * zeta2 * zeta3 * x2 - 4. / 5 * zeta2 * zeta2
                    + 316. / 5 * zeta2 * zeta2 * x
                    - 704. / 5 * zeta2 * zeta2 * x2 - 8 * H0m10 * zeta2
                    - 16 * H0m10 * zeta2 * x - 32 * H0m10 * zeta2 * x2
                    + 8 * Hm1 * zeta2 * zeta2 + 16 * Hm1 * zeta2 * zeta2 * x
                    + 16 * Hm1 * zeta2 * zeta2 * x2 + 16 * Hm1m10 * zeta2
                    + 32 * Hm1m10 * zeta2 * x + 32 * Hm1m10 * zeta2 * x2
                    - 16 * Hm10 * zeta2 - 40 * Hm10 * zeta2 * x
                    - 24 * Hm10 * zeta2 * x2 - 8 * Hm100 * zeta2
                    - 16 * Hm100 * zeta2 * x - 16 * Hm100 * zeta2 * x2 - 26 * H0
                    - 182 * H0 * x + 272 * H0 * x2 - 12 * H0 * zeta3
                    - 80 * H0 * zeta3 * x + 128 * H0 * zeta3 * x2
                    + 29 * H0 * zeta2 + 38 * H0 * zeta2 * x
                    + 122 * H0 * zeta2 * x2 - 228. / 5 * H0 * zeta2 * zeta2
                    + 376. / 5 * H0 * zeta2 * zeta2 * x
                    - 1136. / 5 * H0 * zeta2 * zeta2 * x2 + 63 * H00
                    - 196 * H00 * x + 384 * H00 * x2 - 80. / 3 * H00 * zeta3
                    + 160. / 3 * H00 * zeta3 * x - 352. / 3 * H00 * zeta3 * x2
                    - 8 * H00 * zeta2 + 74 * H00 * zeta2 * x
                    - 104 * H00 * zeta2 * x2 - 10 * H000 - 4 * H000 * x
                    - 72 * H000 * x2 - 6 * H000 * zeta2 + 28 * H000 * zeta2 * x
                    - 16 * H000 * zeta2 * x2 - 12 * H0000 - 64 * H0000 * x
                    + 208 * H0000 * x2 + 8 * H00000 - 16 * H00000 * x
                    + 64 * H00000 * x2 + 39 * H1 - 313 * H1 * x + 272 * H1 * x2
                    - 64. / 3 * H1 * zeta3 - 224. / 3 * H1 * zeta3 * x
                    + 128 * H1 * zeta3 * x2 - 20 * H1 * zeta2
                    - 164 * H1 * zeta2 * x + 122 * H1 * zeta2 * x2
                    - 528. / 5 * H1 * zeta2 * zeta2
                    + 1056. / 5 * H1 * zeta2 * zeta2 * x
                    - 1056. / 5 * H1 * zeta2 * zeta2 * x2 + 50 * H10
                    - 312 * H10 * x + 384 * H10 * x2 - 176. / 3 * H10 * zeta3
                    + 352. / 3 * H10 * zeta3 * x - 352. / 3 * H10 * zeta3 * x2
                    - 54 * H10 * zeta2 + 168 * H10 * zeta2 * x
                    - 128 * H10 * zeta2 * x2 + 28 * H100 + 72 * H100 * x
                    - 72 * H100 * x2 - 16 * H100 * zeta2 + 32 * H100 * zeta2 * x
                    - 32 * H100 * zeta2 * x2 + 48 * H1000 - 272 * H1000 * x
                    + 208 * H1000 * x2 + 32 * H10000 - 64 * H10000 * x
                    + 64 * H10000 * x2 + 148 * H11 - 258 * H11 * x
                    + 272 * H11 * x2 - 16. / 3 * H11 * zeta3
                    + 32. / 3 * H11 * zeta3 * x - 32. / 3 * H11 * zeta3 * x2
                    - 144 * H11 * zeta2 * x + 128 * H11 * zeta2 * x2 + 14 * H110
                    + 172 * H110 * x - 176 * H110 * x2 + 72 * H110 * zeta2
                    - 144 * H110 * zeta2 * x + 144 * H110 * zeta2 * x2
                    + 60 * H1100 - 224 * H1100 * x + 160 * H1100 * x2
                    + 32 * H11000 - 64 * H11000 * x + 64 * H11000 * x2
                    + 124 * H111 + 112 * H111 * x - 248 * H111 * x2
                    + 160 * H111 * zeta2 - 320 * H111 * zeta2 * x
                    + 320 * H111 * zeta2 * x2 + 12 * H1110 + 32 * H1110 * x
                    - 64 * H1110 * x2 - 76 * H1111 + 304 * H1111 * x
                    - 288 * H1111 * x2 - 48 * H11110 + 96 * H11110 * x
                    - 96 * H11110 * x2 - 80 * H11111 + 160 * H11111 * x
                    - 160 * H11111 * x2 - 112 * H11101 + 224 * H11101 * x
                    - 224 * H11101 * x2 + 4 * H1101 + 64 * H1101 * x
                    - 64 * H1101 * x2 - 80 * H11010 + 160 * H11010 * x
                    - 160 * H11010 * x2 - 112 * H11011 + 224 * H11011 * x
                    - 224 * H11011 * x2 - 32 * H11001 + 64 * H11001 * x
                    - 64 * H11001 * x2 + 14 * H101 + 172 * H101 * x
                    - 176 * H101 * x2 + 120 * H101 * zeta2
                    - 240 * H101 * zeta2 * x + 240 * H101 * zeta2 * x2
                    + 24 * H1010 - 128 * H1010 * x + 112 * H1010 * x2
                    - 16 * H10100 + 32 * H10100 * x - 32 * H10100 * x2
                    + 12 * H1011 + 32 * H1011 * x - 64 * H1011 * x2
                    - 48 * H10110 + 96 * H10110 * x - 96 * H10110 * x2
                    - 48 * H10111 + 96 * H10111 * x - 96 * H10111 * x2
                    - 80 * H10101 + 160 * H10101 * x - 160 * H10101 * x2
                    + 60 * H1001 - 224 * H1001 * x + 160 * H1001 * x2
                    - 16 * H10010 + 32 * H10010 * x - 32 * H10010 * x2
                    + 32 * H10001 - 64 * H10001 * x + 64 * H10001 * x2
                    + 28 * H01 - 305 * H01 * x + 272 * H01 * x2
                    + 16. / 3 * H01 * zeta3 - 32. / 3 * H01 * zeta3 * x
                    - 32. / 3 * H01 * zeta3 * x2 + 40 * H01 * zeta2
                    - 164 * H01 * zeta2 * x + 128 * H01 * zeta2 * x2 + 4 * H010
                    + 16 * H010 * x - 176 * H010 * x2 + 100 * H010 * zeta2
                    - 200 * H010 * zeta2 * x + 144 * H010 * zeta2 * x2
                    - 20 * H0100 - 32 * H0100 * x + 160 * H0100 * x2
                    - 16 * H01000 + 32 * H01000 * x + 64 * H01000 * x2
                    + 88 * H011 - 32 * H011 * x - 248 * H011 * x2
                    + 144 * H011 * zeta2 - 288 * H011 * zeta2 * x
                    + 320 * H011 * zeta2 * x2 - 36 * H0110 + 120 * H0110 * x
                    - 64 * H0110 * x2 - 40 * H01100 + 80 * H01100 * x
                    - 44 * H0111 + 208 * H0111 * x - 288 * H0111 * x2
                    - 56 * H01110 + 112 * H01110 * x - 96 * H01110 * x2
                    - 56 * H01111 + 112 * H01111 * x - 160 * H01111 * x2
                    - 104 * H01101 + 208 * H01101 * x - 224 * H01101 * x2
                    - 28 * H0101 + 120 * H0101 * x - 64 * H0101 * x2
                    - 96 * H01010 + 192 * H01010 * x - 160 * H01010 * x2
                    - 120 * H01011 + 240 * H01011 * x - 224 * H01011 * x2
                    - 72 * H01001 + 144 * H01001 * x - 64 * H01001 * x2
                    - 37 * H001 - 10 * H001 * x - 176 * H001 * x2
                    + 68 * H001 * zeta2 - 136 * H001 * zeta2 * x
                    + 240 * H001 * zeta2 * x2 + 12 * H0010 - 48 * H0010 * x
                    + 112 * H0010 * x2 - 24 * H00100 + 48 * H00100 * x
                    - 32 * H00100 * x2 + 24 * H0011 + 32 * H0011 * x
                    - 64 * H0011 * x2 - 32 * H00110 + 64 * H00110 * x
                    - 96 * H00110 * x2 - 24 * H00111 + 48 * H00111 * x
                    - 96 * H00111 * x2 - 48 * H00101 + 96 * H00101 * x
                    - 160 * H00101 * x2 + 12 * H0001 - 52 * H0001 * x
                    + 160 * H0001 * x2 - 8 * H00010 + 16 * H00010 * x
                    - 32 * H00010 * x2 + 12 * H00001 - 24 * H00001 * x
                    + 64 * H00001 * x2)
           + CF * CF * LQm
                 * (-383. / 15 + 32. / 15. / x + 631. / 5 * x - 364. / 5 * x2
                    - 184 * zeta3 + 432 * zeta3 * x - 1024 * zeta3 * x2
                    - 128 * zeta2 - 1064. / 3 * zeta2 * x - 736 * zeta2 * x2
                    + 384. / 5 * zeta2 * x3 + 228. / 5 * zeta2 * zeta2
                    - 1384. / 5 * zeta2 * zeta2 * x
                    + 768. / 5 * zeta2 * zeta2 * x2 - 192 * H00m10
                    - 512 * H00m10 * x2 + 112 * H0m1 * zeta2
                    - 32 * H0m1 * zeta2 * x + 256 * H0m1 * zeta2 * x2
                    + 160 * H0m1m10 - 192 * H0m1m10 * x + 256 * H0m1m10 * x2
                    - 448 * H0m10 - 64 * H0m10 * x - 32 * H0m10 * x2
                    - 160 * H0m100 - 64 * H0m100 * x - 448 * H0m100 * x2
                    - 32 * H0m101 - 64 * H0m101 * x - 128 * H0m101 * x2
                    + 112 * Hm1 * zeta3 + 224 * Hm1 * zeta3 * x
                    + 224 * Hm1 * zeta3 * x2 + 288 * Hm1 * zeta2
                    + 400 * Hm1 * zeta2 * x + 112 * Hm1 * zeta2 * x2
                    + 128 * Hm10m10 + 256 * Hm10m10 * x + 256 * Hm10m10 * x2
                    - 128 * Hm1m1 * zeta2 - 256 * Hm1m1 * zeta2 * x
                    - 256 * Hm1m1 * zeta2 * x2 - 128 * Hm1m1m10
                    - 256 * Hm1m1m10 * x - 256 * Hm1m1m10 * x2 + 448 * Hm1m10
                    + 480 * Hm1m10 * x + 32 * Hm1m10 * x2 + 224 * Hm1m100
                    + 448 * Hm1m100 * x + 448 * Hm1m100 * x2 + 64 * Hm1m101
                    + 128 * Hm1m101 * x + 128 * Hm1m101 * x2 - 976 * Hm10
                    + 32. / 15 * Hm10 / x2 - 3488. / 3 * Hm10 * x
                    - 176 * Hm10 * x2 + 384. / 5 * Hm10 * x3 + 32 * Hm10 * zeta2
                    + 64 * Hm10 * zeta2 * x + 64 * Hm10 * zeta2 * x2
                    - 384 * Hm100 - 640 * Hm100 * x - 256 * Hm100 * x2
                    - 64 * Hm1000 - 128 * Hm1000 * x - 128 * Hm1000 * x2
                    - 64 * Hm101 - 160 * Hm101 * x - 96 * Hm101 * x2
                    - 32 * Hm1001 - 64 * Hm1001 * x - 64 * Hm1001 * x2
                    + 616. / 15 * H0 - 32. / 15 * H0 / x + 1062. / 5 * H0 * x
                    + 976. / 5 * H0 * x2 - 144 * H0 * zeta3
                    + 480 * H0 * zeta3 * x - 704 * H0 * zeta3 * x2
                    - 72 * H0 * zeta2 + 248 * H0 * zeta2 * x
                    - 800 * H0 * zeta2 * x2 + 78 * H00 + 2204. / 3 * H00 * x
                    + 616 * H00 * x2 - 384. / 5 * H00 * x3 - 120 * H00 * zeta2
                    + 176 * H00 * zeta2 * x - 576 * H00 * zeta2 * x2
                    + 232 * H000 * x + 656 * H000 * x2 + 40 * H0000
                    + 48 * H0000 * x + 320 * H0000 * x2 - 85 * H1 - 216 * H1 * x
                    + 272 * H1 * x2 - 240 * H1 * zeta3 + 480 * H1 * zeta3 * x
                    - 480 * H1 * zeta3 * x2 - 40 * H1 * zeta2
                    + 624 * H1 * zeta2 * x - 688 * H1 * zeta2 * x2
                    - 128 * H10m10 + 256 * H10m10 * x - 256 * H10m10 * x2
                    + 232 * H10 - 760 * H10 * x + 440 * H10 * x2
                    - 256 * H10 * zeta2 + 512 * H10 * zeta2 * x
                    - 512 * H10 * zeta2 * x2 - 24 * H100 - 320 * H100 * x
                    + 400 * H100 * x2 + 96 * H1000 - 192 * H1000 * x
                    + 192 * H1000 * x2 + 310 * H11 - 1092 * H11 * x
                    + 736 * H11 * x2 - 256 * H11 * zeta2 + 512 * H11 * zeta2 * x
                    - 512 * H11 * zeta2 * x2 + 256 * H110 - 832 * H110 * x
                    + 704 * H110 * x2 + 160 * H1100 - 320 * H1100 * x
                    + 320 * H1100 * x2 + 288 * H111 - 1152 * H111 * x
                    + 1008 * H111 * x2 + 288 * H1110 - 576 * H1110 * x
                    + 576 * H1110 * x2 + 384 * H1111 - 768 * H1111 * x
                    + 768 * H1111 * x2 + 320 * H1101 - 640 * H1101 * x
                    + 640 * H1101 * x2 + 264 * H101 - 864 * H101 * x
                    + 704 * H101 * x2 + 288 * H1010 - 576 * H1010 * x
                    + 576 * H1010 * x2 + 352 * H1011 - 704 * H1011 * x
                    + 704 * H1011 * x2 + 256 * H1001 - 512 * H1001 * x
                    + 512 * H1001 * x2 + 128 * H01 - 808 * H01 * x
                    + 736 * H01 * x2 - 160 * H01 * zeta2 + 576 * H01 * zeta2 * x
                    - 512 * H01 * zeta2 * x2 + 128 * H010 - 624 * H010 * x
                    + 704 * H010 * x2 + 80 * H0100 - 416 * H0100 * x
                    + 320 * H0100 * x2 + 152 * H011 - 920 * H011 * x
                    + 1008 * H011 * x2 + 224 * H0110 - 448 * H0110 * x
                    + 576 * H0110 * x2 + 288 * H0111 - 576 * H0111 * x
                    + 768 * H0111 * x2 + 240 * H0101 - 480 * H0101 * x
                    + 640 * H0101 * x2 + 72 * H001 - 312 * H001 * x
                    + 800 * H001 * x2 + 160 * H0010 - 320 * H0010 * x
                    + 576 * H0010 * x2 + 200 * H0011 - 400 * H0011 * x
                    + 704 * H0011 * x2 + 120 * H0001 - 176 * H0001 * x
                    + 576 * H0001 * x2)
           + CF * CF * LQm2
                 * (-33. / 2 - 46 * x + 52 * x2 + 28 * zeta3 + 8 * zeta3 * x
                    + 128 * zeta3 * x2 + 24 * zeta2 - 44 * zeta2 * x
                    + 160 * zeta2 * x2 + 16 * H0m10 + 32 * H0m10 * x
                    + 64 * H0m10 * x2 - 16 * Hm1 * zeta2 - 32 * Hm1 * zeta2 * x
                    - 32 * Hm1 * zeta2 * x2 - 32 * Hm1m10 - 64 * Hm1m10 * x
                    - 64 * Hm1m10 * x2 + 32 * Hm10 + 80 * Hm10 * x
                    + 48 * Hm10 * x2 + 16 * Hm100 + 32 * Hm100 * x
                    + 32 * Hm100 * x2 - 22 * H0 + 12 * H0 * x - 36 * H0 * x2
                    + 36 * H0 * zeta2 - 40 * H0 * zeta2 * x
                    + 160 * H0 * zeta2 * x2 - 36 * H00 * x - 128 * H00 * x2
                    - 12 * H000 - 8 * H000 * x - 96 * H000 * x2 - 65 * H1
                    + 126 * H1 * x - 36 * H1 * x2 + 64 * H1 * zeta2
                    - 128 * H1 * zeta2 * x + 128 * H1 * zeta2 * x2 - 40 * H10
                    + 144 * H10 * x - 80 * H10 * x2 - 32 * H100 + 64 * H100 * x
                    - 64 * H100 * x2 - 64 * H11 + 224 * H11 * x - 160 * H11 * x2
                    - 80 * H110 + 160 * H110 * x - 160 * H110 * x2 - 96 * H111
                    + 192 * H111 * x - 192 * H111 * x2 - 80 * H101
                    + 160 * H101 * x - 160 * H101 * x2 - 24 * H01
                    + 124 * H01 * x - 160 * H01 * x2 - 48 * H010 + 96 * H010 * x
                    - 160 * H010 * x2 - 64 * H011 + 128 * H011 * x
                    - 192 * H011 * x2 - 36 * H001 + 72 * H001 * x
                    - 160 * H001 * x2)
           + CF * CF * LQm3
                 * (11. / 3 - 2. / 3 * x - 16. / 3 * zeta2 + 32. / 3 * zeta2 * x
                    - 64. / 3 * zeta2 * x2 - 4 * H0 * x + 4. / 3 * H00
                    - 8. / 3 * H00 * x + 32. / 3 * H00 * x2 + 8. / 3 * H1
                    - 32. / 3 * H1 * x + 16. / 3 * H10 - 32. / 3 * H10 * x
                    + 32. / 3 * H10 * x2 + 32. / 3 * H11 - 64. / 3 * H11 * x
                    + 64. / 3 * H11 * x2 + 16. / 3 * H01 - 32. / 3 * H01 * x
                    + 64. / 3 * H01 * x2)
           + CA
                 * (-5789. / 81 + 13244. / 243. / x + 12266. / 81 * x
                    - 30335. / 243 * x2 + 112. / 27 * zeta3 / x - 16 * zeta3
                    - 200. / 9 * zeta3 * x - 5020. / 27 * zeta3 * x2
                    + 128. / 9 * zeta2 / x - 59. / 9 * zeta2
                    + 1628. / 9 * zeta2 * x - 776. / 9 * zeta2 * x2
                    - 16. / 15 * zeta2 * zeta2 - 616. / 15 * zeta2 * zeta2 * x
                    + 16 * Hm1 * zeta3 + 32 * Hm1 * zeta3 * x
                    + 32 * Hm1 * zeta3 * x2 + 32. / 3 * Hm1 * zeta2 * x
                    + 32. / 3 * Hm1 * zeta2 * x2 - 32. / 3 * Hm1m1 * zeta2
                    - 64. / 3 * Hm1m1 * zeta2 * x - 64. / 3 * Hm1m1 * zeta2 * x2
                    - 64. / 3 * Hm1m1m10 - 128. / 3 * Hm1m1m10 * x
                    - 128. / 3 * Hm1m1m10 * x2 + 64. / 3 * Hm1m10 * x
                    + 64. / 3 * Hm1m10 * x2 + 32. / 3 * Hm1m100
                    + 64. / 3 * Hm1m100 * x + 64. / 3 * Hm1m100 * x2
                    - 32 * Hm10 * x - 32 * Hm10 * x2 - 20. / 3 * Hm10 * zeta2
                    - 40. / 3 * Hm10 * zeta2 * x - 40. / 3 * Hm10 * zeta2 * x2
                    - 32. / 3 * Hm100 * x - 32. / 3 * Hm100 * x2
                    + 16. / 3 * Hm1000 + 32. / 3 * Hm1000 * x
                    + 32. / 3 * Hm1000 * x2 + 32. / 3 * Hm1010
                    + 64. / 3 * Hm1010 * x + 64. / 3 * Hm1010 * x2
                    + 4000. / 81 * H0 - 8780. / 81 * H0 * x
                    + 32606. / 81 * H0 * x2 + 224. / 9 * H0 * zeta3
                    + 656. / 9 * H0 * zeta3 * x + 2. / 3 * H0 * zeta2
                    + 80 * H0 * zeta2 * x + 544. / 9 * H0 * zeta2 * x2
                    - 916. / 27 * H00 - 2368. / 27 * H00 * x
                    - 1096. / 9 * H00 * x2 - 12 * H00 * zeta2
                    - 8. / 3 * H00 * zeta2 * x + 8. / 3 * H000
                    - 16. / 3 * H000 * x + 184. / 9 * H000 * x2
                    - 16. / 3 * H0000 - 32. / 3 * H0000 * x + 212. / 81 * H1
                    - 1256. / 81 * H1 / x - 13756. / 81 * H1 * x
                    + 4546. / 27 * H1 * x2 + 40. / 9 * H1 * zeta3
                    - 80. / 9 * H1 * zeta3 * x + 80. / 9 * H1 * zeta3 * x2
                    - 40. / 9 * H1 * zeta2 + 200. / 9 * H1 * zeta2 * x
                    - 200. / 9 * H1 * zeta2 * x2 + 584. / 27 * H10
                    - 320. / 27 * H10 / x - 4588. / 27 * H10 * x
                    + 1504. / 9 * H10 * x2 + 8 * H100 + 64. / 9 * H100 / x
                    + 128. / 3 * H100 * x - 520. / 9 * H100 * x2
                    + 268. / 27 * H11 - 572. / 27 * H11 * x
                    + 824. / 27 * H11 * x2 - 20. / 3 * H11 * zeta2
                    + 40. / 3 * H11 * zeta2 * x - 40. / 3 * H11 * zeta2 * x2
                    - 32. / 3 * H110 * x + 32. / 3 * H110 * x2 - 16. / 3 * H1100
                    + 32. / 3 * H1100 * x - 32. / 3 * H1100 * x2 + 8. / 3 * H111
                    + 32. / 3 * H111 * x - 40. / 3 * H111 * x2 + 32. / 3 * H1110
                    - 64. / 3 * H1110 * x + 64. / 3 * H1110 * x2
                    - 16. / 3 * H1111 + 32. / 3 * H1111 * x
                    - 32. / 3 * H1111 * x2 + 32. / 3 * H1101
                    - 64. / 3 * H1101 * x + 64. / 3 * H1101 * x2
                    - 32. / 3 * H101 * x + 32. / 3 * H101 * x2 + 32. / 3 * H1010
                    - 64. / 3 * H1010 * x + 64. / 3 * H1010 * x2
                    + 16. / 3 * H1011 - 32. / 3 * H1011 * x
                    + 32. / 3 * H1011 * x2 - 8. / 9 * H01 - 844. / 9 * H01 * x
                    - 32 * H01 * x2 - 16. / 3 * H010 - 136. / 3 * H010 * x
                    - 704. / 9 * H010 * x2 + 32. / 3 * H0100
                    + 128. / 3 * H0100 * x - 8 * H011 * x - 8. / 3 * H011 * x2
                    + 76. / 9 * H001 - 224. / 9 * H001 * x + 32. / 3 * H0010
                    + 64. / 3 * H0010 * x + 8. / 3 * H0001
                    - 16. / 3 * H0001 * x)
           + CA * Lmmu
                 * (752. / 27 + 160. / 27. / x + 1960. / 27 * x
                    - 3196. / 27 * x2 + 80. / 3 * zeta3 + 32. / 3 * zeta3 * x
                    + 32 * zeta3 * x2 - 64. / 9 * zeta2 / x + 32. / 3 * zeta2
                    - 544. / 3 * zeta2 * x + 584. / 3 * zeta2 * x2
                    + 64. / 3 * H0m10 * x2 - 32. / 3 * Hm1 * zeta2
                    - 64. / 3 * Hm1 * zeta2 * x - 32 * Hm1 * zeta2 * x2
                    - 64. / 3 * Hm1m10 * x2 - 32 * Hm10 - 64. / 9 * Hm10 / x
                    + 32. / 3 * Hm10 * x + 416. / 9 * Hm10 * x2 + 16 * Hm100
                    + 32 * Hm100 * x + 128. / 3 * Hm100 * x2 + 32. / 3 * Hm101
                    + 64. / 3 * Hm101 * x + 64. / 3 * Hm101 * x2 + 524. / 9 * H0
                    + 704. / 3 * H0 * x - 11920. / 27 * H0 * x2
                    - 32. / 3 * H0 * zeta2 - 256. / 3 * H0 * zeta2 * x
                    + 64. / 3 * H0 * zeta2 * x2 + 736. / 3 * H00 * x
                    - 152 * H00 * x2 + 64. / 3 * H000 + 64 * H000 * x
                    + 164. / 9 * H1 - 416. / 27 * H1 / x + 1840. / 9 * H1 * x
                    - 6352. / 27 * H1 * x2 + 32. / 3 * H1 * zeta2
                    - 64. / 3 * H1 * zeta2 * x + 32. / 3 * H1 * zeta2 * x2
                    + 8. / 3 * H10 + 128. / 9 * H10 / x + 448. / 3 * H10 * x
                    - 1592. / 9 * H10 * x2 - 16 * H100 + 32 * H100 * x
                    - 64. / 3 * H100 * x2 - 8 * H11 + 64. / 9 * H11 / x
                    + 256. / 3 * H11 * x - 856. / 9 * H11 * x2 - 64. / 3 * H110
                    + 128. / 3 * H110 * x - 128. / 3 * H110 * x2
                    - 32. / 3 * H101 + 64. / 3 * H101 * x - 64. / 3 * H101 * x2
                    - 32. / 3 * H01 + 192 * H01 * x - 584. / 3 * H01 * x2
                    + 32. / 3 * H010 + 320. / 3 * H010 * x - 64. / 3 * H010 * x2
                    + 64 * H011 * x - 64. / 3 * H011 * x2 + 32. / 3 * H001
                    + 256. / 3 * H001 * x - 64. / 3 * H001 * x2)
           + CA * Lmmu2
                 * (-344. / 9 + 32. / 9. / x - 1936. / 9 * x + 2248. / 9 * x2
                    + 64 * zeta2 * x - 64. / 3 * zeta2 * x2 - 16. / 3 * H0
                    - 512. / 3 * H0 * x + 496. / 9 * H0 * x2 - 32. / 3 * H00
                    - 128. / 3 * H00 * x + 16. / 3 * H1 - 64. / 9 * H1 / x
                    - 128 * H1 * x + 1264. / 9 * H1 * x2 + 32. / 3 * H10
                    - 64. / 3 * H10 * x + 64. / 3 * H10 * x2 + 64. / 3 * H11
                    - 128. / 3 * H11 * x + 128. / 3 * H11 * x2 - 64 * H01 * x
                    + 64. / 3 * H01 * x2)
           + CA * LQm
                 * (646. / 27 + 824. / 81. / x + 15704. / 27 * x
                    - 51314. / 81 * x2 + 8 * zeta3 - 176. / 3 * zeta3 * x
                    + 16 * zeta3 * x2 - 32. / 9 * zeta2 / x + 8. / 3 * zeta2
                    - 1472. / 9 * zeta2 * x + 280. / 3 * zeta2 * x2
                    + 32. / 3 * H0m10 * x2 - 8. / 3 * Hm1 * zeta2
                    - 16. / 3 * Hm1 * zeta2 * x - 32. / 3 * Hm1 * zeta2 * x2
                    + 16. / 3 * Hm1m10 + 32. / 3 * Hm1m10 * x - 64. / 9 * Hm10
                    - 32. / 9 * Hm10 / x + 112. / 9 * Hm10 * x
                    + 272. / 9 * Hm10 * x2 + 32. / 3 * Hm100
                    + 64. / 3 * Hm100 * x + 80. / 3 * Hm100 * x2
                    + 16. / 3 * Hm101 + 32. / 3 * Hm101 * x
                    + 32. / 3 * Hm101 * x2 + 496. / 27 * H0
                    + 13312. / 27 * H0 * x - 5200. / 27 * H0 * x2
                    - 160. / 3 * H0 * zeta2 * x + 32. / 3 * H0 * zeta2 * x2
                    + 80. / 9 * H00 + 2368. / 9 * H00 * x - 56 * H00 * x2
                    + 224. / 3 * H000 * x - 520. / 27 * H1 + 160. / 27 * H1 / x
                    + 7712. / 27 * H1 * x - 2752. / 9 * H1 * x2
                    + 8. / 3 * H1 * zeta2 - 16. / 3 * H1 * zeta2 * x
                    - 104. / 9 * H10 + 32. / 9 * H10 / x + 640. / 9 * H10 * x
                    - 232. / 3 * H10 * x2 - 32. / 3 * H100 + 64. / 3 * H100 * x
                    - 16 * H100 * x2 - 104. / 9 * H11 + 32. / 9 * H11 / x
                    + 640. / 9 * H11 * x - 232. / 3 * H11 * x2 - 32. / 3 * H110
                    + 64. / 3 * H110 * x - 64. / 3 * H110 * x2 - 16. / 3 * H111
                    + 32. / 3 * H111 * x - 32. / 3 * H111 * x2 - 8. / 3 * H01
                    + 176 * H01 * x - 280. / 3 * H01 * x2 + 32 * H010 * x
                    - 32. / 3 * H010 * x2 + 32 * H011 * x - 32. / 3 * H011 * x2
                    + 160. / 3 * H001 * x - 32. / 3 * H001 * x2)
           + CA * LQm * Lmmu
                 * (-380. / 9 + 416. / 27. / x - 856. / 9 * x + 3616. / 27 * x2
                    + 128. / 3 * zeta2 * x - 64. / 3 * zeta2 * x2
                    - 32. / 3 * Hm10 - 64. / 3 * Hm10 * x - 64. / 3 * Hm10 * x2
                    - 128 * H0 * x + 400. / 3 * H0 * x2 - 64. / 3 * H00
                    - 64 * H00 * x + 16. / 3 * H1 - 64. / 9 * H1 / x
                    - 320. / 3 * H1 * x + 1072. / 9 * H1 * x2 + 32. / 3 * H10
                    - 64. / 3 * H10 * x + 64. / 3 * H10 * x2 + 32. / 3 * H11
                    - 64. / 3 * H11 * x + 64. / 3 * H11 * x2 - 64 * H01 * x
                    + 64. / 3 * H01 * x2)
           + CA * LQm * Lmmu2
                 * (16. / 3 + 64. / 9. / x + 128. / 3 * x - 496. / 9 * x2
                    + 32. / 3 * H0 + 128. / 3 * H0 * x - 32. / 3 * H1
                    + 64. / 3 * H1 * x - 64. / 3 * H1 * x2)
           + CA * LQm2
                 * (-58. / 9 - 80. / 27. / x - 596. / 9 * x + 2096. / 27 * x2
                    + 32. / 3 * zeta2 * x - 16. / 3 * zeta2 * x2 - 8. / 3 * Hm10
                    - 16. / 3 * Hm10 * x - 16. / 3 * Hm10 * x2 - 52. / 9 * H0
                    - 544. / 9 * H0 * x + 196. / 9 * H0 * x2 - 8. / 3 * H00
                    - 64. / 3 * H00 * x + 52. / 9 * H1 - 16. / 9 * H1 / x
                    - 320. / 9 * H1 * x + 116. / 3 * H1 * x2 + 8. / 3 * H10
                    - 16. / 3 * H10 * x + 16. / 3 * H10 * x2 + 8. / 3 * H11
                    - 16. / 3 * H11 * x + 16. / 3 * H11 * x2 - 16 * H01 * x
                    + 16. / 3 * H01 * x2)
           + CA * LQm2 * Lmmu
                 * (8. / 3 + 32. / 9. / x + 64. / 3 * x - 248. / 9 * x2
                    + 16. / 3 * H0 + 64. / 3 * H0 * x - 16. / 3 * H1
                    + 32. / 3 * H1 * x - 32. / 3 * H1 * x2)
           + CA * LQm3
                 * (4. / 9 + 16. / 27. / x + 32. / 9 * x - 124. / 27 * x2
                    + 8. / 9 * H0 + 32. / 9 * H0 * x - 8. / 9 * H1
                    + 16. / 9 * H1 * x - 16. / 9 * H1 * x2)
           + CA * nf
                 * (-3778. / 243 + 3272. / 729. / x - 18572. / 243 * x
                    + 46084. / 729 * x2 + 128. / 27 * zeta3 / x
                    - 128. / 9 * zeta3 - 356. / 9 * zeta3 * x
                    - 1712. / 27 * zeta3 * x2 - 8. / 9 * zeta2
                    + 1328. / 81 * zeta2 * x + 16. / 9 * zeta2 * x2
                    + 64. / 5 * zeta2 * zeta2 + 224. / 9 * zeta2 * zeta2 * x
                    + 32. / 9 * H0m10 * x - 16 * H0m10 * x2
                    - 16. / 3 * Hm1 * zeta3 - 32. / 3 * Hm1 * zeta3 * x
                    - 32. / 3 * Hm1 * zeta3 * x2 + 40. / 9 * Hm1 * zeta2
                    + 32. / 9 * Hm1 * zeta2 * x + 32. / 9 * Hm1 * zeta2 * x2
                    + 16. / 3 * Hm10m10 + 32. / 3 * Hm10m10 * x
                    + 32. / 3 * Hm10m10 * x2 + 8. / 3 * Hm1m1 * zeta2
                    + 16. / 3 * Hm1m1 * zeta2 * x + 16. / 3 * Hm1m1 * zeta2 * x2
                    + 16. / 3 * Hm1m1m10 + 32. / 3 * Hm1m1m10 * x
                    + 32. / 3 * Hm1m1m10 * x2 + 80. / 9 * Hm1m10
                    + 64. / 9 * Hm1m10 * x + 64. / 9 * Hm1m10 * x2
                    + 8. / 3 * Hm1m100 + 16. / 3 * Hm1m100 * x
                    + 16. / 3 * Hm1m100 * x2 + 448. / 81 * Hm10
                    + 1880. / 81 * Hm10 * x + 1808. / 81 * Hm10 * x2
                    + 8. / 3 * Hm10 * zeta2 + 16. / 3 * Hm10 * zeta2 * x
                    + 16. / 3 * Hm10 * zeta2 * x2 - 40. / 27 * Hm100
                    - 128. / 27 * Hm100 * x - 128. / 27 * Hm100 * x2
                    - 56. / 9 * Hm1000 - 112. / 9 * Hm1000 * x
                    - 112. / 9 * Hm1000 * x2 - 16. / 3 * Hm1010
                    - 32. / 3 * Hm1010 * x - 32. / 3 * Hm1010 * x2
                    + 1012. / 243 * H0 - 3560. / 243 * H0 * x
                    - 7312. / 81 * H0 * x2 + 64. / 9 * H0 * zeta3
                    - 128. / 9 * H0 * zeta3 * x - 8 * H0 * zeta2 * x2
                    + 344. / 81 * H00 + 616. / 27 * H00 * x
                    + 4588. / 81 * H00 * x2 + 64. / 27 * H000
                    + 356. / 27 * H000 * x + 352. / 27 * H000 * x2
                    - 16. / 9 * H0000 + 64. / 3 * H0000 * x + 1100. / 243 * H1
                    + 16. / 27 * H1 / x + 3680. / 243 * H1 * x
                    - 3608. / 243 * H1 * x2 - 16 * H1 * zeta3
                    + 32 * H1 * zeta3 * x - 32 * H1 * zeta3 * x2
                    - 40. / 27 * H1 * zeta2 - 16. / 27 * H1 * zeta2 * x
                    + 16. / 27 * H1 * zeta2 * x2 + 4. / 9 * H10
                    - 208. / 27 * H10 / x - 220. / 9 * H10 * x
                    + 868. / 27 * H10 * x2 + 8. / 9 * H10 * zeta2
                    - 16. / 9 * H10 * zeta2 * x + 16. / 9 * H10 * zeta2 * x2
                    - 256. / 27 * H100 - 64. / 9 * H100 / x
                    - 1168. / 27 * H100 * x + 1576. / 27 * H100 * x2
                    - 8. / 3 * H1000 + 16. / 3 * H1000 * x
                    - 16. / 3 * H1000 * x2 + 520. / 81 * H11
                    + 28. / 81 * H11 * x - 352. / 81 * H11 * x2
                    + 40. / 9 * H11 * zeta2 - 80. / 9 * H11 * zeta2 * x
                    + 80. / 9 * H11 * zeta2 * x2 + 80. / 27 * H110
                    - 112. / 27 * H110 * x + 112. / 27 * H110 * x2
                    + 56. / 9 * H1100 - 112. / 9 * H1100 * x
                    + 112. / 9 * H1100 * x2 - 232. / 27 * H111
                    + 80. / 27 * H111 * x - 8. / 27 * H111 * x2 - 8. / 9 * H1110
                    + 16. / 9 * H1110 * x - 16. / 9 * H1110 * x2
                    + 40. / 9 * H1111 - 80. / 9 * H1111 * x
                    + 80. / 9 * H1111 * x2 - 64. / 9 * H1101
                    + 128. / 9 * H1101 * x - 128. / 9 * H1101 * x2
                    + 160. / 27 * H101 - 80. / 27 * H101 * x
                    + 80. / 27 * H101 * x2 - 40. / 9 * H1010
                    + 80. / 9 * H1010 * x - 80. / 9 * H1010 * x2
                    - 40. / 9 * H1011 + 80. / 9 * H1011 * x
                    - 80. / 9 * H1011 * x2 + 16. / 9 * H1001
                    - 32. / 9 * H1001 * x + 32. / 9 * H1001 * x2 + 8. / 9 * H01
                    + 184. / 27 * H01 * x - 16. / 9 * H01 * x2 - 80. / 9 * H010
                    - 112. / 3 * H010 * x + 8 * H010 * x2 - 32. / 3 * H0100
                    - 128. / 3 * H0100 * x + 4. / 9 * H011 * x
                    + 32. / 3 * H011 * x2 + 32. / 9 * H001 * x + 8 * H001 * x2
                    - 64. / 3 * H0010 * x)
           + CA * nf * Lmmu
                 * (848. / 27 - 104. / 27. / x + 11800. / 27 * x
                    - 12652. / 27 * x2 + 56. / 3 * zeta3 - 16. / 3 * zeta3 * x
                    + 16 * zeta3 * x2 - 32. / 9 * zeta2 / x + 8. / 3 * zeta2
                    - 496. / 3 * zeta2 * x + 92 * zeta2 * x2
                    + 32. / 3 * H0m10 * x2 - 16. / 3 * Hm1 * zeta2
                    - 32. / 3 * Hm1 * zeta2 * x - 16 * Hm1 * zeta2 * x2
                    - 32. / 3 * Hm1m10 * x2 - 16 * Hm10 - 32. / 9 * Hm10 / x
                    + 16. / 3 * Hm10 * x + 208. / 9 * Hm10 * x2 + 8 * Hm100
                    + 16 * Hm100 * x + 64. / 3 * Hm100 * x2 + 16. / 3 * Hm101
                    + 32. / 3 * Hm101 * x + 32. / 3 * Hm101 * x2 + 116. / 9 * H0
                    + 3404. / 9 * H0 * x - 5848. / 27 * H0 * x2
                    - 160. / 3 * H0 * zeta2 * x + 32. / 3 * H0 * zeta2 * x2
                    + 104. / 9 * H00 + 1904. / 9 * H00 * x - 476. / 9 * H00 * x2
                    + 16. / 3 * H000 + 128. / 3 * H000 * x - 8 * H1
                    + 160. / 27 * H1 / x + 2324. / 9 * H1 * x
                    - 7480. / 27 * H1 * x2 + 16. / 3 * H1 * zeta2
                    - 32. / 3 * H1 * zeta2 * x + 16. / 3 * H1 * zeta2 * x2
                    - 68. / 9 * H10 + 64. / 9 * H10 / x + 832. / 9 * H10 * x
                    - 956. / 9 * H10 * x2 - 8 * H100 + 16 * H100 * x
                    - 32. / 3 * H100 * x2 - 196. / 9 * H11 + 32. / 9 * H11 / x
                    + 704. / 9 * H11 * x - 748. / 9 * H11 * x2 - 32. / 3 * H110
                    + 64. / 3 * H110 * x - 64. / 3 * H110 * x2 - 16. / 3 * H101
                    + 32. / 3 * H101 * x - 32. / 3 * H101 * x2 - 8. / 3 * H01
                    + 512. / 3 * H01 * x - 92 * H01 * x2 + 16. / 3 * H010
                    + 160. / 3 * H010 * x - 32. / 3 * H010 * x2 + 32 * H011 * x
                    - 32. / 3 * H011 * x2 + 160. / 3 * H001 * x
                    - 32. / 3 * H001 * x2)
           + CA * nf * Lmmu2
                 * (-86. / 9 + 8. / 9. / x - 484. / 9 * x + 562. / 9 * x2
                    + 16 * zeta2 * x - 16. / 3 * zeta2 * x2 - 4. / 3 * H0
                    - 128. / 3 * H0 * x + 124. / 9 * H0 * x2 - 8. / 3 * H00
                    - 32. / 3 * H00 * x + 4. / 3 * H1 - 16. / 9 * H1 / x
                    - 32 * H1 * x + 316. / 9 * H1 * x2 + 8. / 3 * H10
                    - 16. / 3 * H10 * x + 16. / 3 * H10 * x2 + 16. / 3 * H11
                    - 32. / 3 * H11 * x + 32. / 3 * H11 * x2 - 16 * H01 * x
                    + 16. / 3 * H01 * x2)
           + CA * nf * LQm
                 * (362. / 9 - 32. / 9. / x + 4112. / 9 * x - 4612. / 9 * x2
                    + 8 * zeta3 - 176. / 3 * zeta3 * x + 16 * zeta3 * x2
                    - 32. / 9 * zeta2 / x + 8. / 3 * zeta2
                    - 1448. / 9 * zeta2 * x + 280. / 3 * zeta2 * x2
                    + 32. / 3 * H0m10 * x2 - 8. / 3 * Hm1 * zeta2
                    - 16. / 3 * Hm1 * zeta2 * x - 32. / 3 * Hm1 * zeta2 * x2
                    + 16. / 3 * Hm1m10 + 32. / 3 * Hm1m10 * x - 64. / 9 * Hm10
                    - 32. / 9 * Hm10 / x + 112. / 9 * Hm10 * x
                    + 272. / 9 * Hm10 * x2 + 32. / 3 * Hm100
                    + 64. / 3 * Hm100 * x + 80. / 3 * Hm100 * x2
                    + 16. / 3 * Hm101 + 32. / 3 * Hm101 * x
                    + 32. / 3 * Hm101 * x2 + 356. / 27 * H0
                    + 11288. / 27 * H0 * x - 6080. / 27 * H0 * x2
                    - 160. / 3 * H0 * zeta2 * x + 32. / 3 * H0 * zeta2 * x2
                    + 52. / 3 * H00 + 2168. / 9 * H00 * x - 56 * H00 * x2
                    + 8. / 3 * H000 + 208. / 3 * H000 * x - 332. / 27 * H1
                    + 160. / 27 * H1 / x + 7228. / 27 * H1 * x
                    - 7736. / 27 * H1 * x2 + 8. / 3 * H1 * zeta2
                    - 16. / 3 * H1 * zeta2 * x - 104. / 9 * H10
                    + 32. / 9 * H10 / x + 640. / 9 * H10 * x
                    - 232. / 3 * H10 * x2 - 32. / 3 * H100 + 64. / 3 * H100 * x
                    - 16 * H100 * x2 - 104. / 9 * H11 + 32. / 9 * H11 / x
                    + 640. / 9 * H11 * x - 232. / 3 * H11 * x2 - 32. / 3 * H110
                    + 64. / 3 * H110 * x - 64. / 3 * H110 * x2 - 16. / 3 * H111
                    + 32. / 3 * H111 * x - 32. / 3 * H111 * x2 - 8. / 3 * H01
                    + 520. / 3 * H01 * x - 280. / 3 * H01 * x2 + 32 * H010 * x
                    - 32. / 3 * H010 * x2 + 32 * H011 * x - 32. / 3 * H011 * x2
                    + 160. / 3 * H001 * x - 32. / 3 * H001 * x2)
           + CA * nf * LQm * Lmmu
                 * (-116. / 9 - 160. / 27. / x - 1192. / 9 * x + 4192. / 27 * x2
                    + 64. / 3 * zeta2 * x - 32. / 3 * zeta2 * x2
                    - 16. / 3 * Hm10 - 32. / 3 * Hm10 * x - 32. / 3 * Hm10 * x2
                    - 104. / 9 * H0 - 1088. / 9 * H0 * x + 392. / 9 * H0 * x2
                    - 16. / 3 * H00 - 128. / 3 * H00 * x + 104. / 9 * H1
                    - 32. / 9 * H1 / x - 640. / 9 * H1 * x + 232. / 3 * H1 * x2
                    + 16. / 3 * H10 - 32. / 3 * H10 * x + 32. / 3 * H10 * x2
                    + 16. / 3 * H11 - 32. / 3 * H11 * x + 32. / 3 * H11 * x2
                    - 32 * H01 * x + 32. / 3 * H01 * x2)
           + CA * nf * LQm * Lmmu2
                 * (4. / 3 + 16. / 9. / x + 32. / 3 * x - 124. / 9 * x2
                    + 8. / 3 * H0 + 32. / 3 * H0 * x - 8. / 3 * H1
                    + 16. / 3 * H1 * x - 16. / 3 * H1 * x2)
           + CA * nf * LQm2
                 * (-58. / 9 - 80. / 27. / x - 596. / 9 * x + 2096. / 27 * x2
                    + 32. / 3 * zeta2 * x - 16. / 3 * zeta2 * x2 - 8. / 3 * Hm10
                    - 16. / 3 * Hm10 * x - 16. / 3 * Hm10 * x2 - 52. / 9 * H0
                    - 544. / 9 * H0 * x + 196. / 9 * H0 * x2 - 8. / 3 * H00
                    - 64. / 3 * H00 * x + 52. / 9 * H1 - 16. / 9 * H1 / x
                    - 320. / 9 * H1 * x + 116. / 3 * H1 * x2 + 8. / 3 * H10
                    - 16. / 3 * H10 * x + 16. / 3 * H10 * x2 + 8. / 3 * H11
                    - 16. / 3 * H11 * x + 16. / 3 * H11 * x2 - 16 * H01 * x
                    + 16. / 3 * H01 * x2)
           + CA * nf * LQm2 * Lmmu
                 * (4. / 3 + 16. / 9. / x + 32. / 3 * x - 124. / 9 * x2
                    + 8. / 3 * H0 + 32. / 3 * H0 * x - 8. / 3 * H1
                    + 16. / 3 * H1 * x - 16. / 3 * H1 * x2)
           + CA * nf * LQm3
                 * (4. / 9 + 16. / 27. / x + 32. / 9 * x - 124. / 27 * x2
                    + 8. / 9 * H0 + 32. / 9 * H0 * x - 8. / 9 * H1
                    + 16. / 9 * H1 * x - 16. / 9 * H1 * x2)
           + CA * CF
                 * (2299. / 162 - 466. / 9. / x + 49333. / 81 * x
                    - 9098. / 27 * x2 - 2 * zeta5 + 348 * zeta5 * x
                    - 424 * zeta5 * x2 + 11 * zeta3 - 1660 * zeta3 * x
                    + 1360. / 3 * zeta3 * x2 - 40 * zeta2 / x
                    + 1283. / 36 * zeta2 - 5515. / 18 * zeta2 * x
                    - 2750. / 9 * zeta2 * x2 + 96 * zeta2 * ln2
                    - 192 * zeta2 * ln2 * x + 192 * zeta2 * ln2 * x2
                    - 250. / 3 * zeta2 * zeta3 - 1276. / 3 * zeta2 * zeta3 * x
                    + 472. / 3 * zeta2 * zeta3 * x2 + 202. / 15 * zeta2 * zeta2
                    - 1324. / 3 * zeta2 * zeta2 * x
                    + 4708. / 15 * zeta2 * zeta2 * x2 - 8 * H00m1 * zeta2
                    - 16 * H00m1 * zeta2 * x - 32 * H00m1 * zeta2 * x2
                    - 16 * H00m1m10 - 32 * H00m1m10 * x - 64 * H00m1m10 * x2
                    + 32 * H00m10 * x + 32 * H00m10 * x2 + 8 * H00m100
                    + 16 * H00m100 * x + 32 * H00m100 * x2 - 4 * H0m1 * zeta3
                    - 8 * H0m1 * zeta3 * x - 16 * H0m1 * zeta3 * x2
                    + 4 * H0m1 * zeta2 - 48 * H0m1 * zeta2 * x
                    - 64 * H0m1 * zeta2 * x2 - 16 * H0m10m10 - 32 * H0m10m10 * x
                    - 64 * H0m10m10 * x2 + 16 * H0m1m1 * zeta2
                    + 32 * H0m1m1 * zeta2 * x + 64 * H0m1m1 * zeta2 * x2
                    + 8 * H0m1m10 - 32 * H0m1m10 * x - 64 * H0m1m10 * x2
                    - 8 * H0m1m100 - 16 * H0m1m100 * x - 32 * H0m1m100 * x2
                    - 16 * H0m1m101 - 32 * H0m1m101 * x - 64 * H0m1m101 * x2
                    + 16 * H0m10 + 56 * H0m10 * x - 8 * H0m10 * x2
                    - 20 * H0m10 * zeta2 - 40 * H0m10 * zeta2 * x
                    - 80 * H0m10 * zeta2 * x2 - 4 * H0m100 + 32 * H0m100 * x
                    + 48 * H0m100 * x2 + 16 * H0m1000 + 32 * H0m1000 * x
                    + 64 * H0m1000 * x2 + 32 * H0m101 * x + 32 * H0m101 * x2
                    + 16 * H0m1010 + 32 * H0m1010 * x + 64 * H0m1010 * x2
                    + 8 * H0m1001 + 16 * H0m1001 * x + 32 * H0m1001 * x2
                    - 86 * Hm1 * zeta3 - 112 * Hm1 * zeta3 * x
                    - 104 * Hm1 * zeta3 * x2 - 36 * Hm1 * zeta2
                    - 4 * Hm1 * zeta2 * x + 20 * Hm1 * zeta2 * x2
                    + 36 * Hm1 * zeta2 * zeta2 + 72 * Hm1 * zeta2 * zeta2 * x
                    + 72 * Hm1 * zeta2 * zeta2 * x2 - 32 * Hm100m10
                    - 64 * Hm100m10 * x - 64 * Hm100m10 * x2
                    + 32 * Hm10m1 * zeta2 + 64 * Hm10m1 * zeta2 * x
                    + 64 * Hm10m1 * zeta2 * x2 - 56 * Hm10m10 - 96 * Hm10m10 * x
                    - 64 * Hm10m10 * x2 - 16 * Hm10m100 - 32 * Hm10m100 * x
                    - 32 * Hm10m100 * x2 - 32 * Hm10m101 - 64 * Hm10m101 * x
                    - 64 * Hm10m101 * x2 - 56 * Hm1m1 * zeta3
                    - 112 * Hm1m1 * zeta3 * x - 112 * Hm1m1 * zeta3 * x2
                    + 104 * Hm1m1 * zeta2 + 144 * Hm1m1 * zeta2 * x
                    + 112 * Hm1m1 * zeta2 * x2 + 48 * Hm1m1m1 * zeta2
                    + 96 * Hm1m1m1 * zeta2 * x + 96 * Hm1m1m1 * zeta2 * x2
                    + 96 * Hm1m1m1m10 + 192 * Hm1m1m1m10 * x
                    + 192 * Hm1m1m1m10 * x2 + 96 * Hm1m1m10 + 96 * Hm1m1m10 * x
                    + 96 * Hm1m1m10 * x2 - 48 * Hm1m1m100 - 96 * Hm1m1m100 * x
                    - 96 * Hm1m1m100 * x2 - 40 * Hm1m10 + 8 * Hm1m10 * x
                    + 24 * Hm1m10 * x2 + 24 * Hm1m10 * zeta2
                    + 48 * Hm1m10 * zeta2 * x + 48 * Hm1m10 * zeta2 * x2
                    - 76 * Hm1m100 - 96 * Hm1m100 * x - 80 * Hm1m100 * x2
                    - 32 * Hm1m1000 - 64 * Hm1m1000 * x - 64 * Hm1m1000 * x2
                    - 56 * Hm1m101 - 96 * Hm1m101 * x - 64 * Hm1m101 * x2
                    - 64 * Hm1m1010 - 128 * Hm1m1010 * x - 128 * Hm1m1010 * x2
                    - 32 * Hm1m1011 - 64 * Hm1m1011 * x - 64 * Hm1m1011 * x2
                    - 16 * Hm1m1001 - 32 * Hm1m1001 * x - 32 * Hm1m1001 * x2
                    - 84 * Hm10 - 160 * Hm10 * x - 140 * Hm10 * x2
                    - 8 * Hm10 * zeta3 - 16 * Hm10 * zeta3 * x
                    - 16 * Hm10 * zeta3 * x2 - 52 * Hm10 * zeta2
                    - 108 * Hm10 * zeta2 * x - 68 * Hm10 * zeta2 * x2
                    + 28 * Hm100 - 16 * Hm100 * x2 - 60 * Hm100 * zeta2
                    - 120 * Hm100 * zeta2 * x - 120 * Hm100 * zeta2 * x2
                    + 32 * Hm1000 + 80 * Hm1000 * x + 48 * Hm1000 * x2
                    + 32 * Hm10000 + 64 * Hm10000 * x + 64 * Hm10000 * x2
                    + 16 * Hm101 + 8 * Hm101 * x - 8 * Hm101 * x2
                    - 48 * Hm101 * zeta2 - 96 * Hm101 * zeta2 * x
                    - 96 * Hm101 * zeta2 * x2 + 8 * Hm1010 + 64 * Hm1010 * x
                    + 32 * Hm1010 * x2 + 32 * Hm10100 + 64 * Hm10100 * x
                    + 64 * Hm10100 * x2 + 32 * Hm1011 * x + 32 * Hm1011 * x2
                    + 32 * Hm10110 + 64 * Hm10110 * x + 64 * Hm10110 * x2
                    + 32 * Hm10101 + 64 * Hm10101 * x + 64 * Hm10101 * x2
                    + 28 * Hm1001 + 64 * Hm1001 * x + 48 * Hm1001 * x2
                    + 48 * Hm10010 + 96 * Hm10010 * x + 96 * Hm10010 * x2
                    + 16 * Hm10011 + 32 * Hm10011 * x + 32 * Hm10011 * x2
                    + 32 * Hm10001 + 64 * Hm10001 * x + 64 * Hm10001 * x2
                    - 5833. / 27 * H0 + 11437. / 27 * H0 * x
                    - 9514. / 81 * H0 * x2 - 262. / 9 * H0 * zeta3
                    - 268. / 9 * H0 * zeta3 * x - 3656. / 9 * H0 * zeta3 * x2
                    - 347. / 9 * H0 * zeta2 - 1409. / 9 * H0 * zeta2 * x
                    + 1180. / 9 * H0 * zeta2 * x2 - 88. / 5 * H0 * zeta2 * zeta2
                    - 824. / 5 * H0 * zeta2 * zeta2 * x - 346. / 9 * H00
                    - 2987. / 9 * H00 * x + 10622. / 9 * H00 * x2
                    + 160. / 3 * H00 * zeta3 + 208. / 3 * H00 * zeta3 * x
                    + 32 * H00 * zeta3 * x2 + 49. / 3 * H00 * zeta2
                    + 106. / 3 * H00 * zeta2 * x - 128. / 3 * H00 * zeta2 * x2
                    - 377. / 3 * H000 - 584. / 3 * H000 * x
                    - 7624. / 9 * H000 * x2 + 32 * H000 * zeta2
                    + 40 * H000 * zeta2 * x + 32 * H000 * zeta2 * x2
                    + 26. / 3 * H0000 + 356. / 3 * H0000 * x
                    + 32. / 3 * H0000 * x2 - 8 * H00000 - 7516. / 27 * H1
                    + 730. / 81 * H1 / x + 29350. / 27 * H1 * x
                    - 60778. / 81 * H1 * x2 - 128. / 9 * H1 * zeta3 / x
                    + 130. / 9 * H1 * zeta3 - 2072. / 9 * H1 * zeta3 * x
                    + 2176. / 9 * H1 * zeta3 * x2 + 16 * H1 * zeta2 / x
                    - 1186. / 9 * H1 * zeta2 + 7292. / 9 * H1 * zeta2 * x
                    - 5504. / 9 * H1 * zeta2 * x2
                    + 868. / 5 * H1 * zeta2 * zeta2
                    - 1736. / 5 * H1 * zeta2 * zeta2 * x
                    + 1736. / 5 * H1 * zeta2 * zeta2 * x2 - 934. / 9 * H10
                    + 1480. / 27 * H10 / x - 8336. / 9 * H10 * x
                    + 25538. / 27 * H10 * x2 + 232. / 3 * H10 * zeta3
                    - 464. / 3 * H10 * zeta3 * x + 464. / 3 * H10 * zeta3 * x2
                    - 32 * H10 * zeta2 / x - 94. / 3 * H10 * zeta2
                    - 484. / 3 * H10 * zeta2 * x + 628. / 3 * H10 * zeta2 * x2
                    + 282 * H100 - 80. / 9 * H100 / x - 928. / 3 * H100 * x
                    + 248. / 9 * H100 * x2 + 4 * H100 * zeta2
                    - 8 * H100 * zeta2 * x + 8 * H100 * zeta2 * x2
                    + 100. / 3 * H1000 + 160. / 3 * H1000 / x
                    + 1336. / 3 * H1000 * x - 1616. / 3 * H1000 * x2
                    - 16 * H10000 + 32 * H10000 * x - 32 * H10000 * x2
                    - 166 * H11 - 284. / 27 * H11 / x + 688. / 3 * H11 * x
                    - 4738. / 27 * H11 * x2 + 56. / 3 * H11 * zeta3
                    - 112. / 3 * H11 * zeta3 * x + 112. / 3 * H11 * zeta3 * x2
                    + 32. / 3 * H11 * zeta2 / x - 76 * H11 * zeta2
                    + 408 * H11 * zeta2 * x - 1160. / 3 * H11 * zeta2 * x2
                    + 404. / 3 * H110 - 32. / 9 * H110 / x
                    - 1424. / 3 * H110 * x + 2912. / 9 * H110 * x2
                    - 80 * H110 * zeta2 + 160 * H110 * zeta2 * x
                    - 160 * H110 * zeta2 * x2 + 128. / 3 * H1100 / x
                    + 368 * H1100 * x - 1208. / 3 * H1100 * x2 - 48 * H11000
                    + 96 * H11000 * x - 96 * H11000 * x2 - 96 * H111
                    - 16 * H111 / x + 376 * H111 * x - 228 * H111 * x2
                    - 192 * H111 * zeta2 + 384 * H111 * zeta2 * x
                    - 384 * H111 * zeta2 * x2 + 4 * H1110 + 64. / 3 * H1110 / x
                    + 64 * H1110 * x - 376. / 3 * H1110 * x2 - 16 * H11100
                    + 32 * H11100 * x - 32 * H11100 * x2 + 448. / 3 * H1111
                    + 32. / 3 * H1111 / x - 248. / 3 * H1111 * x
                    - 32 * H1111 * x2 + 80 * H11110 - 160 * H11110 * x
                    + 160 * H11110 * x2 + 144 * H11101 - 288 * H11101 * x
                    + 288 * H11101 * x2 + 52. / 3 * H1101 - 368. / 3 * H1101 * x
                    + 296. / 3 * H1101 * x2 + 112 * H11010 - 224 * H11010 * x
                    + 224 * H11010 * x2 + 128 * H11011 - 256 * H11011 * x
                    + 256 * H11011 * x2 + 32 * H11001 - 64 * H11001 * x
                    + 64 * H11001 * x2 + 132 * H101 - 32. / 9 * H101 / x
                    - 512 * H101 * x + 3272. / 9 * H101 * x2
                    - 144 * H101 * zeta2 + 288 * H101 * zeta2 * x
                    - 288 * H101 * zeta2 * x2 + 16. / 3 * H1010
                    + 64. / 3 * H1010 / x + 208. / 3 * H1010 * x
                    - 344. / 3 * H1010 * x2 + 16 * H10100 - 32 * H10100 * x
                    + 32 * H10100 * x2 + 64. / 3 * H1011 / x + 80 * H1011 * x
                    - 376. / 3 * H1011 * x2 + 64 * H10110 - 128 * H10110 * x
                    + 128 * H10110 * x2 + 32 * H10111 - 64 * H10111 * x
                    + 64 * H10111 * x2 + 96 * H10101 - 192 * H10101 * x
                    + 192 * H10101 * x2 + 4 * H1001 + 128. / 3 * H1001 / x
                    + 352 * H1001 * x - 1208. / 3 * H1001 * x2 + 32 * H10010
                    - 64 * H10010 * x + 64 * H10010 * x2 - 32 * H10001
                    + 64 * H10001 * x - 64 * H10001 * x2 - 124. / 9 * H01
                    + 1690. / 9 * H01 * x + 2318. / 9 * H01 * x2
                    - 164. / 3 * H01 * zeta3 - 56. / 3 * H01 * zeta3 * x
                    - 112. / 3 * H01 * zeta3 * x2 - 144 * H01 * zeta2
                    + 64 * H01 * zeta2 * x + 328. / 3 * H01 * zeta2 * x2
                    - 94. / 3 * H010 - 1444. / 3 * H010 * x - 320 * H010 * x2
                    - 112 * H010 * zeta2 - 64 * H010 * zeta2 * x
                    - 96 * H010 * zeta2 * x2 + 284. / 3 * H0100
                    + 968. / 3 * H0100 * x - 1720. / 3 * H0100 * x2
                    + 96 * H01000 + 288 * H01000 * x - 224. / 3 * H011
                    + 1756. / 3 * H011 * x - 6764. / 9 * H011 * x2
                    - 48 * H011 * zeta2 + 192 * H011 * zeta2 * x
                    - 160 * H011 * zeta2 * x2 + 76 * H0110 + 208 * H0110 * x
                    - 424 * H0110 * x2 + 64 * H01100 + 256 * H01100 * x
                    - 64 * H01100 * x2 - 8. / 3 * H0111 + 904. / 3 * H0111 * x
                    - 128 * H0111 * x2 + 24 * H01110 + 144 * H01110 * x
                    - 24 * H01111 + 144 * H01111 * x - 32 * H01111 * x2
                    + 24 * H01101 - 48 * H01101 * x + 64 * H01101 * x2
                    + 364. / 3 * H0101 + 448. / 3 * H0101 * x - 200 * H0101 * x2
                    + 64 * H01010 + 64 * H01010 * x + 64 * H01010 * x2
                    + 80 * H01011 + 32 * H01011 * x + 96 * H01011 * x2
                    + 104 * H01001 + 176 * H01001 * x + 32 * H01001 * x2
                    + 34. / 3 * H001 + 776. / 3 * H001 * x
                    - 3736. / 9 * H001 * x2 - 40 * H001 * zeta2
                    + 80 * H001 * zeta2 * x - 96 * H001 * zeta2 * x2
                    - 64. / 3 * H0010 + 248. / 3 * H0010 * x
                    - 344. / 3 * H0010 * x2 + 64 * H00100 + 96 * H00100 * x
                    + 64 * H00100 * x2 - 12 * H0011 + 208 * H0011 * x
                    + 328. / 3 * H0011 * x2 + 56 * H00110 + 112 * H00110 * x
                    + 32 * H00110 * x2 + 24 * H00111 + 48 * H00111 * x
                    + 32 * H00111 * x2 + 24 * H00101 - 16 * H00101 * x
                    + 32 * H00101 * x2 - 4 * H0001 + 56 * H0001 * x
                    + 152. / 3 * H0001 * x2 - 16 * H00011 - 32 * H00011 * x
                    - 32 * H00001 - 32 * H00001 * x - 32 * H00001 * x2)
           + CA * CF * Lmmu
                 * (-7174. / 45 + 392. / 9. / x + 1966. / 5 * x - 314 * x2
                    + 64. / 3 * zeta3 / x - 404. / 3 * zeta3
                    + 1336. / 3 * zeta3 * x - 1328 * zeta3 * x2
                    - 192 * zeta3 * x3 - 64. / 15 * zeta2 / x
                    - 884. / 15 * zeta2 + 65416. / 45 * zeta2 * x
                    - 7964. / 5 * zeta2 * x2 - 704. / 5 * zeta2 * x3
                    + 472. / 5 * zeta2 * zeta2 - 288. / 5 * zeta2 * zeta2 * x
                    + 704. / 5 * zeta2 * zeta2 * x2 + 128 * H00m10
                    + 128 * H0m1 * zeta2 + 256 * H0m1 * zeta2 * x
                    + 192 * H0m1 * zeta2 * x2 + 512 * H0m1m10 * x
                    + 128 * H0m1m10 * x2 + 32. / 3 * H0m10
                    - 128. / 3 * H0m10 * x - 544 * H0m10 * x2
                    - 384. / 5 * H0m10 * x3 - 64 * H0m100 - 256 * H0m100 * x
                    - 128 * H0m100 * x2 - 128 * H0m101 - 128 * H0m101 * x2
                    + 224 * Hm1 * zeta3 + 448 * Hm1 * zeta3 * x
                    + 224 * Hm1 * zeta3 * x2 + 16. / 5 * Hm1 * zeta2 / x2
                    + 64. / 3 * Hm1 * zeta2 / x + 848. / 3 * Hm1 * zeta2
                    + 1120. / 3 * Hm1 * zeta2 * x + 272 * Hm1 * zeta2 * x2
                    + 576. / 5 * Hm1 * zeta2 * x3 + 128 * Hm10m10
                    + 256 * Hm10m10 * x + 128 * Hm10m10 * x2
                    - 256 * Hm1m1 * zeta2 - 512 * Hm1m1 * zeta2 * x
                    - 256 * Hm1m1 * zeta2 * x2 - 256 * Hm1m1m10
                    - 512 * Hm1m1m10 * x - 256 * Hm1m1m10 * x2
                    + 544. / 3 * Hm1m10 + 32. / 15 * Hm1m10 / x2
                    + 128. / 3 * Hm1m10 / x + 576 * Hm1m10 * x
                    + 544 * Hm1m10 * x2 + 384. / 5 * Hm1m10 * x3 + 192 * Hm1m100
                    + 384 * Hm1m100 * x + 192 * Hm1m100 * x2 + 128 * Hm1m101
                    + 256 * Hm1m101 * x + 128 * Hm1m101 * x2 + 1072. / 5 * Hm10
                    - 176. / 45 * Hm10 / x2 - 32. / 15 * Hm10 / x
                    + 22048. / 45 * Hm10 * x + 536. / 5 * Hm10 * x2
                    - 704. / 5 * Hm10 * x3 + 128 * Hm10 * zeta2
                    + 256 * Hm10 * zeta2 * x + 128 * Hm10 * zeta2 * x2
                    - 560. / 3 * Hm100 - 32. / 15 * Hm100 / x2
                    - 64. / 3 * Hm100 / x - 992. / 3 * Hm100 * x
                    - 272 * Hm100 * x2 - 384. / 5 * Hm100 * x3 - 64 * Hm1000
                    - 128 * Hm1000 * x - 64 * Hm1000 * x2 - 192 * Hm101
                    - 32. / 15 * Hm101 / x2 - 256. / 3 * Hm101 * x
                    - 384. / 5 * Hm101 * x3 - 64 * Hm1001 - 128 * Hm1001 * x
                    - 64 * Hm1001 * x2 + 104. / 45 * H0 + 16. / 9 * H0 / x
                    - 180 * H0 * x + 6104. / 5 * H0 * x2 + 80 * H0 * zeta3
                    - 128 * H0 * zeta3 * x2 - 176. / 3 * H0 * zeta2
                    + 2080. / 3 * H0 * zeta2 * x - 880 * H0 * zeta2 * x2
                    - 768. / 5 * H0 * zeta2 * x3 - 1316. / 15 * H00
                    + 32. / 15 * H00 / x - 23608. / 45 * H00 * x
                    + 3604. / 5 * H00 * x2 + 704. / 5 * H00 * x3
                    + 32 * H00 * zeta2 + 320 * H00 * zeta2 * x
                    - 128 * H00 * zeta2 * x2 + 88. / 3 * H000 - 176 * H000 * x
                    + 448 * H000 * x2 + 384. / 5 * H000 * x3 - 32 * H0000
                    - 128 * H0000 * x + 1378. / 15 * H1 - 392. / 15 * H1 / x
                    - 5346. / 5 * H1 * x + 4744. / 5 * H1 * x2
                    - 224 * H1 * zeta3 * x2 + 16. / 15 * H1 * zeta2 / x2
                    + 128. / 3 * H1 * zeta2 / x - 400. / 3 * H1 * zeta2
                    + 848 * H1 * zeta2 * x - 944 * H1 * zeta2 * x2
                    - 192. / 5 * H1 * zeta2 * x3 + 200. / 3 * H10
                    - 160. / 3 * H10 / x - 976 * H10 * x + 1020 * H10 * x2
                    - 32 * H10 * zeta2 + 64 * H10 * zeta2 * x
                    - 192 * H10 * zeta2 * x2 + 128 * H100 - 64. / 3 * H100 / x
                    - 352 * H100 * x + 384 * H100 * x2 + 64 * H1000 * x2
                    + 122. / 3 * H11 - 32 * H11 / x - 5048. / 3 * H11 * x
                    + 1708 * H11 * x2 - 128 * H11 * zeta2
                    + 256 * H11 * zeta2 * x - 384 * H11 * zeta2 * x2
                    + 632. / 3 * H110 - 128. / 3 * H110 / x
                    - 2848. / 3 * H110 * x + 992 * H110 * x2 + 64 * H1100
                    - 128 * H1100 * x + 192 * H1100 * x2 + 304 * H111
                    - 64 * H111 / x - 1616 * H111 * x + 1680 * H111 * x2
                    + 224 * H1110 - 448 * H1110 * x + 448 * H1110 * x2
                    + 384 * H1111 - 768 * H1111 * x + 768 * H1111 * x2
                    + 256 * H1101 - 512 * H1101 * x + 512 * H1101 * x2
                    + 224 * H101 - 64 * H101 / x - 1136 * H101 * x
                    + 1216 * H101 * x2 + 160 * H1010 - 320 * H1010 * x
                    + 320 * H1010 * x2 + 288 * H1011 - 576 * H1011 * x
                    + 576 * H1011 * x2 + 96 * H1001 - 192 * H1001 * x
                    + 256 * H1001 * x2 + 884. / 15 * H01 + 32. / 15 * H01 / x
                    - 14456. / 15 * H01 * x + 7964. / 5 * H01 * x2
                    - 48 * H01 * zeta2 + 416 * H01 * zeta2 * x
                    - 256 * H01 * zeta2 * x2 + 128. / 3 * H010
                    - 2656. / 3 * H010 * x + 784 * H010 * x2 + 64 * H0100
                    - 128 * H0100 * x + 128 * H0100 * x2 + 364. / 3 * H011
                    - 4424. / 3 * H011 * x + 1344 * H011 * x2 + 80 * H0110
                    - 544 * H0110 * x + 320 * H0110 * x2 + 144 * H0111
                    - 864 * H0111 * x + 576 * H0111 * x2 + 48 * H0101
                    - 672 * H0101 * x + 320 * H0101 * x2 + 176. / 3 * H001
                    - 736 * H001 * x + 880 * H001 * x2 + 384. / 5 * H001 * x3
                    - 384 * H0010 * x + 192 * H0010 * x2 + 48 * H0011
                    - 576 * H0011 * x + 384 * H0011 * x2 - 32 * H0001
                    - 320 * H0001 * x + 128 * H0001 * x2)
           + CA * CF * LQm
                 * (6041. / 15 + 2512. / 45. / x + 39737. / 45 * x
                    - 67432. / 45 * x2 - 128. / 3 * zeta3 / x - 196. / 3 * zeta3
                    - 2584. / 3 * zeta3 * x + 8 * zeta3 * x2 - 192 * zeta3 * x3
                    - 64. / 15 * zeta2 / x - 1154. / 15 * zeta2
                    + 45256. / 45 * zeta2 * x - 59936. / 45 * zeta2 * x2
                    - 896. / 5 * zeta2 * x3 + 156. / 5 * zeta2 * zeta2
                    - 936. / 5 * zeta2 * zeta2 * x
                    + 752. / 5 * zeta2 * zeta2 * x2 + 224 * H00m10
                    - 64 * H00m10 * x + 128 * H00m10 * x2 + 232 * H0m1 * zeta2
                    + 464 * H0m1 * zeta2 * x + 480 * H0m1 * zeta2 * x2
                    + 16 * H0m1m10 + 800 * H0m1m10 * x + 320 * H0m1m10 * x2
                    + 1088. / 3 * H0m10 - 560. / 3 * H0m10 * x
                    - 640 * H0m10 * x2 - 384. / 5 * H0m10 * x3 - 104 * H0m100
                    - 400 * H0m100 * x - 256 * H0m100 * x2 - 224 * H0m101
                    - 64 * H0m101 * x - 320 * H0m101 * x2 + 408 * Hm1 * zeta3
                    + 816 * Hm1 * zeta3 * x + 592 * Hm1 * zeta3 * x2
                    + 16. / 5 * Hm1 * zeta2 / x2 + 64. / 3 * Hm1 * zeta2 / x
                    + 896. / 3 * Hm1 * zeta2 + 1720. / 3 * Hm1 * zeta2 * x
                    + 480 * Hm1 * zeta2 * x2 + 576. / 5 * Hm1 * zeta2 * x3
                    + 192 * Hm10m10 + 384 * Hm10m10 * x + 256 * Hm10m10 * x2
                    - 448 * Hm1m1 * zeta2 - 896 * Hm1m1 * zeta2 * x
                    - 640 * Hm1m1 * zeta2 * x2 - 384 * Hm1m1m10
                    - 768 * Hm1m1m10 * x - 512 * Hm1m1m10 * x2
                    - 176. / 3 * Hm1m10 + 32. / 15 * Hm1m10 / x2
                    + 128. / 3 * Hm1m10 / x + 432 * Hm1m10 * x
                    + 640 * Hm1m10 * x2 + 384. / 5 * Hm1m10 * x3 + 288 * Hm1m100
                    + 576 * Hm1m100 * x + 384 * Hm1m100 * x2 + 256 * Hm1m101
                    + 512 * Hm1m101 * x + 384 * Hm1m101 * x2 + 4312. / 5 * Hm10
                    - 224. / 45 * Hm10 / x2 - 32. / 15 * Hm10 / x
                    + 48568. / 45 * Hm10 * x + 496. / 5 * Hm10 * x2
                    - 896. / 5 * Hm10 * x3 + 336 * Hm10 * zeta2
                    + 672 * Hm10 * zeta2 * x + 544 * Hm10 * zeta2 * x2
                    - 236. / 3 * Hm100 - 32. / 15 * Hm100 / x2
                    - 64. / 3 * Hm100 / x - 800. / 3 * Hm100 * x
                    - 304 * Hm100 * x2 - 384. / 5 * Hm100 * x3 - 176 * Hm1000
                    - 352 * Hm1000 * x - 288 * Hm1000 * x2 - 328 * Hm101
                    - 32. / 15 * Hm101 / x2 - 1072. / 3 * Hm101 * x
                    - 160 * Hm101 * x2 - 384. / 5 * Hm101 * x3 - 64 * Hm1010
                    - 128 * Hm1010 * x - 128 * Hm1010 * x2 - 64 * Hm1011
                    - 128 * Hm1011 * x - 128 * Hm1011 * x2 - 192 * Hm1001
                    - 384 * Hm1001 * x - 320 * Hm1001 * x2 + 1378. / 9 * H0
                    + 128. / 45 * H0 / x + 23516. / 45 * H0 * x
                    + 484. / 45 * H0 * x2 - 48 * H0 * zeta3
                    - 928 * H0 * zeta3 * x - 64 * H0 * zeta3 * x2
                    - 68 * H0 * zeta2 + 256 * H0 * zeta2 * x
                    - 1616 * H0 * zeta2 * x2 - 768. / 5 * H0 * zeta2 * x3
                    - 2408. / 45 * H00 + 32. / 15 * H00 / x
                    - 14636. / 15 * H00 * x + 65096. / 45 * H00 * x2
                    + 896. / 5 * H00 * x3 + 120 * H00 * zeta2
                    + 432 * H00 * zeta2 * x - 96 * H00 * zeta2 * x2 + 120 * H000
                    - 976. / 3 * H000 * x + 3688. / 3 * H000 * x2
                    + 384. / 5 * H000 * x3 - 80 * H0000 - 288 * H0000 * x
                    + 10349. / 45 * H1 + 544. / 45 * H1 / x
                    + 22856. / 45 * H1 * x - 32084. / 45 * H1 * x2
                    + 392 * H1 * zeta3 - 784 * H1 * zeta3 * x
                    + 560 * H1 * zeta3 * x2 + 16. / 15 * H1 * zeta2 / x2
                    + 64 * H1 * zeta2 / x - 264 * H1 * zeta2
                    + 2536. / 3 * H1 * zeta2 * x - 2576. / 3 * H1 * zeta2 * x2
                    - 192. / 5 * H1 * zeta2 * x3 + 96 * H10m10
                    - 192 * H10m10 * x + 192 * H10m10 * x2 + 892. / 9 * H10
                    + 256. / 9 * H10 / x - 6080. / 9 * H10 * x
                    + 7336. / 9 * H10 * x2 + 128 * H10 * zeta2
                    - 256 * H10 * zeta2 * x + 128 * H10 * zeta2 * x2
                    + 1288. / 3 * H100 - 128. / 3 * H100 / x
                    - 2912. / 3 * H100 * x + 808 * H100 * x2 + 16 * H1000
                    - 32 * H1000 * x + 96 * H1000 * x2 + 520. / 9 * H11
                    + 256. / 9 * H11 / x - 2768. / 9 * H11 * x
                    + 4144. / 9 * H11 * x2 + 64 * H11 * zeta2
                    - 128 * H11 * zeta2 * x + 184 * H110 - 256. / 3 * H110 / x
                    - 976 * H110 * x + 3184. / 3 * H110 * x2 + 112 * H1100
                    - 224 * H1100 * x + 288 * H1100 * x2 + 208 * H111
                    - 64 * H111 / x - 656 * H111 * x + 744 * H111 * x2
                    + 160 * H1110 - 320 * H1110 * x + 320 * H1110 * x2
                    + 128 * H1101 - 256 * H1101 * x + 256 * H1101 * x2
                    + 704. / 3 * H101 - 256. / 3 * H101 / x
                    - 3184. / 3 * H101 * x + 3536. / 3 * H101 * x2 + 160 * H1010
                    - 320 * H1010 * x + 320 * H1010 * x2 + 96 * H1011
                    - 192 * H1011 * x + 192 * H1011 * x2 + 16 * H1001
                    - 32 * H1001 * x + 96 * H1001 * x2 + 1154. / 15 * H01
                    + 32. / 15 * H01 / x + 368. / 5 * H01 * x
                    + 59936. / 45 * H01 * x2 + 24 * H01 * zeta2
                    + 336 * H01 * zeta2 * x - 96 * H01 * zeta2 * x2
                    + 260. / 3 * H010 - 2224. / 3 * H010 * x
                    + 4496. / 3 * H010 * x2 + 112 * H0100 - 224 * H0100 * x
                    + 288 * H0100 * x2 + 428. / 3 * H011 - 2224. / 3 * H011 * x
                    + 4024. / 3 * H011 * x2 - 768 * H0110 * x + 320 * H0110 * x2
                    - 576 * H0111 * x + 192 * H0111 * x2 - 16 * H0101
                    - 736 * H0101 * x + 256 * H0101 * x2 + 68 * H001
                    - 1328. / 3 * H001 * x + 1616 * H001 * x2
                    + 384. / 5 * H001 * x3 - 96 * H0010 - 640 * H0010 * x
                    + 192 * H0010 * x2 - 48 * H0011 - 640 * H0011 * x
                    + 256 * H0011 * x2 - 120 * H0001 - 496 * H0001 * x
                    + 96 * H0001 * x2)
           + CA * CF * LQm * Lmmu
                 * (142. / 3 + 16. / x + 234 * x - 268 * x2 + 48 * zeta3
                    + 128 * zeta3 * x2 + 28 * zeta2 - 552 * zeta2 * x
                    + 608 * zeta2 * x2 + 188. / 3 * H0 + 992. / 3 * H0 * x
                    - 656 * H0 * x2 - 16 * H0 * zeta2 - 256 * H0 * zeta2 * x
                    + 128 * H0 * zeta2 * x2 - 88. / 3 * H00 + 656. / 3 * H00 * x
                    - 448 * H00 * x2 + 32 * H000 + 128 * H000 * x - 30 * H1
                    + 64. / 3 * H1 / x + 2120. / 3 * H1 * x - 688 * H1 * x2
                    + 128 * H1 * zeta2 - 256 * H1 * zeta2 * x
                    + 256 * H1 * zeta2 * x2 - 248. / 3 * H10
                    + 128. / 3 * H10 / x + 1696. / 3 * H10 * x - 608 * H10 * x2
                    - 64 * H100 + 128 * H100 * x - 128 * H100 * x2
                    - 416. / 3 * H11 + 128. / 3 * H11 / x + 2272. / 3 * H11 * x
                    - 768 * H11 * x2 - 128 * H110 + 256 * H110 * x
                    - 256 * H110 * x2 - 192 * H111 + 384 * H111 * x
                    - 384 * H111 * x2 - 128 * H101 + 256 * H101 * x
                    - 256 * H101 * x2 - 28 * H01 + 552 * H01 * x
                    - 608 * H01 * x2 + 16 * H010 + 352 * H010 * x
                    - 128 * H010 * x2 - 32 * H011 + 448 * H011 * x
                    - 256 * H011 * x2 + 16 * H001 + 256 * H001 * x
                    - 128 * H001 * x2)
           + CA * CF * LQm2
                 * (-193. / 6 + 8. / x + 268. / 3 * x - 26. / 3 * x2
                    + 32 * zeta3 * x2 + 88. / 3 * zeta2 - 680. / 3 * zeta2 * x
                    + 1256. / 3 * zeta2 * x2 - 16 * Hm1 * zeta2
                    - 32 * Hm1 * zeta2 * x - 32 * Hm1 * zeta2 * x2 - 12 * Hm10
                    - 24 * Hm10 * x - 24 * Hm10 * x2 + 8 * Hm100
                    + 16 * Hm100 * x + 16 * Hm100 * x2 + 16 * Hm101
                    + 32 * Hm101 * x + 32 * Hm101 * x2 + 109. / 9 * H0
                    + 1126. / 9 * H0 * x - 2320. / 9 * H0 * x2 - 24 * H0 * zeta2
                    - 160 * H0 * zeta2 * x + 64 * H0 * zeta2 * x2
                    - 112. / 3 * H00 + 392. / 3 * H00 * x - 376 * H00 * x2
                    + 24 * H000 + 96 * H000 * x - 382. / 9 * H1
                    - 64. / 9 * H1 / x + 2156. / 9 * H1 * x
                    - 2104. / 9 * H1 * x2 + 48 * H1 * zeta2
                    - 96 * H1 * zeta2 * x + 96 * H1 * zeta2 * x2
                    - 212. / 3 * H10 + 64. / 3 * H10 / x + 928. / 3 * H10 * x
                    - 992. / 3 * H10 * x2 - 40 * H100 + 80 * H100 * x
                    - 80 * H100 * x2 - 80 * H11 + 64. / 3 * H11 / x
                    + 328 * H11 * x - 1048. / 3 * H11 * x2 - 48 * H110
                    + 96 * H110 * x - 96 * H110 * x2 - 48 * H111 + 96 * H111 * x
                    - 96 * H111 * x2 - 48 * H101 + 96 * H101 * x
                    - 96 * H101 * x2 - 88. / 3 * H01 + 608. / 3 * H01 * x
                    - 1256. / 3 * H01 * x2 + 8 * H010 + 176 * H010 * x
                    - 64 * H010 * x2 - 8 * H011 + 208 * H011 * x
                    - 96 * H011 * x2 + 24 * H001 + 160 * H001 * x
                    - 64 * H001 * x2)
           + CA * CF * LQm2 * Lmmu
                 * (-41. / 3 - 52. / 3 * x + 20 * x2 + 8 * zeta2
                    + 80 * zeta2 * x - 32 * zeta2 * x2 + 22. / 3 * H0
                    - 116. / 3 * H0 * x + 112 * H0 * x2 - 8 * H00 - 32 * H00 * x
                    + 32. / 3 * H1 - 32. / 3 * H1 / x - 328. / 3 * H1 * x
                    + 112 * H1 * x2 + 16 * H10 - 32 * H10 * x + 32 * H10 * x2
                    + 32 * H11 - 64 * H11 * x + 64 * H11 * x2 - 8 * H01
                    - 80 * H01 * x + 32 * H01 * x2)
           + CA * CF * LQm3
                 * (-10. / 3 - 32. / 3 * x + 20. / 3 * x2 + 8. / 3 * zeta2
                    + 80. / 3 * zeta2 * x - 32. / 3 * zeta2 * x2 + 44. / 9 * H0
                    - 160. / 9 * H0 * x + 424. / 9 * H0 * x2 - 8. / 3 * H00
                    - 32. / 3 * H00 * x + 76. / 9 * H1 - 32. / 9 * H1 / x
                    - 416. / 9 * H1 * x + 424. / 9 * H1 * x2 + 16. / 3 * H10
                    - 32. / 3 * H10 * x + 32. / 3 * H10 * x2 + 32. / 3 * H11
                    - 64. / 3 * H11 * x + 64. / 3 * H11 * x2 - 8. / 3 * H01
                    - 80. / 3 * H01 * x + 32. / 3 * H01 * x2)
           + CA * CA
                 * (6154. / 9 - 133214. / 243. / x - 261964. / 27 * x
                    + 2317604. / 243 * x2 - 136 * zeta5 + 1088 * zeta5 * x
                    + 896. / 27 * zeta3 / x + 262 * zeta3
                    + 5584. / 3 * zeta3 * x + 32728. / 27 * zeta3 * x2
                    - 2176. / 27 * zeta2 / x + 170 * zeta2
                    - 3403. / 9 * zeta2 * x + 57548. / 27 * zeta2 * x2
                    + 332. / 3 * zeta2 * zeta3 + 1112. / 3 * zeta2 * zeta3 * x
                    + 608. / 15 * zeta2 * zeta2 / x - 48 * zeta2 * zeta2
                    + 13384. / 15 * zeta2 * zeta2 * x
                    + 688. / 15 * zeta2 * zeta2 * x2 + 48 * H0m1 * zeta3
                    - 192 * H0m1 * zeta3 * x - 32 * H0m1 * zeta2 * x
                    - 48 * H0m1 * zeta2 * x2 - 32 * H0m1m1 * zeta2
                    + 128 * H0m1m1 * zeta2 * x - 64 * H0m1m1m10
                    + 256 * H0m1m1m10 * x - 64 * H0m1m10 * x - 96 * H0m1m10 * x2
                    + 32 * H0m1m100 - 128 * H0m1m100 * x + 192 * H0m10 * x
                    + 32 * H0m10 * x2 - 24 * H0m10 * zeta2
                    + 96 * H0m10 * zeta2 * x + 32 * H0m100 * x
                    + 48 * H0m100 * x2 + 16 * H0m1000 - 64 * H0m1000 * x
                    + 32 * H0m1010 - 128 * H0m1010 * x - 32 * Hm1 * zeta3 / x
                    - 20 * Hm1 * zeta3 - 392 * Hm1 * zeta3 * x
                    - 448 * Hm1 * zeta3 * x2 - 160. / 9 * Hm1 * zeta2 / x
                    + 32. / 3 * Hm1 * zeta2 - 1016. / 3 * Hm1 * zeta2 * x
                    - 3304. / 9 * Hm1 * zeta2 * x2
                    - 188. / 5 * Hm1 * zeta2 * zeta2
                    - 376. / 5 * Hm1 * zeta2 * zeta2 * x
                    - 376. / 5 * Hm1 * zeta2 * zeta2 * x2 + 32 * Hm10m1 * zeta2
                    + 64 * Hm10m1 * zeta2 * x + 64 * Hm10m1 * zeta2 * x2
                    + 64 * Hm10m1m10 + 128 * Hm10m1m10 * x
                    + 128 * Hm10m1m10 * x2 - 64 * Hm10m10 * x
                    - 64 * Hm10m10 * x2 - 32 * Hm10m100 - 64 * Hm10m100 * x
                    - 64 * Hm10m100 * x2 + 160 * Hm1m1 * zeta3
                    + 320 * Hm1m1 * zeta3 * x + 320 * Hm1m1 * zeta3 * x2
                    + 64. / 3 * Hm1m1 * zeta2 / x + 40. / 3 * Hm1m1 * zeta2
                    + 944. / 3 * Hm1m1 * zeta2 * x + 352 * Hm1m1 * zeta2 * x2
                    + 64 * Hm1m10m10 + 128 * Hm1m10m10 * x
                    + 128 * Hm1m10m10 * x2 - 160 * Hm1m1m1 * zeta2
                    - 320 * Hm1m1m1 * zeta2 * x - 320 * Hm1m1m1 * zeta2 * x2
                    - 192 * Hm1m1m1m10 - 384 * Hm1m1m1m10 * x
                    - 384 * Hm1m1m1m10 * x2 + 80. / 3 * Hm1m1m10
                    + 128. / 3 * Hm1m1m10 / x + 1504. / 3 * Hm1m1m10 * x
                    + 576 * Hm1m1m10 * x2 + 128 * Hm1m1m100
                    + 256 * Hm1m1m100 * x + 256 * Hm1m1m100 * x2
                    + 64 * Hm1m1m101 + 128 * Hm1m1m101 * x
                    + 128 * Hm1m1m101 * x2 + 64. / 3 * Hm1m10
                    - 320. / 9 * Hm1m10 / x - 1456. / 3 * Hm1m10 * x
                    - 4880. / 9 * Hm1m10 * x2 + 40 * Hm1m10 * zeta2
                    + 80 * Hm1m10 * zeta2 * x + 80 * Hm1m10 * zeta2 * x2
                    - 40. / 3 * Hm1m100 - 64. / 3 * Hm1m100 / x
                    - 848. / 3 * Hm1m100 * x - 320 * Hm1m100 * x2
                    - 16 * Hm1m1000 - 32 * Hm1m1000 * x - 32 * Hm1m1000 * x2
                    - 64 * Hm1m101 * x - 64 * Hm1m101 * x2 + 32 * Hm1m1010
                    + 64 * Hm1m1010 * x + 64 * Hm1m1010 * x2 - 32 * Hm1m1001
                    - 64 * Hm1m1001 * x - 64 * Hm1m1001 * x2 - 400. / 9 * Hm10
                    + 752. / 27 * Hm10 / x + 5864. / 9 * Hm10 * x
                    + 19544. / 27 * Hm10 * x2 - 32 * Hm10 * zeta3
                    - 64 * Hm10 * zeta3 * x - 64 * Hm10 * zeta3 * x2
                    + 16 * Hm10 * zeta2 / x + 8. / 3 * Hm10 * zeta2
                    + 184. / 3 * Hm10 * zeta2 * x + 268. / 3 * Hm10 * zeta2 * x2
                    - 32. / 3 * Hm100 + 160. / 9 * Hm100 / x
                    + 872. / 3 * Hm100 * x + 2872. / 9 * Hm100 * x2
                    + 36 * Hm100 * zeta2 + 72 * Hm100 * zeta2 * x
                    + 72 * Hm100 * zeta2 * x2 - 20. / 3 * Hm1000
                    - 32. / 3 * Hm1000 / x - 184. / 3 * Hm1000 * x
                    - 80 * Hm1000 * x2 - 16 * Hm10000 - 32 * Hm10000 * x
                    - 32 * Hm10000 * x2 + 96 * Hm101 * x + 96 * Hm101 * x2
                    + 48 * Hm101 * zeta2 + 96 * Hm101 * zeta2 * x
                    + 96 * Hm101 * zeta2 * x2 - 40. / 3 * Hm1010
                    - 64. / 3 * Hm1010 / x - 560. / 3 * Hm1010 * x
                    - 224 * Hm1010 * x2 - 32 * Hm10100 - 64 * Hm10100 * x
                    - 64 * Hm10100 * x2 - 32 * Hm10110 - 64 * Hm10110 * x
                    - 64 * Hm10110 * x2 - 32 * Hm10101 - 64 * Hm10101 * x
                    - 64 * Hm10101 * x2 + 32 * Hm1001 * x + 32 * Hm1001 * x2
                    - 32 * Hm10010 - 64 * Hm10010 * x - 64 * Hm10010 * x2
                    - 16 * Hm10001 - 32 * Hm10001 * x - 32 * Hm10001 * x2
                    - 10918. / 27 * H0 - 5248. / 81 * H0 / x
                    - 22276. / 9 * H0 * x - 451924. / 81 * H0 * x2
                    - 32. / 9 * H0 * zeta3 / x - 1076. / 9 * H0 * zeta3
                    + 1360. / 9 * H0 * zeta3 * x + 280. / 3 * H0 * zeta3 * x2
                    - 160. / 9 * H0 * zeta2 / x - 1165. / 9 * H0 * zeta2
                    - 1864. / 9 * H0 * zeta2 * x - 2158. / 3 * H0 * zeta2 * x2
                    - 152. / 5 * H0 * zeta2 * zeta2
                    + 992. / 5 * H0 * zeta2 * zeta2 * x + 164. / 3 * H00
                    - 4382. / 9 * H00 * x + 10496. / 27 * H00 * x2
                    + 160. / 3 * H00 * zeta3 - 992. / 3 * H00 * zeta3 * x
                    + 94. / 3 * H00 * zeta2 - 4 * H00 * zeta2 * x
                    + 76. / 3 * H00 * zeta2 * x2 - 70. / 3 * H000 - 8 * H000 * x
                    + 150 * H000 * x2 - 40 * H000 * zeta2
                    + 64 * H000 * zeta2 * x + 92. / 3 * H0000
                    - 104. / 3 * H0000 * x + 24 * H0000 * x2 - 16 * H00000
                    + 64 * H00000 * x + 2572. / 27 * H1 - 2930. / 27 * H1 / x
                    - 36202. / 27 * H1 * x + 36092. / 27 * H1 * x2
                    + 320. / 9 * H1 * zeta3 / x + 164. / 9 * H1 * zeta3
                    + 2792. / 9 * H1 * zeta3 * x - 3496. / 9 * H1 * zeta3 * x2
                    - 484. / 9 * H1 * zeta2 / x + 232. / 9 * H1 * zeta2
                    - 6104. / 9 * H1 * zeta2 * x + 2074. / 3 * H1 * zeta2 * x2
                    - 92 * H1 * zeta2 * zeta2 + 184 * H1 * zeta2 * zeta2 * x
                    - 184 * H1 * zeta2 * zeta2 * x2 - 700. / 3 * H10
                    + 5552. / 27 * H10 / x + 2420 * H10 * x
                    - 64160. / 27 * H10 * x2 - 80. / 3 * H10 * zeta3
                    + 160. / 3 * H10 * zeta3 * x - 160. / 3 * H10 * zeta3 * x2
                    + 80. / 3 * H10 * zeta2 / x + 28 * H10 * zeta2
                    + 128 * H10 * zeta2 * x - 548. / 3 * H10 * zeta2 * x2
                    + 106 * H100 - 616. / 9 * H100 / x - 1192. / 3 * H100 * x
                    + 3238. / 9 * H100 * x2 - 4 * H100 * zeta2
                    + 8 * H100 * zeta2 * x - 8 * H100 * zeta2 * x2 - 24 * H1000
                    - 64. / 3 * H1000 / x - 128 * H1000 * x
                    + 520. / 3 * H1000 * x2 + 26. / 9 * H11
                    - 268. / 27 * H11 / x - 2072. / 9 * H11 * x
                    + 7072. / 27 * H11 * x2 - 208. / 3 * H11 * zeta3
                    + 416. / 3 * H11 * zeta3 * x - 416. / 3 * H11 * zeta3 * x2
                    - 32. / 3 * H11 * zeta2 / x + 44. / 3 * H11 * zeta2
                    - 376. / 3 * H11 * zeta2 * x + 136 * H11 * zeta2 * x2
                    - 88. / 3 * H110 + 160. / 3 * H110 / x
                    + 1832. / 3 * H110 * x - 1928. / 3 * H110 * x2
                    - 8 * H110 * zeta2 + 16 * H110 * zeta2 * x
                    - 16 * H110 * zeta2 * x2 - 52. / 3 * H1100 - 32 * H1100 / x
                    - 568. / 3 * H1100 * x + 760. / 3 * H1100 * x2 + 16 * H11000
                    - 32 * H11000 * x + 32 * H11000 * x2 + 10. / 3 * H111
                    - 8. / 9 * H111 / x - 320 * H111 * x + 2642. / 9 * H111 * x2
                    + 32 * H111 * zeta2 - 64 * H111 * zeta2 * x
                    + 64 * H111 * zeta2 * x2 - 64. / 3 * H1110
                    + 64. / 3 * H1110 / x + 656. / 3 * H1110 * x
                    - 248 * H1110 * x2 - 76. / 3 * H1111 - 32. / 3 * H1111 / x
                    - 664. / 3 * H1111 * x + 272 * H1111 * x2 - 80 * H11110
                    + 160 * H11110 * x - 160 * H11110 * x2 + 80 * H11111
                    - 160 * H11111 * x + 160 * H11111 * x2 - 80 * H11101
                    + 160 * H11101 * x - 160 * H11101 * x2 - 64. / 3 * H1101
                    + 64. / 3 * H1101 / x + 656. / 3 * H1101 * x
                    - 248 * H1101 * x2 - 96 * H11010 + 192 * H11010 * x
                    - 192 * H11010 * x2 - 64 * H11011 + 128 * H11011 * x
                    - 128 * H11011 * x2 - 16 * H11001 + 32 * H11001 * x
                    - 32 * H11001 * x2 - 88. / 3 * H101 + 160. / 3 * H101 / x
                    + 1832. / 3 * H101 * x - 1928. / 3 * H101 * x2
                    + 56 * H101 * zeta2 - 112 * H101 * zeta2 * x
                    + 112 * H101 * zeta2 * x2 - 40. / 3 * H1010
                    + 64. / 3 * H1010 / x + 752. / 3 * H1010 * x
                    - 288 * H1010 * x2 - 16 * H10100 + 32 * H10100 * x
                    - 32 * H10100 * x2 - 44. / 3 * H1011 + 32. / 3 * H1011 / x
                    + 376. / 3 * H1011 * x - 136 * H1011 * x2 - 80 * H10110
                    + 160 * H10110 * x - 160 * H10110 * x2 - 32 * H10111
                    + 64 * H10111 * x - 64 * H10111 * x2 - 80 * H10101
                    + 160 * H10101 * x - 160 * H10101 * x2 - 24 * H1001
                    - 64. / 3 * H1001 / x - 96 * H1001 * x
                    + 424. / 3 * H1001 * x2 - 32 * H10010 + 64 * H10010 * x
                    - 64 * H10010 * x2 - 16 * H10011 + 32 * H10011 * x
                    - 32 * H10011 * x2 - 200. / 3 * H01 - 448. / 3 * H01 * x
                    - 25784. / 27 * H01 * x2 + 160. / 3 * H01 * zeta3
                    + 640. / 3 * H01 * zeta3 * x + 16. / 3 * H01 * zeta2 / x
                    - 16 * H01 * zeta2 - 216 * H01 * zeta2 * x
                    - 348 * H01 * zeta2 * x2 + 380. / 3 * H010
                    + 320. / 9 * H010 / x + 912 * H010 * x
                    + 2080. / 3 * H010 * x2 + 40 * H010 * zeta2
                    + 160 * H010 * zeta2 * x - 16. / 3 * H0100
                    - 64. / 3 * H0100 / x - 880. / 3 * H0100 * x
                    + 776. / 3 * H0100 * x2 - 32 * H01000 - 128 * H01000 * x
                    - 56. / 3 * H011 - 376. / 3 * H011 * x
                    + 1066. / 9 * H011 * x2 - 16 * H011 * zeta2
                    - 64 * H011 * zeta2 * x + 16 * H0110 + 160 * H0110 * x
                    + 584. / 3 * H0110 * x2 - 48 * H01100 - 192 * H01100 * x
                    + 8 * H0111 - 48 * H0111 * x + 48 * H0111 * x2 + 32 * H01110
                    + 128 * H01110 * x - 16 * H01111 - 64 * H01111 * x
                    + 32 * H01101 + 128 * H01101 * x + 16 * H0101
                    + 160 * H0101 * x + 584. / 3 * H0101 * x2 + 32 * H01010
                    + 128 * H01010 * x + 16 * H01011 + 64 * H01011 * x
                    - 32 * H01001 - 128 * H01001 * x + 112. / 3 * H001
                    + 104. / 3 * H001 * x + 2192. / 9 * H001 * x2
                    + 56 * H001 * zeta2 + 128 * H001 * zeta2 * x
                    - 184. / 3 * H0010 + 304. / 3 * H0010 * x - 320 * H00100 * x
                    - 16 * H0011 * x - 24 * H0011 * x2 - 32 * H00110
                    - 64 * H00110 * x - 32 * H00101 - 64 * H00101 * x
                    - 8 * H0001 - 32 * H0001 * x - 184. / 3 * H0001 * x2
                    + 32 * H00010 - 128 * H00010 * x + 16 * H00001
                    + 32 * H00001 * x)
           + CA * CA * Lmmu
                 * (959. / 3 - 1000. / 9. / x + 6854. / 9 * x - 8587. / 9 * x2
                    - 284. / 3 * zeta3 - 1064. / 3 * zeta3 * x
                    - 2128. / 3 * zeta3 * x2 + 416. / 9 * zeta2 / x
                    - 572. / 3 * zeta2 + 8008. / 3 * zeta2 * x
                    - 39466. / 9 * zeta2 * x2 - 176. / 5 * zeta2 * zeta2
                    - 3552. / 5 * zeta2 * zeta2 * x + 96 * zeta2 * zeta2 * x2
                    + 32 * H00m10 - 128 * H00m10 * x - 64 * H0m1 * zeta2
                    + 352 * H0m1 * zeta2 * x + 128 * H0m1 * zeta2 * x2
                    + 192 * H0m1m10 * x + 128 * H0m1m10 * x2 - 80 * H0m10
                    + 64. / 3 * H0m10 / x - 576 * H0m10 * x
                    - 1024. / 3 * H0m10 * x2 + 64 * H0m100 - 352 * H0m100 * x
                    - 96 * H0m100 * x2 + 64 * H0m101 - 256 * H0m101 * x
                    - 64 * H0m101 * x2 + 120 * Hm1 * zeta3
                    + 240 * Hm1 * zeta3 * x + 352 * Hm1 * zeta3 * x2
                    + 32 * Hm1 * zeta2 / x - 368. / 3 * Hm1 * zeta2
                    + 1856. / 3 * Hm1 * zeta2 * x + 896 * Hm1 * zeta2 * x2
                    + 32 * Hm10m10 + 64 * Hm10m10 * x + 128 * Hm10m10 * x2
                    - 96 * Hm1m1 * zeta2 - 192 * Hm1m1 * zeta2 * x
                    - 320 * Hm1m1 * zeta2 * x2 - 64 * Hm1m1m10
                    - 128 * Hm1m1m10 * x - 256 * Hm1m1m10 * x2 - 48 * Hm1m10
                    + 64. / 3 * Hm1m10 / x + 544 * Hm1m10 * x
                    + 736 * Hm1m10 * x2 + 96 * Hm1m100 + 192 * Hm1m100 * x
                    + 288 * Hm1m100 * x2 + 64 * Hm1m101 + 128 * Hm1m101 * x
                    + 192 * Hm1m101 * x2 + 608. / 3 * Hm10 + 928. / 9 * Hm10 / x
                    - 3224. / 3 * Hm10 * x - 10832. / 9 * Hm10 * x2
                    + 112 * Hm10 * zeta2 + 224 * Hm10 * zeta2 * x
                    + 288 * Hm10 * zeta2 * x2 + 60 * Hm100
                    - 128. / 3 * Hm100 / x - 632 * Hm100 * x - 856 * Hm100 * x2
                    - 64 * Hm1000 - 128 * Hm1000 * x - 160 * Hm1000 * x2
                    + 296. / 3 * Hm101 - 64. / 3 * Hm101 / x
                    - 1040. / 3 * Hm101 * x - 528 * Hm101 * x2 - 32 * Hm1010
                    - 64 * Hm1010 * x - 64 * Hm1010 * x2 - 64 * Hm1011
                    - 128 * Hm1011 * x - 128 * Hm1011 * x2 - 96 * Hm1001
                    - 192 * Hm1001 * x - 224 * Hm1001 * x2 - 398. / 3 * H0
                    - 160. / 9 * H0 / x - 27626. / 9 * H0 * x
                    + 22034. / 3 * H0 * x2 + 128 * H0 * zeta3
                    - 32 * H0 * zeta3 * x + 64. / 3 * H0 * zeta2 / x
                    + 80 * H0 * zeta2 + 6832. / 3 * H0 * zeta2 * x
                    - 824 * H0 * zeta2 * x2 + 5092. / 9 * H00
                    - 29600. / 9 * H00 * x + 9590. / 3 * H00 * x2
                    + 32 * H00 * zeta2 + 832 * H00 * zeta2 * x - 184. / 3 * H000
                    - 6752. / 3 * H000 * x + 192 * H000 * x2 + 96 * H0000
                    - 512 * H0000 * x - 2104. / 9 * H1 + 2336. / 9 * H1 / x
                    - 27550. / 9 * H1 * x + 28538. / 9 * H1 * x2
                    - 8 * H1 * zeta3 + 16 * H1 * zeta3 * x
                    + 96 * H1 * zeta3 * x2 + 96 * H1 * zeta2 / x
                    - 208. / 3 * H1 * zeta2 + 3008. / 3 * H1 * zeta2 * x
                    - 3280. / 3 * H1 * zeta2 * x2 + 1262. / 9 * H10
                    - 1616. / 9 * H10 / x - 33712. / 9 * H10 * x
                    + 3938 * H10 * x2 - 112 * H10 * zeta2
                    + 224 * H10 * zeta2 * x - 160 * H10 * zeta2 * x2 - 4 * H100
                    - 256. / 3 * H100 / x - 920 * H100 * x + 1072 * H100 * x2
                    + 64 * H1000 - 128 * H1000 * x + 96 * H1000 * x2
                    + 1126. / 9 * H11 - 400. / 9 * H11 / x
                    - 23480. / 9 * H11 * x + 8414. / 3 * H11 * x2
                    - 128 * H11 * zeta2 + 256 * H11 * zeta2 * x
                    - 192 * H11 * zeta2 * x2 + 176. / 3 * H110 - 128 * H110 / x
                    - 4384. / 3 * H110 * x + 5056. / 3 * H110 * x2 + 160 * H1100
                    - 320 * H1100 * x + 288 * H1100 * x2 + 72 * H111
                    - 64 * H111 / x - 768 * H111 * x + 856 * H111 * x2
                    + 192 * H1110 - 384 * H1110 * x + 384 * H1110 * x2
                    + 160 * H1101 - 320 * H1101 * x + 320 * H1101 * x2
                    + 136. / 3 * H101 - 320. / 3 * H101 / x
                    - 3824. / 3 * H101 * x + 4384. / 3 * H101 * x2 + 160 * H1010
                    - 320 * H1010 * x + 320 * H1010 * x2 + 128 * H1011
                    - 256 * H1011 * x + 256 * H1011 * x2 + 128 * H1001
                    - 256 * H1001 * x + 224 * H1001 * x2 + 572. / 3 * H01
                    + 512. / 9 * H01 / x - 3744 * H01 * x
                    + 39466. / 9 * H01 * x2 + 64 * H01 * zeta2
                    + 736 * H01 * zeta2 * x - 128 * H01 * zeta2 * x2
                    - 160. / 3 * H010 - 64 * H010 / x - 7408. / 3 * H010 * x
                    + 3208. / 3 * H010 * x2 - 96 * H0100 - 672 * H0100 * x
                    + 96 * H0100 * x2 - 24 * H011 - 64 * H011 / x
                    - 2128 * H011 * x + 1440 * H011 * x2 - 96 * H0110
                    - 960 * H0110 * x + 192 * H0110 * x2 - 576 * H0111 * x
                    + 192 * H0111 * x2 - 64 * H0101 - 832 * H0101 * x
                    + 192 * H0101 * x2 - 80 * H001 - 8560. / 3 * H001 * x
                    + 824 * H001 * x2 - 64 * H0010 - 1152 * H0010 * x
                    + 64 * H0010 * x2 - 128 * H0011 - 1152 * H0011 * x
                    + 128 * H0011 * x2 - 32 * H0001 - 960 * H0001 * x)
           + CA * CA * Lmmu2
                 * (-259. / 9 - 284. / 9. / x + 17962. / 9 * x - 17419. / 9 * x2
                    - 16 * zeta3 - 256 * zeta3 * x - 32. / 3 * zeta2 / x
                    - 16 * zeta2 - 744 * zeta2 * x + 656. / 3 * zeta2 * x2
                    - 106 * H0 - 16. / 3 * H0 / x + 1096 * H0 * x
                    - 5026. / 9 * H0 * x2 - 16 * H0 * zeta2
                    - 256 * H0 * zeta2 * x + 44. / 3 * H00 + 1760. / 3 * H00 * x
                    - 24 * H00 * x2 - 16 * H000 + 128 * H000 * x + 18 * H1
                    + 448. / 9 * H1 / x + 3608. / 3 * H1 * x
                    - 11566. / 9 * H1 * x2 + 48 * H1 * zeta2
                    - 96 * H1 * zeta2 * x + 96 * H1 * zeta2 * x2 - 44. / 3 * H10
                    + 64. / 3 * H10 / x + 856. / 3 * H10 * x
                    - 968. / 3 * H10 * x2 - 16 * H100 + 32 * H100 * x
                    - 32 * H100 * x2 - 88. / 3 * H11 + 128. / 3 * H11 / x
                    + 1712. / 3 * H11 * x - 1936. / 3 * H11 * x2 - 48 * H110
                    + 96 * H110 * x - 96 * H110 * x2 - 96 * H111
                    + 192 * H111 * x - 192 * H111 * x2 - 48 * H101
                    + 96 * H101 * x - 96 * H101 * x2 + 16 * H01
                    + 32. / 3 * H01 / x + 744 * H01 * x - 656. / 3 * H01 * x2
                    + 16 * H010 + 160 * H010 * x - 32 * H010 * x2 + 32 * H011
                    + 320 * H011 * x - 64 * H011 * x2 + 16 * H001
                    + 256 * H001 * x)
           + CA * CA * LQm
                 * (29927. / 27 - 57320. / 81. / x - 176456. / 27 * x
                    + 503837. / 81 * x2 + 128 * zeta3 / x - 20 * zeta3
                    + 4424. / 3 * zeta3 * x - 5008. / 3 * zeta3 * x2
                    + 1088. / 9 * zeta2 / x - 360 * zeta2
                    + 29744. / 9 * zeta2 * x - 36800. / 9 * zeta2 * x2
                    + 314. / 5 * zeta2 * zeta2 - 36. / 5 * zeta2 * zeta2 * x
                    + 704. / 5 * zeta2 * zeta2 * x2 + 104 * H00m10
                    - 144 * H00m10 * x + 8 * H0m1 * zeta2
                    + 272 * H0m1 * zeta2 * x + 192 * H0m1 * zeta2 * x2
                    + 80 * H0m1m10 - 96 * H0m1m10 * x + 128 * H0m1m10 * x2
                    - 152 * H0m10 + 64. / 3 * H0m10 / x - 328 * H0m10 * x
                    - 2080. / 3 * H0m10 * x2 + 72 * H0m100 - 496 * H0m100 * x
                    - 128 * H0m100 * x2 + 32 * H0m101 - 320 * H0m101 * x
                    - 128 * H0m101 * x2 + 24 * Hm1 * zeta3
                    + 48 * Hm1 * zeta3 * x + 160 * Hm1 * zeta3 * x2
                    + 32 * Hm1 * zeta2 / x - 28. / 3 * Hm1 * zeta2
                    + 1696. / 3 * Hm1 * zeta2 * x + 2144. / 3 * Hm1 * zeta2 * x2
                    - 16 * Hm10m10 - 32 * Hm10m10 * x + 32 * Hm10m10 * x2
                    - 16 * Hm1m1 * zeta2 - 32 * Hm1m1 * zeta2 * x
                    - 160 * Hm1m1 * zeta2 * x2 + 96 * Hm1m1m10
                    + 192 * Hm1m1m10 * x + 64 * Hm1m1m10 * x2
                    + 152. / 3 * Hm1m10 - 64. / 3 * Hm1m10 / x
                    + 736. / 3 * Hm1m10 * x + 800. / 3 * Hm1m10 * x2
                    + 80 * Hm1m100 + 160 * Hm1m100 * x + 256 * Hm1m100 * x2
                    + 64 * Hm1m101 + 128 * Hm1m101 * x + 192 * Hm1m101 * x2
                    + 1024. / 9 * Hm10 + 640. / 3 * Hm10 / x
                    - 2056. / 9 * Hm10 * x - 2176. / 9 * Hm10 * x2
                    + 80 * Hm10 * zeta2 + 160 * Hm10 * zeta2 * x
                    + 224 * Hm10 * zeta2 * x2 - 56. / 3 * Hm100
                    - 224. / 3 * Hm100 / x - 2536. / 3 * Hm100 * x
                    - 3112. / 3 * Hm100 * x2 - 72 * Hm1000 - 144 * Hm1000 * x
                    - 176 * Hm1000 * x2 + 104. / 3 * Hm101
                    - 128. / 3 * Hm101 / x - 1328. / 3 * Hm101 * x
                    - 1744. / 3 * Hm101 * x2 - 32 * Hm1010 - 64 * Hm1010 * x
                    - 64 * Hm1010 * x2 - 32 * Hm1011 - 64 * Hm1011 * x
                    - 64 * Hm1011 * x2 - 112 * Hm1001 - 224 * Hm1001 * x
                    - 256 * Hm1001 * x2 - 15674. / 27 * H0 - 2272. / 27 * H0 / x
                    - 136514. / 27 * H0 * x + 24806. / 9 * H0 * x2
                    + 272 * H0 * zeta3 + 1024 * H0 * zeta3 * x
                    + 64. / 3 * H0 * zeta2 / x + 28 * H0 * zeta2
                    + 7204. / 3 * H0 * zeta2 * x - 944 * H0 * zeta2 * x2
                    + 2566. / 3 * H00 - 29168. / 9 * H00 * x
                    + 37376. / 9 * H00 * x2 + 36 * H00 * zeta2
                    + 920 * H00 * zeta2 * x - 96 * H00 * zeta2 * x2
                    - 332. / 3 * H000 - 6592. / 3 * H000 * x + 264 * H000 * x2
                    + 176 * H0000 - 816 * H0000 * x + 938. / 27 * H1
                    + 296 * H1 / x - 109468. / 27 * H1 * x
                    + 105382. / 27 * H1 * x2 - 280 * H1 * zeta3
                    + 560 * H1 * zeta3 * x - 448 * H1 * zeta3 * x2
                    + 256. / 3 * H1 * zeta2 / x - 44. / 3 * H1 * zeta2
                    + 2704. / 3 * H1 * zeta2 * x - 3064. / 3 * H1 * zeta2 * x2
                    - 48 * H10m10 + 96 * H10m10 * x - 96 * H10m10 * x2
                    + 80. / 9 * H10 - 464. / 9 * H10 / x - 24184. / 9 * H10 * x
                    + 8648. / 3 * H10 * x2 - 192 * H10 * zeta2
                    + 384 * H10 * zeta2 * x - 320 * H10 * zeta2 * x2
                    - 244. / 3 * H100 - 96 * H100 / x - 2524. / 3 * H100 * x
                    + 1096 * H100 * x2 + 104 * H1000 - 208 * H1000 * x
                    + 176 * H1000 * x2 + 632. / 9 * H11 - 568. / 9 * H11 / x
                    - 27544. / 9 * H11 * x + 29288. / 9 * H11 * x2
                    - 160 * H11 * zeta2 + 320 * H11 * zeta2 * x
                    - 256 * H11 * zeta2 * x2 + 296. / 3 * H110
                    - 224. / 3 * H110 / x - 3424. / 3 * H110 * x
                    + 1272 * H110 * x2 + 128 * H1100 - 256 * H1100 * x
                    + 224 * H1100 * x2 + 232. / 3 * H111 - 64 * H111 / x
                    - 3056. / 3 * H111 * x + 3392. / 3 * H111 * x2 + 112 * H1110
                    - 224 * H1110 * x + 224 * H1110 * x2 + 96 * H1111
                    - 192 * H1111 * x + 192 * H1111 * x2 + 112 * H1101
                    - 224 * H1101 * x + 224 * H1101 * x2 + 40 * H101
                    - 224. / 3 * H101 / x - 1024 * H101 * x
                    + 3464. / 3 * H101 * x2 + 96 * H1010 - 192 * H1010 * x
                    + 192 * H1010 * x2 + 112 * H1011 - 224 * H1011 * x
                    + 224 * H1011 * x2 + 160 * H1001 - 320 * H1001 * x
                    + 288 * H1001 * x2 + 360 * H01 + 832. / 9 * H01 / x
                    - 10600. / 3 * H01 * x + 36800. / 9 * H01 * x2
                    + 56 * H01 * zeta2 + 688 * H01 * zeta2 * x
                    - 128 * H01 * zeta2 * x2 - 16 * H010 - 128. / 3 * H010 / x
                    - 1952 * H010 * x + 1040 * H010 * x2 - 148 * H0100
                    - 744 * H0100 * x + 32 * H0100 * x2 - 16 * H011
                    - 160. / 3 * H011 / x - 2208 * H011 * x
                    + 3592. / 3 * H011 * x2 - 16 * H0110 - 640 * H0110 * x
                    + 192 * H0110 * x2 - 576 * H0111 * x + 192 * H0111 * x2
                    - 16 * H0101 - 640 * H0101 * x + 192 * H0101 * x2
                    - 28 * H001 - 8188. / 3 * H001 * x + 944 * H001 * x2
                    - 96 * H0010 - 896 * H0010 * x + 64 * H0010 * x2
                    - 80 * H0011 - 1024 * H0011 * x + 128 * H0011 * x2
                    - 36 * H0001 - 1064 * H0001 * x + 96 * H0001 * x2)
           + CA * CA * LQm * Lmmu
                 * (1018. / 9 - 280. / x + 14636. / 9 * x - 4426. / 3 * x2
                    - 96 * zeta3 + 32 * zeta3 * x - 48 * zeta2
                    - 3136. / 3 * zeta2 * x + 672 * zeta2 * x2 - 64 * H0m10
                    + 256 * H0m10 * x - 64 * Hm1 * zeta2 - 128 * Hm1 * zeta2 * x
                    - 128 * Hm1 * zeta2 * x2 - 64 * Hm1m10 - 128 * Hm1m10 * x
                    - 128 * Hm1m10 * x2 - 8. / 3 * Hm10 + 128. / 3 * Hm10 / x
                    + 944. / 3 * Hm10 * x + 1168. / 3 * Hm10 * x2 + 48 * Hm100
                    + 96 * Hm100 * x + 96 * Hm100 * x2 + 32 * Hm101
                    + 64 * Hm101 * x + 64 * Hm101 * x2 - 4036. / 9 * H0
                    - 416. / 9 * H0 / x + 14048. / 9 * H0 * x
                    - 20732. / 9 * H0 * x2 - 96 * H0 * zeta2
                    - 448 * H0 * zeta2 * x + 184. / 3 * H00
                    + 3488. / 3 * H00 * x - 144 * H00 * x2 - 96 * H000
                    + 512 * H000 * x - 236. / 9 * H1 + 224. / 3 * H1 / x
                    + 16816. / 9 * H1 * x - 18196. / 9 * H1 * x2
                    + 64 * H1 * zeta2 - 128 * H1 * zeta2 * x
                    + 128 * H1 * zeta2 * x2 - 40. / 3 * H10 + 64 * H10 / x
                    + 1904. / 3 * H10 * x - 2240. / 3 * H10 * x2 - 48 * H100
                    + 96 * H100 * x - 96 * H100 * x2 - 136. / 3 * H11
                    + 64 * H11 / x + 2480. / 3 * H11 * x - 2816. / 3 * H11 * x2
                    - 96 * H110 + 192 * H110 * x - 192 * H110 * x2 - 96 * H111
                    + 192 * H111 * x - 192 * H111 * x2 - 96 * H101
                    + 192 * H101 * x - 192 * H101 * x2 + 48 * H01
                    + 128. / 3 * H01 / x + 1360 * H01 * x - 672 * H01 * x2
                    + 64 * H010 + 448 * H010 * x - 64 * H010 * x2 + 32 * H011
                    + 512 * H011 * x - 128 * H011 * x2 + 96 * H001
                    + 704 * H001 * x)
           + CA * CA * LQm * Lmmu2
                 * (82 - 544. / 9. / x - 440 * x + 3766. / 9 * x2 + 32 * zeta2
                    + 128 * zeta2 * x - 44. / 3 * H0 - 32. / 3 * H0 / x
                    - 608. / 3 * H0 * x + 24 * H0 * x2 + 16 * H00
                    - 128 * H00 * x - 4. / 3 * H1 - 64. / 3 * H1 / x
                    - 472. / 3 * H1 * x + 584. / 3 * H1 * x2 + 16 * H10
                    - 32 * H10 * x + 32 * H10 * x2 + 32 * H11 - 64 * H11 * x
                    + 64 * H11 * x2 - 32 * H01 - 128 * H01 * x)
           + CA * CA * LQm2
                 * (509. / 9 - 140. / x + 7318. / 9 * x - 2213. / 3 * x2
                    - 48 * zeta3 + 16 * zeta3 * x - 24 * zeta2
                    - 1568. / 3 * zeta2 * x + 336 * zeta2 * x2 - 32 * H0m10
                    + 128 * H0m10 * x - 32 * Hm1 * zeta2 - 64 * Hm1 * zeta2 * x
                    - 64 * Hm1 * zeta2 * x2 - 32 * Hm1m10 - 64 * Hm1m10 * x
                    - 64 * Hm1m10 * x2 - 4. / 3 * Hm10 + 64. / 3 * Hm10 / x
                    + 472. / 3 * Hm10 * x + 584. / 3 * Hm10 * x2 + 24 * Hm100
                    + 48 * Hm100 * x + 48 * Hm100 * x2 + 16 * Hm101
                    + 32 * Hm101 * x + 32 * Hm101 * x2 - 2018. / 9 * H0
                    - 208. / 9 * H0 / x + 7024. / 9 * H0 * x
                    - 10366. / 9 * H0 * x2 - 48 * H0 * zeta2
                    - 224 * H0 * zeta2 * x + 92. / 3 * H00 + 1744. / 3 * H00 * x
                    - 72 * H00 * x2 - 48 * H000 + 256 * H000 * x - 118. / 9 * H1
                    + 112. / 3 * H1 / x + 8408. / 9 * H1 * x
                    - 9098. / 9 * H1 * x2 + 32 * H1 * zeta2
                    - 64 * H1 * zeta2 * x + 64 * H1 * zeta2 * x2 - 20. / 3 * H10
                    + 32 * H10 / x + 952. / 3 * H10 * x - 1120. / 3 * H10 * x2
                    - 24 * H100 + 48 * H100 * x - 48 * H100 * x2 - 68. / 3 * H11
                    + 32 * H11 / x + 1240. / 3 * H11 * x - 1408. / 3 * H11 * x2
                    - 48 * H110 + 96 * H110 * x - 96 * H110 * x2 - 48 * H111
                    + 96 * H111 * x - 96 * H111 * x2 - 48 * H101 + 96 * H101 * x
                    - 96 * H101 * x2 + 24 * H01 + 64. / 3 * H01 / x
                    + 680 * H01 * x - 336 * H01 * x2 + 32 * H010
                    + 224 * H010 * x - 32 * H010 * x2 + 16 * H011
                    + 256 * H011 * x - 64 * H011 * x2 + 48 * H001
                    + 352 * H001 * x)
           + CA * CA * LQm2 * Lmmu
                 * (82 - 544. / 9. / x - 440 * x + 3766. / 9 * x2 + 32 * zeta2
                    + 128 * zeta2 * x - 44. / 3 * H0 - 32. / 3 * H0 / x
                    - 608. / 3 * H0 * x + 24 * H0 * x2 + 16 * H00
                    - 128 * H00 * x - 4. / 3 * H1 - 64. / 3 * H1 / x
                    - 472. / 3 * H1 * x + 584. / 3 * H1 * x2 + 16 * H10
                    - 32 * H10 * x + 32 * H10 * x2 + 32 * H11 - 64 * H11 * x
                    + 64 * H11 * x2 - 32 * H01 - 128 * H01 * x)
           + CA * CA * LQm3
                 * (82. / 3 - 544. / 27. / x - 440. / 3 * x + 3766. / 27 * x2
                    + 32. / 3 * zeta2 + 128. / 3 * zeta2 * x - 44. / 9 * H0
                    - 32. / 9 * H0 / x - 608. / 9 * H0 * x + 8 * H0 * x2
                    + 16. / 3 * H00 - 128. / 3 * H00 * x - 4. / 9 * H1
                    - 64. / 9 * H1 / x - 472. / 9 * H1 * x + 584. / 9 * H1 * x2
                    + 16. / 3 * H10 - 32. / 3 * H10 * x + 32. / 3 * H10 * x2
                    + 32. / 3 * H11 - 64. / 3 * H11 * x + 64. / 3 * H11 * x2
                    - 32. / 3 * H01 - 128. / 3 * H01 * x)
           + CF * nf * (1. - 2. * x + 2. * x2) * (69. - 28. * zeta2)
           // from erratum
           + 8. / 9 * zeta3 - 16. / 9 * zeta3 * x + 16. / 9 * zeta3 * x2
           + a_muindep_->MuIndependentNfIndependentTerm(x)
           + massless_as3_->MuIndependentTerms(x, nf + 1) / (nf + 1.);
}

//==========================================================================================//
//  High scale (Q^2 >> m^2) limit of the quark coefficient functions for F2 at
//  O(as^3) expanded in terms of \alpha_s^{[nf+1]}
//
//  Eq. (B.10) of Ref. [arXiv:1205.5727].
//------------------------------------------------------------------------------------------//

Value HighScaleCoefficientFunction::D2_ps3_highscale(
    double x, double m2Q2, double m2mu2, int nf
) const {

    double x2 = x * x;
    double x3 = x2 * x;

    double LQm = log(1. / m2Q2);
    double LQm2 = LQm * LQm;
    double LQm3 = LQm2 * LQm;

    double Lmmu = log(m2mu2);
    double Lmmu2 = Lmmu * Lmmu;

    // Allocate pointers for the harmonic polylogs
    double wx = x;
    int nw = 5;
    int n1 = -1;
    int n2 = 1;
    int sz = n2 - n1 + 1;
    double *Hr1 = new double[sz];
    double *Hr2 = new double[sz * sz];
    double *Hr3 = new double[sz * sz * sz];
    double *Hr4 = new double[sz * sz * sz * sz];
    double *Hr5 = new double[sz * sz * sz * sz * sz];

    // Call polylogs
    apf_hplog_(&wx, &nw, Hr1, Hr2, Hr3, Hr4, Hr5, &n1, &n2);

    // weight 1
    const double Hm1 = Hr1[0];
    const double H0 = Hr1[1];
    const double H1 = Hr1[2];

    // weight 2
    const double Hm1m1 = Hr2[0];
    const double H0m1 = Hr2[1];
    const double Hm10 = Hr2[3];
    const double H00 = Hr2[4];
    const double H10 = Hr2[5];
    const double H01 = Hr2[7];
    const double H11 = Hr2[8];

    // weight 3
    const double H0m1m1 = Hr3[1];
    const double Hm1m10 = Hr3[9];
    const double H0m10 = Hr3[10];
    const double Hm100 = Hr3[12];
    const double H000 = Hr3[13];
    const double H100 = Hr3[14];
    const double H010 = Hr3[16];
    const double H110 = Hr3[17];
    const double Hm101 = Hr3[21];
    const double H001 = Hr3[22];
    const double H101 = Hr3[23];
    const double H011 = Hr3[25];
    const double H111 = Hr3[26];

    // weight 4
    const double Hm1m1m10 = Hr4[27];
    const double H0m1m10 = Hr4[28];
    const double H00m10 = Hr4[31];
    const double Hm1m100 = Hr4[36];
    const double H0m100 = Hr4[37];
    const double Hm1000 = Hr4[39];
    const double H0000 = Hr4[40];
    const double H1000 = Hr4[41];
    const double H0100 = Hr4[43];
    const double H1100 = Hr4[44];
    const double Hm1010 = Hr4[48];
    const double H0010 = Hr4[49];
    const double H1010 = Hr4[50];
    const double H0110 = Hr4[52];
    const double H1110 = Hr4[53];
    const double H0m101 = Hr4[64];
    const double H0001 = Hr4[67];
    const double H1001 = Hr4[68];
    const double H0101 = Hr4[70];
    const double H1101 = Hr4[71];
    const double H0011 = Hr4[76];
    const double H1011 = Hr4[77];
    const double H0111 = Hr4[79];
    const double H1111 = Hr4[80];

    //  weight 5
    const double H0m1m1m10 = Hr5[82];
    const double H0m1m100 = Hr5[109];
    const double H0m1000 = Hr5[118];
    const double H00000 = Hr5[121];
    const double H01000 = Hr5[124];
    const double H00100 = Hr5[130];
    const double H01100 = Hr5[133];
    const double H0m1010 = Hr5[145];
    const double H00010 = Hr5[148];
    const double H01010 = Hr5[151];
    const double H00110 = Hr5[157];
    const double H01110 = Hr5[160];
    const double H00001 = Hr5[202];
    const double H01001 = Hr5[205];
    const double H00101 = Hr5[211];
    const double H01101 = Hr5[214];
    const double H00011 = Hr5[229];
    const double H01011 = Hr5[232];
    const double H00111 = Hr5[238];
    const double H01111 = Hr5[241];

    delete[] Hr1;
    delete[] Hr2;
    delete[] Hr3;
    delete[] Hr4;
    delete[] Hr5;

    return CF
               * (-6400. / 81 + 11840. / 243 / x + 4912. / 81 * x
                  - 7376. / 243 * x2 + 128. / 27 * zeta3 / x - 296. / 9 * zeta3
                  - 248. / 3 * zeta3 * x - 1184. / 27 * zeta3 * x2
                  + 320. / 27 * zeta2 / x + 1232. / 27 * zeta2
                  + 992. / 27 * zeta2 * x - 224. / 27 * zeta2 * x2
                  + 8 * zeta2 * zeta2 + 8 * zeta2 * zeta2 * x + 3704. / 81 * H0
                  - 3784. / 81 * H0 * x + 6272. / 81 * H0 * x2
                  + 304. / 9 * H0 * zeta3 + 304. / 9 * H0 * zeta3 * x
                  + 8. / 3 * H0 * zeta2 + 40. / 3 * H0 * zeta2 * x
                  + 64. / 9 * H0 * zeta2 * x2 - 1024. / 27 * H00
                  - 1600. / 27 * H00 * x - 400. / 9 * H00 * x2
                  - 16. / 3 * H00 * zeta2 - 16. / 3 * H00 * zeta2 * x
                  + 8. / 3 * H000 + 40. / 3 * H000 * x + 64. / 9 * H000 * x2
                  - 16. / 3 * H0000 - 16. / 3 * H0000 * x + 248. / 27 * H1
                  - 2032. / 81 * H1 / x + 472. / 27 * H1 * x
                  - 128. / 81 * H1 * x2 + 56. / 3 * H10 - 160. / 27 * H10 / x
                  - 104. / 3 * H10 * x + 592. / 27 * H10 * x2 + 16. / 3 * H100
                  + 64. / 9 * H100 / x - 16. / 3 * H100 * x
                  - 64. / 9 * H100 * x2 + 56. / 9 * H11 + 176. / 9 * H11 / x
                  - 56. / 9 * H11 * x - 176. / 9 * H11 * x2 - 8. / 3 * H110
                  - 32. / 9 * H110 / x + 8. / 3 * H110 * x + 32. / 9 * H110 * x2
                  - 32. / 3 * H111 - 128. / 9 * H111 / x + 32. / 3 * H111 * x
                  + 128. / 9 * H111 * x2 - 8. / 3 * H101 - 32. / 9 * H101 / x
                  + 8. / 3 * H101 * x + 32. / 9 * H101 * x2 - 1160. / 27 * H01
                  - 632. / 27 * H01 * x - 176. / 9 * H01 * x2 + 80. / 9 * H010
                  - 64. / 9 * H010 * x - 32. / 3 * H010 * x2 + 32. / 3 * H0100
                  + 32. / 3 * H0100 * x + 40 * H011 + 136. / 3 * H011 * x
                  + 128. / 9 * H011 * x2 - 16. / 3 * H0110 - 16. / 3 * H0110 * x
                  - 64. / 3 * H0111 - 64. / 3 * H0111 * x - 16. / 3 * H0101
                  - 16. / 3 * H0101 * x + 128. / 9 * H001 + 176. / 9 * H001 * x
                  + 32. / 9 * H001 * x2 + 32. / 3 * H0010 + 32. / 3 * H0010 * x
                  - 16. / 3 * H0011 - 16. / 3 * H0011 * x)
           + CF * Lmmu
                 * (608. / 27 + 160. / 27 / x - 2432. / 27 * x + 1664. / 27 * x2
                    + 32. / 3 * zeta3 + 32. / 3 * zeta3 * x
                    - 64. / 9 * zeta2 / x - 64. / 3 * zeta2 * x
                    + 64. / 3 * zeta2 * x2 - 64. / 3 * Hm10 - 64. / 9 * Hm10 / x
                    - 64. / 3 * Hm10 * x - 64. / 9 * Hm10 * x2 + 560. / 9 * H0
                    - 176. / 3 * H0 * x - 1408. / 27 * H0 * x2
                    - 64. / 3 * H0 * zeta2 - 64. / 3 * H0 * zeta2 * x
                    + 160. / 3 * H00 * x - 64. / 3 * H00 * x2 + 64. / 3 * H000
                    + 64. / 3 * H000 * x + 416. / 9 * H1 - 416. / 27 * H1 / x
                    - 320. / 9 * H1 * x + 128. / 27 * H1 * x2 + 32. / 3 * H10
                    + 128. / 9 * H10 / x - 32. / 3 * H10 * x
                    - 128. / 9 * H10 * x2 + 16. / 3 * H11 + 64. / 9 * H11 / x
                    - 16. / 3 * H11 * x - 64. / 9 * H11 * x2
                    - 64. / 3 * H01 * x2 + 64. / 3 * H010 + 64. / 3 * H010 * x
                    + 32. / 3 * H011 + 32. / 3 * H011 * x + 64. / 3 * H001
                    + 64. / 3 * H001 * x)
           + CF * Lmmu2
                 * (-320. / 9 + 32. / 9 / x + 32. / 9 * x + 256. / 9 * x2
                    + 32. / 3 * zeta2 + 32. / 3 * zeta2 * x - 16. / 3 * H0
                    - 80. / 3 * H0 * x + 64. / 9 * H0 * x2 - 32. / 3 * H00
                    - 32. / 3 * H00 * x - 16. / 3 * H1 - 64. / 9 * H1 / x
                    + 16. / 3 * H1 * x + 64. / 9 * H1 * x2 - 32. / 3 * H01
                    - 32. / 3 * H01 * x)
           + CF * LQm
                 * (3664. / 27 + 1984. / 81 / x - 2704. / 27 * x
                    - 4864. / 81 * x2 - 64. / 9 * zeta2 / x
                    - 64. / 3 * zeta2 * x - 64. / 3 * Hm10 - 64. / 9 * Hm10 / x
                    - 64. / 3 * Hm10 * x - 64. / 9 * Hm10 * x2 + 3728. / 27 * H0
                    + 2480. / 27 * H0 * x - 320. / 9 * H0 * x2 + 464. / 9 * H00
                    + 944. / 9 * H00 * x - 64. / 3 * H00 * x2 + 32 * H000
                    + 32 * H000 * x)
           + CF * LQm * Lmmu
                 * (-416. / 9 + 416. / 27 / x + 320. / 9 * x - 128. / 27 * x2
                    + 32. / 3 * zeta2 + 32. / 3 * zeta2 * x + 64. / 3 * H0 * x2
                    - 64. / 3 * H00 - 64. / 3 * H00 * x - 16. / 3 * H1
                    - 64. / 9 * H1 / x + 16. / 3 * H1 * x + 64. / 9 * H1 * x2
                    - 32. / 3 * H01 - 32. / 3 * H01 * x)
           + CF * LQm * Lmmu2
                 * (16. / 3 + 64. / 9 / x - 16. / 3 * x - 64. / 9 * x2
                    + 32. / 3 * H0 + 32. / 3 * H0 * x)
           + CF * LQm2
                 * (-280. / 9 + 16. / 9 / x + 184. / 9 * x + 80. / 9 * x2
                    - 128. / 9 * H0 - 176. / 9 * H0 * x + 64. / 9 * H0 * x2
                    - 32. / 3 * H00 - 32. / 3 * H00 * x)
           + CF * LQm2 * Lmmu
                 * (8. / 3 + 32. / 9 / x - 8. / 3 * x - 32. / 9 * x2
                    + 16. / 3 * H0 + 16. / 3 * H0 * x)
           + CF * LQm3
                 * (8. / 9 + 32. / 27 / x - 8. / 9 * x - 32. / 27 * x2
                    + 16. / 9 * H0 + 16. / 9 * H0 * x)
           // from Moch, Vogt
           + CF * nf
                 * (14800. / 243 - 12032. / 729 / x - 19840. / 243 * x
                    + 27152. / 729 * x2 + 256. / 27 * zeta3 / x
                    - 976. / 27 * zeta3 - 928. / 27 * zeta3 * x
                    + 32. / 27 * zeta3 * x2 + 1600. / 81 * zeta2
                    + 1024. / 81 * zeta2 * x + 592. / 27 * zeta2 * x2
                    + 16 * zeta2 * zeta2 + 16 * zeta2 * zeta2 * x
                    + 6152. / 243 * H0 + 7736. / 243 * H0 * x
                    - 800. / 81 * H0 * x2 - 32. / 9 * H0 * zeta3
                    - 32. / 9 * H0 * zeta3 * x + 160. / 27 * H0 * zeta2
                    - 128. / 27 * H0 * zeta2 * x - 64. / 9 * H0 * zeta2 * x2
                    + 2192. / 81 * H00 + 464. / 81 * H00 * x
                    + 32. / 3 * H00 * x2 + 64. / 9 * H00 * zeta2
                    + 64. / 9 * H00 * zeta2 * x + 392. / 27 * H000
                    + 104. / 27 * H000 * x - 64. / 9 * H000 * x2
                    + 80. / 9 * H0000 + 80. / 9 * H0000 * x - 688. / 81 * H1
                    - 320. / 27 * H1 / x - 752. / 81 * H1 * x
                    + 800. / 27 * H1 * x2 - 64. / 9 * H1 * zeta2 / x
                    - 16. / 3 * H1 * zeta2 + 16. / 3 * H1 * zeta2 * x
                    + 64. / 9 * H1 * zeta2 * x2 - 80. / 3 * H10
                    + 112. / 3 * H10 * x - 32. / 3 * H10 * x2 - 16. / 3 * H100
                    - 64. / 9 * H100 / x + 16. / 3 * H100 * x
                    + 64. / 9 * H100 * x2 + 56. / 9 * H11 - 160. / 81 * H11 / x
                    - 104. / 9 * H11 * x + 592. / 81 * H11 * x2 + 16. / 3 * H110
                    + 64. / 9 * H110 / x - 16. / 3 * H110 * x
                    - 64. / 9 * H110 * x2 - 8. / 9 * H111 - 32. / 27 * H111 / x
                    + 8. / 9 * H111 * x + 32. / 27 * H111 * x2 + 16. / 3 * H101
                    + 64. / 9 * H101 / x - 16. / 3 * H101 * x
                    - 64. / 9 * H101 * x2 - 1600. / 81 * H01
                    - 1024. / 81 * H01 * x - 592. / 27 * H01 * x2
                    - 32. / 3 * H01 * zeta2 - 32. / 3 * H01 * zeta2 * x
                    - 208. / 9 * H010 - 112. / 9 * H010 * x
                    + 64. / 9 * H010 * x2 - 32. / 3 * H0100
                    - 32. / 3 * H0100 * x + 80. / 27 * H011
                    - 64. / 27 * H011 * x - 32. / 9 * H011 * x2
                    + 32. / 3 * H0110 + 32. / 3 * H0110 * x - 16. / 9 * H0111
                    - 16. / 9 * H0111 * x + 32. / 3 * H0101
                    + 32. / 3 * H0101 * x - 160. / 27 * H001
                    + 128. / 27 * H001 * x + 64. / 9 * H001 * x2
                    - 32. / 3 * H0010 - 32. / 3 * H0010 * x + 32. / 9 * H0011
                    + 32. / 9 * H0011 * x - 64. / 9 * H0001
                    - 64. / 9 * H0001 * x)
           // from Bluemline
           // + CF * nf * (
           //     - 944./27 - 14 * zeta2 * 6./9 - 2 * zeta4 * 90./27 +
           //     5248./243/x
           //     + 40 * zeta2 * 6./81/x + 2416. * x/27 + 26. * zeta2 * 6 * x /
           //     9
           //     - 2 * zeta4 * 90 * x / 27 - 18496 * x2/243 - 148 * zeta2 * 6 *
           //     x2/81
           //     + 296./27 * H0 - 20./27 * zeta2 * 6 * H0 - 152./27 * x * H0
           //     + 16./27 * zeta2 * 6 * x * H0 + 6400./81 * x2 * H0 + 8./9 *
           //     zeta2 * 6 * x2 * H0
           //     - 56./9 * H0 * H0 - 4./9 * zeta2 * 6 * H0*H0 - 88./9 * x *
           //     H0*H0
           //     - 4./9 * zeta2 * 6 * x * H0 * H0 - 448./27 * x2 * H0*H0 + 4./9
           //     * H0*H0*H0
           //     + 20./9 * x * H0*H0*H0 + 32./27 * x2 * H0*H0*H0 - 2./9 *
           //     H0*H0*H0*H0
           //     - 2./9 * x * H0*H0*H0*H0 + 2./9 * zeta2 * 6 * H1 + 8 * zeta2 *
           //     6 * H1/27/x
           //     - 2./9 * zeta2 * 6 * x * H1 - 8./27 * zeta2 * 6 * x2 * H1
           //     + 32./3 * H0 * H1
           //     - 320 * H0 * H1/27/x - 32 * x * H0 * H1 + 896./27 * x2 * H0 *
           //     H1 + 8./3 * H0*H0 * H1
           //     + 32 * H0*H0 * H1/9/ x - 8./3 * x * H0*H0 * H1 - 32./9 * x2 *
           //     H0*H0 * H1
           //     - 32./3 * H01 + 4./9 * zeta2 * 6 * H01 + 320 * H01/27/x
           //     + 32 * x * H01 + 4./9 * zeta2 * 6 * x * H01 - 896./27 * x2 *
           //     H01 - 32./3 * H0 * H01
           //     - 64 * H0 * H01 / 9 /x - 64./3 * x * H0 * H01 - 64./9 * x2 *
           //     H0 * H01
           //     + 16./3 * H0*H0 * H01 + 16./3 * x * H0*H0 * H01 + 16 * H001
           //     + 64 * H001/9 / x + 48 * x * H001 + 64./3 * x2 * H001 - 32./3
           //     * H0 * H001
           //     - 32./3 * x * H0 * H001 - 88 * zeta3/9 + 32 * zeta3/27/ x
           //     - 488./9 * x * zeta3 - 800./27 * x2 * zeta3 + 208./9 * H0 *
           //     zeta3 + 208./9 * x * H0 * zeta3
           // )
           // the commented term is from Bluemline and is the only term that
           // disagrees with the one from Moch, Vogt Which one should I trust?
           // Anyway, in the plots the difference is very small
           + CF * nf * Lmmu
                 * (880. / 9 - 208. / 3 * x - 256. / 9 * x2 + 64. / 3 * zeta3
                    + 64. / 3 * zeta3 * x - 64. / 9 * zeta2 / x
                    - 208. / 9 * zeta2 - 304. / 9 * zeta2 * x
                    + 64. / 9 * zeta2 * x2 - 64. / 3 * Hm10 - 64. / 9 * Hm10 / x
                    - 64. / 3 * Hm10 * x - 64. / 9 * Hm10 * x2 + 704. / 9 * H0
                    + 160. / 3 * H0 * x - 416. / 9 * H0 * x2
                    - 32. / 3 * H0 * zeta2 - 32. / 3 * H0 * zeta2 * x
                    + 256. / 9 * H00 + 832. / 9 * H00 * x - 128. / 9 * H00 * x2
                    + 64. / 3 * H000 + 64. / 3 * H000 * x + 80. / 3 * H1
                    - 112. / 3 * H1 * x + 32. / 3 * H1 * x2 + 16. / 3 * H10
                    + 64. / 9 * H10 / x - 16. / 3 * H10 * x - 64. / 9 * H10 * x2
                    - 16. / 3 * H11 - 64. / 9 * H11 / x + 16. / 3 * H11 * x
                    + 64. / 9 * H11 * x2 + 208. / 9 * H01 + 112. / 9 * H01 * x
                    - 64. / 9 * H01 * x2 + 32. / 3 * H010 + 32. / 3 * H010 * x
                    - 32. / 3 * H011 - 32. / 3 * H011 * x + 32. / 3 * H001
                    + 32. / 3 * H001 * x)
           + CF * nf * Lmmu2
                 * (-160. / 9 + 16. / 9 / x + 16. / 9 * x + 128. / 9 * x2
                    + 16. / 3 * zeta2 + 16. / 3 * zeta2 * x - 8. / 3 * H0
                    - 40. / 3 * H0 * x + 32. / 9 * H0 * x2 - 16. / 3 * H00
                    - 16. / 3 * H00 * x - 8. / 3 * H1 - 32. / 9 * H1 / x
                    + 8. / 3 * H1 * x + 32. / 9 * H1 * x2 - 16. / 3 * H01
                    - 16. / 3 * H01 * x)
           + CF * nf * LQm
                 * (3280. / 27 + 1088. / 81 / x - 2608. / 27 * x
                    - 3104. / 81 * x2 + 16. / 3 * zeta3 + 16. / 3 * zeta3 * x
                    - 64. / 9 * zeta2 / x - 128. / 9 * zeta2
                    - 368. / 9 * zeta2 * x - 32. / 9 * zeta2 * x2
                    - 64. / 3 * Hm10 - 64. / 9 * Hm10 / x - 64. / 3 * Hm10 * x
                    - 64. / 9 * Hm10 * x2 + 3040. / 27 * H0
                    + 1408. / 27 * H0 * x - 1264. / 27 * H0 * x2
                    + 464. / 9 * H00 + 944. / 9 * H00 * x - 64. / 3 * H00 * x2
                    + 32 * H000 + 32 * H000 * x + 8 * H1 + 160. / 27 * H1 / x
                    - 8. / 3 * H1 * x - 304. / 27 * H1 * x2 - 8. / 3 * H11
                    - 32. / 9 * H11 / x + 8. / 3 * H11 * x + 32. / 9 * H11 * x2
                    + 128. / 9 * H01 + 176. / 9 * H01 * x + 32. / 9 * H01 * x2
                    - 16. / 3 * H011 - 16. / 3 * H011 * x)
           + CF * nf * LQm * Lmmu
                 * (-560. / 9 + 32. / 9 / x + 368. / 9 * x + 160. / 9 * x2
                    - 256. / 9 * H0 - 352. / 9 * H0 * x + 128. / 9 * H0 * x2
                    - 64. / 3 * H00 - 64. / 3 * H00 * x)
           + CF * nf * LQm * Lmmu2
                 * (8. / 3 + 32. / 9 / x - 8. / 3 * x - 32. / 9 * x2
                    + 16. / 3 * H0 + 16. / 3 * H0 * x)
           + CF * nf * LQm2
                 * (-280. / 9 + 16. / 9 / x + 184. / 9 * x + 80. / 9 * x2
                    - 128. / 9 * H0 - 176. / 9 * H0 * x + 64. / 9 * H0 * x2
                    - 32. / 3 * H00 - 32. / 3 * H00 * x)
           + CF * nf * LQm2 * Lmmu
                 * (8. / 3 + 32. / 9 / x - 8. / 3 * x - 32. / 9 * x2
                    + 16. / 3 * H0 + 16. / 3 * H0 * x)
           + CF * nf * LQm3
                 * (8. / 9 + 32. / 27 / x - 8. / 9 * x - 32. / 27 * x2
                    + 16. / 9 * H0 + 16. / 9 * H0 * x)
           + CF * CF
                 * (-5672. / 27 - 466. / 9 / x + 494. / 3 * x + 2624. / 27 * x2
                    + 120 * zeta5 + 120 * zeta5 * x - 980. / 9 * zeta3
                    - 4276. / 9 * zeta3 * x - 80 * zeta3 * x2 - 40 * zeta2 / x
                    + 290. / 3 * zeta2 + 674. / 3 * zeta2 * x
                    + 3064. / 27 * zeta2 * x2 - 184. / 3 * zeta2 * zeta3
                    - 184. / 3 * zeta2 * zeta3 * x - 32 * zeta2 * zeta2
                    + 184. / 5 * zeta2 * zeta2 * x
                    + 592. / 15 * zeta2 * zeta2 * x2 - 4594. / 27 * H0
                    + 4030. / 27 * H0 * x + 1360. / 9 * H0 * x2
                    - 24 * H0 * zeta3 + 256. / 3 * H0 * zeta3 * x
                    - 32. / 9 * H0 * zeta3 * x2 - 139. / 3 * H0 * zeta2
                    - 467. / 3 * H0 * zeta2 * x + 368. / 9 * H0 * zeta2 * x2
                    - 40 * H0 * zeta2 * zeta2 - 40 * H0 * zeta2 * zeta2 * x
                    - 344. / 3 * H00 - 422. / 3 * H00 * x + 176. / 27 * H00 * x2
                    - 16. / 3 * H00 * zeta3 - 16. / 3 * H00 * zeta3 * x
                    + 24 * H00 * zeta2 + 8 * H00 * zeta2 * x
                    - 32 * H00 * zeta2 * x2 - 92 * H000 - 100 * H000 * x
                    - 848. / 9 * H000 * x2 + 20 * H000 * zeta2
                    + 20 * H000 * zeta2 * x + 8 * H0000 + 24 * H0000 * x
                    + 32. / 3 * H0000 * x2 - 5638. / 27 * H1 - 502. / 9 * H1 / x
                    + 1456. / 27 * H1 * x + 632. / 3 * H1 * x2
                    + 32. / 9 * H1 * zeta3 / x + 8. / 3 * H1 * zeta3
                    - 8. / 3 * H1 * zeta3 * x - 32. / 9 * H1 * zeta3 * x2
                    - 88. / 3 * H1 * zeta2 / x - 106. / 3 * H1 * zeta2
                    + 178. / 3 * H1 * zeta2 * x + 16. / 3 * H1 * zeta2 * x2
                    - 524. / 3 * H10 + 1480. / 27 * H10 / x - 20 * H10 * x
                    + 3776. / 27 * H10 * x2 - 32. / 3 * H10 * zeta2 / x
                    - 8 * H10 * zeta2 + 8 * H10 * zeta2 * x
                    + 32. / 3 * H10 * zeta2 * x2 + 344. / 3 * H100
                    + 80. / 3 * H100 / x - 248. / 3 * H100 * x
                    - 176. / 3 * H100 * x2 + 24 * H1000 + 32 * H1000 / x
                    - 24 * H1000 * x - 32 * H1000 * x2 - 128. / 3 * H11
                    - 284. / 27 * H11 / x - 100. / 3 * H11 * x
                    + 2336. / 27 * H11 * x2 + 16. / 3 * H11 * zeta2 / x
                    + 4 * H11 * zeta2 - 4 * H11 * zeta2 * x
                    - 16. / 3 * H11 * zeta2 * x2 + 164. / 3 * H110
                    + 32 * H110 / x - 164. / 3 * H110 * x - 32 * H110 * x2
                    + 16 * H1100 + 64. / 3 * H1100 / x - 16 * H1100 * x
                    - 64. / 3 * H1100 * x2 + 40. / 3 * H111 - 16 * H111 / x
                    + 56. / 3 * H111 * x - 16 * H111 * x2 + 16 * H1110
                    + 64. / 3 * H1110 / x - 16 * H1110 * x
                    - 64. / 3 * H1110 * x2 + 8 * H1111 + 32. / 3 * H1111 / x
                    - 8 * H1111 * x - 32. / 3 * H1111 * x2 + 60 * H101
                    + 32 * H101 / x - 60 * H101 * x - 32 * H101 * x2
                    + 16 * H1010 + 64. / 3 * H1010 / x - 16 * H1010 * x
                    - 64. / 3 * H1010 * x2 + 16 * H1011 + 64. / 3 * H1011 / x
                    - 16 * H1011 * x - 64. / 3 * H1011 * x2 + 16 * H1001
                    + 64. / 3 * H1001 / x - 16 * H1001 * x
                    - 64. / 3 * H1001 * x2 - 236. / 3 * H01 - 680. / 3 * H01 * x
                    - 2416. / 27 * H01 * x2 + 16. / 3 * H01 * zeta3
                    + 16. / 3 * H01 * zeta3 * x - 52 * H01 * zeta2
                    - 60 * H01 * zeta2 * x - 16. / 3 * H01 * zeta2 * x2
                    - 136. / 3 * H010 - 608. / 3 * H010 * x
                    - 1184. / 9 * H010 * x2 - 16 * H010 * zeta2
                    - 16 * H010 * zeta2 * x + 64 * H0100 + 96 * H0100 * x
                    - 64. / 3 * H0100 * x2 + 48 * H01000 + 48 * H01000 * x
                    + 8. / 3 * H011 - 116. / 3 * H011 * x
                    - 1040. / 9 * H011 * x2 + 8 * H011 * zeta2
                    + 8 * H011 * zeta2 * x + 48 * H0110 + 32 * H0110 * x
                    - 64. / 3 * H0110 * x2 + 32 * H01100 + 32 * H01100 * x
                    - 24 * H0111 - 8 * H0111 * x - 32. / 3 * H0111 * x2
                    + 32 * H01110 + 32 * H01110 * x + 16 * H01111
                    + 16 * H01111 * x + 64 * H0101 + 80 * H0101 * x
                    + 32 * H01010 + 32 * H01010 * x + 32 * H01011
                    + 32 * H01011 * x + 32 * H01001 + 32 * H01001 * x
                    + 62. / 3 * H001 + 310. / 3 * H001 * x
                    - 608. / 9 * H001 * x2 + 32 * H001 * zeta2
                    + 32 * H001 * zeta2 * x - 24 * H0010 + 88 * H0010 * x
                    + 64. / 3 * H0010 * x2 - 8 * H0011 + 88 * H0011 * x
                    + 64. / 3 * H0011 * x2 + 16 * H00110 + 16 * H00110 * x
                    + 16 * H00111 + 16 * H00111 * x - 16 * H00101
                    - 16 * H00101 * x - 16 * H0001 + 64. / 3 * H0001 * x2
                    - 16 * H00010 - 16 * H00010 * x - 16 * H00011
                    - 16 * H00011 * x - 8 * H00001 - 8 * H00001 * x)
           + CF * CF * Lmmu
                 * (-8356. / 45 + 2104. / 45 / x + 1724. / 5 * x
                    - 3088. / 15 * x2 + 64. / 3 * zeta3 / x + 80 * zeta3
                    + 512 * zeta3 * x - 96 * zeta3 * x2 + 1100. / 3 * zeta2
                    + 2204. / 9 * zeta2 * x - 2672. / 9 * zeta2 * x2
                    - 64. / 5 * zeta2 * x3 - 16 * zeta2 * zeta2
                    - 48 * zeta2 * zeta2 * x + 128 * H00m10 - 64 * H0m1 * zeta2
                    + 64 * H0m1 * zeta2 * x - 128 * H0m1m10 + 128 * H0m1m10 * x
                    + 160 * H0m10 + 160. / 3 * H0m10 * x - 64. / 3 * H0m10 * x2
                    + 64 * H0m100 - 64 * H0m100 * x - 32. / 3 * Hm1 * zeta2 / x
                    - 160 * Hm1 * zeta2 - 160 * Hm1 * zeta2 * x
                    - 32. / 3 * Hm1 * zeta2 * x2 - 192 * Hm1m10
                    + 64. / 3 * Hm1m10 / x - 192 * Hm1m10 * x
                    + 64. / 3 * Hm1m10 * x2 + 1136. / 3 * Hm10
                    + 64. / 45 * Hm10 / x2 + 3536. / 9 * Hm10 * x
                    - 64. / 5 * Hm10 * x3 + 128 * Hm100 + 128 * Hm100 * x
                    + 64 * Hm101 + 64. / 3 * Hm101 / x + 64 * Hm101 * x
                    + 64. / 3 * Hm101 * x2 - 4688. / 45 * H0 - 64. / 45 * H0 / x
                    + 23452. / 45 * H0 * x + 4976. / 45 * H0 * x2
                    + 144 * H0 * zeta3 + 16 * H0 * zeta3 * x + 80 * H0 * zeta2
                    + 1264. / 3 * H0 * zeta2 * x - 512. / 3 * H0 * zeta2 * x2
                    - 212 * H00 - 1532. / 9 * H00 * x + 2608. / 9 * H00 * x2
                    + 64. / 5 * H00 * x3 + 176 * H00 * zeta2
                    + 176 * H00 * zeta2 * x + 16 * H000 - 544. / 3 * H000 * x
                    + 320. / 3 * H000 * x2 - 80 * H0000 - 80 * H0000 * x
                    - 3196. / 9 * H1 - 232. / 9 * H1 / x + 2548. / 9 * H1 * x
                    + 880. / 9 * H1 * x2 + 96 * H1 * zeta2 / x - 16 * H1 * zeta2
                    + 16 * H1 * zeta2 * x - 96 * H1 * zeta2 * x2
                    - 920. / 3 * H10 - 352. / 9 * H10 / x + 680. / 3 * H10 * x
                    + 1072. / 9 * H10 * x2 - 64 * H100 / x + 64 * H100 * x2
                    - 424 * H11 - 128. / 9 * H11 / x + 312 * H11 * x
                    + 1136. / 9 * H11 * x2 - 64 * H110 - 256. / 3 * H110 / x
                    + 64 * H110 * x + 256. / 3 * H110 * x2 - 48 * H111
                    - 64 * H111 / x + 48 * H111 * x + 64 * H111 * x2 - 80 * H101
                    - 320. / 3 * H101 / x + 80 * H101 * x + 320. / 3 * H101 * x2
                    - 1100. / 3 * H01 + 148 * H01 * x + 2672. / 9 * H01 * x2
                    + 96 * H01 * zeta2 + 96 * H01 * zeta2 * x - 144 * H010
                    - 96 * H010 * x + 128 * H010 * x2 - 64 * H0100
                    - 64 * H0100 * x - 192 * H011 - 176 * H011 * x
                    + 448. / 3 * H011 * x2 - 128 * H0110 - 128 * H0110 * x
                    - 96 * H0111 - 96 * H0111 * x - 160 * H0101
                    - 160 * H0101 * x - 80 * H001 - 368 * H001 * x
                    + 512. / 3 * H001 * x2 - 160 * H0010 - 160 * H0010 * x
                    - 192 * H0011 - 192 * H0011 * x - 176 * H0001
                    - 176 * H0001 * x)
           + CF * CF * Lmmu2
                 * (44 - 44 * x - 24 * zeta3 - 24 * zeta3 * x - 8 * zeta2
                    - 32 * zeta2 * x + 64. / 3 * zeta2 * x2 + 82. / 3 * H0
                    - 106. / 3 * H0 * x - 128. / 3 * H0 * x2 - 24 * H0 * zeta2
                    - 24 * H0 * zeta2 * x + 16 * H00 * x - 32. / 3 * H00 * x2
                    + 8 * H000 + 8 * H000 * x + 206. / 3 * H1 - 16. / 3 * H1 / x
                    - 62. / 3 * H1 * x - 128. / 3 * H1 * x2 + 8 * H10
                    + 32. / 3 * H10 / x - 8 * H10 * x - 32. / 3 * H10 * x2
                    + 16 * H11 + 64. / 3 * H11 / x - 16 * H11 * x
                    - 64. / 3 * H11 * x2 + 8 * H01 + 32 * H01 * x
                    - 64. / 3 * H01 * x2 + 16 * H010 + 16 * H010 * x + 32 * H011
                    + 32 * H011 * x + 24 * H001 + 24 * H001 * x)
           + CF * CF * LQm
                 * (-10318. / 135 + 2704. / 45 / x + 41278. / 135 * x
                    - 13024. / 45 * x2 - 128. / 3 * zeta3 / x + 136 * zeta3
                    + 584 * zeta3 * x + 1928. / 3 * zeta2 + 872. / 9 * zeta2 * x
                    - 1760. / 9 * zeta2 * x2 - 64. / 5 * zeta2 * x3
                    - 368. / 5 * zeta2 * zeta2 - 528. / 5 * zeta2 * zeta2 * x
                    + 128 * H00m10 - 64 * H0m1 * zeta2 + 64 * H0m1 * zeta2 * x
                    - 128 * H0m1m10 + 128 * H0m1m10 * x + 160 * H0m10
                    + 160. / 3 * H0m10 * x - 64. / 3 * H0m10 * x2 + 64 * H0m100
                    - 64 * H0m100 * x - 32. / 3 * Hm1 * zeta2 / x
                    - 160 * Hm1 * zeta2 - 160 * Hm1 * zeta2 * x
                    - 32. / 3 * Hm1 * zeta2 * x2 - 192 * Hm1m10
                    + 64. / 3 * Hm1m10 / x - 192 * Hm1m10 * x
                    + 64. / 3 * Hm1m10 * x2 + 1136. / 3 * Hm10
                    + 64. / 45 * Hm10 / x2 + 3536. / 9 * Hm10 * x
                    - 64. / 5 * Hm10 * x3 + 128 * Hm100 + 128 * Hm100 * x
                    + 64 * Hm101 + 64. / 3 * Hm101 / x + 64 * Hm101 * x
                    + 64. / 3 * Hm101 * x2 - 808. / 45 * H0 - 64. / 45 * H0 / x
                    + 4868. / 5 * H0 * x - 25192. / 135 * H0 * x2
                    + 32 * H0 * zeta3 - 96 * H0 * zeta3 * x + 160 * H0 * zeta2
                    + 976. / 3 * H0 * zeta2 * x - 896. / 3 * H0 * zeta2 * x2
                    - 342 * H00 + 610. / 9 * H00 * x + 1472. / 9 * H00 * x2
                    + 64. / 5 * H00 * x3 + 288 * H00 * zeta2
                    + 288 * H00 * zeta2 * x + 80 * H000 - 352. / 3 * H000 * x
                    + 832. / 3 * H000 * x2 - 168 * H0000 - 168 * H0000 * x
                    - 660 * H1 - 1168. / 27 * H1 / x + 2548. / 3 * H1 * x
                    - 3944. / 27 * H1 * x2 + 96 * H1 * zeta2 / x
                    - 16 * H1 * zeta2 + 16 * H1 * zeta2 * x
                    - 96 * H1 * zeta2 * x2 - 1324. / 3 * H10
                    + 128. / 3 * H10 / x + 1324. / 3 * H10 * x
                    - 128. / 3 * H10 * x2 - 64 * H100 / x + 64 * H100 * x2
                    - 1300. / 3 * H11 + 80. / 3 * H11 / x + 1252. / 3 * H11 * x
                    - 32. / 3 * H11 * x2 - 80 * H110 - 320. / 3 * H110 / x
                    + 80 * H110 * x + 320. / 3 * H110 * x2 - 56 * H111
                    - 224. / 3 * H111 / x + 56 * H111 * x + 224. / 3 * H111 * x2
                    - 80 * H101 - 320. / 3 * H101 / x + 80 * H101 * x
                    + 320. / 3 * H101 * x2 - 1928. / 3 * H01 + 296 * H01 * x
                    + 1760. / 9 * H01 * x2 + 96 * H01 * zeta2
                    + 96 * H01 * zeta2 * x - 112 * H010 + 32 * H010 * x
                    + 704. / 3 * H010 * x2 - 64 * H0100 - 64 * H0100 * x
                    - 152 * H011 - 40 * H011 * x + 608. / 3 * H011 * x2
                    - 160 * H0110 - 160 * H0110 * x - 112 * H0111
                    - 112 * H0111 * x - 160 * H0101 - 160 * H0101 * x
                    - 160 * H001 - 272 * H001 * x + 896. / 3 * H001 * x2
                    - 240 * H0010 - 240 * H0010 * x - 224 * H0011
                    - 224 * H0011 * x - 288 * H0001 - 288 * H0001 * x)
           + CF * CF * LQm * Lmmu
                 * (1196. / 9 + 16 / x - 1580. / 9 * x + 80. / 3 * x2
                    - 32 * zeta3 - 32 * zeta3 * x - 112 * zeta2
                    - 112 * zeta2 * x + 320. / 3 * zeta2 * x2 + 140 * H0
                    - 236 * H0 * x - 928. / 9 * H0 * x2 - 128 * H0 * zeta2
                    - 128 * H0 * zeta2 * x - 16 * H00 + 32 * H00 * x
                    - 320. / 3 * H00 * x2 + 80 * H000 + 80 * H000 * x + 296 * H1
                    + 64. / 9 * H1 / x - 200 * H1 * x - 928. / 9 * H1 * x2
                    + 48 * H10 + 64 * H10 / x - 48 * H10 * x - 64 * H10 * x2
                    + 48 * H11 + 64 * H11 / x - 48 * H11 * x - 64 * H11 * x2
                    + 112 * H01 + 112 * H01 * x - 320. / 3 * H01 * x2
                    + 96 * H010 + 96 * H010 * x + 96 * H011 + 96 * H011 * x
                    + 128 * H001 + 128 * H001 * x)
           + CF * CF * LQm * Lmmu2
                 * (-46. / 3 + 46. / 3 * x + 16 * zeta2 + 16 * zeta2 * x
                    + 8 * H0 * x + 32. / 3 * H0 * x2 - 8 * H00 - 8 * H00 * x
                    - 8 * H1 - 32. / 3 * H1 / x + 8 * H1 * x + 32. / 3 * H1 * x2
                    - 16 * H01 - 16 * H01 * x)
           + CF * CF * LQm2
                 * (506. / 9 + 8 / x - 818. / 9 * x + 80. / 3 * x2 - 32 * zeta3
                    - 32 * zeta3 * x - 48 * zeta2 - 16 * zeta2 * x
                    + 224. / 3 * zeta2 * x2 + 88 * H0 - 140 * H0 * x
                    - 16. / 9 * H0 * x2 - 80 * H0 * zeta2 - 80 * H0 * zeta2 * x
                    - 16 * H00 - 8 * H00 * x - 224. / 3 * H00 * x2 + 48 * H000
                    + 48 * H000 * x + 164 * H1 - 128. / 9 * H1 / x
                    - 148 * H1 * x - 16. / 9 * H1 * x2 + 24 * H10 + 32 * H10 / x
                    - 24 * H10 * x - 32 * H10 * x2 + 24 * H11 + 32 * H11 / x
                    - 24 * H11 * x - 32 * H11 * x2 + 48 * H01 + 16 * H01 * x
                    - 224. / 3 * H01 * x2 + 48 * H010 + 48 * H010 * x
                    + 48 * H011 + 48 * H011 * x + 80 * H001 + 80 * H001 * x)
           + CF * CF * LQm2 * Lmmu
                 * (-92. / 3 + 92. / 3 * x + 32 * zeta2 + 32 * zeta2 * x
                    + 16 * H0 * x + 64. / 3 * H0 * x2 - 16 * H00 - 16 * H00 * x
                    - 16 * H1 - 64. / 3 * H1 / x + 16 * H1 * x
                    + 64. / 3 * H1 * x2 - 32 * H01 - 32 * H01 * x)
           + CF * CF * LQm3
                 * (-92. / 9 + 92. / 9 * x + 32. / 3 * zeta2
                    + 32. / 3 * zeta2 * x + 16. / 3 * H0 * x + 64. / 9 * H0 * x2
                    - 16. / 3 * H00 - 16. / 3 * H00 * x - 16. / 3 * H1
                    - 64. / 9 * H1 / x + 16. / 3 * H1 * x + 64. / 9 * H1 * x2
                    - 32. / 3 * H01 - 32. / 3 * H01 * x)
           + CA * CF
                 * (58838. / 81 - 14510. / 27 / x - 166340. / 81 * x
                    + 50344. / 27 * x2 - 120 * zeta5 + 280 * zeta5 * x
                    + 304. / 9 * zeta3 / x + 916. / 3 * zeta3
                    + 1544. / 3 * zeta3 * x + 1024. / 3 * zeta3 * x2
                    - 2060. / 27 * zeta2 / x + 1640. / 9 * zeta2
                    - 764. / 3 * zeta2 * x + 7364. / 27 * zeta2 * x2
                    + 148. / 3 * zeta2 * zeta3 + 172. / 3 * zeta2 * zeta3 * x
                    + 608. / 15 * zeta2 * zeta2 / x - 548. / 15 * zeta2 * zeta2
                    + 272. / 3 * zeta2 * zeta2 * x
                    + 16. / 3 * zeta2 * zeta2 * x2 + 48 * H0m1 * zeta3
                    - 48 * H0m1 * zeta3 * x - 32 * H0m1m1 * zeta2
                    + 32 * H0m1m1 * zeta2 * x - 64 * H0m1m1m10
                    + 64 * H0m1m1m10 * x + 32 * H0m1m100 - 32 * H0m1m100 * x
                    + 48 * H0m10 * x - 24 * H0m10 * zeta2
                    + 24 * H0m10 * zeta2 * x + 16 * H0m1000 - 16 * H0m1000 * x
                    + 32 * H0m1010 - 32 * H0m1010 * x - 32 * Hm1 * zeta3 / x
                    + 24 * Hm1 * zeta3 + 24 * Hm1 * zeta3 * x
                    - 32 * Hm1 * zeta3 * x2 - 160. / 9 * Hm1 * zeta2 / x
                    + 32. / 3 * Hm1 * zeta2 - 16. / 3 * Hm1 * zeta2 * x
                    - 304. / 9 * Hm1 * zeta2 * x2 + 64. / 3 * Hm1m1 * zeta2 / x
                    - 16 * Hm1m1 * zeta2 - 16 * Hm1m1 * zeta2 * x
                    + 64. / 3 * Hm1m1 * zeta2 * x2 - 32 * Hm1m1m10
                    + 128. / 3 * Hm1m1m10 / x - 32 * Hm1m1m10 * x
                    + 128. / 3 * Hm1m1m10 * x2 + 64. / 3 * Hm1m10
                    - 320. / 9 * Hm1m10 / x - 32. / 3 * Hm1m10 * x
                    - 608. / 9 * Hm1m10 * x2 + 16 * Hm1m100
                    - 64. / 3 * Hm1m100 / x + 16 * Hm1m100 * x
                    - 64. / 3 * Hm1m100 * x2 - 400. / 9 * Hm10
                    + 752. / 27 * Hm10 / x + 320. / 9 * Hm10 * x
                    + 2912. / 27 * Hm10 * x2 + 16 * Hm10 * zeta2 / x
                    - 12 * Hm10 * zeta2 - 12 * Hm10 * zeta2 * x
                    + 16 * Hm10 * zeta2 * x2 - 32. / 3 * Hm100
                    + 160. / 9 * Hm100 / x + 16. / 3 * Hm100 * x
                    + 304. / 9 * Hm100 * x2 + 8 * Hm1000 - 32. / 3 * Hm1000 / x
                    + 8 * Hm1000 * x - 32. / 3 * Hm1000 * x2 + 16 * Hm1010
                    - 64. / 3 * Hm1010 / x + 16 * Hm1010 * x
                    - 64. / 3 * Hm1010 * x2 - 3796. / 9 * H0
                    - 5248. / 81 * H0 / x - 6226. / 27 * H0 * x
                    - 36136. / 27 * H0 * x2 - 32. / 9 * H0 * zeta3 / x
                    - 1720. / 9 * H0 * zeta3 + 104. / 9 * H0 * zeta3 * x
                    - 160. / 9 * H0 * zeta2 / x - 590. / 9 * H0 * zeta2
                    - 242. / 9 * H0 * zeta2 * x - 316. / 3 * H0 * zeta2 * x2
                    - 8 * H0 * zeta2 * zeta2 + 304. / 5 * H0 * zeta2 * zeta2 * x
                    + 380. / 3 * H00 - 908. / 9 * H00 * x
                    + 8528. / 27 * H00 * x2 + 208. / 3 * H00 * zeta3
                    - 224. / 3 * H00 * zeta3 * x + 136. / 3 * H00 * zeta2
                    - 80. / 3 * H00 * zeta2 * x - 48 * H000
                    - 188. / 3 * H000 * x - 48 * H000 * x2 - 24 * H000 * zeta2
                    + 24 * H000 * zeta2 * x + 136. / 3 * H0000
                    - 32. / 3 * H0000 * x - 16 * H00000 + 16 * H00000 * x
                    - 656. / 27 * H1 - 3542. / 81 * H1 / x + 170. / 27 * H1 * x
                    + 5000. / 81 * H1 * x2 + 160. / 9 * H1 * zeta3 / x
                    + 40. / 3 * H1 * zeta3 - 40. / 3 * H1 * zeta3 * x
                    - 160. / 9 * H1 * zeta3 * x2 - 92. / 9 * H1 * zeta2 / x
                    - 26 * H1 * zeta2 + 2 * H1 * zeta2 * x
                    + 308. / 9 * H1 * zeta2 * x2 - 2128. / 9 * H10
                    + 5392. / 27 * H10 / x + 4432. / 9 * H10 * x
                    - 12304. / 27 * H10 * x2 + 16. / 3 * H10 * zeta2 / x
                    + 4 * H10 * zeta2 - 4 * H10 * zeta2 * x
                    - 16. / 3 * H10 * zeta2 * x2 + 352. / 3 * H100
                    - 904. / 9 * H100 / x - 376. / 3 * H100 * x
                    + 976. / 9 * H100 * x2 - 328. / 9 * H11
                    - 268. / 27 * H11 / x - 20. / 9 * H11 * x
                    + 1312. / 27 * H11 * x2 - 16. / 3 * H11 * zeta2 / x
                    - 4 * H11 * zeta2 + 4 * H11 * zeta2 * x
                    + 16. / 3 * H11 * zeta2 * x2 + 32. / 3 * H110
                    + 160. / 9 * H110 / x + 16. / 3 * H110 * x
                    - 304. / 9 * H110 * x2 - 8 * H1100 - 32. / 3 * H1100 / x
                    + 8 * H1100 * x + 32. / 3 * H1100 * x2 - 40. / 3 * H111
                    - 8. / 9 * H111 / x - 32. / 3 * H111 * x
                    + 224. / 9 * H111 * x2 + 16 * H1110 + 64. / 3 * H1110 / x
                    - 16 * H1110 * x - 64. / 3 * H1110 * x2 - 8 * H1111
                    - 32. / 3 * H1111 / x + 8 * H1111 * x + 32. / 3 * H1111 * x2
                    + 16 * H1101 + 64. / 3 * H1101 / x - 16 * H1101 * x
                    - 64. / 3 * H1101 * x2 + 32. / 3 * H101
                    + 160. / 9 * H101 / x + 16. / 3 * H101 * x
                    - 304. / 9 * H101 * x2 + 16 * H1010 + 64. / 3 * H1010 / x
                    - 16 * H1010 * x - 64. / 3 * H1010 * x2 + 8 * H1011
                    + 32. / 3 * H1011 / x - 8 * H1011 * x - 32. / 3 * H1011 * x2
                    - 304. / 9 * H01 - 92. / 9 * H01 * x - 448. / 27 * H01 * x2
                    + 80. / 3 * H01 * zeta3 + 80. / 3 * H01 * zeta3 * x
                    + 16. / 3 * H01 * zeta2 / x - 44. / 3 * H01 * zeta2
                    - 104. / 3 * H01 * zeta2 * x - 16. / 3 * H01 * zeta2 * x2
                    + 400. / 3 * H010 + 320. / 9 * H010 / x
                    + 640. / 3 * H010 * x + 1408. / 9 * H010 * x2
                    + 8 * H010 * zeta2 + 8 * H010 * zeta2 * x - 152. / 3 * H0100
                    - 64. / 3 * H0100 / x - 152. / 3 * H0100 * x
                    - 56. / 3 * H011 - 124. / 3 * H011 * x - 80. / 9 * H011 * x2
                    - 8 * H011 * zeta2 - 8 * H011 * zeta2 * x - 16 * H01100
                    - 16 * H01100 * x + 8 * H0111 - 8 * H0111 * x + 32 * H01110
                    + 32 * H01110 * x - 16 * H01111 - 16 * H01111 * x
                    + 32 * H01101 + 32 * H01101 * x + 32 * H01010
                    + 32 * H01010 * x + 16 * H01011 + 16 * H01011 * x
                    + 8 * H001 * zeta2 + 8 * H001 * zeta2 * x - 272. / 3 * H0010
                    + 16. / 3 * H0010 * x + 32 * H00100 - 64 * H00100 * x
                    - 8 * H0011 * x + 32 * H00010 - 32 * H00010 * x)
           + CA * CF * Lmmu
                 * (-6788. / 27 - 2848. / 27 / x + 42464. / 27 * x
                    - 32828. / 27 * x2 - 376. / 3 * zeta3 - 904. / 3 * zeta3 * x
                    - 64 * zeta3 * x2 + 128. / 3 * zeta2 / x - 1712. / 9 * zeta2
                    + 232. / 9 * zeta2 * x - 1304. / 3 * zeta2 * x2
                    - 208. / 5 * zeta2 * zeta2 - 672. / 5 * zeta2 * zeta2 * x
                    + 32 * H00m10 - 32 * H00m10 * x - 80 * H0m1 * zeta2
                    + 80 * H0m1 * zeta2 * x - 32 * H0m1m10 + 32 * H0m1m10 * x
                    - 80 * H0m10 + 64. / 3 * H0m10 / x - 96 * H0m10 * x
                    - 128. / 3 * H0m10 * x2 + 80 * H0m100 - 80 * H0m100 * x
                    + 64 * H0m101 - 64 * H0m101 * x + 64 * Hm1 * zeta2 / x
                    - 24 * Hm1 * zeta2 - 24 * Hm1 * zeta2 * x
                    + 64 * Hm1 * zeta2 * x2 + 16 * Hm1m10
                    + 128. / 3 * Hm1m10 / x + 16 * Hm1m10 * x
                    + 128. / 3 * Hm1m10 * x2 + 608. / 3 * Hm10
                    + 896. / 9 * Hm10 / x + 32. / 3 * Hm10 * x
                    - 832. / 9 * Hm10 * x2 + 24 * Hm100 - 64 * Hm100 / x
                    + 24 * Hm100 * x - 64 * Hm100 * x2 + 32 * Hm101
                    - 128. / 3 * Hm101 / x + 32 * Hm101 * x
                    - 128. / 3 * Hm101 * x2 - 4472. / 9 * H0 - 160. / 9 * H0 / x
                    + 352. / 9 * H0 * x + 26104. / 27 * H0 * x2
                    + 96 * H0 * zeta3 - 32 * H0 * zeta3 * x
                    + 64. / 3 * H0 * zeta2 / x + 320. / 3 * H0 * zeta2
                    + 320. / 3 * H0 * zeta2 * x - 128. / 3 * H0 * zeta2 * x2
                    + 3680. / 9 * H00 - 2776. / 9 * H00 * x
                    + 3608. / 9 * H00 * x2 - 32 * H00 * zeta2
                    + 160 * H00 * zeta2 * x - 448. / 3 * H000
                    - 1264. / 3 * H000 * x + 96 * H0000 - 128 * H0000 * x
                    - 2704. / 9 * H1 + 6608. / 27 * H1 / x + 2632. / 9 * H1 * x
                    - 6392. / 27 * H1 * x2 + 128. / 3 * H1 * zeta2 / x
                    + 56 * H1 * zeta2 - 56 * H1 * zeta2 * x
                    - 128. / 3 * H1 * zeta2 * x2 + 104 * H10
                    - 1616. / 9 * H10 / x - 240 * H10 * x + 2840. / 9 * H10 * x2
                    - 56 * H100 - 128. / 3 * H100 / x + 56 * H100 * x
                    + 128. / 3 * H100 * x2 + 184. / 3 * H11 - 400. / 9 * H11 / x
                    - 352. / 3 * H11 * x + 904. / 9 * H11 * x2 - 64 * H110
                    - 256. / 3 * H110 / x + 64 * H110 * x + 256. / 3 * H110 * x2
                    - 48 * H111 - 64 * H111 / x + 48 * H111 * x + 64 * H111 * x2
                    - 48 * H101 - 64 * H101 / x + 48 * H101 * x + 64 * H101 * x2
                    + 1712. / 9 * H01 + 512. / 9 * H01 / x - 136. / 9 * H01 * x
                    + 1304. / 3 * H01 * x2 + 80 * H01 * zeta2
                    + 80 * H01 * zeta2 * x - 248. / 3 * H010 - 64 * H010 / x
                    - 368. / 3 * H010 * x + 64 * H010 * x2 - 80 * H0100
                    - 80 * H0100 * x + 56. / 3 * H011 - 64 * H011 / x
                    + 128. / 3 * H011 * x + 256. / 3 * H011 * x2 - 128 * H0110
                    - 128 * H0110 * x - 96 * H0111 - 96 * H0111 * x - 96 * H0101
                    - 96 * H0101 * x - 320. / 3 * H001 - 608. / 3 * H001 * x
                    + 128. / 3 * H001 * x2 - 32 * H0010 - 224 * H0010 * x
                    - 96 * H0011 - 192 * H0011 * x + 32 * H0001
                    - 192 * H0001 * x)
           + CA * CF * Lmmu2
                 * (20 - 92. / 3 / x + 260 * x - 748. / 3 * x2 - 48 * zeta3 * x
                    - 32. / 3 * zeta2 / x - 112. / 3 * zeta2
                    - 184. / 3 * zeta2 * x + 32. / 3 * zeta2 * x2 - 84 * H0
                    - 16. / 3 * H0 / x + 140 * H0 * x - 176. / 3 * H0 * x2
                    - 48 * H0 * zeta2 * x + 88. / 3 * H00 + 352. / 3 * H00 * x
                    - 16 * H000 + 32 * H000 * x - 20. / 3 * H1
                    + 160. / 3 * H1 / x + 164. / 3 * H1 * x - 304. / 3 * H1 * x2
                    + 8 * H10 + 32. / 3 * H10 / x - 8 * H10 * x
                    - 32. / 3 * H10 * x2 + 16 * H11 + 64. / 3 * H11 / x
                    - 16 * H11 * x - 64. / 3 * H11 * x2 + 112. / 3 * H01
                    + 32. / 3 * H01 / x + 184. / 3 * H01 * x
                    - 32. / 3 * H01 * x2 + 16 * H010 + 16 * H010 * x + 32 * H011
                    + 32 * H011 * x + 48 * H001 * x)
           + CA * CF * LQm
                 * (7184. / 27 - 18416. / 27 / x - 2108. / 27 * x
                    + 13340. / 27 * x2 + 128 * zeta3 / x - 112. / 3 * zeta3
                    - 280. / 3 * zeta3 * x - 608. / 3 * zeta3 * x2
                    + 352. / 3 * zeta2 / x - 3604. / 9 * zeta2
                    + 1484. / 9 * zeta2 * x - 376 * zeta2 * x2
                    + 28 * zeta2 * zeta2 + 164. / 5 * zeta2 * zeta2 * x
                    + 112 * H00m10 - 48 * H00m10 * x - 48 * H0m1 * zeta2
                    + 48 * H0m1 * zeta2 * x + 32 * H0m1m10 - 32 * H0m1m10 * x
                    - 104 * H0m10 + 64. / 3 * H0m10 / x + 8 * H0m10 * x
                    - 320. / 3 * H0m10 * x2 + 112 * H0m100 - 112 * H0m100 * x
                    + 64 * H0m101 - 64 * H0m101 * x + 64 * Hm1 * zeta2 / x
                    + 24 * Hm1 * zeta2 + 24 * Hm1 * zeta2 * x
                    + 64 * Hm1 * zeta2 * x2 + 48 * Hm1m10 + 48 * Hm1m10 * x
                    + 832. / 3 * Hm10 + 1888. / 9 * Hm10 / x
                    + 592. / 3 * Hm10 * x + 1168. / 9 * Hm10 * x2 + 24 * Hm100
                    - 96 * Hm100 / x + 24 * Hm100 * x - 96 * Hm100 * x2
                    - 64 * Hm101 / x - 64 * Hm101 * x2 - 32056. / 27 * H0
                    - 2272. / 27 * H0 / x - 15688. / 27 * H0 * x
                    - 3728. / 27 * H0 * x2 + 240 * H0 * zeta3
                    + 224 * H0 * zeta3 * x + 64. / 3 * H0 * zeta2 / x
                    + 28 * H0 * zeta2 + 124 * H0 * zeta2 * x
                    - 160. / 3 * H0 * zeta2 * x2 + 5188. / 9 * H00
                    - 4160. / 9 * H00 * x + 6032. / 9 * H00 * x2
                    - 24 * H00 * zeta2 + 152 * H00 * zeta2 * x - 272 * H000
                    - 328 * H000 * x + 176 * H0000 - 192 * H0000 * x
                    - 404. / 3 * H1 + 8840. / 27 * H1 / x - 148. / 3 * H1 * x
                    - 3872. / 27 * H1 * x2 + 160. / 3 * H1 * zeta2 / x
                    + 64 * H1 * zeta2 - 64 * H1 * zeta2 * x
                    - 160. / 3 * H1 * zeta2 * x2 + 64. / 3 * H10
                    - 496. / 9 * H10 / x - 304. / 3 * H10 * x
                    + 1216. / 9 * H10 * x2 - 92 * H100 - 224. / 3 * H100 / x
                    + 92 * H100 * x + 224. / 3 * H100 * x2 + 36 * H11
                    - 424. / 9 * H11 / x - 116 * H11 * x + 1144. / 9 * H11 * x2
                    - 40 * H110 - 160. / 3 * H110 / x + 40 * H110 * x
                    + 160. / 3 * H110 * x2 - 40 * H111 - 160. / 3 * H111 / x
                    + 40 * H111 * x + 160. / 3 * H111 * x2 - 40 * H101
                    - 160. / 3 * H101 / x + 40 * H101 * x + 160. / 3 * H101 * x2
                    + 3604. / 9 * H01 + 832. / 9 * H01 / x + 292. / 9 * H01 * x
                    + 376 * H01 * x2 + 96 * H01 * zeta2 + 96 * H01 * zeta2 * x
                    - 16 * H010 - 128. / 3 * H010 / x - 32 * H010 * x
                    + 224. / 3 * H010 * x2 - 136 * H0100 - 136 * H0100 * x
                    - 8. / 3 * H011 - 160. / 3 * H011 / x - 32. / 3 * H011 * x
                    + 224. / 3 * H011 * x2 - 80 * H0110 - 80 * H0110 * x
                    - 80 * H0111 - 80 * H0111 * x - 80 * H0101 - 80 * H0101 * x
                    - 28 * H001 - 116 * H001 * x + 160. / 3 * H001 * x2
                    - 64 * H0010 - 160 * H0010 * x - 80 * H0011
                    - 176 * H0011 * x + 24 * H0001 - 200 * H0001 * x)
           + CA * CF * LQm * Lmmu
                 * (1168. / 3 - 7208. / 27 / x - 176 * x + 1448. / 27 * x2
                    - 64 * zeta3 + 32 * zeta3 * x - 16 * zeta2 - 32 * zeta2 * x
                    + 128. / 3 * zeta2 * x2 - 64 * H0m10 + 64 * H0m10 * x
                    - 32 * Hm10 + 128. / 3 * Hm10 / x - 32 * Hm10 * x
                    + 128. / 3 * Hm10 * x2 - 2624. / 9 * H0 - 416. / 9 * H0 / x
                    + 1696. / 9 * H0 * x - 3232. / 9 * H0 * x2 - 32 * H0 * zeta2
                    - 64 * H0 * zeta2 * x + 448. / 3 * H00 + 448. / 3 * H00 * x
                    - 96 * H000 + 128 * H000 * x - 272. / 3 * H1
                    + 704. / 9 * H1 / x + 368. / 3 * H1 * x - 992. / 9 * H1 * x2
                    + 32 * H10 + 128. / 3 * H10 / x - 32 * H10 * x
                    - 128. / 3 * H10 * x2 + 32 * H11 + 128. / 3 * H11 / x
                    - 32 * H11 * x - 128. / 3 * H11 * x2 + 16 * H01
                    + 128. / 3 * H01 / x - 128. / 3 * H01 * x2 + 64 * H010
                    + 64 * H010 * x + 64 * H011 + 64 * H011 * x + 32 * H001
                    + 128 * H001 * x)
           + CA * CF * LQm * Lmmu2
                 * (60 - 176. / 3 / x - 60 * x + 176. / 3 * x2 + 16 * zeta2
                    + 16 * zeta2 * x - 88. / 3 * H0 - 32. / 3 * H0 / x
                    - 64. / 3 * H0 * x + 16 * H00 - 32 * H00 * x - 8 * H1
                    - 32. / 3 * H1 / x + 8 * H1 * x + 32. / 3 * H1 * x2
                    - 16 * H01 - 16 * H01 * x)
           + CA * CF * LQm2
                 * (584. / 3 - 3604. / 27 / x - 88 * x + 724. / 27 * x2
                    - 32 * zeta3 + 16 * zeta3 * x - 8 * zeta2 - 16 * zeta2 * x
                    + 64. / 3 * zeta2 * x2 - 32 * H0m10 + 32 * H0m10 * x
                    - 16 * Hm10 + 64. / 3 * Hm10 / x - 16 * Hm10 * x
                    + 64. / 3 * Hm10 * x2 - 1312. / 9 * H0 - 208. / 9 * H0 / x
                    + 848. / 9 * H0 * x - 1616. / 9 * H0 * x2 - 16 * H0 * zeta2
                    - 32 * H0 * zeta2 * x + 224. / 3 * H00 + 224. / 3 * H00 * x
                    - 48 * H000 + 64 * H000 * x - 136. / 3 * H1
                    + 352. / 9 * H1 / x + 184. / 3 * H1 * x - 496. / 9 * H1 * x2
                    + 16 * H10 + 64. / 3 * H10 / x - 16 * H10 * x
                    - 64. / 3 * H10 * x2 + 16 * H11 + 64. / 3 * H11 / x
                    - 16 * H11 * x - 64. / 3 * H11 * x2 + 8 * H01
                    + 64. / 3 * H01 / x - 64. / 3 * H01 * x2 + 32 * H010
                    + 32 * H010 * x + 32 * H011 + 32 * H011 * x + 16 * H001
                    + 64 * H001 * x)
           + CA * CF * LQm2 * Lmmu
                 * (60 - 176. / 3 / x - 60 * x + 176. / 3 * x2 + 16 * zeta2
                    + 16 * zeta2 * x - 88. / 3 * H0 - 32. / 3 * H0 / x
                    - 64. / 3 * H0 * x + 16 * H00 - 32 * H00 * x - 8 * H1
                    - 32. / 3 * H1 / x + 8 * H1 * x + 32. / 3 * H1 * x2
                    - 16 * H01 - 16 * H01 * x)
           + CA * CF * LQm3
                 * (20 - 176. / 9 / x - 20 * x + 176. / 9 * x2 + 16. / 3 * zeta2
                    + 16. / 3 * zeta2 * x - 88. / 9 * H0 - 32. / 9 * H0 / x
                    - 64. / 9 * H0 * x + 16. / 3 * H00 - 32. / 3 * H00 * x
                    - 8. / 3 * H1 - 32. / 9 * H1 / x + 8. / 3 * H1 * x
                    + 32. / 9 * H1 * x2 - 16. / 3 * H01 - 16. / 3 * H01 * x)
           + a_muindep_->MuIndependentNfIndependentTerm(x)
           + 1. / (1 + nf) * massless_as3_->MuIndependentTerms(x, 1 + nf);
}
