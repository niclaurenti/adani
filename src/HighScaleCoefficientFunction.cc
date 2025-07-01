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
    double z, double m2Q2, double m2mu2, int nf
) const {

    double z2 = z * z;
    double z3 = z2 * z;
    double z4 = z3 * z;
    double z5 = z4 * z;

    double gammaqg = -4 * (z2 + (1 -2 * z + z2));

    double LM = log(m2mu2);
    double LM2 = LM * LM;
    double LM3 = LM2 * LM;
    double LQ = log(1. / m2Q2) + LM;
    double LQ2 = LQ * LQ;
    double LQ3 = LQ2 * LQ;

    // Allocate pointers for the harmonic polylogs
    double wx = z;
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
//     const double Hm1m1 = Hr2[0];
    const double H0m1 = Hr2[1];
//     const double Hm10 = Hr2[3];
//     const double H00 = Hr2[4];
//     const double H10 = Hr2[5];
    const double H01 = Hr2[7];
//     const double H11 = Hr2[8];

    // weight 3
    const double Hm1m1m1 = Hr3[0];
    const double H0m1m1 = Hr3[1];
    const double Hm10m1 = Hr3[3];
    const double H00m1 = Hr3[4];
    const double H01m1 = Hr3[7];
    const double Hm1m10 = Hr3[9];
    const double H0m10 = Hr3[10];
    const double Hm100 = Hr3[12];
    const double H000 = Hr3[13];
    const double H100 = Hr3[14];
    const double H010 = Hr3[16];
    const double H110 = Hr3[17];
    const double H0m11 = Hr3[19];
    const double Hm101 = Hr3[21];
    const double H001 = Hr3[22];
    const double H101 = Hr3[23];
    const double H011 = Hr3[25];
    const double H111 = Hr3[26];

    // weight 4
    const double H0m1m1m1 = Hr4[1];
    const double H00m1m1 = Hr4[4];
    const double H01m1m1 = Hr4[7];
    const double H0m11m1 = Hr4[19];
    const double H000m1 = Hr4[13];
    const double H001m1 = Hr4[22];
    const double H011m1 = Hr4[25];
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
    const double H0m1m11 = Hr4[55];
    const double H00m11 = Hr4[58];
    const double H01m11 = Hr4[61];
    const double Hm1m101 = Hr4[63];
    const double H0m101 = Hr4[64];
    const double Hm1001 = Hr4[66];
    const double H0001 = Hr4[67];
    const double H1001 = Hr4[68];
    const double H0101 = Hr4[70];
    const double H1101 = Hr4[71];
    const double H0m111 = Hr4[73];
    const double Hm1011 = Hr4[75];
    const double H0011 = Hr4[76];
    const double H1011 = Hr4[77];
    const double H0111 = Hr4[79];
    const double H1111 = Hr4[80];

    //  weight 5
    const double H0m1m1m1m1 = Hr5[1];
    const double H00m1m1m1 = Hr5[4];
    const double H01m1m1m1 = Hr5[7];
    const double H0m10m1m1 = Hr5[10];
    const double H000m1m1 = Hr5[13];
    const double H0m11m1m1 = Hr5[19];
    const double H001m1m1 = Hr5[22];
    const double H011m1m1 = Hr5[25];
    const double H00m10m1 = Hr5[31];
    const double H0000m1 = Hr5[40];
    const double H0010m1 = Hr5[49];
    const double H0m1m11m1 = Hr5[55];
    const double H00m11m1 = Hr5[58];
    const double H01m11m1 = Hr5[61];
    const double H0m101m1 = Hr5[64];
    const double H0001m1 = Hr5[67];
    const double H0m111m1 = Hr5[73];
    const double H0011m1 = Hr5[76];
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
    const double H0m1m1m11 = Hr5[163];
    const double H00m1m11 = Hr5[166];
    const double H01m1m11 = Hr5[169];
    const double H0m10m11 = Hr5[172];
    const double H000m11 = Hr5[175];
    const double H11110 = Hr5[161];
    const double H0m11m11 = Hr5[181];
    const double H001m11 = Hr5[184];
    const double Hm1m1m101 = Hr5[189];
    const double H0m1m101 = Hr5[190];
    const double Hm10m101 = Hr5[192];
    const double H00m101 = Hr5[193];
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
    const double H0m1m111 = Hr5[217];
    const double H00m111 = Hr5[220];
    const double Hm1m1011 = Hr5[225];
    const double H0m1011 = Hr5[226];
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

    double H02 = H0 * H0;
    double H03 = H02 * H0;
    double H04 = H03 * H0;
    double H05 = H04 * H0;

    double H12 = H1 * H1;
    double H13 = H12 * H1;
    double H14 = H13 * H1;
    double H15 = H14 * H1;

    double Hm12 = Hm1 * Hm1;
    double Hm13 = Hm12 * Hm1;
    double Hm14 = Hm13 * Hm1;
    double Hm15 = Hm14 * Hm1;

    double H012 = H01 * H01;
    double H0m12 = H0m1 * H0m1;

    double zeta22 = zeta2 * zeta2;

    return TR*TR*TR*((16*gammaqg*LM3)/9. - (16*gammaqg*LM2*LQ)/9. + LM2*(gammaqg*((16*H0)/9. + (16*H1)/9.) - (64*(1 - 8*z + 8*z2))/9.) - (16*gammaqg*zeta3)/9.) +
         CA*nf*TR*TR*((32*H0*H12*(-1 + z)*z)/3. - (4*H04*(1 + 2*z))/9. + (32*H01*H02*(1 + 4*z))/3. - (16*H011*z*(-4 + 5*z))/3. - (8*H13*(-1 + z)*(1 + 5*z))/9. +
            (32*H1*(-1 - 2*z + 2*z2))/3. - (8*H12*(1 - 6*z + 4*z2))/3. + (8*H03*(3 + 12*z + 23*z2))/27. - (8*H02*H1*(-1 + z)*(8 + 17*z + 65*z2))/(9.*z) +
            (32*H0*H1*(-1 + z)*(20 - 7*z + 254*z2))/(27.*z) - (8*H02*(42 + 39*z + 292*z2))/27. + (16*H0*(111 + 645*z + 3392*z2))/81. +
            LQ3*((8*gammaqg*H1)/9. + (32*H0*(1 + 4*z))/9. - (16*(-1 + z)*(4 + 7*z + 31*z2))/(27.*z)) +
            LM3*((-8*gammaqg*H1)/9. - (32*H0*(1 + 4*z))/9. + (16*(-1 + z)*(4 + 7*z + 31*z2))/(27.*z)) - (16*H0*H01*(8 + 15*z + 96*z2 + 23*z3))/(9.*z) +
            (16*H001*(8 + 21*z + 144*z2 + 111*z3))/(9.*z) - (32*H01*(-20 + 27*z - 261*z2 + 290*z3))/(27.*z) - (32*(-328 + 531*z - 6219*z2 + 5854*z3))/(243.*z) -
            (16*H1*(-5 - 2*z + 2*z2)*zeta2)/9. + (8*H0*(-7 - 16*z + 62*z2)*zeta2)/9. + (4*(-4 + 2*z - 50*z2 + 19*z3)*zeta2)/(9.*z) +
            (1 + 6*z)*((-64*H0*H001)/3. - (8*H02*zeta2)/3.) + z*(1 + z)*
             ((-64*H00m1)/3. + 64*H0m1 + (64*H0*H0m1)/3. + (128*H0m1m1)/3. + (-64*H0 - (32*H02)/3. - (128*H0m1)/3.)*Hm1 + (64*H0*Hm12)/3. + (64*Hm1*zeta2)/3.) +
            LQ2*(gammaqg*((-8*H0*H1)/3. - (4*H12)/3.) - (32*H01*(1 + 4*z))/3. - (16*H02*(1 + 8*z))/3. + ((32*H0m1)/3. - (32*H0*Hm1)/3.)*(1 + 2*z + 2*z2) +
               (16*H0*(-13 - 136*z + 49*z2))/9. + (16*H1*(-4 + 13*z - 80*z2 + 87*z3))/(9.*z) + (8*(-40 - 87*z - 894*z2 + 1048*z3))/(27.*z) - (64*(-2 + z)*z*zeta2)/3.) +
            LM2*((-4*gammaqg*H12)/3. + ((-32*H0m1)/3. + (32*H0*Hm1)/3.)*(1 + 2*z + 2*z2) - (32*H1*(5 - 4*z + 4*z2))/9. - (32*H0*(-5 - 20*z + 9*z2))/9. -
               (8*(-52 + 42*z - 168*z2 + 205*z3))/(27.*z) + z*((64*H02)/3. + (64*zeta2)/3.)) + z*(128*H0001 - (1216*zeta22)/15.) + (32*H0*(13 + 28*z)*zeta3)/9. -
            (16*(-4 + 33*z + 228*z2 + 550*z3)*zeta3)/(27.*z) + gammaqg*
             ((-40*H0011)/3. - (8*H0111)/3. - (8*H012)/3. + ((40*H001)/3. + (8*H011)/3.)*H1 + (2*H02*H12)/3. + H0*((40*H011)/3. - 8*H01*H1 - (8*H13)/9.) + H14/9. +
               (4*H12*zeta2)/3. - (40*H1*zeta3)/9.) + (1 + 2*z + 2*z2)*
             ((-32*H000m1)/3. + (128*H001m1)/3. + (128*H00m11)/3. - (64*H00m1m1)/3. - (16*H02*H0m1)/3. + (64*H0m101)/3. +
               H0*((32*H00m1)/3. - (64*H01m1)/3. - (64*H0m11)/3. + (64*H0m1m1)/3.) + (128*H0m1m1m1)/3. +
               ((-128*H001)/3. + (64*H00m1)/3. + (16*H03)/9. + H0*((64*H01)/3. - (64*H0m1)/3.) - (128*H0m1m1)/3.)*Hm1 + ((16*H02)/3. + (64*H0m1)/3.)*Hm12 -
               (64*H0*Hm13)/9. + ((32*H0m1)/3. - (32*H0*Hm1)/3. - (32*Hm12)/3.)*zeta2 + 32*Hm1*zeta3) +
            LM*((-16*H03*(-1 + 10*z))/9. + (32*H001*(-3 - 18*z + 2*z2))/3. - (32*H0*H01*(-1 - 10*z + 2*z2))/3. + ((64*H0m1)/9. - (64*H0*Hm1)/9.)*(5 + 4*z + 4*z2) +
               (8*H02*(-13 - 66*z + 7*z2))/9. - (8*H12*(23 - 16*z + 13*z2))/9. + (16*H0*(-2 - 269*z + 58*z2))/27. + (16*H1*(29 - 64*z + 64*z2))/27. -
               (16*H0*H1*(-1 + z)*(8 + 17*z + 65*z2))/(9.*z) + (16*H01*(-8 - 9*z - 54*z2 + 68*z3))/(9.*z) + (8*(-4 - 119*z - 268*z2 + 592*z3))/(27.*z) -
               (16*z*(10 + 3*z)*zeta2)/9. + gammaqg*((-32*H011)/3. + (16*H01*H1)/3. - (4*H02*H1)/3. - (8*H13)/9. - (8*H1*zeta2)/3.) +
               (1 + 2*z + 2*z2)*((-32*H00m1)/3. + (32*H0*H0m1)/3. - (64*H0m1m1)/3. + ((-16*H02)/3. + (64*H0m1)/3.)*Hm1 - (32*H0*Hm12)/3. - (32*Hm1*zeta2)/3.) +
               (128*(1 + 5*z)*zeta3)/3.) + LQ*(gammaqg*((-32*H01*H1)/3. + (16*H0*H12)/3. + (8*H13)/9.) +
               ((64*H0m1m1)/3. - (64*H0m1*Hm1)/3. + (32*H0*Hm12)/3.)*(1 + 2*z) + (16*H03*(1 + 26*z))/9. - (64*H001*(2 - 2*z + z2))/3. + (64*H0*H01*(2 + 2*z + z2))/3. +
               (64*H00m1*(2 + 4*z + z2))/3. + ((-64*H01m1)/3. - (64*H0m11)/3. + (64*H01*Hm1)/3.)*(1 + 2*z + 2*z2) - (128*H011*(1 - 5*z + 3*z2))/3. -
               (32*H02*H1*(2 - 4*z + 3*z2))/3. - (64*H0*H0m1*(2 + 4*z + 3*z2))/3. + (32*H02*Hm1*(2 + 4*z + 5*z2))/3. - (8*H02*(-39 - 542*z + 126*z2))/9. -
               (16*H0*(-89 - 2822*z + 1520*z2))/27. + (((-64*H0m1)/9. + (64*H0*Hm1)/9.)*(-2 - 4*z + 7*z2 + 17*z3))/z - (32*H01*(4 - 10*z - 115*z2 + 18*z3))/(9.*z) +
               (((-32*H0*H1)/9. - (16*H12)/9.)*(-4 + 13*z - 80*z2 + 87*z3))/z - (16*H1*(-40 + 83*z - 1807*z2 + 1934*z3))/(27.*z) -
               (8*(16 - 181*z - 2056*z2 + 2306*z3))/(9.*z) + (128*H0*(-5 + z)*z*zeta2)/3. - (32*H1*(-1 + 2*z)*zeta2)/3. - (32*Hm1*(1 + 2*z + 4*z2)*zeta2)/3. +
               (32*(-4 + 3*z - 181*z2 + 105*z3)*zeta2)/(9.*z) + (32*(3 - 22*z + 6*z2)*zeta3)/3.)) +
         CA*TR*TR*((64*H0*H12*(-1 + z)*z)/3. - (8*H04*(1 + 2*z))/9. + (64*H01*H02*(1 + 4*z))/3. - (32*H011*z*(-1 + 5*z))/3. - (16*H13*(-1 + z)*(1 + 5*z))/9. -
            (128*H0*H001*(1 + 6*z))/3. + (32*H0001*(1 + 22*z))/3. + (16*H03*(3 - 6*z + 23*z2))/27. - (16*H02*H1*(-1 + z)*(8 + 17*z + 65*z2))/(9.*z) +
            (8*H12*(67 - 143*z + 206*z2))/27. - (8*H02*(229 + 592*z + 822*z2))/27. + (8*H0*(2000 - 4390*z + 16303*z2))/81. +
            LQ3*((8*gammaqg*H1)/9. + (32*H0*(1 + 4*z))/9. - (16*(-1 + z)*(4 + 7*z + 31*z2))/(27.*z)) +
            LM3*((-56*gammaqg*H1)/9. - (224*H0*(1 + 4*z))/9. + (112*(-1 + z)*(4 + 7*z + 31*z2))/(27.*z)) - (32*H0*H01*(8 + 15*z + 99*z2 + 23*z3))/(9.*z) +
            (16*H001*(16 + 61*z + 244*z2 + 222*z3))/(9.*z) - (32*H01*(-40 + 76*z - 257*z2 + 672*z3))/(27.*z) + (16*H0*H1*(-80 + 146*z - 1147*z2 + 1128*z3))/(27.*z) +
            (8*H1*(-628 + 106*z - 6878*z2 + 6819*z3))/(81.*z) - (4*(-13244 + 17367*z - 36798*z2 + 30335*z3))/(243.*z) - (8*H02*(9 + 2*z)*zeta2)/3. -
            (160*H1*(1 - 5*z + 5*z2)*zeta2)/9. + (8*H0*(3 + 360*z + 272*z2)*zeta2)/9. - (4*(-128 + 59*z - 1628*z2 + 776*z3)*zeta2)/(9.*z) +
            z*(1 + z)*((-128*H00m1)/3. + 128*H0m1 + (128*H0*H0m1)/3. + (256*H0m1m1)/3. + (-128*H0 - (64*H02)/3. - (256*H0m1)/3.)*Hm1 + (128*H0*Hm12)/3. +
               (128*Hm1*zeta2)/3.) + LM2*(16*H02 + gammaqg*((-8*H0*H1)/3. - (20*H12)/3.) - (32*H01*(1 + 4*z))/3. + (-32*H0m1 + 32*H0*Hm1)*(1 + 2*z + 2*z2) -
               (16*H0*(25 + 232*z + 127*z2))/9. + (16*H1*(-4 + 13*z - 128*z2 + 135*z3))/(9.*z) + (8*(-200 - 33*z - 2514*z2 + 2612*z3))/(27.*z) - (64*(-6 + z)*z*zeta2)/3.
               ) + LQ2*(gammaqg*((-8*H0*H1)/3. - (4*H12)/3.) - (32*H01*(1 + 4*z))/3. - (16*H02*(1 + 8*z))/3. + ((32*H0m1)/3. - (32*H0*Hm1)/3.)*(1 + 2*z + 2*z2) +
               (16*H0*(-13 - 136*z + 49*z2))/9. + LM*((8*gammaqg*H1)/3. + (32*H0*(1 + 4*z))/3. - (16*(-1 + z)*(4 + 7*z + 31*z2))/(9.*z)) +
               (16*H1*(-4 + 13*z - 80*z2 + 87*z3))/(9.*z) + (8*(-40 - 87*z - 894*z2 + 1048*z3))/(27.*z) - (64*(-2 + z)*z*zeta2)/3.) - (32*(2 + 77*z)*zeta22)/15. +
            (64*H0*(14 + 41*z)*zeta3)/9. - (16*(-28 + 108*z + 150*z2 + 1255*z3)*zeta3)/(27.*z) +
            gammaqg*((-80*H0011)/3. - (16*H0111)/3. - (16*H012)/3. + ((80*H001)/3. + (16*H011)/3.)*H1 + (4*H02*H12)/3. + H0*((80*H011)/3. - 16*H01*H1 - (16*H13)/9.) +
               (2*H14)/9. + (10*H12*zeta2)/3. - (40*H1*zeta3)/9.) +
            (1 + 2*z + 2*z2)*((-64*H000m1)/3. + (256*H001m1)/3. + (256*H00m11)/3. - (128*H00m1m1)/3. - (32*H02*H0m1)/3. + (128*H0m101)/3. +
               H0*((64*H00m1)/3. - (128*H01m1)/3. - (128*H0m11)/3. + (128*H0m1m1)/3.) + (256*H0m1m1m1)/3. +
               ((-256*H001)/3. + (128*H00m1)/3. + (32*H03)/9. + H0*((128*H01)/3. - (128*H0m1)/3.) - (256*H0m1m1)/3.)*Hm1 + ((32*H02)/3. + (128*H0m1)/3.)*Hm12 -
               (128*H0*Hm13)/9. + ((80*H0m1)/3. - (80*H0*Hm1)/3. - (64*Hm12)/3.)*zeta2 + 64*Hm1*zeta3) +
            LQ*(gammaqg*((-32*H01*H1)/3. + (16*H0*H12)/3. + (8*H13)/9.) + (448*H03*z)/9. + ((64*H0m1m1)/3. - (64*H0m1*Hm1)/3. + (32*H0*Hm12)/3.)*(1 + 2*z) -
               (64*H001*(2 - 2*z + z2))/3. + (64*H0*H01*(2 + 2*z + z2))/3. + (64*H00m1*(2 + 4*z + z2))/3. +
               ((-64*H01m1)/3. - (64*H0m11)/3. + (64*H01*Hm1)/3.)*(1 + 2*z + 2*z2) - (128*H011*(1 - 5*z + 3*z2))/3. - (32*H02*H1*(2 - 4*z + 3*z2))/3. -
               (64*H0*H0m1*(2 + 4*z + 3*z2))/3. + (32*H02*Hm1*(2 + 4*z + 5*z2))/3. - (16*H02*(-10 - 296*z + 63*z2))/9. - (64*H0*(-31 - 832*z + 325*z2))/27. +
               LM2*((8*gammaqg*H1)/3. + (32*H0*(1 + 4*z))/3. - (16*(-1 + z)*(4 + 7*z + 31*z2))/(9.*z)) - (64*H01*(2 - 5*z - 59*z2 + 9*z3))/(9.*z) +
               (((-64*H0m1)/9. + (64*H0*Hm1)/9.)*(-2 - 4*z + 7*z2 + 17*z3))/z + (((-32*H0*H1)/9. - (16*H12)/9.)*(-4 + 13*z - 80*z2 + 87*z3))/z -
               (32*H1*(-20 + 65*z - 964*z2 + 1032*z3))/(27.*z) - (8*(-412 - 969*z - 23556*z2 + 25657*z3))/(81.*z) + (128*H0*(-5 + z)*z*zeta2)/3. -
               (32*H1*(-1 + 2*z)*zeta2)/3. - (32*Hm1*(1 + 2*z + 4*z2)*zeta2)/3. + (32*(-4 + 3*z - 184*z2 + 105*z3)*zeta2)/(9.*z) +
               LM*(gammaqg*((-16*H0*H1)/3. - (8*H12)/3.) - (64*H01*(1 + 4*z))/3. - (32*H02*(3 + 4*z))/3. + ((64*H0m1)/3. - (64*H0*Hm1)/3.)*(1 + 2*z + 2*z2) +
                  (32*H0*(13 - 8*z + 101*z2))/9. - (32*(-8 + 11*z - 14*z2 + 8*z3))/(3.*z) + (32*H1*(-4 - 7*z - 40*z2 + 47*z3))/(9.*z) - (128*(-2 + z)*z*zeta2)/3.) +
               (32*(3 - 22*z + 6*z2)*zeta3)/3.) + LM*(gammaqg*((16*H0*H12)/3. - (8*H13)/9.) - (64*H03*(-2 + z))/9. - (32*H02*H1*(1 + z2 - 2*z))/3. -
               (64*H0*H0m1*(1 + z2 + 2*z))/3. + (64*H001*(-3 - 20*z + z2))/3. - (64*H0*H01*(-3 - 12*z + z2))/3. - (64*H00m1*(-1 - 2*z + z2))/3. +
               (128*H011*(1 + z + z2))/3. + ((-64*H01m1)/3. - (64*H0m11)/3. + (64*H01*Hm1)/3.)*(1 + 2*z + 2*z2) + (32*H02*Hm1*(1 + 2*z + 3*z2))/3. +
               ((-64*H0m1m1)/3. + (64*H0m1*Hm1)/3. - (32*H0*Hm12)/3.)*(1 + 2*z + 4*z2) - (32*H02*(5 + 10*z + 54*z2))/9. - (16*H0*(-269 + 1744*z + 1680*z2))/27. -
               (32*H01*(12 + 25*z + 70*z2 + 2*z3))/(9.*z) - (64*H12*(-1 - z - 4*z2 + 5*z3))/(9.*z) + (((-64*H0m1)/9. + (64*H0*Hm1)/9.)*(-2 - 14*z - z2 + 9*z3))/z -
               (128*H0*H1*(-3 - 4*z - 22*z2 + 28*z3))/(9.*z) + (16*H1*(-144 + 253*z - 548*z2 + 476*z3))/(27.*z) + (8*(-172 + 159*z - 20616*z2 + 20863*z3))/(81.*z) +
               (128*H0*(-1 - 3*z + z2)*zeta2)/3. + (32*H1*(3 - 6*z + 4*z2)*zeta2)/3. - (32*Hm1*(3 + 6*z + 8*z2)*zeta2)/3. +
               (32*(-4 + 9*z - 20*z2 + 114*z3)*zeta2)/(9.*z) + (32*(7 + 26*z + 6*z2)*zeta3)/3.)) +
         CF*nf*TR*TR*((16*H13*(-1 + z)*(1 + 3*z))/9. - (4*H04*(-2 + z)*(1 + 6*z))/9. - (32*H0*H001*(-17 + 28*z + 2*z2))/3. + (32*H011*(-1 - 6*z + 3*z2))/3. +
            (16*H01*H02*(-11 + 4*z + 6*z2))/3. + (8*H12*(4 - 17*z + 12*z2))/3. - (32*H0001*(19 - 74*z + 14*z2))/3. + (16*H1*(-1 - 17*z + 20*z2))/3. +
            (4*H03*(69 - 54*z + 404*z2))/27. - (4*H02*(-210 - 1329*z + 3112*z2))/27. + (8*H0*(2631 + 4929*z + 27824*z2))/81. -
            (8*H02*H1*(-16 + 159*z - 234*z2 + 94*z3))/(9.*z) - (32*H0*H01*(8 - 39*z + 222*z2 + 101*z3))/(9.*z) + (16*H001*(16 + 3*z + 654*z2 + 498*z3))/(9.*z) -
            (16*H01*(-80 + 18*z - 1341*z2 + 1448*z3))/(27.*z) + (16*H0*H1*(-80 + 72*z - 1539*z2 + 1556*z3))/(27.*z) -
            (4*(-5248 + 17163*z - 226026*z2 + 221158*z3))/(243.*z) +
            LQ3*((16*gammaqg*H1)/9. + (16*H02*(-1 + 2*z))/3. + (16*H0*(-11 + 4*z))/9. - (16*(-8 + 84*z - 147*z2 + 62*z3))/(27.*z)) +
            LM3*((8*gammaqg*H1)/9. - (16*H02*(-1 + 2*z))/3. - (32*H0*(-4 - z + 6*z2))/9. + (8*(-16 + 159*z - 258*z2 + 124*z3))/(27.*z)) +
            LM2*((32*H03*(-1 + 2*z))/3. + (16*H02*(3 + 2*z)*(-3 + 4*z))/3. + (32*H1*(5 - 4*z + 4*z2))/9. - (8*H0*(305 + 146*z + 160*z2))/9. +
               (4*(-208 - 2247*z + 1356*z2 + 1000*z3))/(27.*z) + gammaqg*((8*H01)/3. + (4*H12)/3. - (8*zeta2)/3.)) + (16*H1*(-5 - 2*z + 2*z2)*zeta2)/9. -
            (8*H02*(-11 - 5*z + 14*z2)*zeta2)/3. + (8*H0*(193 + 166*z + 240*z2)*zeta2)/9. - (2*(-368 - 555*z - 5646*z2 + 7280*z3)*zeta2)/(27.*z) -
            (64*(-11 + 58*z + 2*z2)*zeta22)/15. + (32*H0*(-28 - 133*z + 78*z2)*zeta3)/9. - (8*(-16 + 1095*z + 2046*z2 + 3784*z3)*zeta3)/(27.*z) +
            gammaqg*((40*H0111)/3. + (16*H012)/3. + ((-32*H001)/3. - (32*H011)/3.)*H1 - (4*H03*H1)/9. + H0*((-32*H011)/3. + (16*H01*H1)/3.) + (8*H01*H12)/3. - H14/9. +
               ((-20*H01)/3. - (4*H0*H1)/3. - 4*H12)*zeta2 + (88*H1*zeta3)/9.) +
            LQ2*(gammaqg*((-32*H0*H1)/3. - 4*H12) - (8*H02*(-31 + 56*z + 8*z2))/3. - (16*H01*(-5 - 8*z + 12*z2))/3. + (8*H0*(571 - 236*z + 244*z2))/9. +
               LM*((-8*gammaqg*H1)/3. - (8*(-1 + 4*z))/3. + (16*H0*(1 - 2*z + 4*z2))/3.) + (16*H1*(-8 + 112*z - 215*z2 + 130*z3))/(9.*z) +
               (4*(256 + 4977*z - 7632*z2 + 2156*z3))/(27.*z) - (16*(13 - 8*z + 4*z2)*zeta2)/3. + (-1 + 2*z)*(-32*H001 - 16*H03 + 32*H0*zeta2 + 32*zeta3)) +
            LQ*(gammaqg*((32*H01*H1)/3. + 8*H0*H12 + (8*H13)/3.) + (128*H0*H0m1*(-2 + z)*(-2 + 3*z))/3. - (64*H02*H1*(2 - 4*z + 5*z2))/3. -
               (128*H00m1*(9 - 14*z + 7*z2))/3. - (32*H001*(15 - 60*z + 8*z2))/3. + (32*H0*H01*(-5 - 8*z + 16*z2))/3. + (8*H03*(-35 + 108*z + 16*z2))/3. +
               (32*H011*(1 - 20*z + 24*z2))/3. + (16*H01*(-16 - 329*z - 254*z2 + 76*z3))/(9.*z) - (8*H02*(5320 - 3745*z + 1120*z2 + 144*z3))/45. -
               (16*H12*(-8 + 131*z - 253*z2 + 168*z3))/(9.*z) - (32*H0*H1*(-8 + 140*z - 283*z2 + 198*z3))/(9.*z) - (16*H0*(4 + 7593*z - 5407*z2 + 1984*z3))/(45.*z) -
               (8*H1*(256 + 5247*z - 8448*z2 + 3080*z3))/(27.*z) - (8*(-6952 + 229863*z - 318513*z2 + 91102*z3))/(405.*z) +
               (((-64*H0m1)/45. + (64*H0*Hm1)/45.)*(1 - 20*z + 225*z2 + 40*z3 - 155*z4 + 36*z5))/z2 - (32*H0*(-33 + 84*z + 4*z2)*zeta2)/3. +
               (64*H1*(3 - 6*z + 8*z2)*zeta2)/3. + (16*(-80 + 3045*z - 1400*z2 + 1600*z3 + 144*z4)*zeta2)/(45.*z) +
               (1 + z2 + 2*z)*((-256*H0m1m1)/3. + ((64*H02)/3. + (256*H0m1)/3.)*Hm1 - (128*H0*Hm12)/3. - (128*Hm1*zeta2)/3.) +
               LM*(gammaqg*((32*H0*H1)/3. + (8*H12)/3.) + (32*H01*(3 - 6*z + 4*z2))/3. + (8*(-61 + 82*z + 36*z2))/9. - (16*H1*(41 - 88*z + 76*z2))/9. -
                  (16*H0*(25 - 38*z + 76*z2))/9. + (1 - 2*z + 4*z2)*((-32*H02)/3. + (32*zeta2)/3.)) + (128*(11 - 21*z + 4*z2)*zeta3)/3. +
               (-1 + 2*z)*(64*H0*H001 + 64*H0011 + (44*H04)/3. - 96*H02*zeta2 - (32*zeta22)/5. - 128*H0*zeta3)) +
            LM*((-128*H001*(-7 + 5*z))/3. - (32*H0*H01*(17 - 16*z + 4*z2))/3. - (32*H011*(7 - 14*z + 12*z2))/3. - (8*H03*(-47 + 28*z + 16*z2))/9. +
               (8*H12*(67 - 104*z + 86*z2))/9. + (8*H1*(163 - 258*z + 152*z2))/9. - (80*H0*(-297 - 213*z + 176*z2))/27. + (8*H02*(410 + 137*z + 422*z2))/9. +
               (16*H01*(-16 + 131*z - 88*z2 + 16*z3))/(9.*z) + (32*H0*H1*(8 - 40*z + 29*z2 + 35*z3))/(9.*z) + (8*(-1064 + 16893*z - 19683*z2 + 3224*z3))/(81.*z) -
               (64*H0*(1 - 2*z + 4*z2)*zeta2)/3. - (16*(51 - 30*z + 86*z2)*zeta2)/9. + gammaqg*((-16*H01*H1)/3. - (32*H02*H1)/3. - (8*H0*H12)/3. + (32*H1*zeta2)/3.) +
               (32*(-19 + 2*z + 16*z2)*zeta3)/3. + (-1 + 2*z)*(-192*H0001 + 64*H0*H001 - (20*H04)/3. + (384*zeta22)/5. + 128*H0*zeta3)) +
            (-1 + 2*z)*(128*H00001 - 128*H0*H0001 + 32*H001*H02 - (4*H05)/15. - 8*H03*zeta2 + (208*H02*zeta3)/3. - 128*zeta5)) +
         CA*CA*TR*(-64*H001*H02*(-1 + z) + 96*H011*H02*(-3 + z)*z +
            (-192*H01m1 - 128*H01m1m1 - 192*H0m11 - 128*H0m11m1 - 128*H0m1m11 + (192*H01 + 128*H01m1 + 128*H0m11)*Hm1 - 64*H01*Hm12)*z*(1 + z) +
            ((4*H05)/15. - 64*H0*H0m101)*(-1 + 4*z) - (32*H01*H03*(1 + 4*z))/3. + 128*H0*H0001*(-1 + 7*z) - 32*H00001*(-1 + 54*z) + 32*H0m12*z2 +
            (256*H0001m1 + 64*H0*H000m1 + 256*H000m11)*(2 - 5*z + z2) +
            (-128*H00m10m1 - 1024*H00m1m1m1 + (32*H03*H0m1)/3. - 512*H0m10m1m1 - 64*H0*H0m12 + H0m1*(128*H00m1 + 256*H0m1m1))*(1 - z + z2) - 64*H01111*(1 + z + z2) -
            192*H001*H01*(2 + 5*z + z2) + 32*H0*H012*(3 + 9*z + z2) - 32*H000m1m1*(7 - 22*z + 2*z2) +
            (-32*H0000m1 + 32*H0010m1 + 32*H00m101 - 32*H001*H0m1)*(5 - 14*z + 2*z2) - 16*H00m1*H02*(3 - 6*z + 2*z2) - 128*H0*H0m1m1m1*(3 + 4*z2) +
            64*H00101*(11 + 29*z + 5*z2) - 32*H0*H0111*(1 - 14*z + 6*z2) - 32*H0*H0011*(9 + 14*z + 6*z2) - 64*H01011*(5 - 7*z + 9*z2) - 128*H00111*(7 - 8*z + 12*z2) +
            64*H00011*(32 + 83*z + 13*z2) + (H04*(23 - 26*z + 18*z2))/9. + (8*H01*H02*(7 - 62*z + 32*z2))/3. + (8*H03*H1*(-1 + z)*(8 + 17*z + 65*z2))/(9.*z) +
            (2*H03*(-35 - 12*z + 225*z2))/9. + (((-32*H0m1m1)/9. + (32*H0m1*Hm1)/9. - (16*H0*Hm12)/9.)*(1 + z)*(20 - 32*z + 305*z2))/z +
            (8*H02*Hm1*(1 + z)*(20 - 32*z + 359*z2))/(9.*z) + (2*H02*H1*(-1 + z)*(308 - 169*z + 1619*z2))/(9.*z) +
            (((-16*H0m1)/27. + (16*H0*Hm1)/27.)*(1 + z)*(94 - 244*z + 2443*z2))/z + (2*H02*(738 - 6573*z + 5248*z2))/27. - (8*H0111*(-8 - 7*z - 10*z2 + 30*z3))/(3.*z) -
            (4*H03*Hm1*(8 + 5*z + 46*z2 + 60*z3))/(9.*z) + (((16*H00m1m1)/3. - (16*H00m1*Hm1)/3.)*(8 + 5*z + 58*z2 + 72*z3))/z +
            (8*H011*H1*(-8 + 5*z - 70*z2 + 84*z3))/(3.*z) + ((H0*((16*H01m1)/3. + (16*H0m11)/3.) - (16*H0*H01*Hm1)/3.)*(8 + 5*z + 70*z2 + 84*z3))/z +
            (((-32*H001m1)/3. - (32*H00m11)/3. + (32*H001*Hm1)/3.)*(8 + 5*z + 76*z2 + 90*z3))/z - (8*H0*H13*(-8 + 8*z - 82*z2 + 93*z3))/(9.*z) +
            (4*H02*H0m1*(8 + 5*z + 70*z2 + 96*z3))/(3.*z) + (((-16*H0m101)/3. + (16*H0*H0m1*Hm1)/3.)*(8 + 5*z + 82*z2 + 96*z3))/z -
            (16*H0*H0m1m1*(8 + 5*z + 82*z2 + 108*z3))/(3.*z) + (((-32*H0m1m1m1)/3. + (32*H0m1m1*Hm1)/3. - (16*H0m1*Hm12)/3. + (16*H0*Hm13)/9.)*
               (8 + 5*z + 94*z2 + 108*z3))/z - (4*H02*Hm12*(8 + 5*z + 106*z2 + 120*z3))/(3.*z) - (16*H0*H001*(-8 + 28*z - 210*z2 + 129*z3))/(3.*z) -
            (8*H0*H00m1*(8 + 5*z + 94*z2 + 132*z3))/(3.*z) - (8*H012*(-16 - 4*z - 130*z2 + 161*z3))/(3.*z) + (8*H000m1*(8 + 5*z + 118*z2 + 168*z3))/(3.*z) +
            (2*H02*H12*(-24 - 13*z - 142*z2 + 190*z3))/(3.*z) + (16*H0001*(-16 + 69*z - 408*z2 + 203*z3))/(3.*z) + (H14*(-8 - 19*z - 166*z2 + 204*z3))/(9.*z) -
            (8*H0*H12*(-20 + 11*z - 229*z2 + 241*z3))/(3.*z) + (16*H00m1*(20 - 12*z - 105*z2 + 287*z3))/(9.*z) - (16*H0*H0m1*(20 - 12*z + 111*z2 + 323*z3))/(9.*z) -
            (8*H0*H01*H1*(-40 - 3*z - 330*z2 + 406*z3))/(3.*z) + (8*H001*H1*(-72 - 11*z - 590*z2 + 728*z3))/(3.*z) + (8*H0*H011*(-56 + 19*z - 398*z2 + 768*z3))/(3.*z) -
            (8*H0011*(-56 + 31*z - 266*z2 + 932*z3))/(3.*z) + (2*H13*(-4 + 15*z - 1440*z2 + 1321*z3))/(27.*z) + (4*H0*H01*(468 + 93*z + 5892*z2 + 1501*z3))/(9.*z) +
            (4*H011*(-240 + 48*z - 3312*z2 + 3425*z3))/(9.*z) - (4*H001*(628 + 495*z + 9840*z2 + 3525*z3))/(9.*z) + (2*H12*(-134 + 39*z - 3108*z2 + 3536*z3))/(27.*z) +
            (8*H01*(-1388 + 1125*z - 17343*z2 + 9594*z3))/(27.*z) - (8*H0*H1*(-1388 + 1575*z - 16335*z2 + 16040*z3))/(27.*z) +
            (4*H1*(-1465 + 1286*z - 18101*z2 + 18046*z3))/(27.*z) - (4*H0*(2624 + 16377*z + 100242*z2 + 225962*z3))/(81.*z) +
            (4*(-66607 + 83079*z - 1178838*z2 + 1158802*z3))/(243.*z) + (8*H03*(-5 + 8*z)*zeta2)/3. + 8*H0*H01*(11 + 38*z + 2*z2)*zeta2 -
            8*H001*(7 + 46*z + 2*z2)*zeta2 + (24*H00m1 - 16*H0m1m1)*(7 - 10*z + 6*z2)*zeta2 - 16*H011*(17 - 22*z + 30*z2)*zeta2 + (2*H02*(47 - 6*z + 38*z2)*zeta2)/3. -
            (8*H0*H1*(-1 + z)*(20 + 41*z + 137*z2)*zeta2)/(3.*z) - (16*Hm1*(1 + z)*(20 - 32*z + 413*z2)*zeta2)/(9.*z) +
            (8*H0*Hm1*(12 + 2*z + 46*z2 + 67*z3)*zeta2)/(3.*z) + (4*H12*(-8 + 11*z - 94*z2 + 102*z3)*zeta2)/(3.*z) - (8*H0m1*(12 + 2*z + 70*z2 + 103*z3)*zeta2)/(3.*z) -
            (8*H01*(16 + 33*z + 258*z2 + 124*z3)*zeta2)/(3.*z) + (8*Hm12*(8 + 5*z + 118*z2 + 132*z3)*zeta2)/(3.*z) +
            (4*H1*(-242 + 116*z - 3052*z2 + 3111*z3)*zeta2)/(9.*z) - (2*H0*(160 + 1165*z + 1864*z2 + 6474*z3)*zeta2)/(9.*z) +
            (2*(-2176 + 4590*z - 10209*z2 + 57548*z3)*zeta2)/(27.*z) + (5 - 2*z + 6*z2)*(16*H02*H0m1m1 - 24*H0*H0m1*zeta2) +
            LM3*(gammaqg*((8*H0*H1)/3. + (8*H12)/3.) + (16*H02*(-1 + 8*z))/3. + (32*H01*(3 + 6*z + 2*z2))/3. - (4*(-1 + z)*(272 - 97*z + 1883*z2))/(27.*z) -
               (8*H0*(-8 - 11*z - 152*z2 + 18*z3))/(9.*z) - (8*H1*(-16 - z - 118*z2 + 146*z3))/(9.*z) - (64*(1 + 4*z)*zeta2)/3.) +
            LQ3*(gammaqg*((-8*H0*H1)/3. - (8*H12)/3.) - (16*H02*(-1 + 8*z))/3. - (32*H01*(3 + 6*z + 2*z2))/3. + (4*(-1 + z)*(272 - 97*z + 1883*z2))/(27.*z) +
               (8*H0*(-8 - 11*z - 152*z2 + 18*z3))/(9.*z) + (8*H1*(-16 - z - 118*z2 + 146*z3))/(9.*z) + (64*(1 + 4*z)*zeta2)/3.) + (16*H0*(-19 + 124*z)*zeta22)/5. +
            (16*(76 - 90*z + 1673*z2 + 86*z3)*zeta22)/(15.*z) - (32*H02*(-5 + 31*z)*zeta3)/3. + (160*H01*(3 + 6*z + 2*z2)*zeta3)/3. + 32*H0m1*(5 - 8*z + 4*z2)*zeta3 -
            (8*Hm1*(8 + 5*z + 98*z2 + 112*z3)*zeta3)/z + (8*H0*(-8 - 269*z + 340*z2 + 210*z3)*zeta3)/(9.*z) - (8*H1*(-80 - 41*z - 698*z2 + 874*z3)*zeta3)/(9.*z) +
            (4*(448 + 3537*z + 25128*z2 + 16364*z3)*zeta3)/(27.*z) + (8*(83 + 278*z)*zeta2*zeta3)/3. +
            gammaqg*(-40*H01*H011 + (-16*H0001 + 24*H0011 + 20*H012)*H1 + 8*H01*H02*H1 + (-44*H001 - 4*H011)*H12 - (2*H03*H12)/3. +
               H0*((-8*H001 - 56*H011)*H1 + 24*H01*H12 + (5*H14)/3.) - H15/3. + (-32*H01*H1 + H02*H1 + 2*H0*H12 - (8*H13)/3.)*zeta2 + 46*H1*zeta22 +
               ((40*H0*H1)/3. + (52*H12)/3.)*zeta3) + (1 + 2*z + 2*z2)*
             (-128*H0011m1 - 128*H001m11 - 192*H001m1m1 - 128*H00m111 - 192*H00m11m1 - 192*H00m1m11 - 128*H01m1m1m1 - 64*H0m1011 - 128*H0m101m1 - 128*H0m10m11 +
               H02*(32*H01m1 + 32*H0m11) - 128*H0m11m1m1 - 128*H0m1m101 +
               H0*(-64*H001m1 - 64*H00m11 - 96*H00m1m1 + 64*H011m1 + 64*H01m11 + 64*H01m1m1 + 64*H0m111 + 64*H0m11m1 + 64*H0m1m11) - 128*H0m1m11m1 - 128*H0m1m1m11 -
               384*H0m1m1m1m1 + (-32*H0001 - 160*H000m1 + 128*H0011 + 192*H001m1 + 192*H00m11 + 128*H01m1m1 - (4*H04)/3. + H02*(-32*H01 - 16*H0m1) + 128*H0m101 +
                  128*H0m11m1 + H0*(64*H001 + 96*H00m1 - 64*H011 - 64*H01m1 - 64*H0m11 + 128*H0m1m1) + 128*H0m1m11 + 384*H0m1m1m1)*Hm1 +
               (-96*H001 - 64*H01m1 - (8*H03)/3. + H0*(32*H01 - 64*H0m1) - 64*H0m11 - 192*H0m1m1)*Hm12 + ((64*H01)/3. + (64*H02)/3. + 64*H0m1)*Hm13 - 16*H0*Hm14 +
               (-96*H01m1 - 96*H0m11 + (96*H01 + 36*H02 - 16*H0m1)*Hm1 + 40*H0*Hm12 - (160*Hm13)/3.)*zeta2 - (376*Hm1*zeta22)/5. + (-64*H0*Hm1 + 160*Hm12)*zeta3) +
            LQ2*((16*H03*(-3 + 16*z))/3. + 64*H011*(2 + 5*z + z2) + 16*H001*(-5 - 6*z + 2*z2) + 16*H0*H01*(7 + 22*z + 2*z2) + 16*H00m1*(11 - 26*z + 6*z2) -
               16*H0*H0m1*(7 - 10*z + 6*z2) - (4*H02*(-23 - 436*z + 54*z2))/3. + (8*H01*(-8 + 23*z + 272*z2 + 28*z3))/(3.*z) +
               (((-8*H0m1)/3. + (8*H0*Hm1)/3.)*(16 - z + 118*z2 + 146*z3))/z - (8*H0*H1*(-24 + 5*z - 238*z2 + 280*z3))/(3.*z) -
               (4*H12*(-24 + 17*z - 310*z2 + 352*z3))/(3.*z) - (4*H1*(-168 + 59*z - 4204*z2 + 4549*z3))/(9.*z) - (4*H0*(104 + 1009*z - 3512*z2 + 5183*z3))/(9.*z) -
               (2*(1260 - 509*z - 7318*z2 + 6639*z3))/(9.*z) - 32*H0*(3 + 14*z)*zeta2 + (16*(-9 - 196*z + 126*z2)*zeta2)/3. +
               gammaqg*(6*H02*H1 + 12*H0*H12 + 4*H13 - 16*H1*zeta2) +
               (1 + 2*z + 2*z2)*(-32*H01m1 - 32*H0m11 - 64*H0m1m1 + (32*H01 + 24*H02 + 64*H0m1)*Hm1 - 32*H0*Hm12 - 64*Hm1*zeta2) + 32*(-3 + z)*zeta3) +
            LM2*((-32*H03*(-1 + 4*z))/3. + 32*H011*(1 + 4*z) + (16*H01*z*(36 + 25*z))/3. - 16*H0*H01*(3 + 6*z + 2*z2) - 16*H00m1*(11 - 26*z + 6*z2) +
               16*H0*H0m1*(7 - 10*z + 6*z2) + (16*H02*(-3 + z + 9*z2))/3. + (16*H0*H1*(-1 + z)*(4 + 7*z + 19*z2))/(3.*z) - (4*H12*(-8 + 5*z - 118*z2 + 132*z3))/(3.*z) +
               (((8*H0m1)/3. - (8*H0*Hm1)/3.)*(16 - z + 118*z2 + 146*z3))/z - (8*H1*(-28 - 70*z - 604*z2 + 617*z3))/(9.*z) +
               (8*H0*(40 + 266*z + 710*z2 + 1335*z3))/(9.*z) - (8*(-244 + 192*z - 2661*z2 + 2695*z3))/(9.*z) - 64*H0*(-1 + z)*zeta2 -
               (16*(4 - 3*z + 83*z2 + 44*z3)*zeta2)/(3.*z) + gammaqg*(-2*H02*H1 + 4*H13 - 8*H1*zeta2) +
               (1 + 2*z + 2*z2)*(16*H001 + 32*H01m1 + 32*H0m11 + 64*H0m1m1 + (-32*H01 - 24*H02 - 64*H0m1)*Hm1 + 32*H0*Hm12 + 64*Hm1*zeta2) - 32*(-2 + 17*z)*zeta3) +
            LQ*(gammaqg*(24*H01*H0m1 + (-48*H00m1 + (88*H01)/3.)*H1 + H0*(24*H0m1*H1 - (28*H13)/3.) - 2*H14) - 64*H0*H01*H1*(1 + z2 - 2*z) - (4*H04*(-11 + 51*z))/3. +
               32*H0*H011*(-5 - 32*z + 2*z2) + (64*H011m1 + 64*H01m11 + 64*H0m111 + (-64*H0*H01 - 64*H011)*Hm1)*(1 + 2*z + 2*z2) +
               (-64*H0m1m1m1 + 64*H0m1m1*Hm1 - 32*H0m1*Hm12 + (32*H0*Hm13)/3.)*(3 + 6*z + 2*z2) - 32*H0111*(7 + 22*z + 2*z2) + 16*H02*H0m1*(9 - 22*z + 3*z2) -
               32*H012*(2 - 4*z + 3*z2) - 32*H0m12*(3 - 2*z + 3*z2) +
               (128*H01m1m1 + 128*H0m11m1 + 128*H0m1m11 + (-128*H01m1 - 128*H0m11)*Hm1 + 64*H01*Hm12)*(1 + 2*z + 3*z2) + 64*H001*H1*(3 - 6*z + 4*z2) +
               64*H0m101*(5 - 8*z + 5*z2) + 128*H0*H0m1m1*(3 + 2*z + 5*z2) + 32*H0*H00m1*(-7 + 44*z + 5*z2) - 32*H0011*(-1 - 8*z + 6*z2) + 16*H02*H12*(4 - 8*z + 7*z2) -
               64*H0*H0m1*Hm1*(3 + 6*z + 7*z2) + (32*H001m1 + 32*H00m11 - 32*H001*Hm1)*(3 + 6*z + 8*z2) + H0*(32*H01m1 + 32*H0m11)*(5 - 2*z + 10*z2) +
               32*H0*H001*(19 + 24*z + 11*z2) + (-32*H00m1m1 + 32*H00m1*Hm1)*(7 + 14*z + 12*z2) + 8*H02*Hm12*(5 + 10*z + 16*z2) + (8*H03*H1*(13 - 26*z + 22*z2))/3. -
               (8*H03*Hm1*(9 + 18*z + 22*z2))/3. - 16*H0001*(37 + 50*z + 22*z2) - 16*H000m1*(3 + 114*z + 26*z2) - 4*H01*H02*(63 + 134*z + 36*z2) +
               (4*H03*(-83 - 1648*z + 198*z2))/9. + (2*H02*(3849 - 14584*z + 18688*z2))/9. - (16*H001*(4 + 29*z - 125*z2 + 15*z3))/(3.*z) -
               (8*H0*H01*(-40 - 49*z + 833*z2 + 42*z3))/(3.*z) + (16*H011*(8 + z - 488*z2 + 60*z3))/(3.*z) +
               (((16*H0m1m1)/3. - (16*H0m1*Hm1)/3. + (8*H0*Hm12)/3.)*(-8 + 19*z + 92*z2 + 100*z3))/z + (16*H0*H0m1*(36 - 50*z + 194*z2 + 129*z3))/(3.*z) +
               (16*H00m1*(-44 + 107*z - 71*z2 + 131*z3))/(3.*z) + (((16*H01m1)/3. + (16*H0m11)/3. - (16*H01*Hm1)/3.)*(16 - 13*z + 166*z2 + 218*z3))/z +
               (((16*H0m1)/9. - (16*H0*Hm1)/9.)*(-240 - 128*z + 257*z2 + 272*z3))/z - (8*H02*Hm1*(28 + 7*z + 317*z2 + 389*z3))/(3.*z) +
               (8*H13*(-24 + 29*z - 382*z2 + 424*z3))/(9.*z) + (8*H0*H12*(-28 + 37*z - 428*z2 + 477*z3))/(3.*z) + (4*H02*H1*(-72 - 61*z - 631*z2 + 822*z3))/(3.*z) +
               (16*H01*(162 + 395*z - 952*z2 + 1357*z3))/(9.*z) + (16*H0*H1*(-58 + 10*z - 3023*z2 + 3243*z3))/(9.*z) + (8*H12*(-71 + 79*z - 3443*z2 + 3661*z3))/(9.*z) +
               (4*H0*(-1136 - 7837*z - 68257*z2 + 37209*z3))/(27.*z) + (4*H1*(3996 + 469*z - 54734*z2 + 52691*z3))/(27.*z) +
               (2*(-57320 + 89781*z - 529368*z2 + 503837*z3))/(81.*z) - 16*H0m1*(9 - 14*z + 4*z2)*zeta2 - 128*H0*H1*(3 - 6*z + 5*z2)*zeta2 -
               16*Hm12*(1 + 2*z + 10*z2)*zeta2 + 32*H0*Hm1*(5 + 10*z + 14*z2)*zeta2 - 4*H02*(-9 - 230*z + 24*z2)*zeta2 + 16*H01*(31 + 38*z + 24*z2)*zeta2 +
               (8*Hm1*(24 - 7*z + 424*z2 + 536*z3)*zeta2)/(3.*z) - (8*H0*(-16 - 21*z - 1801*z2 + 708*z3)*zeta2)/(3.*z) -
               (8*H1*(-64 + 11*z - 676*z2 + 766*z3)*zeta2)/(3.*z) - (16*(-136 + 405*z - 3718*z2 + 4600*z3)*zeta2)/(9.*z) + (4*(157 - 18*z + 352*z2)*zeta22)/5. +
               32*H0*(17 + 64*z)*zeta3 + 16*Hm1*(3 + 6*z + 20*z2)*zeta3 - (8*(-96 + 15*z - 1106*z2 + 1252*z3)*zeta3)/(3.*z) +
               (5 - 10*z + 8*z2)*(-32*H12*zeta2 - 112*H1*zeta3)) + LM*
             (-32*H0*H001*(1 + z)*(7 + 13*z) + (4*H04*(-5 + 19*z))/3. + 16*H02*H0m1*(-1 + 8*z + z2) - 64*H0m101*(1 - 4*z + 2*z2) + 32*H0m12*(4 - 6*z + 3*z2) +
               64*H012*(2 - z + 3*z2) - 32*H0*H00m1*(3 + 16*z + 3*z2) - 64*H0*H0m1m1*(5 - 4*z + 5*z2) +
               (((-64*H01m1)/3. - (64*H0m11)/3. + (64*H01*Hm1)/3.)*(1 + z)*(2 + 4*z + 5*z2))/z + (32*H011*(-11 + 49*z + 8*z2))/3. + 16*H000m1*(23 + 46*z + 10*z2) +
               32*H0011*(9 + 4*z + 12*z2) - 32*H0*H011*(11 + 8*z + 12*z2) + 16*H0001*(13 + 126*z + 22*z2) + (40*H0*H12*(-1 + z)*(4 + 7*z + 31*z2))/(3.*z) +
               4*H01*H02*(23 - 2*z + 36*z2) - (4*H03*(-37 + 40*z + 54*z2))/9. - (8*H13*(2 - 94*z + 103*z2))/9. + (4*H0*H1*(-1 + z)*(192 - 5*z + 1583*z2))/(3.*z) -
               (2*H02*(1303 + 216*z + 4303*z2))/9. - (4*H02*H1*(-8 - 58*z + 59*z2 + 18*z3))/(3.*z) + (8*H0*H01*(-24 - 86*z - 329*z2 + 40*z3))/(3.*z) -
               (16*H01*H1*(-8 - 17*z - 26*z2 + 40*z3))/(3.*z) + (8*H0*H0m1*(-24 - 5*z - 346*z2 + 128*z3))/(3.*z) + (4*H02*Hm1*(24 + 59*z + 160*z2 + 136*z3))/(3.*z) -
               (8*H001*(-40 - 75*z - 624*z2 + 152*z3))/(3.*z) + (((16*H0m1m1)/3. - (16*H0m1*Hm1)/3. + (8*H0*Hm12)/3.)*(16 - 37*z + 112*z2 + 176*z3))/z -
               (8*H00m1*(-24 + 49*z - 532*z2 + 392*z3))/(3.*z) + (((32*H0m1)/9. - (32*H0*Hm1)/9.)*(62 - 50*z + 476*z2 + 541*z3))/z -
               (2*H12*(-84 - 247*z - 2032*z2 + 2023*z3))/(9.*z) - (4*H01*(-416 + 1353*z - 3816*z2 + 3416*z3))/(9.*z) -
               (4*H1*(492 + 3625*z - 13409*z2 + 9884*z3))/(27.*z) + (8*H0*(448 + 3023*z + 13409*z2 + 30972*z3))/(27.*z) -
               (4*(-24160 + 31944*z - 295527*z2 + 290560*z3))/(81.*z) - 16*H0m1*(13 - 2*z + 16*z2)*zeta2 - 16*H01*(9 - 26*z + 20*z2)*zeta2 +
               4*H02*(-1 - 22*z + 24*z2)*zeta2 + 8*H0*(13 - 31*z + 30*z2)*zeta2 + (8*Hm1*(-85 + 40*z + 136*z2)*zeta2)/3. -
               (8*H1*(-8 + 41*z - 76*z2 + 54*z3)*zeta2)/(3.*z) - (4*(336 - 762*z + 2860*z2 + 1333*z3)*zeta2)/(9.*z) - (4*(245 + 1758*z + 112*z2)*zeta22)/5. -
               96*H0*(3 + 22*z)*zeta3 + (32*(-24 - 14*z - 343*z2 + 180*z3)*zeta3)/(3.*z) +
               gammaqg*(-8*H0111 - 24*H01*H0m1 + (64*H001 + 48*H00m1)*H1 + (10*H03*H1)/3. + 8*H01*H12 - 4*H02*H12 +
                  H0*(24*H01m1 + 24*H0m11 + (-16*H01 - 24*H0m1)*H1 - (20*H13)/3.) + 2*H14 + (-40*H0*H1 - 8*H12)*zeta2 - 136*H1*zeta3) +
               (1 + 2*z + 2*z2)*(-32*H001m1 - 32*H00m11 + 160*H00m1m1 + 64*H011m1 + 64*H01m11 + 64*H0m111 + 320*H0m1m1m1 +
                  (32*H001 - 160*H00m1 - 64*H011 + (8*H03)/3. + 64*H0*H0m1 - 320*H0m1m1)*Hm1 + (8*H02 + 160*H0m1)*Hm12 - (160*H0*Hm13)/3. +
                  (64*H0*Hm1 - 80*Hm12)*zeta2 + 192*Hm1*zeta3)) + 272*(-1 + 8*z)*zeta5) +
         CF*CF*TR*(-77 + 192*H0*H0011*(1 + z2 - 2*z) - 84*z + 16*H01*H03*(-1 + 2*z) + (48*H001*H1 - 24*H0*H01*H1)*(-3 + 2*z)*(-1 + 2*z) +
            (-32*H011*H1 + 8*H01*H12)*(-1 + 4*z) - (4*H0*H13*(-3 + 4*z)*(1 + 4*z))/3. - 32*H0111*(-1 + 2*z)*(-5 + 7*z) - 16*H0*H011*(-1 + z)*(-3 + 16*z) -
            8*H0*H001*(-25 + 64*z) - 78*z2 + 64*H01*H011*(2 - 4*z + 3*z2) - 96*H01011*(3 - 6*z + 4*z2) + 64*H001*H01*(1 - 2*z + 5*z2) + (2*H05*(1 - 2*z + 8*z2))/15. -
            4*H01*H02*(17 - 60*z + 12*z2) - 4*H012*(7 - 24*z + 12*z2) + 32*H0*H0001*(-1 + 2*z + 12*z2) - 8*H001*H02*(-5 + 10*z + 12*z2) +
            (16*H03*H1*(3 - 17*z + 13*z2))/3. - 16*H0*H01*(3 + 7*z + 13*z2) + 8*H011*H02*(3 - 6*z + 16*z2) - 4*H02*H1*(-7 - 18*z + 18*z2) -
            8*H0*H012*(7 - 14*z + 20*z2) - 32*H00101*(3 - 6*z + 20*z2) + 8*H0001*(-33 + 67*z + 24*z2) + 8*H0011*(17 - 60*z + 32*z2) - 16*H00111*(25 - 50*z + 36*z2) -
            16*H0*H0111*(19 - 38*z + 36*z2) - (2*H03*(5 + 2*z + 36*z2))/3. - 4*H011*(-37 + 102*z + 36*z2) + 2*H02*H12*(15 - 56*z + 40*z2) -
            16*H01111*(17 - 34*z + 40*z2) - 8*H00001*(1 - 2*z + 48*z2) + (H04*(-3 - 16*z + 52*z2))/3. - (4*H13*(-31 - 28*z + 62*z2))/3. - (H14*(19 - 76*z + 72*z2))/3. -
            2*H0*H12*(-7 - 86*z + 88*z2) - 16*H00011*(31 - 62*z + 92*z2) + 2*H001*(-17 + 30*z + 104*z2) - 2*H01*(22 - 7*z + 112*z2) + 2*H12*(74 - 129*z + 136*z2) +
            4*H0*(-13 - 91*z + 136*z2) + 4*H0*H1*(25 - 156*z + 192*z2) + 2*H1*(39 - 313*z + 272*z2) + H02*(63 - 196*z + 384*z2) + 96*gammaqg*ln2*zeta2 +
            (16*H0m1 - 16*H0*Hm1)*(1 + z)*(2 + 3*z)*zeta2 - (2*H03*(-3 + 2*z)*(-1 + 4*z)*zeta2)/3. + 16*H12*z*(-9 + 8*z)*zeta2 - 4*H0*H1*(-3 + 4*z)*(-9 + 16*z)*zeta2 -
            32*H0*H0m1*z2*zeta2 - 16*H011*(3 - 6*z + 2*z2)*zeta2 + 16*H00m1*(1 + 2*z + 6*z2)*zeta2 - 8*H001*(37 - 74*z + 20*z2)*zeta2 +
            8*H0*H01*(29 - 58*z + 44*z2)*zeta2 - 2*H02*(4 - 37*z + 52*z2)*zeta2 + 4*H1*(-10 - 82*z + 61*z2)*zeta2 + 2*H0*(29 + 38*z + 122*z2)*zeta2 +
            4*H01*(47 - 166*z + 128*z2)*zeta2 - ((-375 - 510*z + 232*z2)*zeta2)/2. +
            LQ3*(gammaqg*((-8*H0*H1)/3. - (8*H12)/3.) - 8*H0*z - (2*(-11 + 2*z))/3. - (16*H1*(-1 + 4*z))/3. + (64*H01*z2)/3. + (4*H02*(1 - 2*z + 8*z2))/3. -
               (32*(1 - 2*z + 4*z2)*zeta2)/3.) + LM3*(gammaqg*((8*H0*H1)/3. + (8*H12)/3.) + 8*H0*z + (2*(-11 + 2*z))/3. + (16*H1*(-1 + 4*z))/3. - (64*H01*z2)/3. -
               (4*H02*(1 - 2*z + 8*z2))/3. + (32*(1 - 2*z + 4*z2)*zeta2)/3.) - (8*(1 - 79*z + 176*z2)*zeta22)/5. - (8*H0*(57 - 94*z + 284*z2)*zeta22)/5. +
            (1 + 2*z + 2*z2)*((32*H0m1m1 + (-8*H02 - 32*H0m1)*Hm1 + 16*H0*Hm12)*zeta2 + 16*Hm1*zeta22) + 8*H0*(-3 + 4*z)*(1 + 8*z)*zeta3 +
            (128*H01*(3 - 6*z + 5*z2)*zeta3)/3. + (64*H1*(-2 - 7*z + 12*z2)*zeta3)/3. - (16*H02*(5 - 10*z + 22*z2)*zeta3)/3. + (4*(13 - 34*z + 252*z2)*zeta3)/3. -
            (16*(-1 + 14*z + 20*z2)*zeta2*zeta3)/3. + gammaqg*(20*H012*H1 + 12*H01*H02*H1 - (2*H04*H1)/3. + (-32*H001 - 16*H011)*H12 - (4*H03*H12)/3. +
               (16*H01*H13)/3. + H0*((-24*H001 - 56*H011)*H1 + 20*H01*H12 + H14) + H15/3. + (-24*H01*H1 + 4*H02*H1 - 18*H0*H12 - (40*H13)/3.)*zeta2 +
               (264*H1*zeta22)/5. + ((88*H0*H1)/3. + (4*H12)/3.)*zeta3) +
            LM2*(-33 - 92*z + (-32*H0m1 + 32*H0*Hm1)*(1 + z)*(2 + 3*z) - 32*H12*(-1 + z)*(-2 + 5*z) - 4*H02*z*(9 + 32*z) + 104*z2 + 64*H0*H0m1*z2 -
               32*H011*(-1 + 2*z + 2*z2) - 32*H0*H01*(1 - 2*z + 6*z2) - 32*H00m1*(1 + 2*z + 6*z2) - 16*H0*H1*(5 - 18*z + 10*z2) - 4*H0*(11 - 6*z + 18*z2) -
               8*H01*(-4 + 5*z + 20*z2) + 8*H001*(7 - 14*z + 24*z2) - (4*H03*(3 + 2*z + 24*z2))/3. - 2*H1*(65 - 126*z + 36*z2) + 8*(6 - 11*z + 40*z2)*zeta2 +
               8*H0*(9 - 10*z + 40*z2)*zeta2 + gammaqg*(8*H02*H1 + 20*H0*H12 + 8*H13 - 32*H1*zeta2) +
               (1 + 2*z + 2*z2)*(-64*H0m1m1 + (16*H02 + 64*H0m1)*Hm1 - 32*H0*Hm12 - 32*Hm1*zeta2) + 8*(7 + 2*z + 32*z2)*zeta3) +
            LQ2*(-33 - 92*z + (-32*H0m1 + 32*H0*Hm1)*(1 + z)*(2 + 3*z) - 32*H12*(-1 + z)*(-2 + 5*z) - 4*H02*z*(9 + 32*z) + 104*z2 + 64*H0*H0m1*z2 -
               32*H011*(-1 + 2*z + 2*z2) - 32*H0*H01*(1 - 2*z + 6*z2) - 32*H00m1*(1 + 2*z + 6*z2) - 16*H0*H1*(5 - 18*z + 10*z2) - 4*H0*(11 - 6*z + 18*z2) -
               8*H01*(-4 + 5*z + 20*z2) + 8*H001*(7 - 14*z + 24*z2) - (4*H03*(3 + 2*z + 24*z2))/3. - 2*H1*(65 - 126*z + 36*z2) + 8*(6 - 11*z + 40*z2)*zeta2 +
               8*H0*(9 - 10*z + 40*z2)*zeta2 + gammaqg*(8*H02*H1 + 20*H0*H12 + 8*H13 - 32*H1*zeta2) +
               (1 + 2*z + 2*z2)*(-64*H0m1m1 + (16*H02 + 64*H0m1)*Hm1 - 32*H0*Hm12 - 32*Hm1*zeta2) +
               LM*(gammaqg*(8*H0*H1 + 8*H12) + 24*H0*z + 2*(-11 + 2*z) + 16*H1*(-1 + 4*z) - 64*H01*z2 - 4*H02*(1 - 2*z + 8*z2) + 32*(1 - 2*z + 4*z2)*zeta2) +
               8*(7 + 2*z + 32*z2)*zeta3) + LM*((-64*H0m1m1 + 64*H0m1*Hm1 - 32*H0*Hm12)*(1 + z)*(14 + z) - 192*H0111*(-1 + 2*z) + 128*H02*Hm1*(1 + z)*(3 + 2*z) +
               (-64*H01m1 - 64*H0m11 + 64*H01*Hm1)*(1 + z)*(2 + 3*z) + 16*H01*H1*(-1 + 4*z) - 32*H0m12*(-1 + 14*z) - (8*H03*z*(29 + 82*z))/3. +
               (-128*H0m101 - 128*H0*H0m1m1)*(2 - 4*z + 3*z2) + 128*H00m1*(-8 + 8*z + 3*z2) - 128*H0*H00m1*(1 + 4*z2) + 128*H0*H011*(3 - 6*z + 4*z2) -
               64*H0*H001*(3 + 10*z + 4*z2) - 48*H13*(2 - 8*z + 7*z2) - 64*H0*H0m1*(-2 + 18*z + 7*z2) - 16*H01*H02*(-1 - 14*z + 8*z2) - 16*H012*(3 - 6*z + 8*z2) +
               32*H02*H0m1*(3 - 2*z + 10*z2) - 64*H000m1*(5 - 2*z + 10*z2) - 64*H0*H12*(4 - 13*z + 11*z2) + 16*H001*(26 - 77*z + 26*z2) - 16*H0011*(17 - 34*z + 32*z2) -
               16*H01*(-13 - 6*z + 37*z2) - 16*H011*(-15 - 3*z + 38*z2) - (2*H04*(5 + 6*z + 40*z2))/3. + 16*H0001*(27 + 34*z + 48*z2) - 8*H02*H1*(-3 - 40*z + 50*z2) -
               16*H0*H1*(29 - 95*z + 55*z2) - 2*H1*(-85 - 216*z + 272*z2) - 2*H12*(155 - 546*z + 368*z2) + (2*H02*(-585 - 5510*z - 4620*z2 + 576*z3))/15. +
               (2*(-32 + 383*z - 1893*z2 + 1092*z3))/(15.*z) - (4*H0*(-16 + 308*z + 1593*z2 + 1464*z3))/(15.*z) +
               (((32*H0m1)/15. - (32*H0*Hm1)/15.)*(2 - 915*z2 - 1090*z3 - 165*z4 + 72*z5))/z2 - 64*H01*(3 + 2*z)*zeta2 - 32*Hm1*(1 + z)*(18 + 7*z)*zeta2 -
               32*H0m1*(5 - 6*z + 12*z2)*zeta2 + 8*H02*(15 - 22*z + 72*z2)*zeta2 + 16*H1*(5 - 78*z + 86*z2)*zeta2 + 16*H0*(9 - 31*z + 100*z2)*zeta2 -
               (16*(-240 - 665*z - 1380*z2 + 144*z3)*zeta2)/15. - (8*(57 - 346*z + 192*z2)*zeta22)/5. + 32*H0*(9 - 30*z + 44*z2)*zeta3 + 16*(23 - 54*z + 128*z2)*zeta3 +
               gammaqg*(-64*H01*H0m1 + (-80*H001 + 128*H00m1)*H1 + 8*H03*H1 + 8*H01*H12 + 20*H02*H12 +
                  H0*(76*H01 + 64*H01m1 + 64*H0m11 + (64*H01 - 64*H0m1)*H1 + 24*H13) + 8*H14 + (-128*H0*H1 - 64*H12)*zeta2 - 120*H1*zeta3) +
               (1 + 2*z + 2*z2)*(-64*H001m1 - 64*H00m11 - 64*H00m1m1 - 128*H01m1m1 - 128*H0m11m1 - 128*H0m1m11 - 256*H0m1m1m1 +
                  (64*H001 + 64*H00m1 + 128*H01m1 + (64*H03)/3. + 192*H0*H0m1 + 128*H0m11 + 256*H0m1m1)*Hm1 + (-64*H01 - 112*H02 - 128*H0m1)*Hm12 + (128*H0*Hm13)/3. +
                  (-64*H0*Hm1 + 128*Hm12)*zeta2 - 224*Hm1*zeta3)) +
            LQ*((64*H0m1m1 - 64*H0m1*Hm1 + 32*H0*Hm12)*(1 + z)*(14 + z) + 192*H0111*(-1 + 2*z) - 128*H02*Hm1*(1 + z)*(3 + 2*z) +
               (64*H01m1 + 64*H0m11 - 64*H01*Hm1)*(1 + z)*(2 + 3*z) - 16*H01*H1*(-1 + 4*z) + 32*H0m12*(-1 + 14*z) + (8*H03*z*(29 + 82*z))/3. +
               (128*H0m101 + 128*H0*H0m1m1)*(2 - 4*z + 3*z2) - 128*H00m1*(-8 + 8*z + 3*z2) + 128*H0*H00m1*(1 + 4*z2) - 128*H0*H011*(3 - 6*z + 4*z2) +
               64*H0*H001*(3 + 10*z + 4*z2) + 48*H13*(2 - 8*z + 7*z2) + 64*H0*H0m1*(-2 + 18*z + 7*z2) + 16*H01*H02*(-1 - 14*z + 8*z2) + 16*H012*(3 - 6*z + 8*z2) -
               32*H02*H0m1*(3 - 2*z + 10*z2) + 64*H000m1*(5 - 2*z + 10*z2) + 64*H0*H12*(4 - 13*z + 11*z2) - 16*H001*(26 - 77*z + 26*z2) + 16*H0011*(17 - 34*z + 32*z2) +
               16*H01*(-13 - 6*z + 37*z2) + 16*H011*(-15 - 3*z + 38*z2) + (2*H04*(5 + 6*z + 40*z2))/3. - 16*H0001*(27 + 34*z + 48*z2) + 8*H02*H1*(-3 - 40*z + 50*z2) +
               16*H0*H1*(29 - 95*z + 55*z2) + 2*H1*(-85 - 216*z + 272*z2) + 2*H12*(155 - 546*z + 368*z2) - (2*H02*(-585 - 5510*z - 4620*z2 + 576*z3))/15. -
               (2*(-32 + 383*z - 1893*z2 + 1092*z3))/(15.*z) + (4*H0*(-16 + 308*z + 1593*z2 + 1464*z3))/(15.*z) +
               (((-32*H0m1)/15. + (32*H0*Hm1)/15.)*(2 - 915*z2 - 1090*z3 - 165*z4 + 72*z5))/z2 + 64*H01*(3 + 2*z)*zeta2 + 32*Hm1*(1 + z)*(18 + 7*z)*zeta2 +
               32*H0m1*(5 - 6*z + 12*z2)*zeta2 - 8*H02*(15 - 22*z + 72*z2)*zeta2 - 16*H1*(5 - 78*z + 86*z2)*zeta2 - 16*H0*(9 - 31*z + 100*z2)*zeta2 +
               (16*(-240 - 665*z - 1380*z2 + 144*z3)*zeta2)/15. + LM2*
                (gammaqg*(-8*H0*H1 - 8*H12) - 24*H0*z - 2*(-11 + 2*z) - 16*H1*(-1 + 4*z) + 64*H01*z2 + 4*H02*(1 - 2*z + 8*z2) - 32*(1 - 2*z + 4*z2)*zeta2) +
               (8*(57 - 346*z + 192*z2)*zeta22)/5. - 32*H0*(9 - 30*z + 44*z2)*zeta3 - 16*(23 - 54*z + 128*z2)*zeta3 +
               gammaqg*(64*H01*H0m1 + (80*H001 - 128*H00m1)*H1 - 8*H03*H1 - 8*H01*H12 - 20*H02*H12 +
                  H0*(-76*H01 - 64*H01m1 - 64*H0m11 + (-64*H01 + 64*H0m1)*H1 - 24*H13) - 8*H14 + (128*H0*H1 + 64*H12)*zeta2 + 120*H1*zeta3) +
               (1 + 2*z + 2*z2)*(64*H001m1 + 64*H00m11 + 64*H00m1m1 + 128*H01m1m1 + 128*H0m11m1 + 128*H0m1m11 + 256*H0m1m1m1 +
                  (-64*H001 - 64*H00m1 - 128*H01m1 - (64*H03)/3. - 192*H0*H0m1 - 128*H0m11 - 256*H0m1m1)*Hm1 + (64*H01 + 112*H02 + 128*H0m1)*Hm12 - (128*H0*Hm13)/3. +
                  (64*H0*Hm1 - 128*Hm12)*zeta2 + 224*Hm1*zeta3) + LM*
                ((64*H0m1 - 64*H0*Hm1)*(1 + z)*(2 + 3*z) + 64*H12*(-1 + z)*(-2 + 5*z) + 8*H02*z*(9 + 32*z) - 128*H0*H0m1*z2 + 64*H011*(-1 + 2*z + 2*z2) +
                  64*H0*H01*(1 - 2*z + 6*z2) + 64*H00m1*(1 + 2*z + 6*z2) + 32*H0*H1*(5 - 18*z + 10*z2) + 8*H0*(11 - 6*z + 18*z2) + 16*H01*(-4 + 5*z + 20*z2) -
                  16*H001*(7 - 14*z + 24*z2) + (8*H03*(3 + 2*z + 24*z2))/3. + 4*H1*(65 - 126*z + 36*z2) - 2*(-33 - 92*z + 104*z2) - 16*(6 - 11*z + 40*z2)*zeta2 -
                  16*H0*(9 - 10*z + 40*z2)*zeta2 + gammaqg*(-16*H02*H1 - 40*H0*H12 - 16*H13 + 64*H1*zeta2) +
                  (1 + 2*z + 2*z2)*(128*H0m1m1 + (-32*H02 - 128*H0m1)*Hm1 + 64*H0*Hm12 + 64*Hm1*zeta2) - 16*(7 + 2*z + 32*z2)*zeta3)) - 8*(1 - 2*z + 56*z2)*zeta5) +
         CA*CF*TR*((-2*H05)/15. - 64*H0m12*(1 + z) + (16*H01m1 + 16*H0m11 - 16*H01*Hm1)*(-2 + z)*(1 + z) +
            (-64*H011m1 - 64*H01m11 - 64*H0m111 + 64*H011*Hm1)*z*(1 + z) + (-24*H00m1m1 + 24*H00m1*Hm1)*(1 + 2*z)*(3 + 2*z) + 8*H0m101*(-1 + 2*z)*(5 + 2*z) +
            (16*H03*Hm1*(1 + z)*(2 + 3*z))/3. - 12*H02*H0m1*(3 + 4*z) + 48*H001*H01*(5 + 6*z) + (16*H01*H1*(-1 + z)*(1 + 15*z))/3. - 8*H000m1*(11 + 20*z) -
            8*H00m1*(1 + 28*z) + (64*H00m101 + 64*H0*H0m101)*z2 + (-192*H0m1m1m1 + 192*H0m1m1*Hm1 - 96*H0m1*Hm12 + 32*H0*Hm13)*(1 + z + z2) +
            16*H001*H02*(-9 - 28*z + 2*z2) - 16*H0*H00m1*(-5 - 6*z + 2*z2) + 8*H0*H0m1*(-3 + 14*z + 2*z2) + (16*H01*H03*(7 + 16*z + 2*z2))/3. +
            (16*H0m1m1 - 16*H0m1*Hm1 + 8*H0*Hm12)*(-5 + z + 3*z2) - 4*H02*Hm1*(-7 + 4*z2) + 8*H0*H0m1*Hm1*(5 + 4*z2) - 8*H0*H0m1m1*(-11 - 16*z + 4*z2) +
            (16*H0*H000m1 - 80*H000m1m1 + 16*H0010m1 - 8*H00m1*H02 - 16*H001*H0m1)*(-1 - 2*z + 4*z2) + (32*H0*H00m1m1 - 8*H02*H0m1m1)*(1 + 2*z + 4*z2) +
            (H0*(-16*H01m1 - 16*H0m11) + 16*H0*H01*Hm1)*(1 + 8*z + 4*z2) + (8*H001m1 + 8*H00m11 - 8*H001*Hm1)*(-3 + 16*z + 4*z2) + 64*H0m1m101*(3 + 6*z + 5*z2) +
            16*H0*H012*(-1 - 10*z + 6*z2) + 16*H0*H0011*(1 - 22*z + 8*z2) +
            (-16*H01m1m1 - 16*H0m11m1 - 16*H0m1m11 + (16*H01m1 + 16*H0m11)*Hm1 - 8*H01*Hm12)*(7 + 12*z + 8*z2) - 16*H00101*(37 + 46*z + 8*z2) -
            16*H00011*(75 + 86*z + 8*z2) - 32*H0*H0001*(-7 - 34*z + 10*z2) - 16*H011*H02*(1 - 26*z + 14*z2) - 32*H01*H011*(9 - 6*z + 14*z2) +
            (16*H0001m1 + 16*H000m11)*(-1 - 2*z + 16*z2) + (H04*(13 + 178*z + 16*z2))/18. + 32*H00001*(-3 - 36*z + 20*z2) - 2*H02*Hm12*(19 + 24*z + 20*z2) +
            16*H01111*(15 - 18*z + 32*z2) + (8*H0m1 - 8*H0*Hm1)*(21 + 40*z + 35*z2) + 16*H0*H0111*(23 - 22*z + 40*z2) + 16*H01011*(41 - 10*z + 56*z2) +
            16*H00111*(79 - 50*z + 120*z2) - (8*H0*H011*(-49 - 328*z + 448*z2))/3. + (8*H0011*(-183 - 240*z + 512*z2))/3. - (H03*(1131 + 1752*z + 7624*z2))/27. -
            (2*H0*(17499 - 34311*z + 9514*z2))/81. + (H02*(-346 - 2987*z + 10622*z2))/9. - (2*H14*(-4 - 56*z + 31*z2 + 12*z3))/(9.*z) -
            (8*H01*H02*(20 - 23*z + 46*z2 + 13*z3))/(3.*z) + (32*H0*H01*H1*(-4 + z - 56*z2 + 54*z3))/(3.*z) - (4*H13*(4 + 24*z - 94*z2 + 57*z3))/(3.*z) +
            (8*H01*H12*(-8 + 5*z - 70*z2 + 84*z3))/(3.*z) - (4*H0*H13*(-16 - 3*z - 48*z2 + 94*z3))/(9.*z) + (2*H02*H1*(-40 + 1269*z - 1392*z2 + 124*z3))/(9.*z) -
            (4*H02*H12*(-16 - 138*z2 + 151*z3))/(3.*z) + (16*H0111*(-32 + 17*z - 139*z2 + 251*z3))/(3.*z) - (8*H011*H1*(-32 + 23*z - 292*z2 + 336*z3))/(3.*z) +
            (8*H0*H001*(40 - 133*z - 88*z2 + 370*z3))/(3.*z) + (4*H012*(-16 + 35*z - 256*z2 + 384*z3))/(3.*z) - (4*H03*H1*(-40 - 25*z - 334*z2 + 404*z3))/(9.*z) -
            (8*H001*H1*(-32 + 5*z - 436*z2 + 432*z3))/(3.*z) + (8*H001*(-20 + 801*z + 2052*z2 + 568*z3))/(9.*z) - (8*H0001*(40 - 233*z - 248*z2 + 590*z3))/(3.*z) +
            (4*H0*H12*(-8 + 303*z - 1068*z2 + 728*z3))/(9.*z) - (8*H0*H01*(-20 + 705*z + 387*z2 + 782*z3))/(9.*z) - (2*H12*(142 + 2241*z - 3096*z2 + 2369*z3))/(27.*z) -
            (8*H011*(-8 + 459*z - 2553*z2 + 2599*z3))/(9.*z) - (4*H01*(740 - 1215*z - 15039*z2 + 9292*z3))/(27.*z) +
            (4*H0*H1*(740 - 1401*z - 12504*z2 + 12769*z3))/(27.*z) - (4*H1*(-365 + 11274*z - 44025*z2 + 30389*z3))/(81.*z) -
            (8388 - 2299*z - 98666*z2 + 54588*z3)/(81.*z) - 48*gammaqg*ln2*zeta2 + 80*H0*H0m1*(1 + z2 + 2*z)*zeta2 + 8*H0m1*(1 + z)*(14 + z)*zeta2 +
            8*H00m1*(-7 - 14*z + 2*z2)*zeta2 - 16*H0m1m1*(3 + 6*z + 2*z2)*zeta2 + (8*H03*(4 + 5*z + 4*z2)*zeta2)/3. + 8*Hm1*(-9 - z + 5*z2)*zeta2 +
            64*H011*(5 - 7*z + 8*z2)*zeta2 + 8*Hm12*(13 + 18*z + 14*z2)*zeta2 - 8*H0*Hm1*(13 + 27*z + 17*z2)*zeta2 - 8*H0*H01*(29 + 14*z + 26*z2)*zeta2 +
            8*H001*(47 + 50*z + 26*z2)*zeta2 - (H02*(-49 - 106*z + 128*z2)*zeta2)/3. + (2*H0*(-347 - 1409*z + 1180*z2)*zeta2)/9. -
            (4*H01*(-48 + 169*z - 338*z2 + 150*z3)*zeta2)/(3.*z) - (4*H12*(-8 + 57*z - 306*z2 + 290*z3)*zeta2)/(3.*z) +
            (4*H0*H1*(-48 - 47*z - 242*z2 + 314*z3)*zeta2)/(3.*z) - (4*H1*(-72 + 593*z - 3646*z2 + 2752*z3)*zeta2)/(9.*z) -
            ((1440 - 1283*z + 11030*z2 + 11000*z3)*zeta2)/(18.*z) +
            LM3*(gammaqg*((-16*H0*H1)/3. - (16*H12)/3.) - 32*H01*(1 + 2*z) - (16*H02*(1 + 4*z))/3. + (2*(-31 - 20*z + 40*z2))/3. + (4*H0*(11 - 94*z + 292*z2))/9. +
               (8*H1*(-16 + 5*z - 142*z2 + 146*z3))/(9.*z) - (32*(-1 - 10*z + 4*z2)*zeta2)/3.) +
            LQ3*(gammaqg*((-8*H0*H1)/3. - (8*H12)/3.) - 16*H01*(1 + 2*z) - (8*H02*(1 + 4*z))/3. + (4*(-5 - 16*z + 10*z2))/3. + (8*H0*(11 - 40*z + 106*z2))/9. +
               (8*H1*(-8 + 19*z - 104*z2 + 106*z3))/(9.*z) - (16*(-1 - 10*z + 4*z2)*zeta2)/3.) - (16*H0*(11 + 103*z)*zeta22)/5. +
            (4*(101 - 3310*z + 2354*z2)*zeta22)/15. + (16*H02*(10 + 13*z + 6*z2)*zeta3)/3. - 8*H01*(33 - 34*z + 48*z2)*zeta3 - 4*Hm1*(43 + 56*z + 52*z2)*zeta3 +
            (2*(33 - 4980*z + 1360*z2)*zeta3)/3. - (4*H0*(131 + 134*z + 1828*z2)*zeta3)/9. + (4*H1*(-64 + 65*z - 1036*z2 + 1088*z3)*zeta3)/(9.*z) +
            (4*(-125 - 638*z + 236*z2)*zeta2*zeta3)/3. + (1 + 2*z)*(32*H0000m1 + 80*H00m10m1 - 48*H00m1*H0m1 - (16*H03*H0m1)/3. + 8*H0*H0m12 + 8*H0m1*zeta3) +
            gammaqg*((16*H0001 - 24*H0011 - 28*H012)*H1 + (H04*H1)/3. + (52*H001 + 20*H011)*H12 + 2*H03*H12 + H02*(-16*H01*H1 + (2*H13)/3.) - (16*H01*H13)/3. +
               H0*((24*H001 + 88*H011)*H1 - 32*H01*H12 - (5*H14)/3.) + (32*H01*H1 - H02*H1 + 20*H0*H12 + 16*H13)*zeta2 - (434*H1*zeta22)/5. +
               ((-116*H0*H1)/3. - (14*H12)/3.)*zeta3) + (1 + 2*z + 2*z2)*
             (96*H0011m1 + 96*H001m11 + 224*H001m1m1 + 96*H00m111 + 224*H00m11m1 + 224*H00m1m11 - 96*H00m1m1m1 - 64*H011m1m1 - 64*H01m11m1 - 64*H01m1m11 + 32*H0m1011 +
               160*H0m101m1 + 160*H0m10m11 + H02*(-32*H01m1 - 32*H0m11) - 64*H0m111m1 - 64*H0m11m11 - 64*H0m1m111 +
               H0*(32*H001m1 + 32*H00m11 - 64*H011m1 - 64*H01m11 - 128*H01m1m1 - 64*H0m111 - 128*H0m11m1 - 128*H0m1m11 + 96*H0m1m1m1) + 192*H0m1m1m1m1 +
               (-32*H0001 + 160*H000m1 - 96*H0011 - 224*H001m1 - 224*H00m11 + 96*H00m1m1 + 64*H011m1 + 64*H01m11 + (8*H04)/3. + H02*(32*H01 + 16*H0m1) - 160*H0m101 +
                  64*H0m111 + H0*(-32*H001 - 64*H00m1 + 64*H011 + 128*H01m1 + 128*H0m11 - 96*H0m1m1) - 192*H0m1m1m1)*Hm1 +
               (112*H001 - 48*H00m1 - 32*H011 - (16*H03)/3. + H0*(-64*H01 + 48*H0m1) + 96*H0m1m1)*Hm12 + (-8*H02 - 32*H0m1)*Hm13 + 8*H0*Hm14 +
               (96*H01m1 + 96*H0m11 + (-96*H01 - 60*H02 + 16*H0m1)*Hm1 + 24*H0*Hm12 + 16*Hm13)*zeta2 + 72*Hm1*zeta22 + (-16*H0*Hm1 - 56*Hm12)*zeta3) +
            LQ2*(8*H03*(1 + 4*z) + 16*H011*(5 + 14*z) + 32*H0*H01*(3 + 6*z + z2) - 32*H001*(2 + 7*z + z2) - (8*H02*(14 - 49*z + 141*z2))/3. -
               (2*H0*(-109 - 1126*z + 2320*z2))/9. + LM*((22*gammaqg*H1)/3. + (22*(-1 + 4*z))/3. - (44*H0*(1 - 2*z + 4*z2))/3.) -
               (-48 + 193*z - 536*z2 + 52*z3)/(3.*z) - (8*H01*(16 - 31*z + 80*z2 + 66*z3))/(3.*z) - (8*H12*(-8 + 30*z - 123*z2 + 131*z3))/(3.*z) -
               (8*H0*H1*(-16 + 53*z - 232*z2 + 248*z3))/(3.*z) - (4*H1*(32 + 191*z - 1078*z2 + 1052*z3))/(9.*z) + 16*H0*(-3 - 20*z + 8*z2)*zeta2 +
               (16*(11 - 85*z + 157*z2)*zeta2)/3. + gammaqg*(10*H02*H1 + 12*H0*H12 + 4*H13 - 24*H1*zeta2) +
               (1 + 2*z + 2*z2)*(16*H00m1 - 32*H01m1 + 24*H0m1 - 16*H0*H0m1 - 32*H0m11 + (-24*H0 + 32*H01 + 8*H02)*Hm1 - 32*Hm1*zeta2) + 64*z2*zeta3) +
            LM2*((-8*H03*(1 + 4*z))/3. - 16*H011*(7 + 10*z) - 32*H001*(-3 - 13*z + z2) + 32*H0*H01*(-2 - 8*z + z2) + 8*H02*(-1 - 11*z + 9*z2) +
               (2*H0*(-455 - 1850*z + 3584*z2))/9. - (16*H01*(-8 + 5*z + 35*z2 + 33*z3))/(3.*z) + (8*H12*(-8 + 22*z - 161*z2 + 157*z3))/(3.*z) +
               (8*H0*H1*(-16 + 9*z - 192*z2 + 208*z3))/(3.*z) + (8*H1*(-64 - 28*z - 1051*z2 + 1022*z3))/(9.*z) + (-48 - 477*z - 868*z2 + 1556*z3)/(3.*z) -
               16*H0*(1 - 12*z + 8*z2)*zeta2 - (8*(-1 - 244*z + 142*z2)*zeta2)/3. + gammaqg*(-6*H02*H1 - 20*H0*H12 - 12*H13 + 40*H1*zeta2) +
               (1 + 2*z + 2*z2)*(16*H00m1 - 32*H01m1 + 24*H0m1 - 16*H0*H0m1 - 32*H0m11 + (-24*H0 + 32*H01 + 8*H02)*Hm1 - 32*Hm1*zeta2) - 96*(1 + 2*z2)*zeta3) +
            LM*((-96*H0m1m1 + 96*H0m1*Hm1 - 48*H0*Hm12)*(1 + z)*(-5 + 2*z) + (4*H04*(3 + 10*z))/3. - 96*H012*(1 + z2) + H0*(64*H01m1 + 64*H0m11)*(1 - 10*z + 2*z2) -
               64*H0*H011*(-1 - 10*z + 3*z2) - 32*H01*H02*(1 - 2*z + 4*z2) + 16*H0m12*(-3 + 10*z + 4*z2) - 64*H0*H00m1*(2 + 5*z2) + 32*H0011*(1 - 10*z + 6*z2) +
               4*H02*Hm1*(-27 - 16*z + 8*z2) - 32*H0111*(7 - 14*z + 8*z2) - 32*H0*H0m1m1*(-1 + 14*z + 8*z2) - 8*H02*H0m1*(9 + 10*z + 12*z2) +
               8*H0*H0m1*(-61 + 52*z + 16*z2) + 32*H0*H001*(11 + 6*z + 18*z2) + (-16*H01m1 - 16*H0m11 + 16*H01*Hm1)*(17 + 34*z + 20*z2) + 8*H13*(4 - 40*z + 39*z2) -
               8*H00m1*(-149 + 88*z + 40*z2) - 16*H0001*(41 + 42*z + 52*z2) - (16*H011*(-10 + 209*z + 53*z2))/3. + 16*H000m1*(37 + 2*z + 68*z2) -
               (8*H03*(34 - 56*z + 293*z2))/9. - (16*H0*H12*(-8 - 5*z - 5*z2 + 13*z3))/(3.*z) + (32*H01*H1*(-4 - 7*z - 19*z2 + 20*z3))/(3.*z) -
               (8*H02*H1*(-8 + 113*z - 232*z2 + 159*z3))/(3.*z) + (8*H001*(16 - 167*z + 460*z2 + 202*z3))/(3.*z) - (8*H0*H01*(16 - 193*z + 572*z2 + 218*z3))/(3.*z) +
               (4*H01*(368 + 65*z - 3316*z2 + 252*z3))/(9.*z) - (4*H02*(385 - 5075*z + 8165*z2 + 432*z3))/45. + (8*H0*H1*(-184 - 73*z - 676*z2 + 461*z3))/(9.*z) +
               (2*H12*(-272 - 77*z - 6188*z2 + 5614*z3))/(9.*z) + (2*H1*(-344 - 1243*z - 14194*z2 + 14956*z3))/(9.*z) +
               (4*H0*(-24 - 3393*z - 15808*z2 + 27226*z3))/(45.*z) + (2*(-552 - 25297*z - 22043*z2 + 53302*z3))/(45.*z) +
               (((-16*H0m1)/15. + (16*H0*Hm1)/15.)*(2 - 1215*z2 - 1105*z3 + 15*z4 + 72*z5))/z2 - 8*H02*(11 + 14*z + 4*z2)*zeta2 + 16*H0m1*(13 + 26*z + 16*z2)*zeta2 +
               16*H01*(11 - 30*z + 20*z2)*zeta2 - 16*Hm1*(2 + 25*z + 26*z2)*zeta2 + (8*H0*(7 + 328*z + 552*z2)*zeta2)/3. -
               (16*H1*(8 - 49*z - z2 + 32*z3)*zeta2)/(3.*z) + (4*(405 + 10080*z - 5870*z2 + 864*z3)*zeta2)/45. - (8*(-79 - 162*z + 12*z2)*zeta22)/5. -
               64*H0*(-4 - 29*z + 2*z2)*zeta3 - (16*(-24 + 26*z - 490*z2 + 501*z3)*zeta3)/(3.*z) +
               gammaqg*(48*H01*H0m1 - 64*H0m101 + (-16*H001 - 96*H00m1)*H1 + (4*H03*H1)/3. - 16*H01*H12 + 6*H02*H12 + H0*((-24*H01 + 48*H0m1)*H1 - (16*H13)/3.) -
                  8*H14 + (80*H0*H1 + 48*H12)*zeta2 + 196*H1*zeta3) +
               (1 + 2*z + 2*z2)*(-64*H00m1m1 - 128*H011m1 - 128*H01m11 - 256*H01m1m1 - 128*H0m111 - 256*H0m11m1 - 256*H0m1m11 - 256*H0m1m1m1 +
                  (64*H00m1 + 128*H011 + 256*H01m1 + (112*H03)/3. + H0*(128*H01 + 64*H0m1) + 256*H0m11 + 256*H0m1m1)*Hm1 + (-128*H01 - 48*H02 - 128*H0m1)*Hm12 +
                  (128*H0*Hm13)/3. + (-416*H0*Hm1 + 192*Hm12)*zeta2 - 368*Hm1*zeta3)) +
            LQ*(gammaqg*(-24*H01*H02 - 48*H01*H0m1 + 96*H00m1*H1 + 8*H01*H12 + H0*(-48*H0m1*H1 - (40*H13)/3.)) + (128*H001m1 + 128*H00m11 - 128*H001*Hm1)*(1 + z2 + 2*z) -
               (4*H04*(5 + 18*z))/3. + 128*H0111*(-1 - 7*z + z2) + H0*(-64*H01m1 - 64*H0m11)*(1 - 10*z + 2*z2) + 32*H0*H01*H1*(3 - 6*z + 2*z2) +
               (128*H011m1 + 128*H01m11 + 128*H0m111 + (-128*H0*H01 - 128*H011)*Hm1)*(1 + 2*z + 2*z2) - 32*H0*H011*(13 + 22*z + 2*z2) +
               (256*H01m1m1 + 256*H0m11m1 + 256*H0m1m11 + (-256*H01m1 - 256*H0m11)*Hm1 + 128*H01*Hm12)*(2 + 4*z + 3*z2) + 64*H0011*(6 + 13*z + 3*z2) -
               128*H001*H1*(3 - 6*z + 4*z2) + 8*H02*H0m1*(9 - 6*z + 4*z2) +
               (64*H00m1m1 + 256*H0m1m1m1 + (-64*H00m1 - 64*H0*H0m1 - 256*H0m1m1)*Hm1 + (48*H02 + 128*H0m1)*Hm12 - (128*H0*Hm13)/3.)*(3 + 6*z + 4*z2) -
               16*H0m12*(-11 + 26*z + 4*z2) + 32*H012*(4 - 8*z + 5*z2) - 128*H0m101*(3 - 6*z + 5*z2) + (16*H03*H1*(1 - 2*z + 6*z2))/3. +
               64*H0*H00m1*(8 + 12*z + 11*z2) + 32*H0*H0m1m1*(-5 + 38*z + 12*z2) + 8*H02*H12*(7 - 14*z + 18*z2) - 32*H0*H001*(19 + 14*z + 18*z2) -
               (16*H03*Hm1*(11 + 22*z + 18*z2))/3. + 16*H0001*(61 + 98*z + 36*z2) + (16*H01*H1*(19 - 32*z + 44*z2))/3. - 16*H000m1*(101 + 82*z + 108*z2) +
               (8*H011*(64 - 107*z + 304*z2 + 34*z3))/(3.*z) + (8*H13*(-8 + 26*z - 82*z2 + 93*z3))/(3.*z) + (8*H03*(225 - 610*z + 2305*z2 + 144*z3))/45. +
               (4*H0*(64 + 3445*z + 11758*z2 + 242*z3))/(45.*z) + (8*H02*H1*(-16 + 161*z - 364*z2 + 303*z3))/(3.*z) + (8*H0*H12*(-32 + 69*z - 366*z2 + 398*z3))/(3.*z) +
               (8*H12*(32 + 65*z - 346*z2 + 518*z3))/(9.*z) + (8*H0*H01*(32 - 257*z + 172*z2 + 518*z3))/(3.*z) + (8*H0*H1*(64 + 223*z - 1520*z2 + 1834*z3))/(9.*z) +
               (4*H01*(-592 - 499*z + 16856*z2 + 11628*z3))/(45.*z) - (2*H1*(-544 - 10349*z - 22856*z2 + 32084*z3))/(45.*z) -
               (2*(-2512 - 18123*z - 39737*z2 + 67432*z3))/(45.*z) + (8*H001*(-160 + 1215*z + 260*z2 - 2150*z3 + 288*z4))/(15.*z) -
               (8*H0*H0m1*(-8 - 80*z - 1655*z2 - 300*z3 + 1260*z4))/(15.*z2) + (4*H02*(24 - 602*z - 10977*z2 + 16274*z3 + 2016*z4))/(45.*z) +
               (((32*H0m1m1)/15. - (32*H0m1*Hm1)/15. + (16*H0*Hm12)/15.)*(2 + 40*z - 55*z2 + 405*z3 + 600*z4 + 72*z5))/z2 +
               (((16*H01m1)/15. + (16*H0m11)/15. - (16*H01*Hm1)/15.)*(4 + 615*z2 + 670*z3 + 300*z4 + 144*z5))/z2 -
               (4*H02*Hm1*(8 + 80*z + 295*z2 + 1000*z3 + 1140*z4 + 288*z5))/(15.*z2) + (8*H00m1*(-8 - 80*z - 3015*z2 + 400*z3 + 3660*z4 + 288*z5))/(15.*z2) +
               (((16*H0m1)/45. - (16*H0*Hm1)/45.)*(28 + 12*z - 4851*z2 - 6071*z3 - 558*z4 + 1008*z5))/z2 + 256*H0*H1*(1 + z2 - 2*z)*zeta2 - 64*H12*(-1 + 2*z)*zeta2 -
               24*H02*(-5 - 18*z + 4*z2)*zeta2 - 16*H0m1*(13 + 26*z + 8*z2)*zeta2 - 64*Hm12*(7 + 14*z + 10*z2)*zeta2 - 16*H01*(13 - 74*z + 28*z2)*zeta2 +
               32*H0*Hm1*(21 + 42*z + 34*z2)*zeta2 - (8*H0*(85 - 320*z + 2020*z2 + 192*z3)*zeta2)/5. - (4*(96 + 1731*z - 22628*z2 + 29968*z3 + 4032*z4)*zeta2)/(45.*z) -
               (16*H1*(-2 - 120*z + 495*z2 - 1585*z3 + 1610*z4 + 72*z5)*zeta2)/(15.*z2) + (16*Hm1*(6 + 40*z + 560*z2 + 1075*z3 + 900*z4 + 216*z5)*zeta2)/(15.*z2) +
               LM2*(gammaqg*(8*H0*H1 + 8*H12) + 48*H01*(1 + 2*z) + 8*H02*(1 + 4*z) - (8*(-1 + z)*(13 + 15*z))/3. - (16*H0*z*(-9 + 31*z))/3. -
                  (8*H1*(-8 - 3*z - 60*z2 + 62*z3))/(3.*z) + 16*(-1 - 10*z + 4*z2)*zeta2) + (8*(39 - 234*z + 188*z2)*zeta22)/5. - 32*H0*(3 + 58*z + 4*z2)*zeta3 +
               112*H1*(7 - 14*z + 10*z2)*zeta3 + 16*Hm1*(51 + 102*z + 74*z2)*zeta3 - (8*(32 + 49*z + 646*z2 - 6*z3 + 144*z4)*zeta3)/(3.*z) +
               LM*(-32*H011*(-1 + 2*z) - (16*H03*(1 + 4*z))/3. + 32*H001*(-1 - 6*z + 2*z2) + (32*H0*H1*(11 - 10*z + 10*z2))/3. - (16*H12*(-4 - 19*z + 13*z2))/3. +
                  8*H01*(-7 + 50*z + 44*z2) + (8*H02*(17 - 16*z + 114*z2))/3. - (4*H0*(-173 - 362*z + 632*z2))/9. - (2*(-335 - 166*z + 752*z2))/3. -
                  (4*H1*(-160 - 247*z - 1024*z2 + 992*z3))/(9.*z) + 64*H0*(1 + 2*z)*zeta2 - (8*(23 + 74*z + 172*z2)*zeta2)/3. +
                  gammaqg*(-4*H02*H1 + H0*(8*H01 + 8*H12) + 8*H13 - 16*H1*zeta2) +
                  (1 + 2*z + 2*z2)*(-32*H00m1 + 64*H01m1 - 48*H0m1 + 32*H0*H0m1 + 64*H0m11 + (48*H0 - 64*H01 - 16*H02)*Hm1 + 64*Hm1*zeta2) + 32*(3 + 4*z2)*zeta3)) -
            4*(1 - 174*z + 212*z2)*zeta5) + CF*TR*TR*((-64*H01*H02*(1 + z2 - 2*z))/3. + (32*H13*(-1 + z)*(1 + 3*z))/9. - (128*H0*H001*(-1 + 2*z + z2))/3. +
            (64*H011*(-1 - 6*z + 3*z2))/3. + (4*H04*(1 - 20*z + 4*z2))/9. - (64*H0*H01*(2 - 4*z + 5*z2))/3. + (16*H03*(4 - 63*z + 10*z2))/9. +
            (16*H02*H1*(1 - 12*z + 10*z2))/3. + (16*H12*(4 - 17*z + 12*z2))/3. - (32*H0*H1*(10 - 27*z + 16*z2))/3. + (16*H001*(56 - 113*z + 20*z2))/3. +
            (16*H01*(140 - 113*z + 20*z2))/3. + (16*H0001*(11 - 28*z + 40*z2))/3. - (8*H02*(-22 + 75*z + 68*z2))/3. + (4*H0*(21 - 1578*z + 538*z2))/3. +
            (4*H1*(16 + 793*z - 1558*z2 + 810*z3))/(3.*z) + (2*(-208 + 939*z - 8610*z2 + 6994*z3))/(9.*z) +
            LQ3*((16*gammaqg*H1)/9. + (16*H02*(-1 + 2*z))/3. + (16*H0*(-11 + 4*z))/9. - (16*(-8 + 84*z - 147*z2 + 62*z3))/(27.*z)) +
            LM3*((32*gammaqg*H1)/9. - (16*H02*(-1 + 2*z))/3. - (16*H0*(-5 - 8*z + 24*z2))/9. + (16*(-8 + 75*z - 111*z2 + 62*z3))/(27.*z)) +
            (80*H1*(-1 - 4*z + 4*z2)*zeta2)/9. + (4*H0*(-1063 + 1124*z + 24*z2)*zeta2)/9. + (4*H02*(-59 + 112*z + 32*z2)*zeta2)/3. +
            (2*(-256 - 12243*z + 7716*z2 + 1984*z3)*zeta2)/(27.*z) - (16*(91 - 194*z + 218*z2)*zeta22)/15. - (32*H0*(37 - 74*z + 36*z2)*zeta3)/9. -
            (16*(-8 + 525*z - 1128*z2 + 170*z3)*zeta3)/(27.*z) + gammaqg*
             ((80*H0111)/3. + (32*H012)/3. + ((-64*H001)/3. - (64*H011)/3.)*H1 - (8*H03*H1)/9. + H0*((-64*H011)/3. + (32*H01*H1)/3.) + (16*H01*H12)/3. - (2*H14)/9. +
               (-12*H01 - (16*H0*H1)/3. - (26*H12)/3.)*zeta2 + (160*H1*zeta3)/9.) +
            LQ2*(gammaqg*((-32*H0*H1)/3. - 4*H12) - (8*H02*(-31 + 56*z + 8*z2))/3. - (16*H01*(-5 - 8*z + 12*z2))/3. + (8*H0*(571 - 236*z + 244*z2))/9. +
               LM*(-16*H02*(-1 + 2*z) - (16*H0*(-9 + 8*z2))/3. + (16*(-1 + z)*(8 - 73*z + 62*z2))/(9.*z)) + (16*H1*(-8 + 112*z - 215*z2 + 130*z3))/(9.*z) +
               (4*(256 + 4977*z - 7632*z2 + 2156*z3))/(27.*z) - (16*(13 - 8*z + 4*z2)*zeta2)/3. + (-1 + 2*z)*(-32*H001 - 16*H03 + 32*H0*zeta2 + 32*zeta3)) +
            LM2*(gammaqg*((32*H0*H1)/3. + (20*H12)/3.) - 16*H01*(-3 + 4*z2) + (8*H0*(547 - 92*z + 4*z2))/9. - (8*H02*(-23 + 40*z + 40*z2))/3. +
               (16*H1*(-8 + 70*z - 71*z2 + 10*z3))/(9.*z) + (4*(256 + 4221*z - 6192*z2 + 1652*z3))/(27.*z) + (16*(-1 - 16*z + 28*z2)*zeta2)/3. +
               (-1 + 2*z)*(-32*H001 - 16*H03 + 32*H0*zeta2 + 32*zeta3)) +
            LQ*(gammaqg*((32*H01*H1)/3. + 8*H0*H12 + (8*H13)/3.) + (128*H0*H0m1*(-2 + z)*(-2 + 3*z))/3. + (64*H03*(-5 + 15*z + 2*z2))/3. -
               (64*H02*H1*(2 - 4*z + 5*z2))/3. - (128*H00m1*(9 - 14*z + 7*z2))/3. - (32*H001*(15 - 60*z + 8*z2))/3. + (32*H0*H01*(-5 - 8*z + 16*z2))/3. +
               (32*H011*(1 - 20*z + 24*z2))/3. - (16*H02*(2975 - 2660*z + 560*z2 + 72*z3))/45. + (16*H01*(-16 - 329*z - 254*z2 + 76*z3))/(9.*z) -
               (16*H12*(-8 + 131*z - 253*z2 + 168*z3))/(9.*z) - (32*H0*H1*(-8 + 140*z - 283*z2 + 198*z3))/(9.*z) - (16*H0*(4 + 9213*z - 5632*z2 + 1444*z3))/(45.*z) -
               (8*H1*(256 + 5247*z - 8448*z2 + 3080*z3))/(27.*z) - (4*(-11744 + 567861*z - 828996*z2 + 269954*z3))/(405.*z) +
               LM2*((-16*gammaqg*H1)/3. + 16*H02*(-1 + 2*z) + (16*H0*(-7 - 4*z + 16*z2))/3. - (16*(-8 + 78*z - 123*z2 + 62*z3))/(9.*z)) +
               (((-64*H0m1)/45. + (64*H0*Hm1)/45.)*(1 - 20*z + 225*z2 + 40*z3 - 155*z4 + 36*z5))/z2 - (32*H0*(-33 + 84*z + 4*z2)*zeta2)/3. +
               (64*H1*(3 - 6*z + 8*z2)*zeta2)/3. + (16*(-80 + 3045*z - 1400*z2 + 1600*z3 + 144*z4)*zeta2)/(45.*z) +
               (1 + z2 + 2*z)*((-256*H0m1m1)/3. + ((64*H02)/3. + (256*H0m1)/3.)*Hm1 - (128*H0*Hm12)/3. - (128*Hm1*zeta2)/3.) +
               LM*((-8*gammaqg*H12)/3. + 16*H02*(-9 + 16*z + 8*z2) - (16*H0*(559 - 164*z + 124*z2))/9. - (32*H1*(-8 + 91*z - 143*z2 + 70*z3))/(9.*z) -
                  (8*(256 + 4599*z - 6912*z2 + 1904*z3))/(27.*z) + (-7 - 4*z + 12*z2)*((32*H01)/3. - (32*zeta2)/3.) +
                  (-1 + 2*z)*(64*H001 + 32*H03 - 64*H0*zeta2 - 64*zeta3)) + (128*(11 - 21*z + 4*z2)*zeta3)/3. +
               (-1 + 2*z)*(64*H0*H001 + 64*H0011 + 16*H04 - 96*H02*zeta2 - (32*zeta22)/5. - 128*H0*zeta3)) +
            LM*(gammaqg*((8*H0*H12)/3. + (8*H13)/3.) + (128*H00m1*(3 - 10*z + z2))/3. - (128*H0*H0m1*(2 - 4*z + z2))/3. + (64*H02*H1*(2 - 4*z + 3*z2))/3. +
               (32*H001*(23 - 76*z + 8*z2))/3. - 32*H0*H01*(1 - 8*z + 8*z2) - (64*H03*(-14 + 43*z + 10*z2))/9. - (32*H011*(-5 - 8*z + 16*z2))/3. +
               (16*H12*(-8 + 98*z - 121*z2 + 42*z3))/(9.*z) - (16*H02*(-2945 + 2640*z - 170*z2 + 72*z3))/45. + (32*H0*H1*(-8 + 104*z - 199*z2 + 120*z3))/(9.*z) -
               (16*H01*(-16 - 365*z - 326*z2 + 172*z3))/(9.*z) + (32*H0*(-2 + 4406*z - 2669*z2 + 173*z3))/(45.*z) + (8*H1*(256 + 4941*z - 7152*z2 + 2108*z3))/(27.*z) +
               (4*(-11456 + 549969*z - 812364*z2 + 274706*z3))/(405.*z) + (((-64*H0m1)/45. + (64*H0*Hm1)/45.)*(1 + 20*z - 45*z2 + 40*z3 + 155*z4 + 36*z5))/z2 +
               (32*H0*(-25 + 68*z + 28*z2)*zeta2)/3. + (16*(80 - 2865*z + 520*z2 - 340*z3 + 144*z4)*zeta2)/(45.*z) +
               (1 + z2 + 2*z)*((-256*H0m1m1)/3. + ((64*H02)/3. + (256*H0m1)/3.)*Hm1 - (128*H0*Hm12)/3. - (128*Hm1*zeta2)/3.) + (64*(-15 + 44*z + 12*z2)*zeta3)/3. +
               (-1 + 2*z)*(-64*H0*H001 - 64*H0011 - 16*H04 + (96*H02 + (64*H1)/3.)*zeta2 + (32*zeta22)/5. + 128*H0*zeta3)) +
            (-1 + 2*z)*(-32*H00001 + (40*H03*zeta2)/3. + (64*H0*zeta22)/5. + (64*H02*zeta3)/3. + 32*zeta5))
           + a_muindep_->MuIndependentTerm(z, nf + 1)
           + massless_as3_->MuIndependentTerms(z, nf + 1) / (nf + 1.);
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
