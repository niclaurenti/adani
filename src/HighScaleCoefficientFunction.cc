#include "adani/HighScaleCoefficientFunction.h"
#include "adani/Constants.h"
#include "adani/SpecialFunctions.h"

#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

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

    SetFunctions();
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

    higher = fxBand(x, m2Q2, m2mu2, nf).GetHigher()
             - a_muindep_->MuIndependentNfIndependentTerm(x).GetHigher()
             + (a_muindep_->NotOrdered(x))[1];
    lower = fxBand(x, m2Q2, m2mu2, nf).GetLower()
            - a_muindep_->MuIndependentNfIndependentTerm(x).GetLower()
            + (a_muindep_->NotOrdered(x))[2];

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
            cout << "Error: quark coefficient function is not present at O(as)!"
                 << endl;
            exit(-1);
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
            cout << "Error: something has gone wrong in "
                    "HighScaleCoefficientFunction::SetFunctions, GetOrder() == "
                    "2!"
                 << endl;
            exit(-1);
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
            cout << "Error: something has gone wrong in "
                    "HighScaleCoefficientFunction::SetFunctions, GetOrder() == "
                    "3!"
                 << endl;
            exit(-1);
        }
    } else {
        cout << "Error: something has gone wrong in "
                "HighScaleCoefficientFunction::SetFunctions!"
             << endl;
        exit(-1);
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

    double H02 = H0 * H0;
    double H03 = H02 * H0;

    return (CF * TR
            * (8.88888888888889 + 33.77777777777777 * z
               - 135.11111111111111 * z2 + H02 * (40. - 16. * z) * z2
               + 92.44444444444444 * z3
               + (10.666666666666666 * H0m1 - 10.666666666666666 * H0 * Hm1)
                     * (1. + z) * (1. + z) * (1. + z)
               + H0 * z
                     * (93.33333333333331
                        + (-87.99999999999999 - 78.22222222222221 * z) * z)
               + 7.111111111111111 * H1 * (-1. + z)
                     * (3.25 + z * (-6.5 + 1. * z))
               + H1 * (5.333333333333333 * H0 + 1.3333333333333333 * H1)
                     * (1 - z) * (4 + z * (7 + 4 * z))
               - 5.333333333333333 * H01 * (4 + z * (3 + z * (-3 + 2 * z)))
               - 10.666666666666666 * (1 - 3 * (-1 + z) * z2) * zeta2
               + z * (1 + z)
                     * (5.333333333333333 * H03 + 32 * H0 * (H01 - zeta2)
                        + 16 * (-2 * H001 + H011 + zeta3))
               + 8.
                     * (-2.2222222222222223 + (2. + (-1. + H0) * H0) * z
                        + (-5.999999999999999 + (-5. + H0) * H0) * z2
                        + (6.222222222222221 - 2.6666666666666665 * H0) * z3)
                     * L_M
               - 1.3333333333333333
                     * (4 + z * (3 + 6 * H0 * (1 + z) - z * (3 + 4 * z))) * L_M2
               + 0.8888888888888888
                     * (36 * H0 * z3 - 2 * (-1 + z) * (13 - 26 * z + 4 * z2)
                        + 3 * H1 * (-1 + z) * (4 + z * (7 + 4 * z))
                        - 18 * z * (1 + z) * (H02 + H01 - zeta2))
                     * L_Q
               + 1.3333333333333333
                     * (4 + z * (3 + 6 * H0 * (1 + z) - z * (3 + 4 * z)))
                     * L_Q2))
           / z;
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
    double x3 = x2 * x;

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

    return Lmmu
               * (-1.3333333333333333
                  + (10.666666666666666 - 10.666666666666666 * x) * x
                  + H0
                        * (-1.3333333333333333
                           + (2.6666666666666665 - 2.6666666666666665 * x) * x)
                  + H1
                        * (-1.3333333333333333
                           + (2.6666666666666665 - 2.6666666666666665 * x) * x)
                  + LQm
                        * (1.3333333333333333
                           + x * (-2.6666666666666665 + 2.6666666666666665 * x))
               )
           + CA
                 * (-0.6666666666666666 - 9.333333333333332 * H0
                    + 1.9999999999999998 * H00 - 3.9999999999999996 * H000
                    + 7.999999999999999 * H010 - 1.9999999999999998 * H1
                    - 3.9999999999999996 * H101 - 1.9999999999999998 * H11
                    - 3.9999999999999996 * H110 + 3.9999999999999996 * H111
                    + 3.9999999999999996 * Hm100 - 7.999999999999999 * Hm1m10
                    - 8.296296296296296 / x + (5.333333333333333 * H10) / x
                    + H10
                          * (6.
                             + (31.999999999999996 - 43.33333333333333 * x) * x)
                    + LQm2
                          * (2 + 4 * H0 - 4 * H1 + 2.6666666666666665 / x
                             + 16 * x + 8 * (2 * H0 + H1) * x
                             + (-20.666666666666664 - 8. * H1) * x2)
                    - (48. * LQm
                       * (-0.24074074074074078
                          + (0.7638888888888888 + 0.3333333333333333 * H00
                             - 0.16666666666666666 * H10
                             - 0.16666666666666666 * H11
                             + 0.16666666666666666 * Hm10)
                                * x
                          + H1
                                * (0.1111111111111111
                                   + x
                                         * (-0.08333333333333333
                                            + (1.6666666666666667
                                               - 1.861111111111111 * x)
                                                  * x))
                          + x2
                                * (1.2777777777777777 + 2. * H0 + H00 + H01
                                   + 0.3333333333333333 * H10
                                   + 0.3333333333333333 * H11
                                   + 0.3333333333333333 * Hm10
                                   - 0.6666666666666666 * zeta2)
                          + x3
                                * (-1.8842592592592593 - 2.0833333333333335 * H0
                                   - 0.3333333333333333 * H01
                                   - 0.3333333333333333 * H10
                                   - 0.3333333333333333 * H11
                                   + 0.3333333333333333 * Hm10
                                   + 0.3333333333333333 * zeta2)))
                          / x
                    - 3.9999999999999996 * Hm1 * zeta2
                    + x2
                          * (117.62962962962962 - 88.88888888888889 * H0
                             + 15.333333333333332 * H00
                             + 1.9999999999999998 * H01 + 7.999999999999999 * H1
                             - 7.999999999999999 * H101 + 10. * H11
                             - 7.999999999999999 * H110
                             + 7.999999999999999 * H111
                             + 7.999999999999999 * Hm10
                             + 7.999999999999999 * Hm100
                             - 15.999999999999998 * Hm1m10
                             - 1.9999999999999998 * zeta2
                             - 7.999999999999999 * Hm1 * zeta2)
                    + Lmmu
                          * (-28.666666666666664 - 4. * H0 - 8. * H00 + 8. * H10
                             + 16. * H11 + 2.666666666666666 / x
                             - (5.333333333333332 * H1) / x
                             + H1 * (4. + x * (-96. + 105.33333333333331 * x))
                             + LQm
                                   * (4 + 5.333333333333333 / x + 32 * x
                                      - 41.33333333333333 * x2
                                      + 8 * H0 * (1 + 4 * x)
                                      - 8 * H1 * (1 + 2 * (-1 + x) * x))
                             + x2
                                   * (187.33333333333331
                                      + 41.33333333333333 * H0 + 16. * H01
                                      + 16. * H10 + 32. * H11 - 16. * zeta2)
                             + x
                                   * (-161.33333333333331 - 128. * H0
                                      - 32. * H00 - 48. * H01 - 16. * H10
                                      - 32. * H11 + 48. * zeta2))
                    + 15.999999999999998 * zeta3
                    + x
                          * (-104.66666666666666 - 28.666666666666668 * H0
                             + 7.999999999999999 * H00
                             - 7.999999999999999 * H000
                             + 31.999999999999996 * H010
                             - 7.999999999999999 * H1 + 7.999999999999999 * H101
                             - 7.999999999999999 * H11
                             + 7.999999999999999 * H110
                             - 7.999999999999999 * H111
                             + 7.999999999999999 * Hm10
                             + 7.999999999999999 * Hm100
                             - 15.999999999999998 * Hm1m10
                             + 7.999999999999999 * zeta2
                             - 7.999999999999999 * Hm1 * zeta2 + 56. * zeta3))
           + CF
                 * (13 - H00 + 2 * H000 + 4 * H01 - 4 * H010 - 4 * H011
                    + 2 * H10 + 4 * H100 + 4 * H11 - 4 * H111 - 41 * x
                    + LQm2
                          * (-1 - 2 * H0 - 4 * H1 + 4 * (1 + H0 + 2 * H1) * x
                             - 8 * (H0 + H1) * x2)
                    - H0 * (8 + 3 * x * (3 + 8 * x)) - 4 * zeta2
                    + 2 * LQm
                          * (9 + 2 * H0 + 4 * H00 + 6 * H01 + 7 * H1 + 8 * H10
                             + 8 * H11 - 17 * x + 4 * H0 * x * (-3 + 5 * x)
                             - 6 * zeta2
                             + 4 * x
                                   * (-6 * H1 - 4 * H10 - 4 * H11 + x
                                      + 5 * H1 * x + H01 * (-3 + 4 * x)
                                      + H00 * (-2 + 4 * x)
                                      + 4 * x * (H10 + H11 - zeta2) + 3 * zeta2)
                          )
                    - 4 * zeta3
                    + 2 * x
                          * (-2 * H000 + 12 * H01 + 4 * H010 + 4 * H011
                             + 13 * H1 - 12 * H10 - 4 * H100 + 4 * H11
                             + 4 * H111 + 20 * x + 2 * H00 * (-3 + 5 * x)
                             - 12 * zeta2 + 4 * zeta3
                             + 2 * x
                                   * (2 * H000 - 3 * H01 - 2 * H011 - 6 * H1
                                      + 5 * H10 + 2 * H100 - 3 * H11 - 2 * H111
                                      + 3 * zeta2 + 2 * zeta3)))
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
    double z5 = z4 * z;

    double L_M = log(m2mu2);
    double L_Q = log(1. / m2Q2) + L_M;

    double H0 = H_0(z);
    double H1 = H_1(z);
    double Hm1 = H_m1(z);
    double H01 = H_01(z);
    double H0m1 = H_0m1(z);

    double H02 = H0 * H0;
    double H12 = H1 * H1;

    return -21.333333333333332 * TR * TR * (-1. + z) * z * L_M
           + CF * TR
                 * (-17.066666666666666 - 13.866666666666669 * H0 - 16. * H1
                    + (4.266666666666667 * H0 * Hm1) / z2
                    + (4.266666666666667 - 4.266666666666667 * H0) / z
                    + (-121.60000000000002
                       + (-83.2 - 21.333333333333336 * H0) * H0 - 32. * H01)
                          * z
                    + (H0m1
                       * (-4.266666666666667 + 21.333333333333336 * z3
                          - 25.6 * z5))
                          / z2
                    + z
                          * (-21.333333333333336 * H0 * Hm1 + 134.4 * z
                             + 38.4 * H0 * z + H1 * (-48. + 64. * z)
                             + 10.666666666666668 * zeta2
                             + z2
                                   * (-12.8 * H02 + 25.6 * H0 * Hm1
                                      + 25.6 * zeta2))
                    + (-0.5 + z * (-0.5 - 1. * H0 + 1. * z))
                          * (32. * L_M - 32. * L_Q))
           + (CA * TR
              * (-3.5555555555555554 + (10.666666666666666 + 32. * H0) * z
                 + H12 * (32. - 32. * z) * z2
                 + H1
                       * (-10.666666666666666
                          + z
                                * (32.
                                   + (288. + H0 * (64. - 64. * z)
                                      - 309.3333333333333 * z)
                                         * z))
                 + z2
                       * (181.33333333333331 + 128. * H01 - 64. * H0m1
                          + H0 * (256. + 96. * H0 + 64. * Hm1) - 128. * zeta2)
                 + z3
                       * (-188.44444444444443 - 64. * H0m1
                          + H0 * (-416. + 64. * Hm1) + 64. * zeta2)
                 + (10.666666666666666
                    + z
                          * (-31.999999999999993
                             + z
                                   * (-159.99999999999997 - 128. * H0
                                      + 181.33333333333331 * z
                                      + H1 * (-64. + 64. * z))))
                       * L_Q))
                 / z;
}

//==========================================================================================//
//  High scale (Q^2 >> m^2) limit of the quark coefficient functions for FL at
//  O(as^2) expanded in terms of \alpha_s^{[nf+1]}
//------------------------------------------------------------------------------------------//

double HighScaleCoefficientFunction::DL_ps2_highscale(
    double z, double m2Q2, double m2mu2
) const {

    double z2 = z * z;
    double z3 = z2 * z;

    double L_M = log(m2mu2);
    double L_Q = log(1. / m2Q2) + L_M;

    double H0 = H_0(z);
    double H1 = H_1(z);
    double H01 = H_01(z);

    return (CF * TR
            * (-3.5555555555555554 + (10.666666666666666 + 32. * H0) * z
               + H1 * (-10.666666666666666 + 32. * z - 21.333333333333332 * z3)
               + z2
                     * (-42.66666666666667 + 32. * H01
                        + H0 * (-32. + 32. * H0 - 64. * z)
                        + 35.55555555555556 * z - 32. * zeta2)
               + (10.666666666666666
                  + z * (-32. - 32. * H0 * z + 21.333333333333332 * z2))
                     * L_Q))
           / z;
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

    double H02 = H0 * H0;
    double H03 = H02 * H0;
    double H04 = H03 * H0;
    double H12 = H1 * H1;
    double H13 = H12 * H1;
    double Hm12 = Hm1 * Hm1;
    double H0m12 = H0m1 * H0m1;

    return -28.444444444444443 * TR * TR * TR * (-1. + z) * z * L_M2
           + CA * CA * TR
                 * ((167.1111111111111 + 46.22222222222222 / z
                     - 2631.1111111111113 * z + H00m1 * (1664. - 384. * z) * z
                     - 128. * H001 * (-3. + z) * z + 2417.777777777778 * z2
                     - 256. * H011 * z * (5. + z)
                     - 128. * H0 * H01 * z * (11. + z)
                     + H0 * H0m1 * z * (-640. + 384. * z)
                     + H1
                           * (-437.3333333333333 + 330.66666666666663 / z
                              - 6716.444444444443 * z + 6823.11111111111 * z2)
                     + H0
                           * (-213.33333333333331 + 14.222222222222221 / z
                              - 4995.555555555556 * z + 9646.22222222222 * z2)
                     + (1717.3333333333333 * H12 * (-1. + z)
                        * (-0.037267080745341616
                           + z * (0.07453416149068323 + 1. * z)))
                           / z
                     + (2922.6666666666665 * H0 * H1 * (-1. + z)
                        * (-0.043795620437956206
                           + z * (0.08759124087591241 + 1. * z)))
                           / z
                     + ((21.333333333333332 * H0m1
                         - 21.333333333333332 * H0 * Hm1)
                        * (1. + z) * (-4. + z * (-8. + 79. * z)))
                           / z
                     + H02 * (128. + z * (-1941.3333333333333 + 384. * z))
                     - (42.666666666666664 * H01
                        * (1. + z * (-6. + z * (53. + z))))
                           / z
                     + H1 * (-1. + z) * z
                           * (192. * H02 + 384. * H0 * H1 + 128. * H12
                              - 512. * zeta2)
                     + (128. + (3413.333333333333 - 2880. * z) * z) * zeta2
                     + z * (1. + z)
                           * (256. * H01m1 + 256. * H0m11 + 512. * H0m1m1
                              + Hm1
                                    * (-192. * H02 - 256. * H01 - 512. * H0m1
                                       + 256. * H0 * Hm1 + 512. * zeta2))
                     + z
                           * (-341.3333333333333 * H03 + 1792. * H0 * zeta2
                              - 128. * zeta3))
                        * L_Q
                    + ((-115.55555555555554 - 42.666666666666664 * H1
                        + 69.33333333333331 * z + 256. * H02 * z2
                        + H0
                              * (-21.333333333333332
                                 + z
                                       * (-64.
                                          + (1002.6666666666666
                                             + H1 * (128. - 128. * z) - 128. * z
                                            ) * z))
                        + z
                              * (H12 * (128. - 128. * z) * z
                                 + H1
                                       * (128.
                                          + (757.3333333333333
                                             - 842.6666666666666 * z)
                                                * z)
                                 + z
                                       * (1882.6666666666665
                                          - 1836.4444444444443 * z
                                          + H01 * (384. + 128. * z)
                                          - 512. * zeta2)))
                       * L_Q2)
                          / z)
           + CF * CF * TR
                 * (48. + 272. * z + 48. * H12 * (-1. + z) * z - 320. * z2
                    + 192. * H1 * (-1. + z) * (-0.16666666666666666 + 1. * z)
                    - 160. * H0 * H1 * (-1. + z) * (-0.1 + 1. * z)
                    + (32. * H011 - 16. * H02 * H1 + 5.333333333333333 * H13)
                          * (-1. + z) * (1. + 2. * z)
                    + H02 * (4. + (44. - 80. * z) * z)
                    + H001 * (96. + (-544. - 64. * z) * z)
                    - 10.666666666666666 * H03 * (-0.25 + (-1. + z) * z)
                    + H0 * H01 * (-64. + z * (160. + 64. * z))
                    + H0 * (-56. + z * (296. + 192. * z))
                    + H01 * (48. + z * (-384. + 256. * z))
                    - 96. * (-2. + z) * (-0.16666666666666666 + 1. * z) * zeta2
                    + z
                          * (1.3333333333333333 * H04 + 384. * H0001
                             - 64. * H0011 + 32. * H02 * H01 - 64. * H0111
                             - 121.60000000000001 * zeta2_2
                             + H0 * (-192. * H001 + 192. * zeta2 - 64. * zeta3))
                    + (-32. + (608. - 64. * z) * z) * zeta3
                    + (-3.7333333333333485 - 17.066666666666666 / z + 510.4 * z
                       - 489.6 * z2 + 64. * H1 * (-1. + z) * (-0.75 + 1. * z)
                       + H01 * (32. + 192. * z)
                       + H001 * (128. + (320. - 256. * z) * z)
                       + H0
                             * (-72.53333333333333 + 17.066666666666666 / z
                                + 300.8 * z + 166.4 * z2)
                       + H0 * H01 * (-128. + z * (-384. + 256. * z))
                       + H02
                             * (-32.
                                + z
                                      * (-202.66666666666669
                                         + z * (192. + 51.2 * z)))
                       + (102.4 * (1. * H0m1 - 1. * H0 * Hm1) * (1. + z)
                          * (0.16666666666666666
                             + z
                                   * (-0.16666666666666666
                                      + z
                                            * (-6.083333333333333
                                               + z * (0.25 + 1. * z)))))
                             / z2
                       + (96.
                          + z * (533.3333333333334 + (-256. - 102.4 * z) * z))
                             * zeta2
                       + H1 * (-1. + z) * (1. + 2. * z)
                             * ((128. - 64. * H0) * H0 + 64. * H1 + 128. * zeta2
                             )
                       + (1. + z) * (-1. + 2. * z)
                             * (384. * H00m1 - 128. * H0 * H0m1 + 256. * H0m1m1
                                - 64. * H02 * Hm1 - 256. * H0m1 * Hm1
                                + 128. * H0 * Hm12 + 128. * Hm1 * zeta2)
                       + z
                             * (10.666666666666666 * H03 + 768. * H0001
                                - 768. * H000m1 - 256. * H011
                                + H02 * (128. * H01 - 128. * H0m1)
                                - 256. * H0m12 - 256. * H01 * zeta2
                                + 256. * H0m1 * zeta2
                                + 153.60000000000002 * zeta2_2
                                + H0
                                      * (-512. * H001 + 512. * H00m1
                                         + 512. * H0m1m1 + 64. * zeta2
                                         - 512. * zeta3)
                                - 576. * zeta3))
                          * L_M
                    + (-24. - 16. * H02 * z + H0 * (-16. + 64. * z2)
                       + H1 * (-32. + z * (-32. + 64. * z))
                       + z * (24. - 64. * H01 + 64. * zeta2))
                          * L_M2
                    + (3.7333333333333485 + H01 * (-32. - 192. * z)
                       + 17.066666666666666 / z - 510.4 * z + 489.6 * z2
                       - 64. * H1 * (-1. + z) * (-0.75 + 1. * z)
                       + H0 * H01 * (128. + (384. - 256. * z) * z)
                       + H0
                             * (72.53333333333333 - 17.066666666666666 / z
                                - 300.8 * z - 166.4 * z2)
                       + H001 * (-128. + z * (-320. + 256. * z))
                       + H02
                             * (32.
                                + z
                                      * (202.66666666666669
                                         + (-192. - 51.2 * z) * z))
                       + ((-8.533333333333333 * H0m1
                           + 8.533333333333333 * H0 * Hm1)
                          * (1. + z)
                          * (2. + z * (-2. + z * (-73. + z * (3 + 12. * z)))))
                             / z2
                       + H1 * (-1. + z) * (1. + 2. * z)
                             * (H0 * (-128. + 64. * H0) - 64. * H1
                                - 128. * zeta2)
                       + (-96.
                          + z * (-533.3333333333334 + z * (256. + 102.4 * z)))
                             * zeta2
                       + (1. + z) * (-1. + 2. * z)
                             * (-384. * H00m1 + 128. * H0 * H0m1 - 256. * H0m1m1
                                + 64. * H02 * Hm1 + 256. * H0m1 * Hm1
                                - 128. * H0 * Hm12 - 128. * Hm1 * zeta2)
                       + z
                             * (-10.666666666666666 * H03 - 768. * H0001
                                + 768. * H000m1 + 256. * H011 + 256. * H0m12
                                + H02 * (-128. * H01 + 128. * H0m1)
                                + 256. * H01 * zeta2 - 256. * H0m1 * zeta2
                                - 153.60000000000002 * zeta2_2 + 576. * zeta3
                                + H0
                                      * (512. * H001 - 512. * H00m1
                                         - 512. * H0m1m1 - 64. * zeta2
                                         + 512. * zeta3))
                       + (48. + 32. * H02 * z
                          + H1 * (64. + (64. - 128. * z) * z)
                          + H0 * (32. - 128. * z2)
                          + z * (-48. + 128. * H01 - 128. * zeta2))
                             * L_M)
                          * L_Q
                    + (-24. - 16. * H02 * z + H0 * (-16. + 64. * z2)
                       + H1 * (-32. + z * (-32. + 64. * z))
                       + z * (24. - 64. * H01 + 64. * zeta2))
                          * L_Q2)
           + CF * TR * TR
                 * (1311.9999999999998 + H02 * (96. - 608. * z)
                    + H03 * (10.666666666666664 - 74.66666666666666 * z)
                    - 21.333333333333332 / z - 3120. * z
                    - 5.333333333333333 * H04 * z + 1829.3333333333333 * z2
                    + H0 * (512. + (-448. - 384. * z) * z)
                    + ((5.688888888888888 * H0 * Hm1
                        + (0.9481481481481482 - 5.688888888888888 * H0
                           + 28.444444444444443 * H1)
                              * z
                        + (1125.6888888888889
                           + H0 * (514.8444444444444 + 128. * H0) + 128. * H01
                           + 256. * H1)
                              * z2
                        + (-2007.4666666666667
                           + H0
                                 * (-786.488888888889
                                    + (-433.77777777777777 - 128. * H0) * H0)
                           - 256. * H001 - 128. * H01)
                              * z3
                        + H0m1
                              * (-5.688888888888888 + 28.444444444444446 * z3
                                 - 34.13333333333333 * z5)
                        + z2
                              * (-128. * zeta2
                                 + z
                                       * (H02 * (-256. - 17.066666666666666 * z)
                                              * z
                                          + H1
                                                * (-810.6666666666666
                                                   + 526.2222222222222 * z)
                                          + 99.55555555555554 * zeta2
                                          + H0
                                                * (-28.444444444444446 * Hm1
                                                   + 605.8666666666667 * z
                                                   + 34.13333333333333 * Hm1
                                                         * z2
                                                   + 256. * zeta2)
                                          + z
                                                * (880.8296296296297
                                                   - 170.66666666666666 * H01
                                                   + 170.66666666666669 * zeta2
                                                   + 34.13333333333333 * z
                                                         * zeta2)
                                          + 256. * zeta3)))
                       * L_M)
                          / z2
                    + ((14.222222222222223
                        + z
                              * (106.66666666666669 - 64. * H02 * z
                                 + z
                                       * (-405.33333333333337
                                          + 284.44444444444446 * z)
                                 + H0
                                       * (64.
                                          + (-106.66666666666666
                                             - 85.33333333333333 * z)
                                                * z)))
                       * L_M2)
                          / z
                    + ((5.68888888888889 * H0 * Hm1
                        + (10.42962962962963 - 5.68888888888889 * H0
                           - 28.444444444444446 * H1)
                              * z
                        + (-1171.1999999999998
                           + (-551.8222222222223 - 128. * H0) * H0 - 128. * H01
                           - 298.66666666666663 * H1)
                              * z2
                        + (1747.1999999999998
                           + H0
                                 * (564.6222222222223
                                    + H0 * (376.8888888888889 + 128. * H0))
                           + 256. * H001 + 42.66666666666666 * H01)
                              * z3
                        + H0m1
                              * (-5.68888888888889 + 28.444444444444443 * z3
                                 - 34.13333333333333 * z5)
                        + z2
                              * (128. * zeta2
                                 + z
                                       * (H1
                                              * (682.6666666666666
                                                 - 355.55555555555554 * z)
                                          + H02
                                                * (256. - 17.066666666666666 * z
                                                )
                                                * z
                                          + H0
                                                * (-28.444444444444443 * Hm1
                                                   - 503.4666666666667 * z
                                                   + 34.13333333333333 * Hm1
                                                         * z2
                                                   - 256. * zeta2)
                                          - 71.11111111111111 * zeta2
                                          + z
                                                * (-586.4296296296296
                                                   + 170.66666666666666 * H01
                                                   - 170.66666666666666 * zeta2
                                                   + 34.13333333333333 * z
                                                         * zeta2)
                                          - 256. * zeta3))
                        + z
                              * (-28.444444444444443
                                 + z
                                       * (-256. + 128. * H02 * z
                                          + (768. - 483.55555555555554 * z) * z
                                          + H0
                                                * (-128.
                                                   + z
                                                         * (128.
                                                            + 170.66666666666666
                                                                  * z))))
                              * L_M)
                       * L_Q)
                          / z2
                    + ((14.222222222222221
                        + z
                              * (149.33333333333331 - 64. * H02 * z
                                 + z
                                       * (-362.66666666666663
                                          + 199.1111111111111 * z)
                                 + H0
                                       * (64.
                                          + (-21.33333333333333
                                             - 85.33333333333333 * z)
                                                * z)))
                       * L_Q2)
                          / z
                    + nf
                          * (((5.68888888888889 * H0 * Hm1
                               + (10.42962962962963 - 5.68888888888889 * H0
                                  - 28.444444444444446 * H1)
                                     * z
                               + (-1171.1999999999998
                                  + (-551.8222222222223 - 128. * H0) * H0
                                  - 128. * H01 - 298.66666666666663 * H1)
                                     * z2
                               + (1747.1999999999998
                                  + H0
                                        * (564.6222222222223
                                           + H0
                                                 * (376.8888888888889
                                                    + 128. * H0))
                                  + 256. * H001 + 42.66666666666666 * H01)
                                     * z3
                               + H0m1
                                     * (-5.68888888888889
                                        + 28.444444444444443 * z3
                                        - 34.13333333333333 * z5)
                               + z2
                                     * (128. * zeta2
                                        + z
                                              * (H1
                                                     * (682.6666666666666
                                                        - 355.55555555555554 * z
                                                     )
                                                 + H02
                                                       * (256.
                                                          - 17.066666666666666
                                                                * z)
                                                       * z
                                                 + H0
                                                       * (-28.444444444444443
                                                              * Hm1
                                                          - 503.4666666666667
                                                                * z
                                                          + 34.13333333333333
                                                                * Hm1 * z2
                                                          - 256. * zeta2)
                                                 - 71.11111111111111 * zeta2
                                                 + z
                                                       * (-586.4296296296296
                                                          + 170.66666666666666
                                                                * H01
                                                          - 170.66666666666669
                                                                * zeta2
                                                          + 34.13333333333333
                                                                * z * zeta2)
                                                 - 256. * zeta3)))
                              * L_Q)
                                 / z2
                             + L_M
                                   * (88.88888888888889
                                      + 42.666666666666664 * H02 * z
                                      + H1
                                            * (21.333333333333332
                                               + (21.333333333333332
                                                  - 42.666666666666664 * z)
                                                     * z)
                                      + H0
                                            * (21.333333333333332
                                               + (220.44444444444446
                                                  - 42.666666666666664 * z)
                                                     * z)
                                      + z
                                            * (152.88888888888889
                                               + 42.666666666666664 * H01
                                               - 241.77777777777777 * z
                                               - 42.666666666666664 * zeta2)
                                      + (-21.333333333333332
                                         + z
                                               * (-21.333333333333332
                                                  - 42.666666666666664 * H0
                                                  + 42.666666666666664 * z))
                                            * L_Q)
                             + ((14.222222222222221
                                 + z
                                       * (149.33333333333331 - 64. * H02 * z
                                          + z
                                                * (-362.66666666666663
                                                   + 199.1111111111111 * z)
                                          + H0
                                                * (64.
                                                   + (-21.33333333333333
                                                      - 85.33333333333333 * z)
                                                         * z)))
                                * L_Q2)
                                   / z))
           + CA
                 * (CF * TR
                        * (-74.66666666666666 + H0 * H01 * (64. - 224. * z)
                           + H03 * (-5.333333333333333 - 16. * z)
                           - 33.18518518518518 / z + 2677.333333333333 * z
                           - 2569.4814814814813 * z2
                           + (192. * H0m1 - 192. * H0 * Hm1) * z * (1. + z)
                           - 72. * H12 * (-1. + z)
                                 * (-0.1111111111111111 + 1. * z)
                           - 256. * H1 * (-1. + z) * (-0.0625 + 1. * z)
                           + (16. * H0 - 5.333333333333333 * H1) * H12
                                 * (-1. + z) * (1. + 2. * z)
                           + H011 * (32. + (64. - 64. * z) * z)
                           + H02
                                 * (-7.999999999999999
                                    + (48. - 29.333333333333332 * z) * z)
                           + H01
                                 * (-112. - 21.333333333333332 / z + 624. * z
                                    - 362.66666666666663 * z2)
                           + (410.66666666666663 * H0 * H1 * (-1. + z)
                              * (-0.05194805194805195
                                 + z * (-0.3246753246753247 + 1. * z)))
                                 / z
                           + H0
                                 * (-69.33333333333333
                                    + z
                                          * (853.3333333333333
                                             + 1372.4444444444443 * z))
                           + (-272. - 48. * z) * z * zeta2
                           + (1. + z) * (-1. + 2. * z)
                                 * (-32. * H00m1 + 32. * H0 * H0m1
                                    + 64. * H0m1m1 - 16. * H02 * Hm1
                                    - 64. * H0m1 * Hm1 + 32. * H0 * Hm12
                                    + 32. * Hm1 * zeta2)
                           + z
                                 * (2.6666666666666665 * H04 + 768. * H0001
                                    - 192. * H000m1 + 128. * H0011 + 64. * H0111
                                    - 32. * H02 * H0m1 - 64. * H0m12
                                    + 64. * H0m1 * zeta2
                                    - 294.40000000000003 * zeta2_2
                                    + H0
                                          * (-256. * H001 + 128. * H00m1
                                             - 64. * H011 + 128. * H0m1m1
                                             - 32. * zeta2 - 448. * zeta3))
                           + (-4. + 15. * z) * (32. * H001 - 32. * zeta3)
                           + (-203.91111111111113 - 27.022222222222222 / z
                              + 899.0222222222224 * z - 668.088888888889 * z2
                              - 10.666666666666666 * H1 * (-1. + z)
                                    * (-11.5 + 1. * z)
                              + H00m1 * (192. + (64. - 384. * z) * z)
                              + H0
                                    * (-22.400000000000002
                                       - 8.533333333333333 / z
                                       - 15.28888888888889 * z
                                       + 791.4666666666667 * z2)
                              + H02 * (32. + z * (64. + (-32. - 25.6 * z) * z))
                              - (51.2 * (1. * H0m1 - 1. * H0 * Hm1) * (1. + z)
                                 * (0.16666666666666666
                                    + z
                                          * (-0.16666666666666666
                                             + z
                                                   * (-4.833333333333333
                                                      + z * (-2.25 + 1. * z)))))
                                    / z2
                              + (1. + z) * (-1. + 2. * z)
                                    * (-128. * H0m1m1
                                       + Hm1
                                             * (32. * H02 + 128. * H0m1
                                                - 64. * H0 * Hm1 - 64. * zeta2))
                              + z * (-181.33333333333334 + 51.2 * z2) * zeta2
                              + (-1. + z) * (1. + 2. * z)
                                    * (64. * H001 - 64. * H0 * H01
                                       + 64. * H0 * H0m1 + 32. * H02 * H1
                                       - 32. * H12 - 64. * H1 * zeta2)
                              + z
                                    * (-42.666666666666664 * H03 - 384. * H0001
                                       + 384. * H000m1
                                       - 181.33333333333331 * H01 + 128. * H011
                                       + 128. * H0m12
                                       + H02 * (-64. * H01 + 64. * H0m1)
                                       + 128. * H01 * zeta2
                                       - 128. * H0m1 * zeta2
                                       - 76.80000000000001 * zeta2_2
                                       + H0
                                             * (256. * H001 - 256. * H00m1
                                                - 256. * H0m1m1 + 256. * zeta3))
                             ) * L_M
                           + ((-10.666666666666666
                               + (-48. - 32. * H0 + 32. * H1) * z
                               + (-229.33333333333331 - 64. * H1) * z3
                               + z2
                                     * (288. + H0 * (96. + 64. * H0) + 64. * H01
                                        + 32. * H1 - 64. * zeta2))
                              * L_M2)
                                 / z
                           + (267.2 - 29.866666666666667 / z
                              + 1729.2444444444445 * z - 1966.5777777777778 * z2
                              + 128. * H0 * H01 * (1. + z) * (-0.5 + 1. * z)
                              - 192. * H12 * (-1. + z)
                                    * (0.16666666666666666 + 1. * z)
                              - 256. * H0 * H1 * (-1. + z) * (0.25 + 1. * z)
                              - 64. * H02 * H1 * (-1. + z) * (0.5 + 1. * z)
                              + H0 * H0m1
                                    * (64. + 17.066666666666666 / z2
                                       + (149.33333333333334 - 128. * z) * z)
                              + H0
                                    * (46.22222222222222
                                       + 22.755555555555556 / z
                                       + 1860.6222222222223 * z
                                       - 1659.7333333333333 * z2)
                              + H1
                                    * (217.60000000000002
                                       - 59.733333333333334 / z
                                       + 1425.0666666666666 * z
                                       - 1582.9333333333334 * z2)
                              + H03 * z
                                    * (99.55555555555554
                                       + 17.066666666666666 * z2)
                              + H01
                                    * (-72.53333333333333
                                       + 17.066666666666666 / z
                                       + 887.4666666666667 * z + 102.4 * z2)
                              + H00m1
                                    * (-192. - 17.066666666666666 / z2
                                       - 320. * z + 384. * z2 + 102.4 * z3)
                              + ((17.066666666666666 * H01m1
                                  + 17.066666666666666 * H0m11
                                  - 17.066666666666666 * H01 * Hm1)
                                 * (1. + z)
                                 * (1. + z * (-1 + z - 6. * z2 + 6. * z3)))
                                    / z2
                              + H001
                                    * (64.
                                       + z
                                             * (234.66666666666666
                                                + z * (-128. + 102.4 * z)))
                              + (238.93333333333334
                                 * (1. * H0m1 - 1. * H0 * Hm1) * (1. + z)
                                 * (0.16666666666666666
                                    + z
                                          * (-0.09523809523809523
                                             + z
                                                   * (-0.8333333333333334
                                                      + z
                                                            * (-1.6428571428571428
                                                               + 1. * z)))))
                                    / z2
                              - (51.2 * H02 * Hm1 * (1. + z)
                                 * (0.16666666666666666
                                    + z
                                          * (-0.16666666666666666
                                             + z
                                                   * (-0.4583333333333333
                                                      + z * (0.25 + 1. * z)))))
                                    / z2
                              + (102.4
                                 * (1. * H0m1m1
                                    + Hm1 * (-1. * H0m1 + 0.5 * H0 * Hm1))
                                 * (1. + z)
                                 * (0.16666666666666666
                                    + z
                                          * (-0.16666666666666666
                                             + z
                                                   * (-1.0833333333333333
                                                      + z * (1.5 + 1. * z)))))
                                    / z2
                              + (H02
                                 * (8.533333333333333
                                    + z
                                          * (-68.26666666666667
                                             + z
                                                   * (386.84444444444443
                                                      + z
                                                            * (-12.799999999999999
                                                               + 119.46666666666667
                                                                     * z)))))
                                    / z
                              + H0 * z * (-256. - 204.8 * z2) * zeta2
                              + H1
                                    * (-64. + 8.533333333333333 / z2
                                       + z
                                             * (-21.333333333333343
                                                + (128. - 51.2 * z) * z))
                                    * zeta2
                              + Hm1
                                    * (-64. + 25.6 / z2
                                       + z
                                             * (-64.
                                                + z
                                                      * (128.00000000000003
                                                         + 153.60000000000002
                                                               * z)))
                                    * zeta2
                              + ((-34.13333333333333
                                  + z
                                        * (8.533333333333333
                                           + z
                                                 * (-487.8222222222223
                                                    + (153.6
                                                       - 238.93333333333334 * z)
                                                          * z)))
                                 * zeta2)
                                    / z
                              + z
                                    * (384. * H0001 - 384. * H000m1
                                       + 128. * H011
                                       + H02 * (64. * H01 - 64. * H0m1)
                                       - 128. * H0m12 - 128. * H01 * zeta2
                                       + 128. * H0m1 * zeta2
                                       + 76.80000000000001 * zeta2_2
                                       + H0
                                             * (-256. * H001 + 256. * H00m1
                                                + 256. * H0m1m1 - 256. * zeta3))
                              + z * (-170.66666666666666 - 256. * z2) * zeta3
                              + (58.666666666666664
                                 + (58.666666666666664 + 117.33333333333333 * H0
                                    - 117.33333333333333 * z)
                                       * z)
                                    * L_M)
                                 * L_Q
                           + ((10.666666666666666
                               + (-10.666666666666666 + 32. * H0 - 32. * H1) * z
                               + (346.66666666666663 + 64. * H1) * z3
                               + z2
                                     * (-346.66666666666663
                                        + (-213.33333333333331 - 64. * H0) * H0
                                        - 64. * H01 - 32. * H1 + 64. * zeta2))
                              * L_Q2)
                                 / z)
                    + TR * TR
                          * ((54.913580246913575
                              + z
                                    * (-13.037037037037031
                                       + H02
                                             * (10.666666666666666
                                                - 99.55555555555554 * z)
                                       - 7.111111111111111 * H03 * z
                                       + H0
                                             * (78.22222222222221
                                                + (-656.5925925925925
                                                   - 260.7407407407407 * z)
                                                      * z)
                                       + z
                                             * (-1114.074074074074
                                                - 21.333333333333332 * H01
                                                - 154.07407407407408 * H1
                                                + 1072.1975308641975 * z
                                                + 154.07407407407408 * H1 * z
                                                + 21.333333333333332 * zeta2))
                              + (7.111111111111111
                                 + z
                                       * (-21.333333333333332
                                          + z
                                                * (-106.66666666666666
                                                   - 85.33333333333333 * H0
                                                   + 120.88888888888889 * z
                                                   + H1
                                                         * (-42.666666666666664
                                                            + 42.666666666666664
                                                                  * z))))
                                    * L_M2
                              + (-59.25925925925925
                                 + H12
                                       * (42.666666666666664
                                          - 42.666666666666664 * z)
                                       * z2
                                 + H1
                                       * (-14.222222222222221
                                          + z
                                                * (42.666666666666664
                                                   + (526.2222222222222
                                                      + H0
                                                            * (85.33333333333333
                                                               - 85.33333333333333
                                                                     * z)
                                                      - 554.6666666666666 * z)
                                                         * z))
                                 + z
                                       * (85.33333333333333
                                          + z
                                                * (1066.6666666666665
                                                   + 170.66666666666666 * H01
                                                   - 85.33333333333333 * H0m1
                                                   + H0
                                                         * (839.1111111111111
                                                            + 170.66666666666666
                                                                  * H0
                                                            + 85.33333333333333
                                                                  * Hm1)
                                                   - 170.66666666666666 * zeta2)
                                          + z2
                                                * (-1092.7407407407406
                                                   - 369.7777777777778 * H0
                                                   - 85.33333333333333 * H0m1
                                                   + 85.33333333333333 * H0
                                                         * Hm1
                                                   + 85.33333333333333 * zeta2))
                                ) * L_Q
                              + L_M
                                    * (49.77777777777777
                                       + (-56.888888888888886
                                          + 85.33333333333333 * H0)
                                             * z
                                       + H12
                                             * (42.666666666666664
                                                - 42.666666666666664 * z)
                                             * z2
                                       + H1
                                             * (-14.222222222222221
                                                + z
                                                      * (42.666666666666664
                                                         + (241.77777777777777
                                                            + H0
                                                                  * (85.33333333333333
                                                                     - 85.33333333333333
                                                                           * z)
                                                            - 270.22222222222223
                                                                  * z)
                                                               * z))
                                       + z2
                                             * (-476.44444444444446
                                                + 170.66666666666666 * H01
                                                - 85.33333333333333 * H0m1
                                                + H0
                                                      * (-156.44444444444443
                                                         + 85.33333333333333
                                                               * H0
                                                         + 85.33333333333333
                                                               * Hm1)
                                                - 170.66666666666666 * zeta2)
                                       + z3
                                             * (483.55555555555554
                                                - 85.33333333333333 * H0m1
                                                + H0
                                                      * (-739.5555555555557
                                                         + 85.33333333333333
                                                               * Hm1)
                                                + 85.33333333333333 * zeta2)
                                       + (14.222222222222221
                                          + z
                                                * (-42.666666666666664
                                                   + z
                                                         * (-213.33333333333331
                                                            - 170.66666666666666
                                                                  * H0
                                                            + 241.77777777777777
                                                                  * z
                                                            + H1
                                                                  * (-85.33333333333333
                                                                     + 85.33333333333333
                                                                           * z))
                                                ))
                                             * L_Q)
                              + 7.111111111111111 * L_Q2
                              + z
                                    * (-21.333333333333332
                                       + z
                                             * (-106.66666666666666
                                                - 85.33333333333333 * H0
                                                + 120.88888888888889 * z
                                                + H1
                                                      * (-42.666666666666664
                                                         + 42.666666666666664
                                                               * z)))
                                    * L_Q2)
                                 / z
                             + (nf
                                * ((-59.25925925925925
                                    + H12
                                          * (42.666666666666664
                                             - 42.666666666666664 * z)
                                          * z2
                                    + H1
                                          * (-14.222222222222221
                                             + z
                                                   * (42.666666666666664
                                                      + (526.2222222222222
                                                         + H0
                                                               * (85.33333333333333
                                                                  - 85.33333333333333
                                                                        * z)
                                                         - 554.6666666666666 * z
                                                        ) * z))
                                    + z
                                          * (85.33333333333333
                                             + z
                                                   * (1066.6666666666665
                                                      + 170.66666666666666 * H01
                                                      - 85.33333333333333 * H0m1
                                                      + H0
                                                            * (839.1111111111111
                                                               + 170.66666666666666
                                                                     * H0
                                                               + 85.33333333333333
                                                                     * Hm1)
                                                      - 170.66666666666666
                                                            * zeta2)
                                             + z2
                                                   * (-1092.7407407407406
                                                      - 369.7777777777778 * H0
                                                      - 85.33333333333333 * H0m1
                                                      + 85.33333333333333 * H0
                                                            * Hm1
                                                      + 85.33333333333333
                                                            * zeta2)))
                                       * L_Q
                                   + (7.111111111111111
                                      + z
                                            * (-21.333333333333332
                                               + z
                                                     * (-106.66666666666666
                                                        - 85.33333333333333 * H0
                                                        - 42.666666666666664
                                                              * H1
                                                        + 120.88888888888889 * z
                                                        + 42.666666666666664
                                                              * H1 * z)))
                                         * L_Q2))
                                   / z))
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

    double H02 = H0 * H0;
    double H03 = H02 * H0;
    double H04 = H03 * H0;

    double H12 = H1 * H1;

    return (CF * TR * TR
            * (44.24691358024691
               + H12
                     * (7.111111111111111 - 21.333333333333332 * z
                        + 14.222222222222221 * z3)
               + H1
                     * (-23.703703703703702
                        + z
                              * (99.55555555555554
                                 + (14.222222222222229 - 90.07407407407408 * z)
                                       * z))
               + z
                     * (-180.14814814814812
                        + z2
                              * (173.82716049382714 - 90.07407407407408 * H0
                                 + 28.444444444444443 * H01
                                 - 28.444444444444443 * zeta2)
                        + z
                              * (-37.925925925925924 - 317.6296296296296 * H0
                                 + 156.44444444444443 * H01
                                 - 42.666666666666664 * H011
                                 - 156.44444444444443 * zeta2
                                 + 42.666666666666664 * zeta3))
               + (47.407407407407405
                  + H1
                        * (-28.444444444444443 + 85.33333333333333 * z
                           - 56.888888888888886 * z3)
                  + z
                        * (-199.1111111111111
                           + z
                                 * (-28.444444444444457
                                    - 312.88888888888886 * H0
                                    + 85.33333333333333 * H01
                                    + 180.14814814814815 * z
                                    - 56.888888888888886 * H0 * z
                                    - 85.33333333333333 * zeta2)))
                     * L_M
               + (14.222222222222221
                  + z
                        * (-42.666666666666664 - 42.666666666666664 * H0 * z
                           + 28.444444444444443 * z2))
                     * L_M2
               + (-56.888888888888886 - 56.888888888888886 * nf) * L_Q
               + z
                     * (227.55555555555554 + 85.33333333333333 * H0
                        + 227.55555555555554 * nf + 85.33333333333333 * H0 * nf
                        + (-85.33333333333333 - 85.33333333333333 * nf
                           + H0
                                 * (227.55555555555554 + 227.55555555555554 * nf
                                    + H0
                                          * (85.33333333333333
                                             + 85.33333333333333 * nf)))
                              * z
                        + (-85.33333333333333
                           + H0
                                 * (-113.77777777777777
                                    - 113.77777777777777 * nf)
                           - 85.33333333333333 * nf)
                              * z2)
                     * L_Q
               + 14.222222222222221 * L_Q2
               + (14.222222222222221 * nf - 42.666666666666664 * z
                  - 42.666666666666664 * nf * z
                  + H0 * (-42.666666666666664 - 42.666666666666664 * nf) * z2
                  + (28.444444444444443 + 28.444444444444443 * nf) * z3)
                     * L_Q2))
               / z
           + (CA * CF * TR
              * ((54.51851851851852 - 85.33333333333333 * H0m1
                  + 231.1111111111111 * H1 - 344.8888888888889 * z
                  + 256. * H01 * z - 256. * H0m1 * z
                  + H02 * (128. - 202.66666666666666 * z) * z
                  - 85.33333333333333 * H03 * z2
                  + H12
                        * (42.666666666666664 - 128. * z
                           + 85.33333333333333 * z3)
                  + H1 * z
                        * (-149.33333333333331
                           + z * (-458.66666666666663 + 376.8888888888889 * z))
                  + H0
                        * (14.222222222222221
                           + Hm1
                                 * (85.33333333333333 + 256. * z
                                    - 170.66666666666666 * z3)
                           + H1
                                 * (85.33333333333333 - 256. * z
                                    + 170.66666666666666 * z3)
                           + z
                                 * (-405.3333333333333
                                    + z
                                          * (-487.1111111111111 - 256. * H01
                                             - 256. * H0m1
                                             + 1500.4444444444443 * z
                                             + 256. * zeta2)))
                  + z2
                        * (1084.4444444444443 + 512. * H00m1 - 64. * H01
                           - 256. * H011
                           + z
                                 * (-794.074074074074
                                    + 170.66666666666666 * H0m1
                                    - 170.66666666666666 * zeta2)
                           + 64. * zeta2 - 128. * zeta3))
                     * L_Q
                 + (-112. + 64. * H02 * z2
                    + H1
                          * (-21.333333333333332 + 64. * z
                             - 42.666666666666664 * z3)
                    + H0
                          * (-21.333333333333332
                             + z * (-64. + 149.33333333333331 * z))
                    + z
                          * (106.66666666666669
                             + z
                                   * (250.66666666666666 + 64. * H01
                                      - 245.33333333333331 * z - 64. * zeta2)))
                       * L_Q2))
                 / z
           + CF * CF * TR
                 * (-74.66666666666666
                    + H03 * (-5.333333333333333 - 13.333333333333332 * z)
                    - 33.18518518518518 / z + 629.3333333333333 * z
                    - 521.4814814814814 * z2
                    + H02 * (-8. - 21.333333333333332 * z2)
                    + ((-10.666666666666666 * H01 + 10.666666666666666 * H0 * H1
                       )
                       * (-1. + z) * (-2. + z * (-11. + 4. * z)))
                          / z
                    + H0
                          * (-69.33333333333333
                             + z * (165.33333333333334 + 284.44444444444446 * z)
                          )
                    + (2. + z) * (-64. * H001 + 32. * H0 * H01 + 64. * zeta3)
                    + z
                          * (1.3333333333333333 * H04 + 192. * H0001
                             - 64. * H0 * H001 - 76.80000000000001 * zeta2_2
                             - 128. * H0 * zeta3)
                    + ((-35.55555555555556
                        + (96.00000000000001 + H0 * (32. + 32. * H0)) * z
                        + (224.00000000000003
                           + H0 * (192. + (80. - 10.666666666666666 * H0) * H0))
                              * z2
                        + (-284.44444444444446 + 85.33333333333333 * H0) * z3)
                       * L_M)
                          / z
                    + (-48. + H0 * (-32. - 16. * z) - 10.666666666666666 / z
                       + 80. * z + 16. * H02 * z - 21.333333333333332 * z2)
                          * L_M2
                    + (-27.733333333333334 + 4.266666666666667 / z
                       + 84.62222222222223 * z - 61.15555555555556 * z2
                       + H0
                             * (-100.97777777777779 - 11.377777777777778 / z
                                + 578.1333333333333 * z
                                - 39.82222222222222 * z2)
                       + H1
                             * (-426.66666666666663 + 35.55555555555556 / z
                                + 362.66666666666663 * z
                                + 28.444444444444443 * z2)
                       + H01
                             * (-128. - 42.666666666666664 / z - 256. * z
                                + 170.66666666666666 * z2)
                       + (H1
                          * (42.666666666666664 * H0 + 21.333333333333332 * H1)
                          * (-1. + z) * (-1. + z * (2. + 2. * z)))
                             / z
                       + (17.066666666666666 * (1. * H0m1 - 1. * H0 * Hm1)
                          * (1. + z)
                          * (-0.6666666666666666
                             + z
                                   * (0.6666666666666666
                                      + z
                                            * (1.8333333333333333
                                               + z * (-1. + 1. * z)))))
                             / z2
                       + (-90. + z * (-85. + z * (90. + 6. * z)))
                             * (1.4222222222222223 * H02
                                - 2.8444444444444446 * zeta2)
                       + z
                             * (-7.111111111111111 * H03
                                - 85.33333333333333 * H00m1 - 128. * H011
                                + H0
                                      * (-128. * H01 + 42.666666666666664 * H0m1
                                         + 298.66666666666663 * zeta2)
                                + 192. * zeta3))
                          * L_Q
                    + ((10.666666666666666
                        + H1
                              * (-21.333333333333332 + 64. * z
                                 - 42.666666666666664 * z3)
                        + z
                              * (85.33333333333333
                                 + H0 * (64. - 42.666666666666664 * z2)
                                 + z
                                       * (-117.33333333333333 + 64. * H01
                                          + 21.333333333333332 * z
                                          - 64. * zeta2)))
                       * L_Q2)
                          / z)
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
    double x4 = x3 * x;
    double x5 = x4 * x;

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

    return Lmmu2
               * (-0.8888888888888888
                  + (7.111111111111111 - 7.111111111111111 * x) * x
                  + H0
                        * (-0.8888888888888888
                           + (1.7777777777777777 - 1.7777777777777777 * x) * x)
                  + H1
                        * (-0.8888888888888888
                           + (1.7777777777777777 - 1.7777777777777777 * x) * x)
                  + LQm
                        * (0.8888888888888888
                           + x * (-1.7777777777777777 + 1.7777777777777777 * x))
               )
           + (0.8888888888888888
              + x * (-1.7777777777777777 + 1.7777777777777777 * x))
                 * zeta3
           + CF
                 * (LQm3
                        * (-12.444444444444443 - 2.6666666666666665 * H00
                           - 1.7777777777777777 * H1 + 1.1851851851851851 / x
                           + 21.777777777777775 * x
                           + (5.333333333333333 * H00
                              + H1
                                    * (3.5555555555555554
                                       - 3.5555555555555554 * x)
                              - 9.185185185185185 * x)
                                 * x
                           + H0 * (-4.888888888888888 + 1.7777777777777777 * x))
                    + LQm2
                          * (184.33333333333331 + 126.88888888888887 * H0
                             + 41.33333333333333 * H00 + 24. * H000 + 8. * H001
                             + 17.333333333333332 * H01
                             + 10.666666666666666 * H10 + 8. * H11
                             + 9.481481481481481 / x
                             - (3.555555555555556 * H1) / x
                             + H1
                                   * (49.77777777777777
                                      + x
                                            * (-95.55555555555556
                                               + 57.77777777777777 * x))
                             + x2
                                   * (79.85185185185185 + 54.22222222222222 * H0
                                      - 10.666666666666666 * H00
                                      + 5.333333333333333 * H01
                                      + 21.333333333333332 * H10 + 16. * H11
                                      - 5.333333333333333 * zeta2)
                             - 17.333333333333332 * zeta2 - 8. * H0 * zeta2
                             - 8. * zeta3
                             + x
                                   * (-282.66666666666663
                                      - 74.66666666666666 * H00 - 48. * H000
                                      - 16. * H001 - 10.666666666666666 * H01
                                      - 21.333333333333332 * H10 - 16. * H11
                                      + 10.666666666666666 * zeta2
                                      + H0 * (-52.44444444444444 + 16. * zeta2)
                                      + 16. * zeta3))
                    + Lmmu
                          * (LQm2
                                 * (-1.3333333333333333 + 5.333333333333333 * x
                                    + H0
                                          * (-2.6666666666666665
                                             + (5.333333333333333
                                                - 10.666666666666666 * x)
                                                   * x)
                                    + H1
                                          * (-5.333333333333333
                                             + (10.666666666666666
                                                - 10.666666666666666 * x)
                                                   * x))
                             + LQm
                                   * (28. + 10.666666666666666 * H00 + 16. * H01
                                      + 18.666666666666664 * H1
                                      + 21.333333333333332 * H10
                                      + 21.333333333333332 * H11
                                      - 53.33333333333333 * x
                                      + H0
                                            * (5.333333333333333
                                               + x
                                                     * (-32.
                                                        + 53.33333333333333 * x)
                                            )
                                      - 16. * zeta2
                                      + x
                                            * (-32. * H01 - 64. * H1
                                               - 42.666666666666664 * H10
                                               - 42.666666666666664 * H11
                                               + 18.666666666666664 * x
                                               + H00
                                                     * (-21.333333333333332
                                                        + 42.666666666666664 * x
                                                     )
                                               + x
                                                     * (42.666666666666664 * H01
                                                        + 53.33333333333333 * H1
                                                        + 42.666666666666664
                                                              * H10
                                                        + 42.666666666666664
                                                              * H11
                                                        - 42.666666666666664
                                                              * zeta2)
                                               + 32. * zeta2))
                             + (64.
                                * (Hm10
                                       * (0.011111111111111112 + x2
                                          + 0.4444444444444444 * x3 + 0.4 * x5)
                                   + x
                                         * (0.011111111111111112
                                            + H0
                                                  * (-0.011111111111111112
                                                     + x
                                                           * (-0.5569444444444445
                                                              + x
                                                                    * (0.4083333333333333
                                                                       + x
                                                                             * (-1.5250000000000001
                                                                                + zeta2
                                                                             )
                                                                       - 0.6666666666666666
                                                                             * zeta2
                                                                    )
                                                              + 0.3333333333333333
                                                                    * zeta2))
                                            + x4 * (-0.4 * H00 + 0.4 * zeta2)
                                            + x2
                                                  * (0.6416666666666666
                                                     + 0.05555555555555555 * H00
                                                     + 0.3333333333333333 * H000
                                                     + 0.6666666666666666 * H001
                                                     + 1.6666666666666665 * H01
                                                     + 0.6666666666666666 * H010
                                                     + 0.8333333333333333 * H011
                                                     + 1.5 * H1
                                                     + 1.1666666666666665 * H10
                                                     + H101
                                                     + 1.8333333333333333 * H11
                                                     + 0.6666666666666666 * H110
                                                     + H111
                                                     + 0.6666666666666666
                                                           * Hm100
                                                     - 1.3333333333333333
                                                           * Hm1m10
                                                     - 1.222222222222222 * zeta2
                                                     - 0.3333333333333333 * H1
                                                           * zeta2
                                                     - 0.6666666666666666 * Hm1
                                                           * zeta2
                                                     + 0.16666666666666666
                                                           * zeta3)
                                            + x
                                                  * (-0.6902777777777778
                                                     - 0.08333333333333333 * H00
                                                     - 0.16666666666666666
                                                           * H000
                                                     - 0.3333333333333333 * H001
                                                     - 0.25 * H01
                                                     - 0.3333333333333333 * H010
                                                     - 0.41666666666666663
                                                           * H011
                                                     + 0.6666666666666666
                                                           * H0m10
                                                     - 0.35416666666666663 * H1
                                                     - 0.5 * H10 - 0.5 * H101
                                                     - 0.4583333333333333 * H11
                                                     - 0.3333333333333333 * H110
                                                     - 0.5 * H111
                                                     + 0.3333333333333333
                                                           * Hm100
                                                     - 0.6666666666666666
                                                           * Hm1m10
                                                     + 0.25 * zeta2
                                                     + 0.16666666666666666 * H1
                                                           * zeta2
                                                     - 0.3333333333333333 * Hm1
                                                           * zeta2
                                                     + 0.5833333333333334
                                                           * zeta3)
                                            + x3
                                                  * (0.18333333333333332
                                                     - 1.0833333333333333 * H00
                                                     - 0.6666666666666666 * H000
                                                     - 1. * H001 - 1.75 * H01
                                                     - 0.6666666666666666 * H010
                                                     - 1. * H011
                                                     + 0.6666666666666666
                                                           * H0m10
                                                     - 1.125 * H1
                                                     - 1.0833333333333333 * H10
                                                     - 0.3333333333333333 * H100
                                                     - 1. * H101 - 1.75 * H11
                                                     - 0.6666666666666666 * H110
                                                     - 1. * H111
                                                     + 0.3333333333333333
                                                           * Hm100
                                                     - 0.6666666666666666
                                                           * Hm1m10
                                                     + 1.75 * zeta2
                                                     + 0.6666666666666666 * H1
                                                           * zeta2
                                                     - 0.3333333333333333 * Hm1
                                                           * zeta2
                                                     + 1.6666666666666665
                                                           * zeta3))))
                                   / x2)
                    + (LQm
                       * ((28.99753086419753 - 0.3555555555555556 * H0
                           - 18.962962962962962 * H1 + 7.111111111111111 * H10
                           + 7.111111111111111 * H11)
                              * x
                          + Hm10
                                * (0.35555555555555557
                                   + x
                                         * (-7.111111111111111
                                            + x
                                                  * (80.
                                                     + x
                                                           * (14.222222222222221
                                                              + x
                                                                    * (-55.11111111111111
                                                                       + 12.8
                                                                             * x
                                                                    )))))
                          + x
                                * (-7.111111111111111 * zeta2
                                   + x4 * (-12.8 * H00 + 12.8 * zeta2)
                                   + x2
                                         * (2046.9037037037037
                                            + 472.88888888888886 * H00
                                            + 480. * H000 + 192. * H0000
                                            + 96. * H0001 + 160. * H001
                                            + 32. * H0010 + 32. * H0011
                                            + 138.66666666666666 * H01
                                            + 21.333333333333332 * H010
                                            + 21.333333333333332 * H011
                                            - 64. * H0m10
                                            + 625.7777777777777 * H1
                                            + 251.55555555555554 * H10
                                            + 42.666666666666664 * H100
                                            + 53.33333333333333 * H101
                                            + 224.88888888888889 * H11
                                            + 32. * H110 + 32. * H111
                                            + 21.333333333333332 * Hm100
                                            - 42.666666666666664 * Hm1m10
                                            - 124.44444444444444 * zeta2
                                            + (-96. * H00 - 32. * H1
                                               - 21.333333333333332 * Hm1
                                               - 3.2 * zeta2)
                                                  * zeta2
                                            + H0
                                                  * (500.62222222222226
                                                     - 224. * zeta2
                                                     - 64. * zeta3)
                                            - 224. * zeta3)
                                   + x3
                                         * (-666.5530864197531
                                            - 99.55555555555554 * H00
                                            + 64. * H000
                                            + 10.666666666666666 * H001
                                            - 142.22222222222223 * H01
                                            - 10.666666666666666 * H010
                                            - 10.666666666666666 * H011
                                            + 42.666666666666664 * H0m10
                                            - 228.14814814814812 * H1
                                            - 176. * H10
                                            - 53.33333333333333 * H100
                                            - 53.33333333333333 * H101
                                            - 149.33333333333331 * H11
                                            - 32. * H110 - 32. * H111
                                            + 10.666666666666666 * Hm100
                                            - 21.333333333333332 * Hm1m10
                                            + H0
                                                  * (-128.35555555555555
                                                     - 10.666666666666666
                                                           * zeta2)
                                            + 142.22222222222223 * zeta2
                                            + 42.666666666666664 * H1 * zeta2
                                            - 10.666666666666666 * Hm1 * zeta2
                                            + 42.666666666666664 * zeta3)
                                   + x
                                         * (-1402.125925925926
                                            - 528.8888888888889 * H00
                                            - 160. * H000 - 96. * H0000
                                            - 48. * H0001 - 88. * H001
                                            - 16. * H0010 - 16. * H0011
                                            - 270.66666666666663 * H01
                                            - 34.666666666666664 * H010
                                            - 34.666666666666664 * H011
                                            + 53.33333333333333 * H0m10
                                            - 388.66666666666663 * H1
                                            - 124.44444444444444 * H10
                                            - 21.333333333333332 * H100
                                            - 26.666666666666664 * H101
                                            - 116.44444444444444 * H11
                                            - 16. * H110 - 16. * H111
                                            + 10.666666666666666 * Hm100
                                            - 21.333333333333332 * Hm1m10
                                            + 270.66666666666663 * zeta2
                                            + zeta2
                                                  * (48. * H00 + 16. * H1
                                                     - 10.666666666666666 * Hm1
                                                     + 1.6 * zeta2)
                                            + 117.33333333333333 * zeta3
                                            + H0
                                                  * (-818.9333333333334
                                                     + 88. * zeta2
                                                     + 32. * zeta3)))))
                          / x2
                    + (-11.555555555555557 - 4.7407407407407405 * zeta2
                       + x3
                             * (388.55555555555554 - 90.66666666666666 * H00
                                + 26.666666666666664 * H000
                                + 10.666666666666666 * H0000
                                - 21.333333333333332 * H0010 - 16. * H01
                                - 21.333333333333332 * H0101 + 16. * H011
                                + 10.666666666666666 * H0111
                                - 42.666666666666664 * H10
                                + 26.666666666666664 * H100
                                + 10.666666666666666 * H1000
                                - 21.333333333333332 * H1010 + 32. * H11
                                - 21.333333333333332 * H1101 + 16. * H111
                                + 10.666666666666666 * H1111
                                + 36.74074074074074 * zeta2
                                + (21.333333333333332 * H00
                                   + 34.666666666666664 * H01
                                   + 10.666666666666666 * H10
                                   + 34.666666666666664 * H11
                                   - 58.13333333333333 * zeta2)
                                      * zeta2
                                + H0
                                      * (179.33333333333331
                                         + 2.6666666666666665 * zeta2
                                         - 32. * zeta3)
                                - 25.185185185185183 * zeta3)
                       + 1.1851851851851851 * zeta3
                       + H1
                             * (5.333333333333333
                                + x
                                      * (264.3333333333333
                                         - 2.2222222222222223 * zeta2
                                         - 17.77777777777778 * zeta3
                                         + x
                                               * (-519.3333333333333
                                                  - 8.88888888888889 * zeta2
                                                  + x
                                                        * (270.
                                                           + 8.88888888888889
                                                                 * zeta2
                                                           - 35.55555555555556
                                                                 * zeta3)
                                                  + 35.55555555555556 * zeta3)))
                       + x
                             * (52.166666666666664 + 10.666666666666666 * H000
                                + 2.6666666666666665 * H0000 + 8. * H00001
                                + 20. * H0001 + 56. * H001
                                - 5.333333333333333 * H0010 + 160. * H01
                                - 8. * H010 - 5.333333333333333 * H0100
                                - 10.666666666666666 * H0101
                                - 5.333333333333333 * H011
                                + 5.333333333333333 * H0111
                                - 26.666666666666664 * H10
                                + 2.6666666666666665 * H100
                                + 5.333333333333333 * H1000
                                - 10.666666666666666 * H1010
                                + 10.666666666666666 * H11
                                - 10.666666666666666 * H1101
                                - 5.333333333333333 * H111
                                + 5.333333333333333 * H1111
                                - 226.7222222222222 * zeta2
                                + (-20. * H000 + 17.333333333333332 * H01
                                   + 5.333333333333333 * H10
                                   + 17.333333333333332 * H11
                                   - 24.266666666666666 * zeta2)
                                      * zeta2
                                + H0
                                      * (7.
                                         + (-118.1111111111111 - 3.2 * zeta2)
                                               * zeta2
                                         - 32.888888888888886 * zeta3)
                                + H00
                                      * (29.333333333333332
                                         - 39.33333333333333 * zeta2
                                         - 10.666666666666666 * zeta3)
                                - 77.77777777777777 * zeta3 - 8. * zeta5)
                       + x2
                             * (-478.33333333333326 - 168. * H000
                                - 53.33333333333333 * H0000 - 16. * H00001
                                - 48. * H0001 - 140. * H001
                                + 10.666666666666666 * H0010
                                - 78.66666666666666 * H01
                                - 10.666666666666666 * H010
                                + 10.666666666666666 * H0100
                                + 21.333333333333332 * H0101 - 32. * H011
                                - 10.666666666666666 * H0111 + 72. * H10
                                - 32. * H100 - 10.666666666666666 * H1000
                                + 21.333333333333332 * H1010
                                - 45.33333333333333 * H11
                                + 21.333333333333332 * H1101
                                - 10.666666666666666 * H111
                                - 10.666666666666666 * H1111
                                + 142.88888888888889 * zeta2
                                + zeta2
                                      * (40. * H000 - 34.666666666666664 * H01
                                         - 10.666666666666666 * H10
                                         - 34.666666666666664 * H11
                                         + 51.733333333333334 * zeta2)
                                + 167.11111111111111 * zeta3
                                + H00
                                      * (-100. + 74.66666666666666 * zeta2
                                         + 21.333333333333332 * zeta3)
                                + H0
                                      * (-526.
                                         + zeta2
                                               * (124.88888888888889
                                                  + 6.4 * zeta2)
                                         + 65.77777777777777 * zeta3)
                                + 16. * zeta5))
                          / x)
           + CF * nf
                 * (LQm3
                        * (-12.444444444444443 - 2.6666666666666665 * H00
                           - 1.7777777777777777 * H1 + 1.1851851851851851 / x
                           + 21.777777777777775 * x
                           + (5.333333333333333 * H00
                              + H1
                                    * (3.5555555555555554
                                       - 3.5555555555555554 * x)
                              - 9.185185185185185 * x)
                                 * x
                           + H0 * (-4.888888888888888 + 1.7777777777777777 * x))
                    + LQm2
                          * (184.33333333333331 + 126.88888888888887 * H0
                             + 41.33333333333333 * H00 + 24. * H000 + 8. * H001
                             + 17.333333333333332 * H01
                             + 10.666666666666666 * H10 + 8. * H11
                             + 9.481481481481481 / x
                             - (3.555555555555556 * H1) / x
                             + H1
                                   * (49.77777777777777
                                      + x
                                            * (-95.55555555555556
                                               + 57.77777777777777 * x))
                             + x2
                                   * (79.85185185185185 + 54.22222222222222 * H0
                                      - 10.666666666666666 * H00
                                      + 5.333333333333333 * H01
                                      + 21.333333333333332 * H10 + 16. * H11
                                      - 5.333333333333333 * zeta2)
                             - 17.333333333333332 * zeta2 - 8. * H0 * zeta2
                             - 8. * zeta3
                             + x
                                   * (-282.66666666666663
                                      - 74.66666666666666 * H00 - 48. * H000
                                      - 16. * H001 - 10.666666666666666 * H01
                                      - 21.333333333333332 * H10 - 16. * H11
                                      + 10.666666666666666 * zeta2
                                      + H0 * (-52.44444444444444 + 16. * zeta2)
                                      + 16. * zeta3))
                    + Lmmu2
                          * (LQm
                                 * (-4 * (9 + 3 * H0 + 2 * H00)
                                    + 3.5555555555555554 / x
                                    + 4 * (15 + 4 * H00) * x
                                    + (-27.555555555555557
                                       + 10.666666666666666 * H0)
                                          * x2)
                             + (1.7777777777777777
                                + H1
                                      * (-3.5555555555555554
                                         + x
                                               * (36.
                                                  + x
                                                        * (-60.00000000000001
                                                           + 27.555555555555554
                                                                 * x)))
                                + x
                                      * (87.55555555555556 + 12. * H00
                                         + 8. * H000 + 8. * H001 + 12. * H01
                                         - 214.2222222222222 * x - 12. * zeta2
                                         + H0
                                               * (48. - 8. * zeta2
                                                  + x
                                                        * (-68.
                                                           - 15.11111111111111
                                                                 * x
                                                           + 16. * zeta2))
                                         - 8. * zeta3
                                         + x
                                               * (-16. * H000 - 16. * H001
                                                  + H00
                                                        * (-48.
                                                           - 10.666666666666666
                                                                 * x)
                                                  + x
                                                        * (124.8888888888889
                                                           - 10.666666666666666
                                                                 * H01
                                                           + 10.666666666666666
                                                                 * zeta2)
                                                  + 16. * zeta3)))
                                   / x)
                    + (LQm
                       * ((34.330864197530865 - 0.3555555555555556 * H0
                           - 18.962962962962962 * H1 + 7.111111111111111 * H10
                           + 7.111111111111111 * H11)
                              * x
                          + Hm10
                                * (0.35555555555555557
                                   + x
                                         * (-7.111111111111111
                                            + x
                                                  * (80.
                                                     + x
                                                           * (14.222222222222221
                                                              + x
                                                                    * (-55.11111111111111
                                                                       + 12.8
                                                                             * x
                                                                    )))))
                          + x
                                * (-7.111111111111111 * zeta2
                                   + x4 * (-12.8 * H00 + 12.8 * zeta2)
                                   + x2
                                         * (1572.9037037037037
                                            + 332.88888888888886 * H00
                                            + 432. * H000 + 176. * H0000
                                            + 96. * H0001 + 160. * H001
                                            + 32. * H0010 + 32. * H0011
                                            + 138.66666666666666 * H01
                                            + 21.333333333333332 * H010
                                            + 21.333333333333332 * H011
                                            - 64. * H0m10
                                            + 625.7777777777777 * H1
                                            + 251.55555555555554 * H10
                                            + 42.666666666666664 * H100
                                            + 53.33333333333333 * H101
                                            + 224.88888888888889 * H11
                                            + 32. * H110 + 32. * H111
                                            + 21.333333333333332 * Hm100
                                            - 42.666666666666664 * Hm1m10
                                            - 124.44444444444444 * zeta2
                                            + (-96. * H00 - 32. * H1
                                               - 21.333333333333332 * Hm1
                                               - 3.2 * zeta2)
                                                  * zeta2
                                            + H0
                                                  * (480.62222222222226
                                                     - 224. * zeta2
                                                     - 64. * zeta3)
                                            - 224. * zeta3)
                                   + x3
                                         * (-449.8864197530864
                                            - 99.55555555555554 * H00
                                            + 64. * H000
                                            + 10.666666666666666 * H001
                                            - 142.22222222222223 * H01
                                            - 10.666666666666666 * H010
                                            - 10.666666666666666 * H011
                                            + 42.666666666666664 * H0m10
                                            - 228.14814814814812 * H1
                                            - 176. * H10
                                            - 53.33333333333333 * H100
                                            - 53.33333333333333 * H101
                                            - 149.33333333333331 * H11
                                            - 32. * H110 - 32. * H111
                                            + 10.666666666666666 * Hm100
                                            - 21.333333333333332 * Hm1m10
                                            + H0
                                                  * (-176.35555555555555
                                                     - 10.666666666666666
                                                           * zeta2)
                                            + 142.22222222222223 * zeta2
                                            + 42.666666666666664 * H1 * zeta2
                                            - 10.666666666666666 * Hm1 * zeta2
                                            + 42.666666666666664 * zeta3)
                                   + x
                                         * (-1135.125925925926
                                            - 472.88888888888886 * H00
                                            - 140. * H000 - 88. * H0000
                                            - 48. * H0001 - 88. * H001
                                            - 16. * H0010 - 16. * H0011
                                            - 270.66666666666663 * H01
                                            - 34.666666666666664 * H010
                                            - 34.666666666666664 * H011
                                            + 53.33333333333333 * H0m10
                                            - 388.66666666666663 * H1
                                            - 124.44444444444444 * H10
                                            - 21.333333333333332 * H100
                                            - 26.666666666666664 * H101
                                            - 116.44444444444444 * H11
                                            - 16. * H110 - 16. * H111
                                            + 10.666666666666666 * Hm100
                                            - 21.333333333333332 * Hm1m10
                                            + 270.66666666666663 * zeta2
                                            + zeta2
                                                  * (48. * H00 + 16. * H1
                                                     - 10.666666666666666 * Hm1
                                                     + 1.6 * zeta2)
                                            + 117.33333333333333 * zeta3
                                            + H0
                                                  * (-674.9333333333333
                                                     + 88. * zeta2
                                                     + 32. * zeta3)))))
                          / x2
                    + Lmmu
                          * (LQm2
                                 * (-36.666666666666664
                                    - 13.333333333333332 * H0 - 8. * H00
                                    - 2.6666666666666665 * H1
                                    + 3.5555555555555554 / x
                                    + (62.666666666666664
                                       + 2.6666666666666665 * H0 + 16. * H00
                                       + 5.333333333333333 * H1)
                                          * x
                                    + (-27.555555555555557
                                       + 5.333333333333333 * H0
                                       - 5.333333333333333 * H1)
                                          * x2)
                             + LQm
                                   * (355.1111111111111
                                      + 242.66666666666666 * H0
                                      + 77.33333333333333 * H00 + 48. * H000
                                      + 16. * H001 + 32. * H01
                                      + 10.666666666666666 * H10
                                      + 10.666666666666666 * H11
                                      + 18.962962962962962 / x
                                      - (7.111111111111111 * H1) / x
                                      + (167.7037037037037
                                         + 74.66666666666666 * H0
                                         - 42.666666666666664 * H00
                                         + 21.333333333333332 * H10
                                         + 21.333333333333332 * H11)
                                            * x2
                                      + H1
                                            * (81.33333333333333
                                               + x
                                                     * (-152.
                                                        + 81.77777777777777 * x)
                                            )
                                      - 32. * zeta2 - 16. * H0 * zeta2
                                      - 16. * zeta3
                                      + x
                                            * (-547.1111111111111
                                               - 138.66666666666666 * H00
                                               - 96. * H000 - 32. * H001
                                               - 16. * H01
                                               - 21.333333333333332 * H10
                                               - 21.333333333333332 * H11
                                               + 16. * zeta2
                                               + H0
                                                     * (-88.00000000000001
                                                        + 32. * zeta2)
                                               + 32. * zeta3))
                             + ((8.05925925925926 - 0.3555555555555556 * H0
                                 - 18.962962962962962 * H1
                                 + 14.222222222222221 * H10
                                 + 7.111111111111111 * H11)
                                    * x
                                + Hm10
                                      * (0.3555555555555556
                                         + x
                                               * (-7.111111111111111
                                                  + x
                                                        * (80.
                                                           + x
                                                                 * (14.222222222222221
                                                                    + x
                                                                          * (-55.111111111111114
                                                                             + 12.8
                                                                                   * x
                                                                          )))))
                                + x
                                      * (-7.111111111111111 * zeta2
                                         + x4 * (-12.8 * H00 + 12.8 * zeta2)
                                         + x2
                                               * (1086.9037037037037
                                                  + 393.7777777777778 * H00
                                                  + 394.66666666666663 * H000
                                                  + 96. * H0000 + 96. * H0001
                                                  + 149.33333333333331 * H001
                                                  + 64. * H0010 + 32. * H0011
                                                  + 125.33333333333333 * H01
                                                  + 21.333333333333332 * H010
                                                  + 26.666666666666664 * H011
                                                  - 64. * H0m10
                                                  + 568.4444444444445 * H1
                                                  + 277.3333333333333 * H10
                                                  + 32. * H101
                                                  + 178.66666666666666 * H11
                                                  + 21.333333333333332 * H110
                                                  + 32. * H111
                                                  + 21.333333333333332 * Hm100
                                                  - 42.666666666666664 * Hm1m10
                                                  + H0
                                                        * (638.4000000000001
                                                           - 213.33333333333331
                                                                 * zeta2)
                                                  - 111.1111111111111 * zeta2
                                                  + zeta2
                                                        * (-96. * H00
                                                           - 10.666666666666666
                                                                 * H1
                                                           - 21.333333333333332
                                                                 * Hm1
                                                           + 35.2 * zeta2)
                                                  - 218.66666666666666 * zeta3)
                                         + x
                                               * (-718.0148148148148
                                                  - 290.66666666666663 * H00
                                                  - 77.33333333333333 * H000
                                                  - 48. * H0000 - 48. * H0001
                                                  - 82.66666666666666 * H001
                                                  - 32. * H0010 - 16. * H0011
                                                  - 248. * H01
                                                  - 58.666666666666664 * H010
                                                  - 37.33333333333333 * H011
                                                  + 53.33333333333333 * H0m10
                                                  - 352.4444444444444 * H1
                                                  - 160. * H10 - 16. * H101
                                                  - 86.66666666666666 * H11
                                                  - 10.666666666666666 * H110
                                                  - 16. * H111
                                                  + 10.666666666666666 * Hm100
                                                  - 21.333333333333332 * Hm1m10
                                                  + 248. * zeta2
                                                  + (48. * H00
                                                     + 5.333333333333333 * H1
                                                     - 10.666666666666666 * Hm1
                                                     - 17.6 * zeta2)
                                                        * zeta2
                                                  + H0
                                                        * (-454.9333333333334
                                                           + 82.66666666666666
                                                                 * zeta2)
                                                  + 66.66666666666666 * zeta3)
                                         + x3
                                               * (-370.2814814814815 + 88. * H00
                                                  + 42.666666666666664 * H000
                                                  + 32. * H001 - 104. * H01
                                                  + 21.333333333333332 * H010
                                                  - 10.666666666666666 * H011
                                                  + 42.666666666666664 * H0m10
                                                  - 194.37037037037035 * H1
                                                  - 144.88888888888889 * H10
                                                  - 10.666666666666666 * H100
                                                  - 32. * H101
                                                  - 111.1111111111111 * H11
                                                  - 21.333333333333332 * H110
                                                  - 32. * H111
                                                  + 10.666666666666666 * Hm100
                                                  - 21.333333333333332 * Hm1m10
                                                  + H0
                                                        * (-306.72592592592594
                                                           - 32. * zeta2)
                                                  + 104. * zeta2
                                                  + 21.333333333333332 * H1
                                                        * zeta2
                                                  - 10.666666666666666 * Hm1
                                                        * zeta2
                                                  + 85.33333333333333 * zeta3)))
                                   / x2)
                    + (-41.98628257887517
                       + (-467.4989711934156 - 386.0411522633745 * H0
                          - 172.71604938271605 * H00 - 125.7037037037037 * H000
                          - 41.77777777777778 * H0000 - 32. * H00000
                          + 32. * H00010 + 48. * H0010 + 32. * H00100
                          - 13.975308641975309 * H01 + 144.44444444444443 * H010
                          + 53.33333333333333 * H0100
                          + 5.333333333333333 * H0101 + 11.25925925925926 * H011
                          - 4.444444444444445 * H0111 + 5.843621399176955 * H1)
                             * x
                       + H10
                             * (15.407407407407407
                                + x
                                      * (158.34567901234567
                                         + (-91.06172839506172
                                            - 86.5679012345679 * x)
                                               * x))
                       + H100
                             * (-14.222222222222221
                                + x
                                      * (133.92592592592592
                                         + x
                                               * (-194.96296296296293
                                                  + 70.51851851851852 * x)))
                       + 9.481481481481483 * zeta3
                       + x
                             * (5.333333333333333 * H1101
                                + 11.25925925925926 * H111
                                - 4.444444444444445 * H1111
                                + (208.1748971193416 - 90.41152263374487 * H0
                                   - 150.04938271604937 * H00
                                   - 25.037037037037035 * H000
                                   + 3.5555555555555554 * H0000 + 64. * H00000
                                   - 64. * H00010 - 32. * H0010 - 64. * H00100
                                   + 22.320987654320987 * H01
                                   + 57.77777777777777 * H010
                                   - 10.666666666666666 * H0100
                                   - 10.666666666666666 * H0101
                                   + 18.37037037037037 * H011
                                   + 8.88888888888889 * H0111
                                   + 4.559670781893004 * H1)
                                      * x
                                - 10.666666666666666 * H1101 * x
                                - 2.962962962962963 * H111 * x
                                + 8.88888888888889 * H1111 * x
                                + 288.5418381344307 * x2
                                - 386.37037037037027 * H0 * x2
                                + 285.9753086419753 * H00 * x2
                                - 147.55555555555554 * H000 * x2
                                + 39.111111111111114 * H0000 * x2
                                - 42.666666666666664 * H0010 * x2
                                - 26.32098765432099 * H01 * x2
                                + 78.22222222222223 * H010 * x2
                                - 42.666666666666664 * H0100 * x2
                                + 10.666666666666666 * H0101 * x2
                                - 2.3703703703703707 * H011 * x2
                                - 8.88888888888889 * H0111 * x2
                                - 11.522633744855968 * H1 * x2
                                + 10.666666666666666 * H1101 * x2
                                - 2.3703703703703707 * H111 * x2
                                - 8.88888888888889 * H1111 * x2
                                + H1000
                                      * (-8.88888888888889
                                         + (17.77777777777778
                                            - 17.77777777777778 * x)
                                               * x)
                                - 2.0246913580246932 * zeta2
                                - 5.333333333333333 * H01 * zeta2
                                + 9.679012345679013 * x * zeta2
                                + 10.666666666666666 * H01 * x * zeta2
                                - 5.679012345679013 * x2 * zeta2
                                - 10.666666666666666 * H01 * x2 * zeta2
                                - 0.8888888888888888 * zeta2_2
                                - 36.62222222222223 * x * zeta2_2
                                + 11.022222222222222 * x2 * zeta2_2
                                + H11
                                      * (-11.30864197530864
                                         - 5.333333333333333 * zeta2
                                         + x
                                               * (34.32098765432099
                                                  - 26.32098765432099 * x
                                                  + 10.666666666666666 * zeta2
                                                  - 10.666666666666666 * x
                                                        * zeta2))
                                + 183.4074074074074 * zeta3
                                + 67.55555555555556 * H0 * zeta3
                                + 42.666666666666664 * H00 * zeta3
                                + 7.111111111111111 * H1 * zeta3
                                + 250.07407407407405 * x * zeta3
                                - 71.11111111111111 * H0 * x * zeta3
                                - 85.33333333333333 * H00 * x * zeta3
                                - 14.222222222222221 * H1 * x * zeta3
                                + 85.33333333333333 * x2 * zeta3
                                - 42.666666666666664 * H0 * x2 * zeta3
                                + 14.222222222222221 * H1 * x2 * zeta3
                                - 64. * zeta5 + 128. * x * zeta5))
                          / x)
           + CF * CF
                 * (-38.5 - 10. * H000 - 12. * H0000
                    + 7.999999999999999 * H00000 + 12. * H00001 + 12. * H0001
                    - 7.999999999999999 * H00010 - 37. * H001 + 12. * H0010
                    - 24. * H00100 - 48. * H00101 + 24. * H0011
                    - 31.999999999999996 * H00110 - 24. * H00111 + 28. * H01
                    + 3.9999999999999996 * H010 - 20. * H0100
                    - 15.999999999999998 * H01000 - 72. * H01001 - 28. * H0101
                    - 96. * H01010 - 120. * H01011 + 88. * H011 - 36. * H0110
                    - 40. * H01100 - 104. * H01101 - 44. * H0111 - 56. * H01110
                    - 56. * H01111 + 39. * H1 + 50. * H10 + 28. * H100
                    + 48. * H1000 + 31.999999999999996 * H10000
                    + 31.999999999999996 * H10001 + 60. * H1001
                    - 15.999999999999998 * H10010 + 14. * H101 + 24. * H1010
                    - 15.999999999999998 * H10100 - 80. * H10101 + 12. * H1011
                    - 48. * H10110 - 48. * H10111 + 148. * H11 + 14. * H110
                    + 60. * H1100 + 31.999999999999996 * H11000
                    - 31.999999999999996 * H11001 + 3.9999999999999996 * H1101
                    - 80. * H11010 - 112. * H11011 + 124. * H111 + 12. * H1110
                    - 112. * H11101 - 76. * H1111 - 48. * H11110 - 80. * H11111
                    - 42. * x - 3.9999999999999996 * H000 * x
                    - 63.99999999999999 * H0000 * x
                    - 15.999999999999998 * H00000 * x - 24. * H00001 * x
                    - 52. * H0001 * x + 15.999999999999998 * H00010 * x
                    - 10. * H001 * x - 48. * H0010 * x + 48. * H00100 * x
                    + 96. * H00101 * x + 31.999999999999996 * H0011 * x
                    + 63.99999999999999 * H00110 * x + 48. * H00111 * x
                    - 305. * H01 * x + 15.999999999999998 * H010 * x
                    - 31.999999999999996 * H0100 * x
                    + 31.999999999999996 * H01000 * x + 144. * H01001 * x
                    + 120. * H0101 * x + 192. * H01010 * x + 240. * H01011 * x
                    - 31.999999999999996 * H011 * x + 120. * H0110 * x
                    + 80. * H01100 * x + 208. * H01101 * x + 208. * H0111 * x
                    + 112. * H01110 * x + 112. * H01111 * x - 313. * H1 * x
                    - 312. * H10 * x + 72. * H100 * x - 272. * H1000 * x
                    - 63.99999999999999 * H10000 * x
                    - 63.99999999999999 * H10001 * x - 224. * H1001 * x
                    + 31.999999999999996 * H10010 * x + 172. * H101 * x
                    - 127.99999999999999 * H1010 * x
                    + 31.999999999999996 * H10100 * x + 160. * H10101 * x
                    + 31.999999999999996 * H1011 * x + 96. * H10110 * x
                    + 96. * H10111 * x - 258. * H11 * x + 172. * H110 * x
                    - 224. * H1100 * x - 63.99999999999999 * H11000 * x
                    + 63.99999999999999 * H11001 * x
                    + 63.99999999999999 * H1101 * x + 160. * H11010 * x
                    + 224. * H11011 * x + 112. * H111 * x
                    + 31.999999999999996 * H1110 * x + 224. * H11101 * x
                    + 304. * H1111 * x + 96. * H11110 * x + 160. * H11111 * x
                    - 39. * x2 - 72. * H000 * x2 + 208. * H0000 * x2
                    + 63.99999999999999 * H00000 * x2
                    + 63.99999999999999 * H00001 * x2 + 160. * H0001 * x2
                    - 31.999999999999996 * H00010 * x2 - 176. * H001 * x2
                    + 112. * H0010 * x2 - 31.999999999999996 * H00100 * x2
                    - 160. * H00101 * x2 - 63.99999999999999 * H0011 * x2
                    - 96. * H00110 * x2 - 96. * H00111 * x2 + 272. * H01 * x2
                    - 176. * H010 * x2 + 160. * H0100 * x2
                    + 63.99999999999999 * H01000 * x2
                    - 63.99999999999999 * H01001 * x2
                    - 63.99999999999999 * H0101 * x2 - 160. * H01010 * x2
                    - 224. * H01011 * x2 - 248. * H011 * x2
                    - 63.99999999999999 * H0110 * x2 - 224. * H01101 * x2
                    - 288. * H0111 * x2 - 96. * H01110 * x2 - 160. * H01111 * x2
                    + 272. * H1 * x2 + 384. * H10 * x2 - 72. * H100 * x2
                    + 208. * H1000 * x2 + 63.99999999999999 * H10000 * x2
                    + 63.99999999999999 * H10001 * x2 + 160. * H1001 * x2
                    - 31.999999999999996 * H10010 * x2 - 176. * H101 * x2
                    + 112. * H1010 * x2 - 31.999999999999996 * H10100 * x2
                    - 160. * H10101 * x2 - 63.99999999999999 * H1011 * x2
                    - 96. * H10110 * x2 - 96. * H10111 * x2 + 272. * H11 * x2
                    - 176. * H110 * x2 + 160. * H1100 * x2
                    + 63.99999999999999 * H11000 * x2
                    - 63.99999999999999 * H11001 * x2
                    - 63.99999999999999 * H1101 * x2 - 160. * H11010 * x2
                    - 224. * H11011 * x2 - 248. * H111 * x2
                    - 63.99999999999999 * H1110 * x2 - 224. * H11101 * x2
                    - 288. * H1111 * x2 - 96. * H11110 * x2 - 160. * H11111 * x2
                    + 93.75 * zeta2 - 6. * H000 * zeta2 + 68. * H001 * zeta2
                    + 40. * H01 * zeta2 + 100. * H010 * zeta2
                    + 144. * H011 * zeta2 - 7.999999999999999 * H0m10 * zeta2
                    - 20. * H1 * zeta2 - 54. * H10 * zeta2
                    - 15.999999999999998 * H100 * zeta2 + 120. * H101 * zeta2
                    + 72. * H110 * zeta2 + 160. * H111 * zeta2
                    - 15.999999999999998 * Hm10 * zeta2
                    - 7.999999999999999 * Hm100 * zeta2
                    + 15.999999999999998 * Hm1m10 * zeta2 - 192. * ln2 * zeta2
                    + 127.49999999999999 * x * zeta2 + 28. * H000 * x * zeta2
                    - 136. * H001 * x * zeta2 - 164. * H01 * x * zeta2
                    - 200. * H010 * x * zeta2 - 288. * H011 * x * zeta2
                    - 15.999999999999998 * H0m10 * x * zeta2
                    - 164. * H1 * x * zeta2 + 168. * H10 * x * zeta2
                    + 31.999999999999996 * H100 * x * zeta2
                    - 240. * H101 * x * zeta2 - 144. * H11 * x * zeta2
                    - 144. * H110 * x * zeta2 - 320. * H111 * x * zeta2
                    - 40. * Hm10 * x * zeta2
                    - 15.999999999999998 * Hm100 * x * zeta2
                    + 31.999999999999996 * Hm1m10 * x * zeta2
                    + 384. * ln2 * x * zeta2 - 58. * x2 * zeta2
                    - 15.999999999999998 * H000 * x2 * zeta2
                    + 240. * H001 * x2 * zeta2
                    + 127.99999999999999 * H01 * x2 * zeta2
                    + 144. * H010 * x2 * zeta2 + 320. * H011 * x2 * zeta2
                    - 31.999999999999996 * H0m10 * x2 * zeta2
                    + 122.00000000000001 * H1 * x2 * zeta2
                    - 127.99999999999999 * H10 * x2 * zeta2
                    - 31.999999999999996 * H100 * x2 * zeta2
                    + 240. * H101 * x2 * zeta2
                    + 127.99999999999999 * H11 * x2 * zeta2
                    + 144. * H110 * x2 * zeta2 + 320. * H111 * x2 * zeta2
                    - 24. * Hm10 * x2 * zeta2
                    - 15.999999999999998 * Hm100 * x2 * zeta2
                    + 31.999999999999996 * Hm1m10 * x2 * zeta2
                    - 384. * ln2 * x2 * zeta2 - 0.8 * zeta2_2
                    - 105.60000000000001 * H1 * zeta2_2
                    + 7.999999999999999 * Hm1 * zeta2_2
                    + 63.20000000000001 * x * zeta2_2
                    + 211.20000000000002 * H1 * x * zeta2_2
                    + 15.999999999999998 * Hm1 * x * zeta2_2
                    - 140.8 * x2 * zeta2_2
                    - 211.20000000000002 * H1 * x2 * zeta2_2
                    + 15.999999999999998 * Hm1 * x2 * zeta2_2
                    + LQm3
                          * (3.6666666666666665 + 5.333333333333333 * H01
                             + 2.6666666666666665 * H1 + 5.333333333333333 * H10
                             + 10.666666666666666 * H11 - 0.6666666666666666 * x
                             - 4. * H0 * x
                             + H00
                                   * (1.3333333333333333
                                      + x
                                            * (-2.6666666666666665
                                               + 10.666666666666666 * x))
                             - 5.333333333333333 * zeta2
                             + x
                                   * (-10.666666666666666 * H1
                                      - 10.666666666666666 * H10
                                      - 21.333333333333332 * H11
                                      + 10.666666666666666 * H10 * x
                                      + 21.333333333333332 * H11 * x
                                      + H01
                                            * (-10.666666666666666
                                               + 21.333333333333332 * x)
                                      + 10.666666666666666 * zeta2
                                      - 21.333333333333332 * x * zeta2))
                    + 8.666666666666666 * zeta3
                    + 5.333333333333333 * H01 * zeta3
                    - 21.333333333333332 * H1 * zeta3
                    - 58.66666666666666 * H10 * zeta3
                    - 5.333333333333333 * H11 * zeta3
                    - 22.666666666666664 * x * zeta3
                    - 10.666666666666666 * H01 * x * zeta3
                    - 74.66666666666666 * H1 * x * zeta3
                    + 117.33333333333331 * H10 * x * zeta3
                    + 10.666666666666666 * H11 * x * zeta3 + 168. * x2 * zeta3
                    - 10.666666666666666 * H01 * x2 * zeta3
                    + 127.99999999999999 * H1 * x2 * zeta3
                    - 117.33333333333331 * H10 * x2 * zeta3
                    - 10.666666666666666 * H11 * x2 * zeta3
                    + 2.6666666666666665 * zeta2 * zeta3
                    - 37.33333333333333 * x * zeta2 * zeta3
                    - 53.33333333333333 * x2 * zeta2 * zeta3
                    + LQm2
                          * (-16.5 - 36. * H001 - 24. * H01 - 48. * H010
                             - 64. * H011 + 16. * H0m10 - 65. * H1 - 40. * H10
                             - 32. * H100 - 80. * H101 - 64. * H11 - 80. * H110
                             - 96. * H111 + 32. * Hm10 + 16. * Hm100
                             - 32. * Hm1m10 - 46. * x - 36. * H00 * x
                             + 72. * H001 * x + 124. * H01 * x + 96. * H010 * x
                             + 128. * H011 * x + 32. * H0m10 * x + 126. * H1 * x
                             + 144. * H10 * x + 64. * H100 * x + 160. * H101 * x
                             + 224. * H11 * x + 160. * H110 * x
                             + 192. * H111 * x + 80. * Hm10 * x
                             + 32. * Hm100 * x - 64. * Hm1m10 * x + 52. * x2
                             - 128. * H00 * x2 - 160. * H001 * x2
                             - 160. * H01 * x2 - 160. * H010 * x2
                             - 192. * H011 * x2 + 64. * H0m10 * x2
                             - 36. * H1 * x2 - 80. * H10 * x2 - 64. * H100 * x2
                             - 160. * H101 * x2 - 160. * H11 * x2
                             - 160. * H110 * x2 - 192. * H111 * x2
                             + 48. * Hm10 * x2 + 32. * Hm100 * x2
                             - 64. * Hm1m10 * x2
                             + H000 * (-12. + (-8. - 96. * x) * x) + 24. * zeta2
                             + 64. * H1 * zeta2 - 16. * Hm1 * zeta2
                             - 44. * x * zeta2 - 128. * H1 * x * zeta2
                             - 32. * Hm1 * x * zeta2 + 160. * x2 * zeta2
                             + 128. * H1 * x2 * zeta2 - 32. * Hm1 * x2 * zeta2
                             + H0
                                   * (-22. + 36. * zeta2
                                      + x
                                            * (12. - 40. * zeta2
                                               + x * (-36. + 160. * zeta2)))
                             + 28. * zeta3 + 8. * x * zeta3 + 128. * x2 * zeta3)
                    + H00
                          * (63. - 7.999999999999999 * zeta2
                             - 26.666666666666664 * zeta3
                             + x
                                   * (-196. + 74. * zeta2
                                      + x
                                            * (384. - 104. * zeta2
                                               - 117.33333333333331 * zeta3)
                                      + 53.33333333333333 * zeta3))
                    + H0
                          * (-26. + (29. - 45.6 * zeta2) * zeta2
                             + x
                                   * (-182. + zeta2 * (38. + 75.2 * zeta2)
                                      - 80. * zeta3)
                             - 12. * zeta3
                             + x2
                                   * (272.
                                      + (122.00000000000001
                                         - 227.20000000000002 * zeta2)
                                            * zeta2
                                      + 127.99999999999999 * zeta3))
                    + (LQm
                       * (Hm10
                              * (2.1333333333333333
                                 + x2
                                       * (-976. + 32. * zeta2
                                          + x
                                                * (-1162.6666666666665
                                                   + 64. * zeta2
                                                   + x
                                                         * (-176.00000000000003
                                                            + 76.80000000000001
                                                                  * x
                                                            + 64. * zeta2))))
                          + x
                                * (2.1333333333333333
                                   + x4
                                         * (-76.80000000000001 * H00
                                            + 76.80000000000001 * zeta2)
                                   + x
                                         * (-25.533333333333328 + 78. * H00
                                            + 40. * H0000 + 120. * H0001
                                            + 72. * H001 + 160. * H0010
                                            + 199.99999999999997 * H0011
                                            - 192. * H00m10 + 128. * H01
                                            + 128. * H010 + 80. * H0100
                                            + 240. * H0101 + 152. * H011
                                            + 224.00000000000003 * H0110
                                            + 288. * H0111
                                            - 448.00000000000006 * H0m10
                                            - 160. * H0m100 - 32. * H0m101
                                            + 160. * H0m1m10 - 85. * H1
                                            + 231.99999999999997 * H10
                                            - 24. * H100 + 96. * H1000
                                            + 256. * H1001 + 264. * H101
                                            + 288. * H1010
                                            + 352.00000000000006 * H1011
                                            - 128. * H10m10 + 310. * H11
                                            + 256. * H110 + 160. * H1100
                                            + 320. * H1101 + 288. * H111
                                            + 288. * H1110 + 384. * H1111
                                            - 384. * Hm100 - 64. * Hm1000
                                            - 32. * Hm1001 - 64. * Hm101
                                            + 128. * Hm10m10
                                            + 448.00000000000006 * Hm1m10
                                            + 224.00000000000003 * Hm1m100
                                            + 64. * Hm1m101 - 128. * Hm1m1m10
                                            - 128. * zeta2 - 120. * H00 * zeta2
                                            - 160. * H01 * zeta2
                                            + 112.00000000000001 * H0m1 * zeta2
                                            - 40. * H1 * zeta2
                                            - 256. * H10 * zeta2
                                            - 256. * H11 * zeta2
                                            + 288. * Hm1 * zeta2
                                            - 128. * Hm1m1 * zeta2
                                            + 45.6 * zeta2_2 - 184. * zeta3
                                            - 239.99999999999997 * H1 * zeta3
                                            + 112.00000000000001 * Hm1 * zeta3)
                                   + x3
                                         * (-72.8 + 616. * H00 + 656. * H000
                                            + 320. * H0000 + 576. * H0001
                                            + 799.9999999999999 * H001
                                            + 576. * H0010
                                            + 704.0000000000001 * H0011
                                            - 512. * H00m10 + 736. * H01
                                            + 704.0000000000001 * H010
                                            + 320. * H0100 + 640. * H0101
                                            + 1008. * H011 + 576. * H0110
                                            + 768. * H0111 - 32. * H0m10
                                            - 448.00000000000006 * H0m100
                                            - 128. * H0m101 + 256. * H0m1m10
                                            + 272. * H1
                                            + 440.0000000000001 * H10
                                            + 399.99999999999994 * H100
                                            + 192. * H1000 + 512. * H1001
                                            + 704.0000000000001 * H101
                                            + 576. * H1010
                                            + 704.0000000000001 * H1011
                                            - 256. * H10m10 + 736. * H11
                                            + 704.0000000000001 * H110
                                            + 320. * H1100 + 640. * H1101
                                            + 1008. * H111 + 576. * H1110
                                            + 768. * H1111 - 256. * Hm100
                                            - 128. * Hm1000 - 64. * Hm1001
                                            - 96. * Hm101 + 256. * Hm10m10
                                            + 32. * Hm1m10
                                            + 448.00000000000006 * Hm1m100
                                            + 128. * Hm1m101 - 256. * Hm1m1m10
                                            - 736. * zeta2 - 576. * H00 * zeta2
                                            - 512. * H01 * zeta2
                                            + 256. * H0m1 * zeta2
                                            - 688. * H1 * zeta2
                                            - 512. * H10 * zeta2
                                            - 512. * H11 * zeta2
                                            + 112.00000000000001 * Hm1 * zeta2
                                            - 256. * Hm1m1 * zeta2
                                            + 153.60000000000002 * zeta2_2
                                            - 1024. * zeta3
                                            - 479.99999999999994 * H1 * zeta3
                                            + 224.00000000000003 * Hm1 * zeta3)
                                   + x2
                                         * (126.2 + 734.6666666666666 * H00
                                            + 231.99999999999997 * H000
                                            + 48. * H0000
                                            - 176.00000000000003 * H0001
                                            - 312. * H001 - 320. * H0010
                                            - 399.99999999999994 * H0011
                                            - 808. * H01 - 624. * H010
                                            - 416. * H0100 - 480. * H0101
                                            - 920. * H011
                                            - 448.00000000000006 * H0110
                                            - 576. * H0111 - 64. * H0m10
                                            - 64. * H0m100 - 64. * H0m101
                                            - 192. * H0m1m10 - 216. * H1
                                            - 760. * H10 - 320. * H100
                                            - 192. * H1000 - 512. * H1001
                                            - 864. * H101 - 576. * H1010
                                            - 704.0000000000001 * H1011
                                            + 256. * H10m10 - 1092. * H11
                                            - 832. * H110 - 320. * H1100
                                            - 640. * H1101 - 1152. * H111
                                            - 576. * H1110 - 768. * H1111
                                            - 640. * Hm100 - 128. * Hm1000
                                            - 64. * Hm1001 - 160. * Hm101
                                            + 256. * Hm10m10 + 480. * Hm1m10
                                            + 448.00000000000006 * Hm1m100
                                            + 128. * Hm1m101 - 256. * Hm1m1m10
                                            - 354.66666666666663 * zeta2
                                            + 176.00000000000003 * H00 * zeta2
                                            + 576. * H01 * zeta2
                                            - 32. * H0m1 * zeta2
                                            + 624. * H1 * zeta2
                                            + 512. * H10 * zeta2
                                            + 512. * H11 * zeta2
                                            + 399.99999999999994 * Hm1 * zeta2
                                            - 256. * Hm1m1 * zeta2
                                            - 276.8 * zeta2_2 + 432. * zeta3
                                            + 479.99999999999994 * H1 * zeta3
                                            + 224.00000000000003 * Hm1 * zeta3)
                                   + H0
                                         * (-2.1333333333333333
                                            + x
                                                  * (41.06666666666666
                                                     - 72. * zeta2
                                                     - 144. * zeta3
                                                     + x
                                                           * (212.4
                                                              + 248. * zeta2
                                                              + x
                                                                    * (195.20000000000002
                                                                       - 799.9999999999999
                                                                             * zeta2
                                                                       - 704.0000000000001
                                                                             * zeta3
                                                                    )
                                                              + 479.99999999999994
                                                                    * zeta3)))))
                      ) / x2
                    - 3.9999999999999996 * zeta5 + 7.999999999999999 * x * zeta5
                    - 224. * x2 * zeta5)
           + CA
                 * (LQm3
                        * (0.4444444444444444 + 0.5925925925925926 / x
                           + 3.5555555555555554 * x - 4.592592592592593 * x2
                           + H0 * (0.8888888888888888 + 3.5555555555555554 * x)
                           + H1
                                 * (-0.8888888888888888
                                    + (1.7777777777777777
                                       - 1.7777777777777777 * x)
                                          * x))
                    + LQm2
                          * (-6.444444444444444 - 5.777777777777778 * H0
                             - 2.6666666666666665 * H00
                             + 2.6666666666666665 * H10
                             + 2.6666666666666665 * H11
                             - 2.6666666666666665 * Hm10
                             - 2.9629629629629624 / x
                             - (1.7777777777777777 * H1) / x
                             + H1
                                   * (5.777777777777778
                                      + x
                                            * (-35.55555555555556
                                               + 38.666666666666664 * x))
                             + x2
                                   * (77.62962962962962
                                      + 21.777777777777775 * H0
                                      + 5.333333333333333 * H01
                                      + 5.333333333333333 * H10
                                      + 5.333333333333333 * H11
                                      - 5.333333333333333 * Hm10
                                      - 5.333333333333333 * zeta2)
                             + x
                                   * (-66.22222222222221
                                      - 60.44444444444444 * H0
                                      - 21.333333333333332 * H00
                                      - 15.999999999999998 * H01
                                      - 5.333333333333333 * H10
                                      - 5.333333333333333 * H11
                                      - 5.333333333333333 * Hm10
                                      + 10.666666666666666 * zeta2))
                    + Lmmu2
                          * (-38.22222222222222 - 5.333333333333333 * H0
                             - 10.666666666666666 * H00
                             + 10.666666666666666 * H10
                             + 21.333333333333332 * H11 + 3.555555555555556 / x
                             - (7.111111111111112 * H1) / x
                             + H1
                                   * (5.333333333333333
                                      + x * (-128. + 140.44444444444443 * x))
                             + LQm
                                   * (5.333333333333333 + 7.111111111111111 / x
                                      + 42.666666666666664 * x
                                      - 55.11111111111111 * x2
                                      + H0
                                            * (10.666666666666666
                                               + 42.666666666666664 * x)
                                      + H1
                                            * (-10.666666666666666
                                               + (21.333333333333332
                                                  - 21.333333333333332 * x)
                                                     * x))
                             + x2
                                   * (249.77777777777777
                                      + 55.11111111111111 * H0
                                      + 21.333333333333332 * H01
                                      + 21.333333333333332 * H10
                                      + 42.666666666666664 * H11
                                      - 21.333333333333332 * zeta2)
                             + x
                                   * (-215.1111111111111
                                      - 170.66666666666666 * H0
                                      - 42.666666666666664 * H00 - 64. * H01
                                      - 21.333333333333332 * H10
                                      - 42.666666666666664 * H11 + 64. * zeta2))
                    + Lmmu
                          * (LQm2
                                 * (2.6666666666666665 + 3.5555555555555554 / x
                                    + 21.333333333333332 * x
                                    - 27.555555555555554 * x2
                                    + H0
                                          * (5.333333333333333
                                             + 21.333333333333332 * x)
                                    + H1
                                          * (-5.333333333333333
                                             + (10.666666666666666
                                                - 10.666666666666666 * x)
                                                   * x))
                             - (64. * LQm
                                * (-0.2407407407407407
                                   + (0.6597222222222222
                                      + 0.3333333333333333 * H00
                                      - 0.16666666666666666 * H10
                                      - 0.16666666666666666 * H11
                                      + 0.16666666666666666 * Hm10)
                                         * x
                                   + H1
                                         * (0.1111111111111111
                                            + x
                                                  * (-0.08333333333333333
                                                     + (1.6666666666666665
                                                        - 1.861111111111111 * x)
                                                           * x))
                                   + x2
                                         * (1.486111111111111 + 2. * H0 + H00
                                            + H01 + 0.3333333333333333 * H10
                                            + 0.3333333333333333 * H11
                                            + 0.3333333333333333 * Hm10
                                            - 0.6666666666666666 * zeta2)
                                   + x3
                                         * (-2.0925925925925926
                                            - 2.083333333333333 * H0
                                            - 0.3333333333333333 * H01
                                            - 0.3333333333333333 * H10
                                            - 0.3333333333333333 * H11
                                            + 0.3333333333333333 * Hm10
                                            + 0.3333333333333333 * zeta2)))
                                   / x
                             + (5.925925925925926 + 7.111111111111111 * H11
                                - 7.111111111111111 * Hm10
                                + (27.85185185185185 + 58.22222222222223 * H0
                                   + 21.333333333333332 * H000
                                   + 10.666666666666666 * H001
                                   - 10.666666666666666 * H01
                                   + 10.666666666666666 * H010)
                                      * x
                                + H10
                                      * (14.222222222222221
                                         + x
                                               * (2.6666666666666665
                                                  + (149.33333333333331
                                                     - 176.8888888888889 * x)
                                                        * x))
                                - 7.111111111111111 * zeta2
                                + H1
                                      * (-15.407407407407407
                                         + x
                                               * (18.22222222222222
                                                  + 10.666666666666666 * zeta2
                                                  + x
                                                        * (204.44444444444446
                                                           - 21.333333333333332
                                                                 * zeta2
                                                           + x
                                                                 * (-235.25925925925924
                                                                    + 10.666666666666666
                                                                          * zeta2
                                                                 ))))
                                + x
                                      * (-8. * H11 - 21.333333333333332 * H110
                                         - 32. * Hm10 + 16. * Hm100
                                         + 10.666666666666666 * Hm101
                                         + (72.5925925925926
                                            + 234.66666666666666 * H0
                                            + 245.33333333333331 * H00
                                            + 64. * H000
                                            + 85.33333333333333 * H001
                                            + 192. * H01
                                            + 106.66666666666666 * H010
                                            + 64. * H011)
                                               * x
                                         + 85.33333333333333 * H11 * x
                                         + 42.666666666666664 * H110 * x
                                         + 10.666666666666666 * Hm10 * x
                                         + 32. * Hm100 * x
                                         + 21.333333333333332 * Hm101 * x
                                         - 118.37037037037038 * x2
                                         - 441.4814814814815 * H0 * x2
                                         - 152. * H00 * x2
                                         - 21.333333333333332 * H001 * x2
                                         - 194.66666666666666 * H01 * x2
                                         - 21.333333333333332 * H010 * x2
                                         - 21.333333333333332 * H011 * x2
                                         + 21.333333333333332 * H0m10 * x2
                                         - 95.1111111111111 * H11 * x2
                                         - 42.666666666666664 * H110 * x2
                                         + 46.22222222222223 * Hm10 * x2
                                         + 42.666666666666664 * Hm100 * x2
                                         + 21.333333333333332 * Hm101 * x2
                                         - 21.333333333333332 * Hm1m10 * x2
                                         + H101
                                               * (-10.666666666666666
                                                  + (21.333333333333332
                                                     - 21.333333333333332 * x)
                                                        * x)
                                         + H100
                                               * (-16.
                                                  + (32.
                                                     - 21.333333333333332 * x)
                                                        * x)
                                         + 10.666666666666666 * zeta2
                                         - 10.666666666666666 * H0 * zeta2
                                         - 10.666666666666666 * Hm1 * zeta2
                                         - 181.33333333333331 * x * zeta2
                                         - 85.33333333333333 * H0 * x * zeta2
                                         - 21.333333333333332 * Hm1 * x * zeta2
                                         + 194.66666666666666 * x2 * zeta2
                                         + 21.333333333333332 * H0 * x2 * zeta2
                                         - 32. * Hm1 * x2 * zeta2
                                         + 26.66666666666666 * zeta3
                                         + 10.666666666666664 * x * zeta3
                                         + 31.999999999999993 * x2 * zeta3))
                                   / x)
                    + (LQm
                       * (10.172839506172838 + 3.555555555555555 * H10
                          + 3.555555555555555 * H11 - 3.555555555555555 * Hm10
                          + (23.925925925925924 + 18.37037037037037 * H0
                             + 8.88888888888889 * H00 - 2.6666666666666665 * H01
                            ) * x
                          - 3.555555555555555 * zeta2
                          + H1
                                * (5.925925925925926
                                   + x
                                         * (-19.25925925925926
                                            + x
                                                  * (285.6296296296296
                                                     - 305.77777777777777 * x
                                                     - 5.333333333333333 * zeta2
                                                  )
                                            + 2.6666666666666665 * zeta2))
                          + x
                                * (-10.666666666666666 * H100
                                   - 11.555555555555555 * H11
                                   - 10.666666666666666 * H110
                                   - 5.333333333333333 * H111
                                   - 7.11111111111111 * Hm10
                                   + 10.666666666666666 * Hm100
                                   + 5.333333333333333 * Hm101
                                   + 5.333333333333333 * Hm1m10
                                   + (581.6296296296296 + 493.037037037037 * H0
                                      + 263.1111111111111 * H00
                                      + 74.66666666666666 * H000
                                      + 53.33333333333333 * H001 + 176. * H01
                                      + 32. * H010 + 32. * H011)
                                         * x
                                   + H10
                                         * (-11.555555555555555
                                            + (71.11111111111111
                                               - 77.33333333333333 * x)
                                                  * x)
                                   + 2.6666666666666665 * zeta2
                                   - 2.6666666666666665 * Hm1 * zeta2
                                   + 8. * zeta3
                                   + x
                                         * (71.11111111111111 * H11
                                            + 21.333333333333332 * H110
                                            + 10.666666666666666 * H111
                                            + 12.444444444444443 * Hm10
                                            + 21.333333333333332 * Hm100
                                            + 10.666666666666666 * Hm101
                                            + 10.666666666666666 * Hm1m10
                                            + H100
                                                  * (21.333333333333332
                                                     - 16. * x)
                                            + (-633.5061728395061
                                               - 192.59259259259255 * H0
                                               - 55.99999999999999 * H00
                                               - 10.666666666666666 * H001
                                               - 93.33333333333333 * H01
                                               - 10.666666666666666 * H010
                                               - 10.666666666666666 * H011
                                               + 10.666666666666666 * H0m10)
                                                  * x
                                            + (-163.55555555555554
                                               - 53.33333333333333 * H0
                                               - 5.333333333333333 * Hm1)
                                                  * zeta2
                                            - 58.666666666666664 * zeta3
                                            + x
                                                  * (-77.33333333333333 * H11
                                                     - 21.333333333333332 * H110
                                                     - 10.666666666666666 * H111
                                                     + 30.22222222222222 * Hm10
                                                     + 26.666666666666664
                                                           * Hm100
                                                     + 10.666666666666666
                                                           * Hm101
                                                     + (93.33333333333333
                                                        + 10.666666666666666
                                                              * H0
                                                        - 10.666666666666666
                                                              * Hm1)
                                                           * zeta2
                                                     + 16. * zeta3)))))
                          / x
                    + (54.50205761316874 + 7.111111111111111 * H100
                       + (-71.46913580246913 + 49.382716049382715 * H0
                          - 33.925925925925924 * H00 + 2.6666666666666665 * H000
                          - 5.333333333333333 * H0000
                          + 2.6666666666666665 * H0001
                          + 8.444444444444445 * H001
                          + 10.666666666666666 * H0010
                          - 0.8888888888888888 * H01 - 5.333333333333333 * H010
                          + 10.666666666666666 * H0100)
                             * x
                       + H10
                             * (-11.851851851851851
                                + x
                                      * (21.62962962962963
                                         + x
                                               * (-169.92592592592592
                                                  + 167.11111111111111 * x)))
                       + 14.222222222222221 * zeta2 + 4.148148148148148 * zeta3
                       + x
                             * (5.333333333333333 * H1011
                                + 9.925925925925926 * H11
                                - 5.333333333333333 * H1100
                                + 10.666666666666666 * H1101
                                + 2.6666666666666665 * H111
                                + 10.666666666666666 * H1110
                                - 5.333333333333333 * H1111
                                + 5.333333333333333 * Hm1000
                                + 10.666666666666666 * Hm1010
                                + 10.666666666666666 * Hm1m100
                                - 21.333333333333332 * Hm1m1m10
                                + (151.4320987654321 - 108.39506172839505 * H0
                                   - 87.70370370370368 * H00
                                   - 5.333333333333333 * H000
                                   - 10.666666666666666 * H0000
                                   - 5.333333333333333 * H0001
                                   - 24.888888888888886 * H001
                                   + 21.333333333333332 * H0010
                                   - 93.77777777777777 * H01
                                   - 45.33333333333333 * H010
                                   + 42.666666666666664 * H0100 - 8. * H011
                                   - 10.666666666666666 * H101)
                                      * x
                                - 10.666666666666666 * H1011 * x
                                - 21.185185185185183 * H11 * x
                                - 10.666666666666666 * H110 * x
                                + 10.666666666666666 * H1100 * x
                                - 21.333333333333332 * H1101 * x
                                + 10.666666666666666 * H111 * x
                                - 21.333333333333332 * H1110 * x
                                + 10.666666666666666 * H1111 * x
                                - 32. * Hm10 * x
                                - 10.666666666666666 * Hm100 * x
                                + 10.666666666666666 * Hm1000 * x
                                + 21.333333333333332 * Hm1010 * x
                                + 21.333333333333332 * Hm1m10 * x
                                + 21.333333333333332 * Hm1m100 * x
                                - 42.666666666666664 * Hm1m1m10 * x
                                - 124.83539094650207 * x2
                                + 402.5432098765432 * H0 * x2
                                - 121.77777777777777 * H00 * x2
                                + 20.444444444444443 * H000 * x2
                                - 32. * H01 * x2 - 78.22222222222221 * H010 * x2
                                - 2.6666666666666665 * H011 * x2
                                + 10.666666666666666 * H101 * x2
                                + 10.666666666666666 * H1011 * x2
                                + 30.518518518518515 * H11 * x2
                                + 10.666666666666666 * H110 * x2
                                - 10.666666666666666 * H1100 * x2
                                + 21.333333333333332 * H1101 * x2
                                - 13.333333333333332 * H111 * x2
                                + 21.333333333333332 * H1110 * x2
                                - 10.666666666666666 * H1111 * x2
                                - 32. * Hm10 * x2
                                - 10.666666666666666 * Hm100 * x2
                                + 10.666666666666666 * Hm1000 * x2
                                + 21.333333333333332 * Hm1010 * x2
                                + 21.333333333333332 * Hm1m10 * x2
                                + 21.333333333333332 * Hm1m100 * x2
                                - 42.666666666666664 * Hm1m1m10 * x2
                                + H100
                                      * (8.
                                         + (42.666666666666664
                                            - 57.77777777777778 * x)
                                               * x)
                                + H1010
                                      * (10.666666666666666
                                         + x
                                               * (-21.333333333333332
                                                  + 21.333333333333332 * x))
                                - 6.555555555555555 * zeta2
                                + 0.6666666666666666 * H0 * zeta2
                                - 12. * H00 * zeta2
                                - 6.666666666666666 * H11 * zeta2
                                - 6.666666666666666 * Hm10 * zeta2
                                - 10.666666666666666 * Hm1m1 * zeta2
                                + 180.88888888888886 * x * zeta2
                                + 80. * H0 * x * zeta2
                                - 2.6666666666666665 * H00 * x * zeta2
                                + 13.333333333333332 * H11 * x * zeta2
                                + 10.666666666666666 * Hm1 * x * zeta2
                                - 13.333333333333332 * Hm10 * x * zeta2
                                - 21.333333333333332 * Hm1m1 * x * zeta2
                                - 86.22222222222221 * x2 * zeta2
                                + 60.44444444444444 * H0 * x2 * zeta2
                                - 13.333333333333332 * H11 * x2 * zeta2
                                + 10.666666666666666 * Hm1 * x2 * zeta2
                                - 13.333333333333332 * Hm10 * x2 * zeta2
                                - 21.333333333333332 * Hm1m1 * x2 * zeta2
                                - 1.0666666666666667 * zeta2_2
                                - 41.06666666666666 * x * zeta2_2 - 16. * zeta3
                                + ((-22.222222222222218 - 185.92592592592592 * x
                                   ) * x
                                   + H0
                                         * (24.888888888888886
                                            + 72.88888888888889 * x)
                                   + Hm1 * (16. + x * (32. + 32. * x)))
                                      * zeta3)
                       + H1
                             * (-15.506172839506172
                                + x
                                      * (2.617283950617284
                                         - 4.444444444444445 * zeta2
                                         + 4.444444444444445 * zeta3
                                         + x
                                               * (-169.82716049382717
                                                  + 22.222222222222218 * zeta2
                                                  - 8.88888888888889 * zeta3
                                                  + x
                                                        * (168.37037037037032
                                                           - 22.222222222222218
                                                                 * zeta2
                                                           + 8.88888888888889
                                                                 * zeta3)))))
                          / x)
           + CA * nf
                 * (LQm3
                        * (0.4444444444444444 + 0.5925925925925926 / x
                           + 3.5555555555555554 * x - 4.592592592592593 * x2
                           + H0 * (0.8888888888888888 + 3.5555555555555554 * x)
                           + H1
                                 * (-0.8888888888888888
                                    + (1.7777777777777777
                                       - 1.7777777777777777 * x)
                                          * x))
                    + LQm2
                          * (-6.444444444444444 - 5.777777777777778 * H0
                             - 2.6666666666666665 * H00
                             + 2.6666666666666665 * H10
                             + 2.6666666666666665 * H11
                             - 2.6666666666666665 * Hm10
                             - 2.9629629629629624 / x
                             - (1.7777777777777777 * H1) / x
                             + H1
                                   * (5.777777777777778
                                      + x
                                            * (-35.55555555555556
                                               + 38.666666666666664 * x))
                             + x2
                                   * (77.62962962962962
                                      + 21.777777777777775 * H0
                                      + 5.333333333333333 * H01
                                      + 5.333333333333333 * H10
                                      + 5.333333333333333 * H11
                                      - 5.333333333333333 * Hm10
                                      - 5.333333333333333 * zeta2)
                             + x
                                   * (-66.22222222222221
                                      - 60.44444444444444 * H0
                                      - 21.333333333333332 * H00
                                      - 15.999999999999998 * H01
                                      - 5.333333333333333 * H10
                                      - 5.333333333333333 * H11
                                      - 5.333333333333333 * Hm10
                                      + 10.666666666666666 * zeta2))
                    + Lmmu2
                          * (-9.555555555555555 - 1.3333333333333333 * H0
                             - 2.6666666666666665 * H00
                             + 2.6666666666666665 * H10
                             + 5.333333333333333 * H11 + 0.888888888888889 / x
                             - (1.777777777777778 * H1) / x
                             + H1
                                   * (1.3333333333333333
                                      + x * (-32. + 35.11111111111111 * x))
                             + LQm
                                   * (1.3333333333333333
                                      + 1.7777777777777777 / x
                                      + 10.666666666666666 * x
                                      - 13.777777777777777 * x2
                                      + H0
                                            * (2.6666666666666665
                                               + 10.666666666666666 * x)
                                      + H1
                                            * (-2.6666666666666665
                                               + (5.333333333333333
                                                  - 5.333333333333333 * x)
                                                     * x))
                             + x2
                                   * (62.44444444444444
                                      + 13.777777777777777 * H0
                                      + 5.333333333333333 * H01
                                      + 5.333333333333333 * H10
                                      + 10.666666666666666 * H11
                                      - 5.333333333333333 * zeta2)
                             + x
                                   * (-53.77777777777777
                                      - 42.666666666666664 * H0
                                      - 10.666666666666666 * H00 - 16. * H01
                                      - 5.333333333333333 * H10
                                      - 10.666666666666666 * H11 + 16. * zeta2))
                    + (LQm
                       * (-3.5555555555555554 + 3.5555555555555554 * H10
                          + 3.5555555555555554 * H11 - 3.5555555555555554 * Hm10
                          + (40.22222222222222 + 13.185185185185185 * H0
                             + 17.333333333333332 * H00
                             + 2.6666666666666665 * H000
                             - 2.6666666666666665 * H01)
                                * x
                          - 3.5555555555555554 * zeta2
                          + H1
                                * (5.925925925925926
                                   + x
                                         * (-12.296296296296296
                                            + x
                                                  * (267.7037037037037
                                                     - 286.5185185185185 * x
                                                     - 5.333333333333333 * zeta2
                                                  )
                                            + 2.6666666666666665 * zeta2))
                          + x
                                * (-10.666666666666666 * H100
                                   - 11.555555555555555 * H11
                                   - 10.666666666666666 * H110
                                   - 5.333333333333333 * H111
                                   - 7.111111111111111 * Hm10
                                   + 10.666666666666666 * Hm100
                                   + 5.333333333333333 * Hm101
                                   + 5.333333333333333 * Hm1m10
                                   + (456.88888888888886
                                      + 418.0740740740741 * H0
                                      + 240.88888888888889 * H00
                                      + 69.33333333333333 * H000
                                      + 53.33333333333332 * H001
                                      + 173.33333333333331 * H01 + 32. * H010
                                      + 32. * H011)
                                         * x
                                   + H10
                                         * (-11.555555555555555
                                            + (71.11111111111111
                                               - 77.33333333333333 * x)
                                                  * x)
                                   + 2.6666666666666665 * zeta2
                                   - 2.6666666666666665 * Hm1 * zeta2
                                   + 8. * zeta3
                                   + x
                                         * (71.11111111111111 * H11
                                            + 21.333333333333332 * H110
                                            + 10.666666666666666 * H111
                                            + 12.444444444444443 * Hm10
                                            + 21.333333333333332 * Hm100
                                            + 10.666666666666666 * Hm101
                                            + 10.666666666666666 * Hm1m10
                                            + H100
                                                  * (21.333333333333332
                                                     - 16. * x)
                                            + (-512.4444444444445
                                               - 225.18518518518516 * H0
                                               - 56. * H00
                                               - 10.666666666666666 * H001
                                               - 93.33333333333333 * H01
                                               - 10.666666666666666 * H010
                                               - 10.666666666666666 * H011
                                               + 10.666666666666666 * H0m10)
                                                  * x
                                            + (-160.88888888888889
                                               - 53.33333333333332 * H0
                                               - 5.333333333333333 * Hm1)
                                                  * zeta2
                                            - 58.666666666666664 * zeta3
                                            + x
                                                  * (-77.33333333333333 * H11
                                                     - 21.333333333333332 * H110
                                                     - 10.666666666666666 * H111
                                                     + 30.22222222222222 * Hm10
                                                     + 26.66666666666666 * Hm100
                                                     + 10.666666666666666
                                                           * Hm101
                                                     + (93.33333333333333
                                                        + 10.666666666666666
                                                              * H0
                                                        - 10.666666666666666
                                                              * Hm1)
                                                           * zeta2
                                                     + 16. * zeta3)))))
                          / x
                    + (4.4883401920438954 - 7.11111111111111 * H100
                       + (-15.54732510288066 + 4.1646090534979425 * H0
                          + 4.246913580246913 * H00 + 2.3703703703703702 * H000
                          - 1.7777777777777775 * H0000
                          + 0.8888888888888887 * H01 - 8.88888888888889 * H010
                          - 10.666666666666666 * H0100)
                             * x
                       + H10
                             * (-7.703703703703703
                                + x
                                      * (0.44444444444444436
                                         + 0.8888888888888887 * zeta2
                                         + x
                                               * (-24.444444444444443
                                                  - 1.7777777777777775 * zeta2
                                                  + x
                                                        * (32.148148148148145
                                                           + 1.7777777777777775
                                                                 * zeta2))))
                       + 4.7407407407407405 * zeta3
                       + x
                             * (1.7777777777777775 * H1001
                                + 5.925925925925926 * H101
                                - 4.444444444444445 * H1010
                                - 4.444444444444445 * H1011
                                + 6.419753086419753 * H11
                                + 2.962962962962963 * H110
                                + 6.222222222222221 * H1100
                                - 7.11111111111111 * H1101
                                - 8.592592592592592 * H111
                                - 0.8888888888888887 * H1110
                                + 4.444444444444445 * H1111
                                + 5.530864197530864 * Hm10
                                - 1.4814814814814814 * Hm100
                                - 6.222222222222221 * Hm1000
                                - 5.333333333333333 * Hm1010
                                + 5.333333333333333 * Hm10m10
                                + 8.88888888888889 * Hm1m10
                                + 2.6666666666666665 * Hm1m100
                                + 5.333333333333333 * Hm1m1m10
                                + (-76.42798353909465 - 14.650205761316874 * H0
                                   + 22.814814814814813 * H00
                                   + 13.185185185185185 * H000
                                   + 21.333333333333332 * H0000
                                   + 3.555555555555555 * H001
                                   - 21.333333333333332 * H0010
                                   + 6.814814814814814 * H01
                                   - 37.33333333333333 * H010
                                   - 42.666666666666664 * H0100
                                   + 0.44444444444444436 * H011
                                   + 3.555555555555555 * H0m10)
                                      * x
                                - 3.555555555555555 * H1001 * x
                                - 2.962962962962963 * H101 * x
                                + 8.88888888888889 * H1010 * x
                                + 8.88888888888889 * H1011 * x
                                + 0.345679012345679 * H11 * x
                                - 4.148148148148148 * H110 * x
                                - 12.444444444444443 * H1100 * x
                                + 14.22222222222222 * H1101 * x
                                + 2.962962962962963 * H111 * x
                                + 1.7777777777777775 * H1110 * x
                                - 8.88888888888889 * H1111 * x
                                + 23.209876543209877 * Hm10 * x
                                - 4.7407407407407405 * Hm100 * x
                                - 12.444444444444443 * Hm1000 * x
                                - 10.666666666666666 * Hm1010 * x
                                + 10.666666666666666 * Hm10m10 * x
                                + 7.11111111111111 * Hm1m10 * x
                                + 5.333333333333333 * Hm1m100 * x
                                + 10.666666666666666 * Hm1m1m10 * x
                                + 63.21536351165981 * x2
                                - 90.2716049382716 * H0 * x2
                                + 56.641975308641975 * H00 * x2
                                + 13.037037037037036 * H000 * x2
                                + 8. * H001 * x2 - 1.7777777777777775 * H01 * x2
                                + 8. * H010 * x2
                                + 10.666666666666666 * H011 * x2
                                - 16. * H0m10 * x2
                                + 3.555555555555555 * H1001 * x2
                                + 2.962962962962963 * H101 * x2
                                - 8.88888888888889 * H1010 * x2
                                - 8.88888888888889 * H1011 * x2
                                - 4.345679012345679 * H11 * x2
                                + 4.148148148148148 * H110 * x2
                                + 12.444444444444443 * H1100 * x2
                                - 14.22222222222222 * H1101 * x2
                                - 0.2962962962962963 * H111 * x2
                                - 1.7777777777777775 * H1110 * x2
                                + 8.88888888888889 * H1111 * x2
                                + 22.320987654320987 * Hm10 * x2
                                - 4.7407407407407405 * Hm100 * x2
                                - 12.444444444444443 * Hm1000 * x2
                                - 10.666666666666666 * Hm1010 * x2
                                + 10.666666666666666 * Hm10m10 * x2
                                + 7.11111111111111 * Hm1m10 * x2
                                + 5.333333333333333 * Hm1m100 * x2
                                + 10.666666666666666 * Hm1m1m10 * x2
                                + H1000
                                      * (-2.6666666666666665
                                         + (5.333333333333333
                                            - 5.333333333333333 * x)
                                               * x)
                                + H100
                                      * (-9.481481481481481
                                         + x
                                               * (-43.25925925925926
                                                  + 58.370370370370374 * x))
                                - 0.8888888888888887 * zeta2
                                + 4.444444444444445 * H11 * zeta2
                                + 4.444444444444445 * Hm1 * zeta2
                                + 2.6666666666666665 * Hm10 * zeta2
                                + 2.6666666666666665 * Hm1m1 * zeta2
                                + 16.39506172839506 * x * zeta2
                                - 8.88888888888889 * H11 * x * zeta2
                                + 3.555555555555555 * Hm1 * x * zeta2
                                + 5.333333333333333 * Hm10 * x * zeta2
                                + 5.333333333333333 * Hm1m1 * x * zeta2
                                + 1.7777777777777775 * x2 * zeta2
                                - 8. * H0 * x2 * zeta2
                                + 8.88888888888889 * H11 * x2 * zeta2
                                + 3.555555555555555 * Hm1 * x2 * zeta2
                                + 5.333333333333333 * Hm10 * x2 * zeta2
                                + 5.333333333333333 * Hm1m1 * x2 * zeta2
                                + 12.8 * zeta2_2
                                + 24.888888888888886 * x * zeta2_2
                                - 14.222222222222221 * zeta3
                                + (H0
                                       * (7.111111111111111
                                          - 14.222222222222221 * x)
                                   + (-39.55555555555555 - 63.4074074074074 * x)
                                         * x
                                   + Hm1
                                         * (-5.333333333333333
                                            + (-10.666666666666666
                                               - 10.666666666666666 * x)
                                                  * x))
                                      * zeta3)
                       + H1
                             * (0.5925925925925926
                                + x
                                      * (4.526748971193416
                                         - 1.4814814814814814 * zeta2
                                         - 16. * zeta3
                                         + x
                                               * (15.1440329218107
                                                  - 0.5925925925925926 * zeta2
                                                  + x
                                                        * (-14.847736625514404
                                                           + 0.5925925925925926
                                                                 * zeta2
                                                           - 32. * zeta3)
                                                  + 32. * zeta3))))
                          / x
                    + Lmmu
                          * (LQm2
                                 * (1.3333333333333333 + 1.7777777777777777 / x
                                    + 10.666666666666666 * x
                                    - 13.777777777777777 * x2
                                    + H0
                                          * (2.6666666666666665
                                             + 10.666666666666666 * x)
                                    + H1
                                          * (-2.6666666666666665
                                             + (5.333333333333333
                                                - 5.333333333333333 * x)
                                                   * x))
                             + LQm
                                   * (-12.888888888888888
                                      - 11.555555555555555 * H0
                                      - 5.333333333333333 * H00
                                      + 5.333333333333333 * H10
                                      + 5.333333333333333 * H11
                                      - 5.333333333333333 * Hm10
                                      - 5.925925925925925 / x
                                      - (3.5555555555555554 * H1) / x
                                      + H1
                                            * (11.555555555555555
                                               + x
                                                     * (-71.11111111111111
                                                        + 77.33333333333333 * x)
                                            )
                                      + x2
                                            * (155.25925925925924
                                               + 43.55555555555555 * H0
                                               + 10.666666666666666 * H01
                                               + 10.666666666666666 * H10
                                               + 10.666666666666666 * H11
                                               - 10.666666666666666 * Hm10
                                               - 10.666666666666666 * zeta2)
                                      + x
                                            * (-132.44444444444443
                                               - 120.88888888888889 * H0
                                               - 42.666666666666664 * H00
                                               - 31.999999999999996 * H01
                                               - 10.666666666666666 * H10
                                               - 10.666666666666666 * H11
                                               - 10.666666666666666 * Hm10
                                               + 21.333333333333332 * zeta2))
                             + (-3.8518518518518516 + 7.111111111111111 * H10
                                + 3.5555555555555554 * H11
                                - 3.5555555555555554 * Hm10
                                + (31.407407407407405 + 12.888888888888888 * H0
                                   + 11.555555555555557 * H00
                                   + 5.333333333333333 * H000
                                   - 2.6666666666666665 * H01
                                   + 5.333333333333333 * H010)
                                      * x
                                - 3.5555555555555554 * zeta2
                                + H1
                                      * (5.925925925925926
                                         + x
                                               * (-8.
                                                  + 5.333333333333333 * zeta2
                                                  + x
                                                        * (258.22222222222223
                                                           - 10.666666666666666
                                                                 * zeta2
                                                           + x
                                                                 * (-277.037037037037
                                                                    + 5.333333333333333
                                                                          * zeta2
                                                                 ))))
                                + x
                                      * (-8. * H100 - 5.333333333333333 * H101
                                         - 21.777777777777775 * H11
                                         - 10.666666666666666 * H110
                                         - 16. * Hm10 + 8. * Hm100
                                         + 5.333333333333333 * Hm101
                                         + (437.037037037037
                                            + 378.22222222222223 * H0
                                            + 211.55555555555554 * H00
                                            + 42.666666666666664 * H000
                                            + 53.33333333333333 * H001
                                            + 170.66666666666666 * H01
                                            + 53.33333333333333 * H010
                                            + 32. * H011)
                                               * x
                                         + H10
                                               * (-7.555555555555555
                                                  + (92.44444444444446
                                                     - 106.22222222222223 * x)
                                                        * x)
                                         + 2.6666666666666665 * zeta2
                                         - 5.333333333333333 * Hm1 * zeta2
                                         + 18.66666666666666 * zeta3
                                         + x
                                               * (10.666666666666666 * H101
                                                  + 78.22222222222221 * H11
                                                  + 21.333333333333332 * H110
                                                  + 5.333333333333333 * Hm10
                                                  + 16. * Hm100
                                                  + 10.666666666666666 * Hm101
                                                  + H100
                                                        * (16.
                                                           - 10.666666666666666
                                                                 * x)
                                                  + (-468.5925925925926
                                                     - 216.5925925925926 * H0
                                                     - 52.888888888888886 * H00
                                                     - 10.666666666666666 * H001
                                                     - 92. * H01
                                                     - 10.666666666666666 * H010
                                                     - 10.666666666666666 * H011
                                                     + 10.666666666666666
                                                           * H0m10)
                                                        * x
                                                  + (-165.33333333333331
                                                     - 53.33333333333333 * H0
                                                     - 10.666666666666666 * Hm1)
                                                        * zeta2
                                                  - 5.333333333333333 * zeta3
                                                  + x
                                                        * (-10.666666666666666
                                                               * H101
                                                           - 83.1111111111111
                                                                 * H11
                                                           - 21.333333333333332
                                                                 * H110
                                                           + 23.111111111111114
                                                                 * Hm10
                                                           + 21.333333333333332
                                                                 * Hm100
                                                           + 10.666666666666666
                                                                 * Hm101
                                                           - 10.666666666666666
                                                                 * Hm1m10
                                                           + (92.
                                                              + 10.666666666666666
                                                                    * H0
                                                              - 16. * Hm1)
                                                                 * zeta2
                                                           + 16. * zeta3))))
                                   / x))
           + CA * CF
                 * (14.191358024691358 - 125.66666666666666 * H000
                    + 8.666666666666666 * H0000 - 8. * H00000 - 32. * H00001
                    - 4. * H0001 - 16. * H00011 + 11.333333333333332 * H001
                    - 21.333333333333332 * H0010 + 64. * H00100 + 24. * H00101
                    - 12. * H0011 + 56. * H00110 + 24. * H00111 + 8. * H00m100
                    - 16. * H00m1m10 - 13.777777777777777 * H01
                    - 31.333333333333332 * H010 + 94.66666666666666 * H0100
                    + 96. * H01000 + 104. * H01001 + 121.33333333333333 * H0101
                    + 64. * H01010 + 80. * H01011 - 74.66666666666666 * H011
                    + 76. * H0110 + 64. * H01100 + 24. * H01101
                    - 2.6666666666666665 * H0111 + 24. * H01110 - 24. * H01111
                    + 16. * H0m10 - 4. * H0m100 + 16. * H0m1000 + 8. * H0m1001
                    + 16. * H0m1010 - 16. * H0m10m10 + 8. * H0m1m10
                    - 8. * H0m1m100 - 16. * H0m1m101 - 278.3703703703704 * H1
                    - 103.77777777777777 * H10 + 282. * H100
                    + 33.33333333333333 * H1000 - 16. * H10000 - 32. * H10001
                    + 4. * H1001 + 32. * H10010 + 132. * H101
                    + 5.333333333333333 * H1010 + 16. * H10100 + 96. * H10101
                    + 64. * H10110 + 32. * H10111 - 166.00000000000003 * H11
                    + 134.66666666666666 * H110 - 48. * H11000 + 32. * H11001
                    + 17.333333333333332 * H1101 + 112. * H11010 + 128. * H11011
                    - 96. * H111 + 4. * H1110 - 16. * H11100 + 144. * H11101
                    + 149.33333333333331 * H1111 + 80. * H11110 - 84. * Hm10
                    + 28. * Hm100 + 32. * Hm1000 + 32. * Hm10000 + 32. * Hm10001
                    + 28. * Hm1001 + 48. * Hm10010 + 16. * Hm10011
                    - 32. * Hm100m10 + 16. * Hm101 + 8. * Hm1010 + 32. * Hm10100
                    + 32. * Hm10101 + 32. * Hm10110 - 56. * Hm10m10
                    - 16. * Hm10m100 - 32. * Hm10m101 - 40. * Hm1m10
                    - 76. * Hm1m100 - 32. * Hm1m1000 - 16. * Hm1m1001
                    - 56. * Hm1m101 - 64. * Hm1m1010 - 32. * Hm1m1011
                    + 96. * Hm1m1m10 - 48. * Hm1m1m100 + 96. * Hm1m1m1m10
                    - 51.77777777777778 / x
                    + (9.012345679012343 * H1 + 54.81481481481481 * H10
                       - 8.888888888888884 * H100 + 53.33333333333333 * H1000
                       + 42.666666666666664 * H1001 - 3.5555555555555545 * H101
                       + 21.333333333333332 * H1010 + 21.333333333333332 * H1011
                       - 10.518518518518515 * H11 - 3.5555555555555545 * H110
                       + 42.666666666666664 * H1100 - 15.999999999999993 * H111
                       + 21.333333333333332 * H1110 + 10.666666666666666 * H1111
                      ) / x
                    + 609.0493827160493 * x - 194.66666666666666 * H000 * x
                    + 118.66666666666666 * H0000 * x - 32 * H00001 * x
                    + 56 * H0001 * x - 32 * H00011 * x
                    + 258.66666666666663 * H001 * x
                    + 82.66666666666666 * H0010 * x + 96 * H00100 * x
                    - 16 * H00101 * x + 208 * H0011 * x + 112 * H00110 * x
                    + 48 * H00111 * x + 32 * H00m10 * x + 16 * H00m100 * x
                    - 32 * H00m1m10 * x + 187.77777777777777 * H01 * x
                    - 481.3333333333333 * H010 * x
                    + 322.66666666666663 * H0100 * x + 288 * H01000 * x
                    + 176 * H01001 * x + 149.33333333333331 * H0101 * x
                    + 64 * H01010 * x + 32 * H01011 * x
                    + 585.3333333333333 * H011 * x + 208 * H0110 * x
                    + 256 * H01100 * x - 48 * H01101 * x
                    + 301.3333333333333 * H0111 * x + 144 * H01110 * x
                    + 144 * H01111 * x + 56 * H0m10 * x + 32 * H0m100 * x
                    + 32 * H0m1000 * x + 16 * H0m1001 * x + 32 * H0m101 * x
                    + 32 * H0m1010 * x - 32 * H0m10m10 * x - 32 * H0m1m10 * x
                    - 16 * H0m1m100 * x - 32 * H0m1m101 * x
                    + 1087.037037037037 * H1 * x - 926.2222222222222 * H10 * x
                    - 309.3333333333333 * H100 * x
                    + 445.3333333333333 * H1000 * x + 32 * H10000 * x
                    + 64 * H10001 * x + 352 * H1001 * x - 64 * H10010 * x
                    - 512 * H101 * x + 69.33333333333333 * H1010 * x
                    - 32 * H10100 * x - 192 * H10101 * x + 80 * H1011 * x
                    - 128 * H10110 * x - 64 * H10111 * x
                    + 229.33333333333331 * H11 * x
                    - 474.66666666666663 * H110 * x + 368 * H1100 * x
                    + 96 * H11000 * x - 64 * H11001 * x
                    - 122.66666666666666 * H1101 * x - 224 * H11010 * x
                    - 256 * H11011 * x + 376 * H111 * x + 64 * H1110 * x
                    + 32 * H11100 * x - 288 * H11101 * x
                    - 82.66666666666666 * H1111 * x - 160 * H11110 * x
                    - 160 * Hm10 * x + 80 * Hm1000 * x + 64 * Hm10000 * x
                    + 64 * Hm10001 * x + 64 * Hm1001 * x + 96 * Hm10010 * x
                    + 32 * Hm10011 * x - 64 * Hm100m10 * x + 8 * Hm101 * x
                    + 64 * Hm1010 * x + 64 * Hm10100 * x + 64 * Hm10101 * x
                    + 32 * Hm1011 * x + 64 * Hm10110 * x - 96 * Hm10m10 * x
                    - 32 * Hm10m100 * x - 64 * Hm10m101 * x + 8 * Hm1m10 * x
                    - 96 * Hm1m100 * x - 64 * Hm1m1000 * x - 32 * Hm1m1001 * x
                    - 96 * Hm1m101 * x - 128 * Hm1m1010 * x - 64 * Hm1m1011 * x
                    + 96 * Hm1m1m10 * x - 96 * Hm1m1m100 * x
                    + 192 * Hm1m1m1m10 * x - 336.96296296296293 * x2
                    - 847.1111111111111 * H000 * x2
                    + 10.666666666666666 * H0000 * x2 - 32 * H00001 * x2
                    + 50.666666666666664 * H0001 * x2
                    - 415.1111111111111 * H001 * x2
                    - 114.66666666666666 * H0010 * x2 + 64 * H00100 * x2
                    + 32 * H00101 * x2 + 109.33333333333333 * H0011 * x2
                    + 32 * H00110 * x2 + 32 * H00111 * x2 + 32 * H00m10 * x2
                    + 32 * H00m100 * x2 - 64 * H00m1m10 * x2
                    + 257.55555555555554 * H01 * x2 - 320 * H010 * x2
                    - 573.3333333333333 * H0100 * x2 + 32 * H01001 * x2
                    - 200 * H0101 * x2 + 64 * H01010 * x2 + 96 * H01011 * x2
                    - 751.5555555555555 * H011 * x2 - 424 * H0110 * x2
                    - 64 * H01100 * x2 + 64 * H01101 * x2 - 128 * H0111 * x2
                    - 32 * H01111 * x2 - 8 * H0m10 * x2 + 48 * H0m100 * x2
                    + 64 * H0m1000 * x2 + 32 * H0m1001 * x2 + 32 * H0m101 * x2
                    + 64 * H0m1010 * x2 - 64 * H0m10m10 * x2 - 64 * H0m1m10 * x2
                    - 32 * H0m1m100 * x2 - 64 * H0m1m101 * x2
                    - 750.3456790123456 * H1 * x2 + 945.8518518518518 * H10 * x2
                    + 27.555555555555554 * H100 * x2
                    - 538.6666666666666 * H1000 * x2 - 32 * H10000 * x2
                    - 64 * H10001 * x2 - 402.66666666666663 * H1001 * x2
                    + 64 * H10010 * x2 + 363.55555555555554 * H101 * x2
                    - 114.66666666666666 * H1010 * x2 + 32 * H10100 * x2
                    + 192 * H10101 * x2 - 125.33333333333333 * H1011 * x2
                    + 128 * H10110 * x2 + 64 * H10111 * x2
                    - 175.48148148148147 * H11 * x2
                    + 323.55555555555554 * H110 * x2
                    - 402.66666666666663 * H1100 * x2 - 96 * H11000 * x2
                    + 64 * H11001 * x2 + 98.66666666666666 * H1101 * x2
                    + 224 * H11010 * x2 + 256 * H11011 * x2 - 228 * H111 * x2
                    - 125.33333333333333 * H1110 * x2 - 32 * H11100 * x2
                    + 288 * H11101 * x2 - 32 * H1111 * x2 + 160 * H11110 * x2
                    - 140 * Hm10 * x2 - 16 * Hm100 * x2 + 48 * Hm1000 * x2
                    + 64 * Hm10000 * x2 + 64 * Hm10001 * x2 + 48 * Hm1001 * x2
                    + 96 * Hm10010 * x2 + 32 * Hm10011 * x2 - 64 * Hm100m10 * x2
                    - 8 * Hm101 * x2 + 32 * Hm1010 * x2 + 64 * Hm10100 * x2
                    + 64 * Hm10101 * x2 + 32 * Hm1011 * x2 + 64 * Hm10110 * x2
                    - 64 * Hm10m10 * x2 - 32 * Hm10m100 * x2
                    - 64 * Hm10m101 * x2 + 24 * Hm1m10 * x2 - 80 * Hm1m100 * x2
                    - 64 * Hm1m1000 * x2 - 32 * Hm1m1001 * x2
                    - 64 * Hm1m101 * x2 - 128 * Hm1m1010 * x2
                    - 64 * Hm1m1011 * x2 + 96 * Hm1m1m10 * x2
                    - 96 * Hm1m1m100 * x2 + 192 * Hm1m1m1m10 * x2
                    + 35.638888888888886 * zeta2 + 32 * H000 * zeta2
                    - 40 * H001 * zeta2 - 8 * H00m1 * zeta2 - 144 * H01 * zeta2
                    - 112 * H010 * zeta2 - 48 * H011 * zeta2 + 4 * H0m1 * zeta2
                    - 20 * H0m10 * zeta2 + 16 * H0m1m1 * zeta2
                    - 131.77777777777777 * H1 * zeta2
                    - 31.333333333333332 * H10 * zeta2 + 4 * H100 * zeta2
                    - 144 * H101 * zeta2 - 76 * H11 * zeta2 - 80 * H110 * zeta2
                    - 192 * H111 * zeta2 - 36 * Hm1 * zeta2 - 52 * Hm10 * zeta2
                    - 60 * Hm100 * zeta2 - 48 * Hm101 * zeta2
                    + 32 * Hm10m1 * zeta2 + 104 * Hm1m1 * zeta2
                    + 24 * Hm1m10 * zeta2 + 48 * Hm1m1m1 * zeta2
                    + 96 * ln2 * zeta2 - (40 * zeta2) / x
                    + (16 * H1 * zeta2) / x - (32 * H10 * zeta2) / x
                    + (10.666666666666666 * H11 * zeta2) / x
                    - 306.38888888888886 * x * zeta2 + 40 * H000 * x * zeta2
                    + 80 * H001 * x * zeta2 - 16 * H00m1 * x * zeta2
                    + 64 * H01 * x * zeta2 - 64 * H010 * x * zeta2
                    + 192 * H011 * x * zeta2 - 48 * H0m1 * x * zeta2
                    - 40 * H0m10 * x * zeta2 + 32 * H0m1m1 * x * zeta2
                    + 810.2222222222222 * H1 * x * zeta2
                    - 161.33333333333331 * H10 * x * zeta2
                    - 8 * H100 * x * zeta2 + 288 * H101 * x * zeta2
                    + 408 * H11 * x * zeta2 + 160 * H110 * x * zeta2
                    + 384 * H111 * x * zeta2 - 4 * Hm1 * x * zeta2
                    - 108 * Hm10 * x * zeta2 - 120 * Hm100 * x * zeta2
                    - 96 * Hm101 * x * zeta2 + 64 * Hm10m1 * x * zeta2
                    + 144 * Hm1m1 * x * zeta2 + 48 * Hm1m10 * x * zeta2
                    + 96 * Hm1m1m1 * x * zeta2 - 192 * ln2 * x * zeta2
                    - 305.55555555555554 * x2 * zeta2 + 32 * H000 * x2 * zeta2
                    - 96 * H001 * x2 * zeta2 - 32 * H00m1 * x2 * zeta2
                    + 109.33333333333333 * H01 * x2 * zeta2
                    - 96 * H010 * x2 * zeta2 - 160 * H011 * x2 * zeta2
                    - 64 * H0m1 * x2 * zeta2 - 80 * H0m10 * x2 * zeta2
                    + 64 * H0m1m1 * x2 * zeta2
                    - 611.5555555555555 * H1 * x2 * zeta2
                    + 209.33333333333331 * H10 * x2 * zeta2
                    + 8 * H100 * x2 * zeta2 - 288 * H101 * x2 * zeta2
                    - 386.66666666666663 * H11 * x2 * zeta2
                    - 160 * H110 * x2 * zeta2 - 384 * H111 * x2 * zeta2
                    + 20 * Hm1 * x2 * zeta2 - 68 * Hm10 * x2 * zeta2
                    - 120 * Hm100 * x2 * zeta2 - 96 * Hm101 * x2 * zeta2
                    + 64 * Hm10m1 * x2 * zeta2 + 112 * Hm1m1 * x2 * zeta2
                    + 48 * Hm1m10 * x2 * zeta2 + 96 * Hm1m1m1 * x2 * zeta2
                    + 192 * ln2 * x2 * zeta2 + 13.466666666666667 * zeta2_2
                    + 173.60000000000002 * H1 * zeta2_2 + 36 * Hm1 * zeta2_2
                    - 441.3333333333333 * x * zeta2_2
                    - 347.20000000000005 * H1 * x * zeta2_2
                    + 72 * Hm1 * x * zeta2_2 + 313.8666666666667 * x2 * zeta2_2
                    + 347.20000000000005 * H1 * x2 * zeta2_2
                    + 72 * Hm1 * x2 * zeta2_2
                    + LQm3
                          * (-3.333333333333333 - 2.6666666666666665 * H00
                             - 2.6666666666666665 * H01 + 8.444444444444445 * H1
                             + 5.333333333333333 * H10
                             + 10.666666666666666 * H11
                             - (3.5555555555555554 * H1) / x
                             - 10.666666666666666 * x
                             + H0
                                   * (4.888888888888888
                                      + x
                                            * (-17.77777777777778
                                               + 47.111111111111114 * x))
                             + 2.6666666666666665 * zeta2
                             + x
                                   * (-10.666666666666666 * H00
                                      - 46.22222222222223 * H1
                                      - 10.666666666666666 * H10
                                      - 21.333333333333332 * H11
                                      + 6.666666666666666 * x
                                      + H01 * (-26.666666666666664 + 10.666666666666666 * x)
                                      + x
                                            * (47.111111111111114 * H1
                                               + 10.666666666666666 * H10
                                               + 21.333333333333332 * H11
                                               - 10.666666666666666 * zeta2)
                                      + 26.666666666666664 * zeta2))
                    + H0
                          * (-216.03703703703704
                             + (-38.55555555555555 - 17.6 * zeta2) * zeta2
                             + x2
                                   * (-117.45679012345678
                                      + 131.11111111111111 * zeta2
                                      - 406.2222222222222 * zeta3)
                             + x
                                   * (423.59259259259255
                                      + (-156.55555555555554 - 164.8 * zeta2)
                                            * zeta2
                                      - 29.777777777777775 * zeta3)
                             - 29.11111111111111 * zeta3)
                    + 11. * zeta3 - 54.666666666666664 * H01 * zeta3
                    - 4 * H0m1 * zeta3 + 14.444444444444443 * H1 * zeta3
                    + 77.33333333333333 * H10 * zeta3
                    + 18.666666666666664 * H11 * zeta3 - 86 * Hm1 * zeta3
                    - 8 * Hm10 * zeta3 - 56 * Hm1m1 * zeta3
                    - (14.222222222222221 * H1 * zeta3) / x
                    - 1660 * x * zeta3 - 18.666666666666664 * H01 * x * zeta3 - 8 * H0m1 * x * zeta3 - 230.2222222222222 * H1 * x * zeta3 - 154.66666666666666 * H10 * x * zeta3 - 37.33333333333333 * H11 * x * zeta3 - 112 * Hm1 * x * zeta3 - 16 * Hm10 * x * zeta3 - 112 * Hm1m1 * x * zeta3 + 453.3333333333333 * x2 * zeta3 - 37.33333333333333 * H01 * x2 * zeta3 - 16 * H0m1 * x2 * zeta3 + 241.77777777777777 * H1 * x2 * zeta3 + 154.66666666666666 * H10 * x2 * zeta3 + 37.33333333333333 * H11 * x2 * zeta3 - 104 * Hm1 * x2 * zeta3 - 16 * Hm10 * x2 * zeta3 - 112 * Hm1m1 * x2 * zeta3 - 83.33333333333333 * zeta2 * zeta3 - 425.3333333333333 * x * zeta2 * zeta3 + 157.33333333333331 * x2 * zeta2 * zeta3 + LQm * (402.73333333333335 + 120.00000000000001 * H000 - 80. * H0000 - 120.00000000000001 * H0001 + 68. * H001 - 96. * H0010 - 48. * H0011 + 224. * H00m10 + 76.93333333333334 * H01 + 86.66666666666664 * H010 + 112. * H0100 - 16. * H0101 + 142.66666666666666 * H011 + 362.6666666666667 * H0m10 - 104. * H0m100 - 224. * H0m101 + 16. * H0m1m10 + 229.9777777777778 * H1 + 99.1111111111111 * H10 + 429.3333333333333 * H100 + 16. * H1000 + 16. * H1001 + 234.66666666666669 * H101 + 160. * H1010 + 96. * H1011 + 96. * H10m10 + 57.77777777777777 * H11 + 183.99999999999997 * H110 + 112. * H1100 + 128. * H1101 + 208. * H111 + 160. * H1110 - 78.66666666666666 * Hm100 - 176.00000000000003 * Hm1000 - 192. * Hm1001 - 328. * Hm101 - 64. * Hm1010 - 64. * Hm1011 + 192. * Hm10m10 - 58.666666666666664 * Hm1m10 + 288. * Hm1m100 + 256. * Hm1m101 - 384. * Hm1m1m10 + (2.1333333333333333 * H00) / x + (2.133333333333333 * H01 + 12.088888888888889 * H1 + 28.44444444444444 * H10 - 42.666666666666664 * H100 - 85.33333333333333 * H101 + 28.44444444444444 * H11 - 85.33333333333333 * H110 - 64. * H111 - 2.133333333333333 * Hm10 - 21.333333333333332 * Hm100 + 42.666666666666664 * Hm1m10) / x + 883.0444444444445 * x - 325.33333333333326 * H000 * x - 288. * H0000 * x - 496. * H0001 * x - 442.66666666666663 * H001 * x - 640. * H0010 * x - 640. * H0011 * x - 64. * H00m10 * x + 73.60000000000001 * H01 * x - 741.3333333333333 * H010 * x - 224. * H0100 * x - 736. * H0101 * x - 741.3333333333333 * H011 * x - 768. * H0110 * x - 576. * H0111 * x - 186.66666666666666 * H0m10 * x - 399.99999999999994 * H0m100 * x - 64. * H0m101 * x + 799.9999999999999 * H0m1m10 * x + 507.91111111111115 * H1 * x - 675.5555555555554 * H10 * x - 970.6666666666665 * H100 * x - 32. * H1000 * x - 32. * H1001 * x - 1061.3333333333333 * H101 * x - 320. * H1010 * x - 192. * H1011 * x - 192. * H10m10 * x - 307.55555555555554 * H11 * x - 976. * H110 * x - 224. * H1100 * x - 256. * H1101 * x - 656. * H111 * x - 320. * H1110 * x + 1079.288888888889 * Hm10 * x - 266.66666666666663 * Hm100 * x - 352.00000000000006 * Hm1000 * x - 384. * Hm1001 * x - 357.3333333333333 * Hm101 * x - 128. * Hm1010 * x - 128. * Hm1011 * x + 384. * Hm10m10 * x + 432. * Hm1m10 * x + 576. * Hm1m100 * x + 512. * Hm1m101 * x - 768. * Hm1m1m10 * x - 1498.488888888889 * x2 + 1229.333333333333 * H000 * x2 + 96. * H0001 * x2 + 1616. * H001 * x2 + 192. * H0010 * x2 + 256. * H0011 * x2 + 128. * H00m10 * x2 + 1331.9111111111113 * H01 * x2 + 1498.6666666666667 * H010 * x2 + 288. * H0100 * x2 + 256. * H0101 * x2 + 1341.3333333333333 * H011 * x2 + 320. * H0110 * x2 + 192. * H0111 * x2 - 640. * H0m10 * x2 - 256. * H0m100 * x2 - 320. * H0m101 * x2 + 320. * H0m1m10 * x2 - 712.9777777777778 * H1 * x2 + 815.1111111111111 * H10 * x2 + 808. * H100 * x2 + 96. * H1000 * x2 + 96. * H1001 * x2 + 1178.6666666666665 * H101 * x2 + 320. * H1010 * x2 + 192. * H1011 * x2 + 192. * H10m10 * x2 + 460.44444444444446 * H11 * x2 + 1061.3333333333333 * H110 * x2 + 288. * H1100 * x2 + 256. * H1101 * x2 + 744. * H111 * x2 + 320. * H1110 * x2 + 99.2 * Hm10 * x2 - 304. * Hm100 * x2 - 288. * Hm1000 * x2 - 320. * Hm1001 * x2 - 160. * Hm101 * x2 - 128. * Hm1010 * x2 - 128. * Hm1011 * x2 + 256. * Hm10m10 * x2 + 640. * Hm1m10 * x2 + 384. * Hm1m100 * x2 + 384. * Hm1m101 * x2 - 512. * Hm1m1m10 * x2 + 76.80000000000001 * H000 * x3 + 76.80000000000001 * H001 * x3 - 76.80000000000001 * H0m10 * x3 - 179.20000000000002 * Hm10 * x3 - 76.80000000000001 * Hm100 * x3 - 76.80000000000001 * Hm101 * x3 + 76.80000000000001 * Hm1m10 * x3 + (-2.1333333333333333 * Hm100 - 2.1333333333333333 * Hm101 + 2.1333333333333333 * Hm1m10 + 55.82222222222222 * x + Hm10 * (-4.977777777777778 + 862.4000000000001 * x2)) / x2 - 76.93333333333334 * zeta2 + 24. * H01 * zeta2 + 232. * H0m1 * zeta2 - 264. * H1 * zeta2 + 128. * H10 * zeta2 + 64. * H11 * zeta2 + 298.66666666666663 * Hm1 * zeta2 + 336. * Hm10 * zeta2 - 448. * Hm1m1 * zeta2 + (1.0666666666666667 * H1 * zeta2) / x2 + (3.2 * Hm1 * zeta2) / x2 - (4.266666666666667 * zeta2) / x + (64. * H1 * zeta2) / x + (21.333333333333332 * Hm1 * zeta2) / x + 1005.688888888889 * x * zeta2 + 336. * H01 * x * zeta2 + 464. * H0m1 * x * zeta2 + 845.3333333333333 * H1 * x * zeta2 - 256. * H10 * x * zeta2 - 128. * H11 * x * zeta2 + 573.3333333333333 * Hm1 * x * zeta2 + 672. * Hm10 * x * zeta2 - 896. * Hm1m1 * x * zeta2 - 1331.9111111111113 * x2 * zeta2 - 96. * H01 * x2 * zeta2 + 480. * H0m1 * x2 * zeta2 - 858.6666666666666 * H1 * x2 * zeta2 + 128. * H10 * x2 * zeta2 + 480. * Hm1 * x2 * zeta2 + 544. * Hm10 * x2 * zeta2 - 640. * Hm1m1 * x2 * zeta2 - 179.20000000000002 * x3 * zeta2 - 38.400000000000006 * H1 * x3 * zeta2 + 115.2 * Hm1 * x3 * zeta2 + 31.200000000000003 * zeta2_2 - 187.20000000000002 * x * zeta2_2 + 150.4 * x2 * zeta2_2 + H00 * (-53.51111111111111 + 120. * zeta2 + x * (-975.7333333333335 + x * (1446.5777777777778 + 179.20000000000002 * x - 96. * zeta2) + 432. * zeta2)) + (H0 * (2.8444444444444454 + x * (153.11111111111111 - 68. * zeta2 + x * (522.5777777777778 + 256. * zeta2 + x * (10.755555555555556 + (-1616. - 153.60000000000002 * x) * zeta2 - 64. * zeta3) - 928. * zeta3) - 48. * zeta3))) / x - 65.33333333333333 * zeta3 + 392. * H1 * zeta3 + 408. * Hm1 * zeta3 - (42.666666666666664 * zeta3) / x - 861.3333333333333 * x * zeta3 - 784. * H1 * x * zeta3 + 816. * Hm1 * x * zeta3 + 8. * x2 * zeta3 + 560. * H1 * x2 * zeta3 + 592. * Hm1 * x2 * zeta3 - 192. * x3 * zeta3)
                    + (96. * LQm2
                       * (0.08333333333333333 + 0.2222222222222222 * H11
                          + (-0.3350694444444444 + 0.1261574074074074 * H0
                             - 0.38888888888888884 * H00 + 0.25 * H000
                             + 0.25 * H001 - 0.3055555555555555 * H01
                             + 0.08333333333333333 * H010
                             - 0.08333333333333333 * H011)
                                * x
                          + H10 * (0.2222222222222222 + x * (-0.736111111111111 + (3.222222222222222 - 3.444444444444444 * x) * x))
                          + H1
                                * (-0.07407407407407407
                                   + x
                                         * (-0.4421296296296296
                                            + x
                                                  * (2.4953703703703702
                                                     + x
                                                           * (-2.435185185185185
                                                              + zeta2)
                                                     - 1. * zeta2)
                                            + 0.5 * zeta2))
                          + x
                                * (-0.8333333333333334 * H11 - 0.5 * H110
                                   - 0.5 * H111 - 0.125 * Hm10
                                   + 0.08333333333333333 * Hm100
                                   + 0.16666666666666666 * Hm101
                                   + (0.9305555555555555
                                      + 1.3032407407407407 * H0
                                      + 1.3611111111111112 * H00 + 1. * H000
                                      + 1.666666666666667 * H001
                                      + 2.111111111111111 * H01
                                      + 1.8333333333333335 * H010
                                      + 2.1666666666666665 * H011)
                                         * x
                                   + 3.4166666666666665 * H11 * x + H110 * x
                                   + H111 * x - 0.25 * Hm10 * x
                                   + 0.16666666666666666 * Hm100 * x
                                   + 0.3333333333333333 * Hm101 * x
                                   - 0.09027777777777778 * x2
                                   - 2.685185185185185 * H0 * x2
                                   - 3.9166666666666665 * H00 * x2
                                   - 0.6666666666666666 * H001 * x2
                                   - 4.361111111111111 * H01 * x2
                                   - 0.6666666666666666 * H010 * x2
                                   - 1. * H011 * x2
                                   - 3.638888888888889 * H11 * x2
                                   - 1. * H110 * x2 - 1. * H111 * x2
                                   - 0.25 * Hm10 * x2
                                   + 0.16666666666666666 * Hm100 * x2
                                   + 0.3333333333333333 * Hm101 * x2
                                   + H100
                                         * (-0.4166666666666667
                                            + (0.8333333333333334
                                               - 0.8333333333333334 * x)
                                                  * x)
                                   + H101 * (-0.5 + x - 1. * x2)
                                   + 0.3055555555555555 * zeta2
                                   - 0.25 * H0 * zeta2
                                   - 0.16666666666666666 * Hm1 * zeta2
                                   - 2.361111111111111 * x * zeta2
                                   - 1.6666666666666667 * H0 * x * zeta2
                                   - 0.3333333333333333 * Hm1 * x * zeta2
                                   + 4.361111111111111 * x2 * zeta2
                                   + 0.6666666666666666 * H0 * x2 * zeta2
                                   - 0.3333333333333333 * Hm1 * x2 * zeta2
                                   + 0.3333333333333333 * x2 * zeta3)))
                          / x
                    + H00 * (-38.44444444444444 + 16.333333333333332 * zeta2 + 53.33333333333333 * zeta3 + x * (-331.88888888888886 + 35.33333333333333 * zeta2 + 69.33333333333333 * zeta3 + x * (1180.2222222222222 - 42.666666666666664 * zeta2 + 32 * zeta3)))
                    + Lmmu
                          * (-159.42222222222222 + 29.333333333333332 * H000
                             - 32. * H0000 - 32. * H0001
                             + 58.666666666666664 * H001 + 48. * H0011
                             + 128. * H00m10 + 58.93333333333332 * H01
                             + 42.666666666666664 * H010 + 64. * H0100
                             + 48. * H0101 + 121.33333333333331 * H011
                             + 80. * H0110 + 144. * H0111
                             + 10.666666666666666 * H0m10 - 64. * H0m100
                             - 128. * H0m101 + 91.86666666666666 * H1
                             + 66.66666666666666 * H10 + 128. * H100
                             + 96. * H1001 + 224.00000000000003 * H101
                             + 160. * H1010 + 288. * H1011
                             + 40.66666666666666 * H11
                             + 210.66666666666666 * H110 + 64. * H1100
                             + 256. * H1101 + 304. * H111
                             + 224.00000000000003 * H1110 + 384. * H1111
                             - 186.66666666666666 * Hm100 - 64. * Hm1000
                             - 64. * Hm1001 - 192. * Hm101 + 128. * Hm10m10
                             + 181.33333333333334 * Hm1m10 + 192. * Hm1m100
                             + 128. * Hm1m101 - 256. * Hm1m1m10
                             + (2.1333333333333333 * H00) / x
                             + (2.1333333333333333 * H01
                                - 26.133333333333336 * H1
                                - 53.33333333333332 * H10
                                - 21.33333333333333 * H100
                                - 63.99999999999999 * H101
                                - 31.999999999999996 * H11
                                - 42.66666666666666 * H110
                                - 63.99999999999999 * H111
                                - 2.1333333333333333 * Hm10
                                - 21.33333333333333 * Hm100
                                + 42.66666666666666 * Hm1m10)
                                   / x
                             + 393.2 * x - 176.00000000000003 * H000 * x
                             - 128. * H0000 * x - 320. * H0001 * x
                             - 736. * H001 * x - 384. * H0010 * x
                             - 576. * H0011 * x - 963.7333333333333 * H01 * x
                             - 885.3333333333333 * H010 * x - 128. * H0100 * x
                             - 672. * H0101 * x - 1474.6666666666667 * H011 * x
                             - 544. * H0110 * x - 864. * H0111 * x
                             - 42.666666666666664 * H0m10 * x
                             - 256. * H0m100 * x + 512. * H0m1m10 * x
                             - 1069.2 * H1 * x - 976. * H10 * x
                             - 352.00000000000006 * H100 * x - 192. * H1001 * x
                             - 1136. * H101 * x - 320. * H1010 * x
                             - 576. * H1011 * x - 1682.6666666666665 * H11 * x
                             - 949.3333333333333 * H110 * x - 128. * H1100 * x
                             - 512. * H1101 * x - 1616. * H111 * x
                             - 448. * H1110 * x - 768. * H1111 * x
                             + 489.9555555555556 * Hm10 * x
                             - 330.66666666666663 * Hm100 * x
                             - 128. * Hm1000 * x - 128. * Hm1001 * x
                             - 85.33333333333333 * Hm101 * x
                             + 256. * Hm10m10 * x + 576. * Hm1m10 * x
                             + 384. * Hm1m100 * x + 256. * Hm1m101 * x
                             - 512. * Hm1m1m10 * x - 314. * x2
                             + 448. * H000 * x2 + 128. * H0001 * x2
                             + 880. * H001 * x2 + 192. * H0010 * x2
                             + 384. * H0011 * x2 + 1592.8000000000002 * H01 * x2
                             + 784. * H010 * x2 + 128. * H0100 * x2
                             + 320. * H0101 * x2 + 1344. * H011 * x2
                             + 320. * H0110 * x2 + 576. * H0111 * x2
                             - 544. * H0m10 * x2 - 128. * H0m100 * x2
                             - 128. * H0m101 * x2 + 128. * H0m1m10 * x2
                             + 948.8000000000001 * H1 * x2 + 1020. * H10 * x2
                             + 384. * H100 * x2 + 64. * H1000 * x2
                             + 256. * H1001 * x2 + 1216. * H101 * x2
                             + 320. * H1010 * x2 + 576. * H1011 * x2
                             + 1708. * H11 * x2 + 992. * H110 * x2
                             + 192. * H1100 * x2 + 512. * H1101 * x2
                             + 1680. * H111 * x2 + 448. * H1110 * x2
                             + 768. * H1111 * x2
                             + 107.20000000000002 * Hm10 * x2
                             - 272. * Hm100 * x2 - 64. * Hm1000 * x2
                             - 64. * Hm1001 * x2 + 128. * Hm10m10 * x2
                             + 544. * Hm1m10 * x2 + 192. * Hm1m100 * x2
                             + 128. * Hm1m101 * x2 - 256. * Hm1m1m10 * x2
                             + 76.80000000000001 * H000 * x3
                             + 76.80000000000001 * H001 * x3
                             - 76.80000000000001 * H0m10 * x3
                             - 140.8 * Hm10 * x3
                             - 76.80000000000001 * Hm100 * x3
                             - 76.80000000000001 * Hm101 * x3
                             + 76.80000000000001 * Hm1m10 * x3
                             + (-2.1333333333333333 * Hm100
                                - 2.1333333333333333 * Hm101
                                + 2.1333333333333333 * Hm1m10
                                + 43.55555555555555 * x
                                + Hm10
                                      * (-3.911111111111112
                                         + 214.40000000000003 * x2))
                                   / x2
                             - 58.93333333333333 * zeta2 - 48. * H01 * zeta2
                             + 128. * H0m1 * zeta2
                             - 133.33333333333331 * H1 * zeta2
                             - 32. * H10 * zeta2 - 128. * H11 * zeta2
                             + 282.66666666666663 * Hm1 * zeta2
                             + 128. * Hm10 * zeta2 - 256. * Hm1m1 * zeta2
                             + (1.0666666666666667 * H1 * zeta2) / x2
                             + (3.2 * Hm1 * zeta2) / x2
                             - (4.266666666666667 * zeta2) / x
                             + (42.666666666666664 * H1 * zeta2) / x
                             + (21.333333333333332 * Hm1 * zeta2) / x
                             + 1453.6888888888889 * x * zeta2
                             + 416. * H01 * x * zeta2 + 256. * H0m1 * x * zeta2
                             + 848. * H1 * x * zeta2 + 64. * H10 * x * zeta2
                             + 256. * H11 * x * zeta2
                             + 373.3333333333333 * Hm1 * x * zeta2
                             + 256. * Hm10 * x * zeta2
                             - 512. * Hm1m1 * x * zeta2
                             - 1592.8000000000002 * x2 * zeta2
                             - 256. * H01 * x2 * zeta2
                             + 192. * H0m1 * x2 * zeta2 - 944. * H1 * x2 * zeta2
                             - 192. * H10 * x2 * zeta2 - 384. * H11 * x2 * zeta2
                             + 272. * Hm1 * x2 * zeta2
                             + 128. * Hm10 * x2 * zeta2
                             - 256. * Hm1m1 * x2 * zeta2 - 140.8 * x3 * zeta2
                             - 38.400000000000006 * H1 * x3 * zeta2
                             + 115.2 * Hm1 * x3 * zeta2 + 94.4 * zeta2_2
                             - 57.6 * x * zeta2_2 + 140.8 * x2 * zeta2_2
                             + LQm2
                                   * (-13.666666666666666 - 8. * H00 - 8. * H01
                                      + 10.666666666666666 * H1 + 16. * H10
                                      + 32. * H11
                                      - (10.666666666666666 * H1) / x
                                      - 17.333333333333332 * x
                                      + H0
                                            * (7.333333333333334
                                               + x
                                                     * (-38.666666666666664
                                                        + 112. * x))
                                      + 8. * zeta2
                                      + x
                                            * (-32. * H00
                                               - 109.33333333333333 * H1
                                               - 32. * H10 - 64. * H11 + 20. * x
                                               + H01 * (-80. + 32. * x)
                                               + x
                                                     * (112. * H1 + 32. * H10
                                                        + 64. * H11
                                                        - 32. * zeta2)
                                               + 80. * zeta2))
                             + H00 * (-87.73333333333333 + 32. * zeta2 + x * (-524.6222222222223 + x * (720.8000000000001 + 140.8 * x - 128. * zeta2) + 320. * zeta2))
                             - 134.66666666666666 * zeta3 + 224. * Hm1 * zeta3
                             + (21.333333333333332 * zeta3) / x
                             + 445.33333333333326 * x * zeta3
                             + 448. * Hm1 * x * zeta3 - 1328. * x2 * zeta3
                             - 224. * H1 * x2 * zeta3 + 224. * Hm1 * x2 * zeta3
                             - 192. * x3 * zeta3
                             + (H0
                                * (1.7777777777777777
                                   + x
                                         * (2.3111111111111113
                                            - 58.666666666666664 * zeta2
                                            + x
                                                  * (-180.
                                                     + 693.3333333333333 * zeta2
                                                     + x
                                                           * (1220.8
                                                              + (-880.
                                                                 - 153.60000000000002
                                                                       * x)
                                                                    * zeta2
                                                              - 128. * zeta3))
                                            + 80. * zeta3)))
                                   / x
                             + (LQm
                                * (16. + 42.666666666666664 * H11
                                   + (47.33333333333333
                                      + 62.666666666666664 * H0
                                      - 29.333333333333332 * H00 + 32. * H000
                                      + 16. * H001 - 28. * H01 + 16. * H010
                                      - 32. * H011)
                                         * x
                                   + H10 * (42.666666666666664 + x * (-82.66666666666666 + (565.3333333333333 - 608. * x) * x))
                                   + H1
                                         * (21.333333333333332
                                            + x
                                                  * (-30. + 128. * zeta2
                                                     + x
                                                           * (706.6666666666666
                                                              - 256. * zeta2
                                                              + x
                                                                    * (-688.
                                                                       + 256.
                                                                             * zeta2
                                                                    ))))
                                   + x
                                         * (-138.66666666666666 * H11
                                            - 128. * H110 - 192. * H111
                                            + (234. + 330.66666666666663 * H0
                                               + 218.66666666666666 * H00
                                               + 128. * H000 + 256. * H001
                                               + 552. * H01 + 352. * H010
                                               + 448. * H011)
                                                  * x
                                            + 757.3333333333333 * H11 * x
                                            + 256. * H110 * x + 384. * H111 * x
                                            - 268. * x2 - 656. * H0 * x2
                                            - 448. * H00 * x2 - 128. * H001 * x2
                                            - 608. * H01 * x2 - 128. * H010 * x2
                                            - 256. * H011 * x2 - 768. * H11 * x2
                                            - 256. * H110 * x2
                                            - 384. * H111 * x2
                                            + H101
                                                  * (-128.
                                                     + (256. - 256. * x) * x)
                                            + H100
                                                  * (-64.
                                                     + (128. - 128. * x) * x)
                                            + 28. * zeta2 - 16. * H0 * zeta2
                                            - 552. * x * zeta2
                                            - 256. * H0 * x * zeta2
                                            + 608. * x2 * zeta2
                                            + 128. * H0 * x2 * zeta2
                                            + 48. * zeta3 + 128. * x2 * zeta3)))
                                   / x)
                    - 2.0000000000000004 * zeta5
                    + 348.00000000000006 * x * zeta5
                    - 424.00000000000006 * x2 * zeta5)
           + CA * CA
                 * (683.7777777777777 + 54.666666666666664 * H00
                    - 23.333333333333332 * H000 + 30.66666666666666 * H0000
                    - 15.999999999999996 * H00000 + 15.999999999999996 * H00001
                    - 7.999999999999998 * H0001 + 31.999999999999993 * H00010
                    + 37.33333333333332 * H001 - 61.33333333333332 * H0010
                    - 31.999999999999993 * H00101 - 31.999999999999993 * H00110
                    - 66.66666666666666 * H01 + 126.66666666666663 * H010
                    - 5.333333333333332 * H0100 - 31.999999999999993 * H01000
                    - 31.999999999999993 * H01001 + 15.999999999999996 * H0101
                    + 31.999999999999993 * H01010 + 15.999999999999996 * H01011
                    - 18.66666666666666 * H011 + 15.999999999999996 * H0110 - 47.99999999999999 * H01100 + 31.999999999999993 * H01101 + 7.999999999999998 * H0111 + 31.999999999999993 * H01110 - 15.999999999999996 * H01111 + 15.999999999999996 * H0m1000 + 31.999999999999993 * H0m1010 + 31.999999999999993 * H0m1m100 - 63.999999999999986 * H0m1m1m10 + 95.25925925925924 * H1 - 233.3333333333333 * H10 + 105.99999999999999 * H100 - 23.999999999999996 * H1000 - 23.999999999999996 * H1001 - 31.999999999999993 * H10010 - 15.999999999999996 * H10011 - 29.33333333333333 * H101 - 13.33333333333333 * H1010 - 15.999999999999996 * H10100 - 79.99999999999999 * H10101 - 14.666666666666664 * H1011 - 79.99999999999999 * H10110 - 31.999999999999993 * H10111 + 2.888888888888889 * H11 - 29.33333333333333 * H110 - 17.33333333333333 * H1100 + 15.999999999999996 * H11000 - 15.999999999999996 * H11001 - 21.33333333333333 * H1101 - 95.99999999999999 * H11010 - 63.999999999999986 * H11011 + 3.3333333333333326 * H111 - 21.33333333333333 * H1110 - 79.99999999999999 * H11101 - 25.33333333333333 * H1111 - 79.99999999999999 * H11110 + 79.99999999999999 * H11111 - 44.444444444444436 * Hm10 - 10.666666666666664 * Hm100 - 6.666666666666665 * Hm1000 - 15.999999999999996 * Hm10000 - 15.999999999999996 * Hm10001 - 31.999999999999993 * Hm10010 - 13.33333333333333 * Hm1010 - 31.999999999999993 * Hm10100 - 31.999999999999993 * Hm10101 - 31.999999999999993 * Hm10110 - 31.999999999999993 * Hm10m100 + 63.999999999999986 * Hm10m1m10 + 21.33333333333333 * Hm1m10 - 13.33333333333333 * Hm1m100 - 15.999999999999996 * Hm1m1000 - 31.999999999999993 * Hm1m1001 + 31.999999999999993 * Hm1m1010 + 63.999999999999986 * Hm1m10m10 + 26.66666666666666 * Hm1m1m10 + 127.99999999999997 * Hm1m1m100 + 63.999999999999986 * Hm1m1m101 - 191.99999999999997 * Hm1m1m1m10 - 548.2057613168724 / x + 170. * zeta2 + (LQm3 * (-20.148148148148145 + H0 * (-3.5555555555555554 + x * (-4.888888888888888 + x * (-67.55555555555556 + 8. * x))) + H1 * (-7.111111111111111 + x * (-0.4444444444444444 + x * (-52.44444444444444 + 64.88888888888889 * x))) + x * (27.333333333333332 - 10.666666666666666 * H01 + 5.333333333333333 * H10 + 10.666666666666666 * H11 + H00 * (5.333333333333333 - 42.666666666666664 * x) - 146.66666666666666 * x + 10.666666666666666 * zeta2 + x * (-42.666666666666664 * H01 - 10.666666666666666 * H10 - 21.333333333333332 * H11 + 139.48148148148147 * x + 10.666666666666666 * H10 * x + 21.333333333333332 * H11 * x + 42.666666666666664 * zeta2)))) / x + 262. * zeta3 + (256. * LQm2 * (-0.546875 + 0.14583333333333331 * H1 + 0.125 * H10 + 0.125 * H11 + 0.08333333333333333 * Hm10 + (0.22092013888888887 + 0.11979166666666666 * H00 - 0.1875 * H000 + 0.1875 * H001) * x + H01 * (0.08333333333333333 + x * (0.09375 + (2.65625 - 1.3125 * x) * x)) + H0 * (-0.09027777777777778 + x * (-0.8758680555555555 + x * (3.048611111111111 - 4.499131944444444 * x - 0.875 * zeta2) - 0.1875 * zeta2)) + x * (-0.12499999999999999 * H0m10 - 0.051215277777777776 * H1 - 0.026041666666666664 * H10 - 0.09375 * H100 - 0.1875 * H101 - 0.08854166666666666 * H11 - 0.1875 * H110 - 0.1875 * H111 - 0.005208333333333333 * Hm10 + 0.09375 * Hm100 + 0.0625 * Hm101 - 0.125 * Hm1m10 + (3.1762152777777777 + 2.270833333333333 * H00 + H000 + 1.375 * H001) * x + 0.5 * H0m10 * x + 3.6493055555555554 * H1 * x + 1.2395833333333333 * H10 * x + 0.1875 * H100 * x + 0.375 * H101 * x + 1.6145833333333333 * H11 * x + 0.375 * H110 * x + 0.375 * H111 * x + 0.6145833333333333 * Hm10 * x + 0.1875 * Hm100 * x + 0.125 * Hm101 * x - 0.25 * Hm1m10 * x - 2.8815104166666665 * x2 - 0.28125 * H00 * x2 - 3.948784722222222 * H1 * x2 - 1.4583333333333333 * H10 * x2 - 0.1875 * H100 * x2 - 0.375 * H101 * x2 - 1.8333333333333333 * H11 * x2 - 0.375 * H110 * x2 - 0.375 * H111 * x2 + 0.7604166666666666 * Hm10 * x2 + 0.1875 * Hm100 * x2 + 0.125 * Hm101 * x2 - 0.25 * Hm1m10 * x2 + H010 * (0.12499999999999999 + (0.875 - 0.125 * x) * x) + H011 * (0.06249999999999999 + x - 0.25 * x2) - 0.09375 * zeta2 + 0.125 * H1 * zeta2 - 0.125 * Hm1 * zeta2 - 2.0416666666666665 * x * zeta2 - 0.25 * H1 * x * zeta2 - 0.25 * Hm1 * x * zeta2 + 1.3125 * x2 * zeta2 + 0.25 * H1 * x2 * zeta2 - 0.25 * Hm1 * x2 * zeta2 - 0.1875 * zeta3 + 0.0625 * x * zeta3))) / x + Lmmu2 * ((-128. * LQm * (0.4722222222222222 + H1 * (0.16666666666666666 + x * (0.010416666666666666 + (1.2291666666666665 - 1.5208333333333333 * x) * x)) + H0 * (0.08333333333333333 + x * (0.11458333333333331 + (1.5833333333333335 - 0.1875 * x) * x)) + x * (-0.640625 + 0.25 * H01 - 0.125 * H10 - 0.25 * H11 + H00 * (-0.125 + x) + 3.4375 * x + x * (H01 + 0.25 * H10 + 0.5 * H11 - 3.269097222222222 * x - 0.25 * H10 * x - 0.5 * H11 * x - 1. * zeta2) - 0.25 * zeta2))) / x + (-31.555555555555557 + 49.77777777777777 * H1 + 21.333333333333332 * H10 + 42.666666666666664 * H11 + (-28.777777777777775 + 14.666666666666666 * H00 - 16. * H000 + 16. * H001) * x + H01 * (10.666666666666666 + x * (16. + (744. - 218.66666666666666 * x) * x)) + H0 * (-5.333333333333333 + x * (-106. + x * (1096. - 558.4444444444445 * x - 256. * zeta2) - 16. * zeta2)) - 10.666666666666666 * zeta2 + x * (18. * H1 - 14.666666666666666 * H10 - 16. * H100 - 48. * H101 - 29.333333333333332 * H11 - 48. * H110 - 96. * H111 + (1995.7777777777776 + 586.6666666666666 * H00 + 128. * H000 + 256. * H001) * x + 1202.6666666666665 * H1 * x + 285.3333333333333 * H10 * x + 32. * H100 * x + 96. * H101 * x + 570.6666666666666 * H11 * x + 96. * H110 * x + 192. * H111 * x - 1935.4444444444443 * x2 - 24. * H00 * x2 - 1285.111111111111 * H1 * x2 - 322.66666666666663 * H10 * x2 - 32. * H100 * x2 - 96. * H101 * x2 - 645.3333333333333 * H11 * x2 - 96. * H110 * x2 - 192. * H111 * x2 + H011 * (32. + (320. - 64. * x) * x) + H010 * (16. + (160. - 32. * x) * x) - 16. * zeta2 + 48. * H1 * zeta2 - 744. * x * zeta2 - 96. * H1 * x * zeta2 + 218.66666666666666 * x2 * zeta2 + 96. * H1 * x2 * zeta2 - 16. * zeta3 - 256. * x * zeta3)) / x) + Lmmu * (319.66666666666663 + 565.7777777777777 * H00 - 61.33333333333333 * H000 + 95.99999999999999 * H0000 - 32. * H0001 - 80. * H001 - 64. * H0010 - 128. * H0011 + 32. * H00m10 + 190.66666666666666 * H01 - 53.33333333333332 * H010 - 95.99999999999999 * H0100 - 64. * H0101 - 23.999999999999996 * H011 - 95.99999999999999 * H0110 - 80. * H0m10 + 64. * H0m100 + 64. * H0m101 - 233.7777777777778 * H1 + 140.22222222222223 * H10 - 4. * H100 + 64. * H1000 + 128. * H1001 + 45.33333333333333 * H101 + 160. * H1010 + 128. * H1011 + 125.1111111111111 * H11 + 58.66666666666666 * H110 + 160. * H1100 + 160. * H1101 + 72. * H111 + 191.99999999999997 * H1110 + 202.66666666666666 * Hm10 + 59.999999999999986 * Hm100 - 64. * Hm1000 - 95.99999999999999 * Hm1001 + 98.66666666666664 * Hm101 - 32. * Hm1010 - 64. * Hm1011 + 32. * Hm10m10 - 47.99999999999999 * Hm1m10 + 95.99999999999999 * Hm1m100 + 64. * Hm1m101 - 64. * Hm1m1m10 - 111.11111111111111 / x - (128. * LQm2 * (0.4722222222222222 + H1 * (0.16666666666666666 + x * (0.010416666666666666 + (1.2291666666666665 - 1.5208333333333333 * x) * x)) + H0 * (0.08333333333333333 + x * (0.11458333333333331 + (1.5833333333333335 - 0.1875 * x) * x)) + x * (-0.640625 + 0.25 * H01 - 0.125 * H10 - 0.25 * H11 + H00 * (-0.125 + x) + 3.4375 * x + x * (H01 + 0.25 * H10 + 0.5 * H11 - 3.269097222222222 * x - 0.25 * H10 * x - 0.5 * H11 * x - 1. * zeta2) - 0.25 * zeta2))) / x - 190.66666666666666 * zeta2 - 94.66666666666666 * zeta3 + (LQm * (-280. + 74.66666666666666 * H1 + 64. * H10 + 64. * H11 + 42.666666666666664 * Hm10 + (113.1111111111111 + 61.33333333333333 * H00 - 96. * H000 + 96. * H001) * x + H01 * (42.666666666666664 + x * (48. + (1360. - 672. * x) * x)) + H0 * (-46.22222222222222 + x * (-448.4444444444444 + x * (1560.888888888889 - 2303.555555555555 * x - 448. * zeta2) - 96. * zeta2)) + x * (-63.99999999999999 * H0m10 - 26.22222222222222 * H1 - 13.333333333333332 * H10 - 48. * H100 - 96. * H101 - 45.33333333333333 * H11 - 96. * H110 - 96. * H111 - 2.6666666666666665 * Hm10 + 48. * Hm100 + 32. * Hm101 - 64. * Hm1m10 + (1626.2222222222222 + 1162.6666666666665 * H00 + 512. * H000 + 704. * H001) * x + 256. * H0m10 * x + 1868.4444444444443 * H1 * x + 634.6666666666666 * H10 * x + 96. * H100 * x + 192. * H101 * x + 826.6666666666666 * H11 * x + 192. * H110 * x + 192. * H111 * x + 314.66666666666663 * Hm10 * x + 96. * Hm100 * x + 64. * Hm101 * x - 128. * Hm1m10 * x - 1475.3333333333333 * x2 - 144. * H00 * x2 - 2021.7777777777776 * H1 * x2 - 746.6666666666666 * H10 * x2 - 96. * H100 * x2 - 192. * H101 * x2 - 938.6666666666666 * H11 * x2 - 192. * H110 * x2 - 192. * H111 * x2 + 389.3333333333333 * Hm10 * x2 + 96. * Hm100 * x2 + 64. * Hm101 * x2 - 128. * Hm1m10 * x2 + H011 * (31.999999999999996 + (512. - 128. * x) * x) + H010 * (63.99999999999999 + (448. - 64. * x) * x) - 48. * zeta2 + 64. * H1 * zeta2 - 64. * Hm1 * zeta2 - 1045.3333333333333 * x * zeta2 - 128. * H1 * x * zeta2 - 128. * Hm1 * x * zeta2 + 672. * x2 * zeta2 + 128. * H1 * x2 * zeta2 - 128. * Hm1 * x2 * zeta2 - 96. * zeta3 + 32. * x * zeta3))) / x + (56.888888888888886 * H01 - 63.99999999999999 * H010 - 63.99999999999999 * H011 + 21.333333333333332 * H0m10 + 259.55555555555554 * H1 - 179.55555555555554 * H10 - 85.33333333333333 * H100 - 106.66666666666663 * H101 - 44.44444444444444 * H11 - 127.99999999999999 * H110 - 63.99999999999999 * H111 + 103.11111111111109 * Hm10 - 42.666666666666664 * Hm100 - 21.333333333333332 * Hm101 + 21.333333333333332 * Hm1m10 + 46.22222222222223 * zeta2 + 96. * H1 * zeta2 + 31.999999999999996 * Hm1 * zeta2 + x * (31.999999999999996 * H00 * zeta2 + 63.99999999999999 * H01 * zeta2 - 63.99999999999999 * H0m1 * zeta2 - 69.33333333333333 * H1 * zeta2 - 112. * H10 * zeta2 - 127.99999999999999 * H11 * zeta2 - 122.66666666666666 * Hm1 * zeta2 + 112. * Hm10 * zeta2 - 96. * Hm1m1 * zeta2 - 35.2 * zeta2_2 - 8. * H1 * zeta3 + 119.99999999999997 * Hm1 * zeta3) + x2 * (761.5555555555555 - 2250.6666666666665 * H000 - 511.99999999999994 * H0000 - 959.9999999999998 * H0001 - 2853.333333333333 * H001 - 1152. * H0010 - 1152. * H0011 - 127.99999999999999 * H00m10 - 3744. * H01 - 2469.333333333333 * H010 - 672. * H0100 - 832. * H0101 - 2128. * H011 - 959.9999999999999 * H0110 - 576. * H0111 - 576. * H0m10 - 352. * H0m100 - 255.99999999999997 * H0m101 + 192. * H0m1m10 - 3061.111111111111 * H1 - 3745.7777777777774 * H10 - 920. * H100 - 127.99999999999999 * H1000 - 255.99999999999997 * H1001 - 1274.6666666666665 * H101 - 320. * H1010 - 255.99999999999997 * H1011 - 2608.8888888888887 * H11 - 1461.3333333333333 * H110 - 320. * H1100 - 320. * H1101 - 768. * H111 - 384. * H1110 - 1074.6666666666665 * Hm10 - 632. * Hm100 - 127.99999999999999 * Hm1000 - 192. * Hm1001 - 346.66666666666663 * Hm101 - 63.99999999999999 * Hm1010 - 127.99999999999999 * Hm1011 + 63.99999999999999 * Hm10m10 + 544. * Hm1m10 + 192. * Hm1m100 + 127.99999999999999 * Hm1m101 - 127.99999999999999 * Hm1m1m10 + 2669.333333333333 * zeta2 + 736. * H01 * zeta2 + 352. * H0m1 * zeta2 + 1002.6666666666666 * H1 * zeta2 + 224. * H10 * zeta2 + 255.99999999999997 * H11 * zeta2 + 618.6666666666666 * Hm1 * zeta2 + 224. * Hm10 * zeta2 - 192. * Hm1m1 * zeta2 - 710.4000000000001 * zeta2_2 + H00 * (-3288.8888888888887 + 832. * zeta2) - 354.66666666666663 * zeta3 + 16. * H1 * zeta3 + 239.99999999999994 * Hm1 * zeta3) + x3 * (-954.1111111111111 + 3196.6666666666665 * H00 + 192. * H000 + 824. * H001 + 63.99999999999999 * H0010 + 127.99999999999999 * H0011 + 1069.3333333333333 * H010 + 96. * H0100 + 192. * H0101 + 1440. * H011 + 192. * H0110 + 192. * H0111 - 341.3333333333333 * H0m10 - 96. * H0m100 - 63.99999999999999 * H0m101 + 127.99999999999999 * H0m1m10 + 3938.0000000000005 * H10 + 1072. * H100 + 96. * H1000 + 224. * H1001 + 1461.3333333333333 * H101 + 320. * H1010 + 255.99999999999997 * H1011 + 2804.6666666666665 * H11 + 1685.3333333333335 * H110 + 288. * H1100 + 320. * H1101 + 856. * H111 + 384. * H1110 - 1203.5555555555554 * Hm10 - 856. * Hm100 - 160. * Hm1000 - 224. * Hm1001 - 528. * Hm101 - 63.99999999999999 * Hm1010 - 127.99999999999999 * Hm1011 + 127.99999999999999 * Hm10m10 + 736. * Hm1m10 + 288. * Hm1m100 + 192. * Hm1m101 - 255.99999999999997 * Hm1m1m10 + H01 * (4385.111111111111 - 127.99999999999999 * zeta2) - 4385.111111111111 * zeta2 + 127.99999999999999 * H0m1 * zeta2 - 160. * H10 * zeta2 - 192. * H11 * zeta2 + 896. * Hm1 * zeta2 + 288. * Hm10 * zeta2 - 320. * Hm1m1 * zeta2 + 96. * zeta2_2 - 709.3333333333333 * zeta3 + 352. * Hm1 * zeta3 + H1 * (3170.8888888888887 - 1093.3333333333333 * zeta2 + 96. * zeta3)) + H0 * (-17.77777777777778 + 21.333333333333332 * zeta2 + x * (-132.66666666666666 + 80. * zeta2 + x * (-3069.555555555555 + x * (7344.666666666666 - 824. * zeta2) + 2277.333333333333 * zeta2 - 31.999999999999996 * zeta3) + 127.99999999999999 * zeta3))) / x) + LQm * (1108.4074074074072 + 855.3333333333333 * H00 - 110.66666666666666 * H000 + 176. * H0000 - 36. * H0001 - 28. * H001 - 96. * H0010 - 80. * H0011 + 103.99999999999999 * H00m10 + 360. * H01 - 15.999999999999998 * H010 - 148.00000000000003 * H0100 - 15.999999999999998 * H0101 - 15.999999999999998 * H011 - 15.999999999999998 * H0110 - 152. * H0m10 + 72. * H0m100 + 31.999999999999996 * H0m101 + 80. * H0m1m10 + 34.74074074074074 * H1 + 8.888888888888888 * H10 - 81.33333333333333 * H100 + 104.00000000000001 * H1000 + 160. * H1001 + 40. * H101 + 96. * H1010 + 112. * H1011 - 48. * H10m10 + 70.22222222222221 * H11 + 98.66666666666666 * H110 + 127.99999999999999 * H1100 + 112. * H1101 + 77.33333333333333 * H111 + 112. * H1110 + 96. * H1111 + 113.77777777777777 * Hm10 - 18.666666666666664 * Hm100 - 72. * Hm1000 - 112. * Hm1001 + 34.666666666666664 * Hm101 - 31.999999999999996 * Hm1010 - 31.999999999999996 * Hm1011 - 15.999999999999998 * Hm10m10 + 50.66666666666666 * Hm1m10 + 80. * Hm1m100 + 63.99999999999999 * Hm1m101 + 96. * Hm1m1m10 - 707.6543209876543 / x - 360. * zeta2 - 20. * zeta3 + (92.44444444444446 * H01 - 42.666666666666664 * H010 - 53.33333333333333 * H011 + 21.333333333333332 * H0m10 + 296.00000000000006 * H1 - 51.55555555555556 * H10 - 96. * H100 - 74.66666666666666 * H101 - 63.11111111111113 * H11 - 74.66666666666666 * H110 - 64. * H111 + 213.33333333333331 * Hm10 - 74.66666666666666 * Hm100 - 42.666666666666664 * Hm101 - 21.333333333333332 * Hm1m10 + 120.88888888888889 * zeta2 + 85.33333333333333 * H1 * zeta2 + 32. * Hm1 * zeta2 + 128. * zeta3 + x * (36. * H00 * zeta2 + 56. * H01 * zeta2 + 8. * H0m1 * zeta2 - 14.666666666666666 * H1 * zeta2 - 192. * H10 * zeta2 - 160. * H11 * zeta2 - 9.333333333333332 * Hm1 * zeta2 + 80. * Hm10 * zeta2 - 16. * Hm1m1 * zeta2 + 62.800000000000004 * zeta2_2 - 280. * H1 * zeta3 + 24. * Hm1 * zeta3) + x2 * (-6535.407407407407 - 2197.333333333333 * H000 - 816. * H0000 - 1064. * H0001 - 2729.333333333333 * H001 - 896. * H0010 - 1024. * H0011 - 144. * H00m10 - 3533.333333333333 * H01 - 1951.9999999999998 * H010 - 744. * H0100 - 640. * H0101 - 2208. * H011 - 640. * H0110 - 576. * H0111 - 328. * H0m10 - 496.00000000000006 * H0m100 - 320. * H0m101 - 96. * H0m1m10 - 4054.37037037037 * H1 - 2687.111111111111 * H10 - 841.3333333333333 * H100 - 208. * H1000 - 320. * H1001 - 1024. * H101 - 192. * H1010 - 224. * H1011 + 96. * H10m10 - 3060.4444444444443 * H11 - 1141.3333333333333 * H110 - 256. * H1100 - 224. * H1101 - 1018.6666666666666 * H111 - 224. * H1110 - 192. * H1111 - 228.44444444444443 * Hm10 - 845.3333333333333 * Hm100 - 144. * Hm1000 - 224. * Hm1001 - 442.66666666666663 * Hm101 - 64. * Hm1010 - 64. * Hm1011 - 32. * Hm10m10 + 245.33333333333331 * Hm1m10 + 160. * Hm1m100 + 128. * Hm1m101 + 192. * Hm1m1m10 + 3304.8888888888887 * zeta2 + 688. * H01 * zeta2 + 272. * H0m1 * zeta2 + 901.3333333333333 * H1 * zeta2 + 384. * H10 * zeta2 + 320. * H11 * zeta2 + 565.3333333333333 * Hm1 * zeta2 + 160. * Hm10 * zeta2 - 32. * Hm1m1 * zeta2 - 7.2 * zeta2_2 + H00 * (-3240.8888888888887 + 920. * zeta2) + 1474.6666666666665 * zeta3 + 560. * H1 * zeta3 + 48. * Hm1 * zeta3) + x3 * (6220.209876543209 + 264. * H000 + 96. * H0001 + 944. * H001 + 64. * H0010 + 128. * H0011 + 4088.8888888888887 * H01 + 1040. * H010 + 32. * H0100 + 192. * H0101 + 1197.3333333333333 * H011 + 192. * H0110 + 192. * H0111 - 693.3333333333333 * H0m10 - 128. * H0m100 - 128. * H0m101 + 128. * H0m1m10 + 3903.0370370370365 * H1 + 2882.6666666666665 * H10 + 1096. * H100 + 176. * H1000 + 288. * H1001 + 1154.6666666666665 * H101 + 192. * H1010 + 224. * H1011 - 96. * H10m10 + 3254.222222222222 * H11 + 1272. * H110 + 224. * H1100 + 224. * H1101 + 1130.6666666666665 * H111 + 224. * H1110 + 192. * H1111 - 241.77777777777777 * Hm10 - 1037.3333333333333 * Hm100 - 176. * Hm1000 - 256. * Hm1001 - 581.3333333333333 * Hm101 - 64. * Hm1010 - 64. * Hm1011 + 32. * Hm10m10 + 266.66666666666663 * Hm1m10 + 256. * Hm1m100 + 192. * Hm1m101 + 64. * Hm1m1m10 + H00 * (4152.888888888889 - 96. * zeta2) - 4088.8888888888887 * zeta2 - 128. * H01 * zeta2 + 192. * H0m1 * zeta2 - 1021.3333333333331 * H1 * zeta2 - 320. * H10 * zeta2 - 256. * H11 * zeta2 + 714.6666666666666 * Hm1 * zeta2 + 224. * Hm10 * zeta2 - 160. * Hm1m1 * zeta2 + 140.8 * zeta2_2 - 1669.3333333333333 * zeta3 - 448. * H1 * zeta3 + 160. * Hm1 * zeta3) + H0 * (-84.14814814814814 + 21.333333333333332 * zeta2 + x * (-580.5185185185185 + 28. * zeta2 + 272. * zeta3 + x * (-5056.074074074074 + x * (2756.222222222222 - 944. * zeta2) + 2401.333333333333 * zeta2 + 1024. * zeta3)))) / x)
                    - 136. * zeta5
                    + (35.55555555555555 * H010 - 21.33333333333333 * H0100
                       - 108.51851851851849 * H1 + 205.62962962962965 * H10
                       - 68.44444444444443 * H100 - 21.33333333333333 * H1000
                       - 21.33333333333333 * H1001 + 53.33333333333332 * H101
                       + 21.33333333333333 * H1010 + 10.666666666666664 * H1011
                       - 9.925925925925922 * H11 + 53.33333333333332 * H110
                       - 31.999999999999993 * H1100 + 21.33333333333333 * H1101
                       - 0.8888888888888888 * H111 + 21.33333333333333 * H1110
                       - 10.666666666666664 * H1111 + 27.851851851851848 * Hm10
                       + 17.777777777777775 * Hm100
                       - 10.666666666666664 * Hm1000
                       - 21.33333333333333 * Hm1010 - 35.55555555555555 * Hm1m10
                       - 21.33333333333333 * Hm1m100
                       + 42.66666666666666 * Hm1m1m10
                       - 80.59259259259258 * zeta2
                       + 5.333333333333333 * H01 * zeta2
                       - 53.777777777777764 * H1 * zeta2
                       + 26.666666666666668 * H10 * zeta2
                       - 10.666666666666666 * H11 * zeta2
                       - 17.77777777777778 * Hm1 * zeta2 + 16. * Hm10 * zeta2
                       + 21.333333333333332 * Hm1m1 * zeta2
                       + 40.53333333333333 * zeta2_2
                       + x * zeta2
                             * (31.333333333333332 * H00 - 40. * H000
                                + 55.99999999999999 * H001 - 16. * H01
                                + 40. * H010 - 16. * H011 - 24. * H0m10
                                - 32. * H0m1m1 + 25.777777777777775 * H1
                                + 27.999999999999996 * H10 - 4. * H100
                                + 55.99999999999999 * H101
                                + 14.666666666666668 * H11 - 8. * H110
                                + 32. * H111 + 10.666666666666666 * Hm1
                                + 2.6666666666666665 * Hm10 + 36. * Hm100
                                + 48. * Hm101 + 32. * Hm10m1
                                + 13.333333333333334 * Hm1m1 + 40. * Hm1m10
                                - 160. * Hm1m1m1 - 48. * zeta2
                                - 92. * H1 * zeta2 - 37.6 * Hm1 * zeta2)
                       + 33.18518518518518 * zeta3
                       + 35.55555555555556 * H1 * zeta3 - 32. * Hm1 * zeta3
                       + x
                             * (53.333333333333336 * H00
                                + 53.333333333333336 * H01 + 48. * H0m1
                                + 18.22222222222222 * H1
                                - 26.666666666666668 * H10
                                - 69.33333333333333 * H11 - 20. * Hm1
                                - 32. * Hm10 + 160. * Hm1m1
                                + 110.66666666666666 * zeta2)
                             * zeta3
                       + x3
                             * (9537.465020576132 + 150. * H000 + 24. * H0000
                                - 61.33333333333333 * H0001
                                + 243.55555555555554 * H001 - 24. * H0011
                                - 954.9629629629628 * H01
                                + 693.3333333333333 * H010
                                + 258.66666666666663 * H0100
                                + 194.66666666666666 * H0101
                                + 118.44444444444444 * H011
                                + 194.66666666666666 * H0110 + 48. * H0111
                                + 32. * H0m10 + 48. * H0m100 - 96. * H0m1m10
                                + 1336.7407407407406 * H1
                                - 2376.296296296296 * H10
                                + 359.77777777777777 * H100
                                + 173.33333333333331 * H1000
                                + 141.33333333333331 * H1001 - 64. * H10010
                                - 32. * H10011 - 642.6666666666666 * H101
                                - 288. * H1010 - 32. * H10100 - 160. * H10101
                                - 136. * H1011 - 160. * H10110 - 64. * H10111
                                + 261.9259259259259 * H11
                                - 642.6666666666666 * H110
                                + 253.33333333333331 * H1100 + 32. * H11000
                                - 32. * H11001 - 248. * H1101 - 192. * H11010
                                - 128. * H11011 + 293.55555555555554 * H111
                                - 248. * H1110 - 160. * H11101 + 272. * H1111
                                - 160. * H11110 + 160. * H11111
                                + 723.8518518518517 * Hm10
                                + 319.1111111111111 * Hm100 - 80. * Hm1000
                                - 32. * Hm10000 - 32. * Hm10001 + 32. * Hm1001
                                - 64. * Hm10010 + 96. * Hm101
                                - 223.99999999999997 * Hm1010 - 64. * Hm10100
                                - 64. * Hm10101 - 64. * Hm10110 - 64. * Hm10m10
                                - 64. * Hm10m100 + 128. * Hm10m1m10
                                - 542.2222222222222 * Hm1m10 - 320. * Hm1m100
                                - 32. * Hm1m1000 - 64. * Hm1m1001
                                - 64. * Hm1m101 + 64. * Hm1m1010
                                + 128. * Hm1m10m10 + 576. * Hm1m1m10
                                + 256. * Hm1m1m100 + 128. * Hm1m1m101
                                - 384. * Hm1m1m1m10 + 2131.4074074074074 * zeta2
                                - 348. * H01 * zeta2 - 48. * H0m1 * zeta2
                                + 691.3333333333333 * H1 * zeta2
                                - 182.66666666666669 * H10 * zeta2
                                - 8. * H100 * zeta2
                                + 111.99999999999999 * H101 * zeta2
                                + 136. * H11 * zeta2 - 16. * H110 * zeta2
                                + 64. * H111 * zeta2
                                - 367.11111111111114 * Hm1 * zeta2
                                + 89.33333333333333 * Hm10 * zeta2
                                + 72. * Hm100 * zeta2 + 96. * Hm101 * zeta2
                                + 64. * Hm10m1 * zeta2 + 352. * Hm1m1 * zeta2
                                + 80. * Hm1m10 * zeta2 - 320. * Hm1m1m1 * zeta2
                                + 45.86666666666667 * zeta2_2
                                - 184. * H1 * zeta2_2 - 75.2 * Hm1 * zeta2_2
                                + H00 * (388.7407407407407 + 25.333333333333332 * zeta2)
                                + (1212.148148148148 - 388.4444444444444 * H1
                                   - 53.333333333333336 * H10
                                   - 138.66666666666666 * H11 - 448. * Hm1
                                   - 64. * Hm10 + 320. * Hm1m1)
                                      * zeta3)
                       + H0
                             * (-64.79012345679013 - 17.77777777777778 * zeta2
                                - 3.5555555555555554 * zeta3
                                + x
                                      * (-404.3703703703703
                                         + (-129.44444444444443
                                            - 30.400000000000002 * zeta2)
                                               * zeta2
                                         - 119.55555555555554 * zeta3
                                         + x2
                                               * (-5579.308641975308
                                                  - 719.3333333333333 * zeta2
                                                  + 93.33333333333333 * zeta3)
                                         + x
                                               * (-2475.111111111111
                                                  + zeta2
                                                        * (-207.1111111111111
                                                           + 198.4 * zeta2)
                                                  + 151.11111111111111 * zeta3))
                             )
                       + x2
                             * (-9702.37037037037 - 486.88888888888886 * H00
                                - 8. * H000 - 34.666666666666664 * H0000
                                + 64. * H00000 + 32. * H00001 - 32. * H0001
                                - 128. * H00010 + 34.666666666666664 * H001
                                + 101.33333333333333 * H0010 - 320. * H00100
                                - 64. * H00101 - 16. * H0011 - 64. * H00110
                                - 149.33333333333331 * H01
                                + 912.0000000000001 * H010
                                - 293.3333333333333 * H0100 - 128. * H01000
                                - 128. * H01001 + 160. * H0101 + 128. * H01010
                                + 64. * H01011 - 125.33333333333333 * H011
                                + 160. * H0110 - 192. * H01100 + 128. * H01101
                                - 48. * H0111 + 128. * H01110 - 64. * H01111
                                + 192. * H0m10 + 32. * H0m100 - 64. * H0m1000
                                - 128. * H0m1010 - 64. * H0m1m10
                                - 128. * H0m1m100 + 256. * H0m1m1m10
                                - 1340.8148148148148 * H1 + 2420. * H10
                                - 397.3333333333333 * H100 - 128. * H1000
                                - 96. * H1001 + 64. * H10010 + 32. * H10011
                                + 610.6666666666666 * H101
                                + 250.66666666666666 * H1010 + 32. * H10100
                                + 160. * H10101 + 125.33333333333333 * H1011
                                + 160. * H10110 + 64. * H10111
                                - 230.2222222222222 * H11
                                + 610.6666666666666 * H110
                                - 189.33333333333334 * H1100 - 32. * H11000
                                + 32. * H11001 + 218.66666666666666 * H1101
                                + 192. * H11010 + 128. * H11011 - 320. * H111
                                + 218.66666666666666 * H1110 + 160. * H11101
                                - 221.33333333333331 * H1111 + 160. * H11110
                                - 160. * H11111 + 651.5555555555555 * Hm10
                                + 290.66666666666663 * Hm100
                                - 61.33333333333333 * Hm1000 - 32. * Hm10000
                                - 32. * Hm10001 + 32. * Hm1001 - 64. * Hm10010
                                + 96. * Hm101 - 186.66666666666666 * Hm1010
                                - 64. * Hm10100 - 64. * Hm10101 - 64. * Hm10110
                                - 64. * Hm10m10 - 64. * Hm10m100
                                + 128. * Hm10m1m10 - 485.3333333333333 * Hm1m10
                                - 282.66666666666663 * Hm1m100 - 32. * Hm1m1000
                                - 64. * Hm1m1001 - 64. * Hm1m101
                                + 64. * Hm1m1010 + 128. * Hm1m10m10
                                + 501.3333333333333 * Hm1m1m10
                                + 256. * Hm1m1m100 + 128. * Hm1m1m101
                                - 384. * Hm1m1m1m10 - 378.1111111111111 * zeta2
                                + (-4. * H00 + 64. * H000 + 128. * H001
                                   - 216. * H01 + 160. * H010 - 64. * H011
                                   - 32. * H0m1 + 96. * H0m10 + 128. * H0m1m1
                                   - 678.2222222222222 * H1 + 128. * H10
                                   + 8. * H100 - 111.99999999999999 * H101
                                   - 125.33333333333333 * H11 + 16. * H110
                                   - 64. * H111 - 338.66666666666663 * Hm1
                                   + 61.33333333333333 * Hm10 + 72. * Hm100
                                   + 96. * Hm101 + 64. * Hm10m1
                                   + 314.66666666666663 * Hm1m1 + 80. * Hm1m10
                                   - 320. * Hm1m1m1)
                                      * zeta2
                                + (892.2666666666667 + 184. * H1 - 75.2 * Hm1)
                                      * zeta2_2
                                + 1861.333333333333 * zeta3
                                + (-330.66666666666663 * H00
                                   + 213.33333333333334 * H01 - 192. * H0m1
                                   + 310.22222222222223 * H1
                                   + 53.333333333333336 * H10
                                   + 138.66666666666666 * H11 - 392. * Hm1
                                   - 64. * Hm10 + 320. * Hm1m1
                                   + 370.6666666666666 * zeta2)
                                      * zeta3
                                + 1088. * zeta5))
                          / x)
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
               * (LQm3
                      * (0.8888888888888888 + 1.1851851851851851 / x
                         - 0.8888888888888888 * x - 1.1851851851851851 * x2
                         + H0 * (1.7777777777777777 + 1.7777777777777777 * x))
                  + LQm2
                        * (-31.11111111111111
                           + H00
                                 * (-10.666666666666666 - 10.666666666666666 * x
                                 )
                           + 1.7777777777777777 / x + 20.444444444444443 * x
                           + 8.88888888888889 * x2
                           + H0
                                 * (-14.222222222222221
                                    + x
                                          * (-19.555555555555554
                                             + 7.111111111111111 * x)))
                  + (LQm
                     * (24.493827160493826
                        + (135.7037037037037 + 138.07407407407408 * H0
                           + 51.55555555555555 * H00 + 32. * H000)
                              * x
                        - 7.111111111111111 * Hm10 * (1. + x) * (1. + x)
                              * (1. + x)
                        + x2
                              * (-100.14814814814814 + 32. * H000
                                 + H0
                                       * (91.85185185185185
                                          - 35.55555555555556 * x)
                                 + H00
                                       * (104.88888888888889
                                          - 21.333333333333332 * x)
                                 - 60.04938271604938 * x
                                 - 21.333333333333332 * zeta2)
                        - 7.111111111111111 * zeta2))
                        / x
                  + Lmmu2
                        * (LQm
                               * (5.333333333333333 + 7.111111111111111 / x
                                  - 5.333333333333333 * x
                                  - 7.111111111111111 * x2
                                  + H0
                                        * (10.666666666666666
                                           + 10.666666666666666 * x))
                           + (3.5555555555555554
                              + H1
                                    * (-7.111111111111111
                                       + x
                                             * (-5.333333333333333
                                                + x
                                                      * (5.333333333333333
                                                         + 7.111111111111111 * x
                                                      )))
                              + x
                                    * (-35.55555555555556
                                       - 10.666666666666666 * H00
                                       - 10.666666666666666 * H01
                                       + 3.5555555555555554 * x
                                       + H0
                                             * (-5.333333333333333
                                                + x
                                                      * (-26.666666666666664
                                                         + 7.111111111111111 * x
                                                      ))
                                       + 10.666666666666666 * zeta2
                                       + x
                                             * (-10.666666666666666 * H00
                                                - 10.666666666666666 * H01
                                                + 28.444444444444443 * x
                                                + 10.666666666666666 * zeta2)))
                                 / x)
                  + (48.72427983539095 + 7.111111111111111 * H100
                     - 3.5555555555555554 * H101 + 19.555555555555554 * H11
                     - 3.5555555555555554 * H110 - 14.222222222222221 * H111
                     + (-79.01234567901234 + 45.72839506172839 * H0
                        - 37.925925925925924 * H00 + 2.6666666666666665 * H000
                        - 5.333333333333333 * H0000 + 14.222222222222221 * H001
                        + 10.666666666666666 * H0010 - 5.333333333333333 * H0011
                        - 42.96296296296296 * H01 + 8.88888888888889 * H010
                        + 10.666666666666666 * H0100 - 5.333333333333333 * H0101
                        + 40. * H011 - 5.333333333333333 * H0110
                        - 21.333333333333332 * H0111)
                           * x
                     + H1
                           * (-25.086419753086417
                              + x
                                    * (9.185185185185185
                                       + (17.48148148148148
                                          - 1.5802469135802468 * x)
                                             * x))
                     + H10
                           * (-5.925925925925926
                              + x
                                    * (18.666666666666664
                                       + x
                                             * (-34.666666666666664
                                                + 21.925925925925924 * x)))
                     + 11.851851851851851 * zeta2 + 4.7407407407407405 * zeta3
                     + x
                           * (6.222222222222221 * H11
                              - 2.6666666666666665 * H110
                              - 10.666666666666666 * H111
                              + (60.641975308641975 - 46.716049382716044 * H0
                                 - 59.259259259259245 * H00
                                 + 13.333333333333332 * H000
                                 - 5.333333333333333 * H0000
                                 + 19.555555555555554 * H001
                                 + 10.666666666666666 * H0010
                                 - 5.333333333333333 * H0011
                                 - 23.407407407407405 * H01
                                 - 7.111111111111111 * H010
                                 + 10.666666666666666 * H0100
                                 - 5.333333333333333 * H0101
                                 + 45.33333333333333 * H011
                                 - 5.333333333333333 * H0110
                                 - 21.333333333333332 * H0111)
                                    * x
                              - 6.222222222222221 * H11 * x
                              + 2.6666666666666665 * H110 * x
                              + 10.666666666666666 * H111 * x
                              - 30.35390946502058 * x2
                              + 77.4320987654321 * H0 * x2
                              - 44.44444444444444 * H00 * x2
                              + 7.111111111111111 * H000 * x2
                              + 3.5555555555555554 * H001 * x2
                              - 19.555555555555554 * H01 * x2
                              - 10.666666666666666 * H010 * x2
                              + 14.222222222222221 * H011 * x2
                              - 19.555555555555554 * H11 * x2
                              + 3.5555555555555554 * H110 * x2
                              + 14.222222222222221 * H111 * x2
                              + H100
                                    * (5.333333333333333
                                       + (-5.333333333333333
                                          - 7.111111111111111 * x)
                                             * x)
                              + H101
                                    * (-2.6666666666666665
                                       + x
                                             * (2.6666666666666665
                                                + 3.5555555555555554 * x))
                              + 45.629629629629626 * zeta2
                              + 2.6666666666666665 * H0 * zeta2
                              - 5.333333333333333 * H00 * zeta2
                              + 36.74074074074074 * x * zeta2
                              + 13.333333333333332 * H0 * x * zeta2
                              - 5.333333333333333 * H00 * x * zeta2
                              - 8.296296296296296 * x2 * zeta2
                              + 7.111111111111111 * H0 * x2 * zeta2
                              + 8. * zeta2_2 + 8. * x * zeta2_2
                              - 32.888888888888886 * zeta3
                              + 33.77777777777778 * H0 * zeta3
                              - 82.66666666666666 * x * zeta3
                              + 33.77777777777778 * H0 * x * zeta3
                              - 43.85185185185185 * x2 * zeta3))
                        / x
                  + Lmmu
                        * (LQm2
                               * (2.6666666666666665 + 3.5555555555555554 / x
                                  - 2.6666666666666665 * x
                                  - 3.5555555555555554 * x2
                                  + H0
                                        * (5.333333333333333
                                           + 5.333333333333333 * x))
                           + (LQm
                              * (15.407407407407405
                                 + H1
                                       * (-7.111111111111111
                                          + x
                                                * (-5.333333333333333
                                                   + x
                                                         * (5.333333333333333
                                                            + 7.111111111111111
                                                                  * x)))
                                 + x
                                       * (-46.22222222222223
                                          + H00
                                                * (-21.333333333333332
                                                   - 21.333333333333332 * x)
                                          + H01
                                                * (-10.666666666666666
                                                   - 10.666666666666666 * x)
                                          + 10.666666666666666 * zeta2
                                          + x
                                                * (35.55555555555556
                                                   - 4.7407407407407405 * x
                                                   + 21.333333333333332 * H0 * x
                                                   + 10.666666666666666 * zeta2)
                                       )))
                                 / x
                           + (5.925925925925926 + 7.111111111111111 * H11
                              - 7.111111111111111 * Hm10
                              + (22.51851851851852 + 62.22222222222223 * H0
                                 + 21.333333333333332 * H000
                                 + 21.333333333333332 * H001
                                 + 21.333333333333332 * H010
                                 + 10.666666666666666 * H011)
                                    * x
                              + H10
                                    * (14.222222222222221
                                       + x
                                             * (10.666666666666666
                                                + (-10.666666666666666
                                                   - 14.222222222222221 * x)
                                                      * x))
                              + H1
                                    * (-15.407407407407407
                                       + x
                                             * (46.22222222222223
                                                + x
                                                      * (-35.55555555555556
                                                         + 4.7407407407407405
                                                               * x)))
                              - 7.111111111111111 * zeta2
                              + x
                                    * ((-90.07407407407408
                                        - 58.666666666666664 * H0
                                        + 53.33333333333333 * H00
                                        + 21.333333333333332 * H000
                                        + 21.333333333333332 * H001
                                        + 21.333333333333332 * H010
                                        + 10.666666666666666 * H011)
                                           * x
                                       + 61.629629629629626 * x2
                                       - 52.14814814814815 * H0 * x2
                                       - 21.333333333333332 * H00 * x2
                                       - 21.333333333333332 * H01 * x2
                                       + Hm10
                                             * (-21.333333333333332
                                                + (-21.333333333333332
                                                   - 7.111111111111111 * x)
                                                      * x)
                                       + H11
                                             * (5.333333333333333
                                                + (-5.333333333333333
                                                   - 7.111111111111111 * x)
                                                      * x)
                                       - 21.333333333333332 * H0 * zeta2
                                       - 21.333333333333332 * x * zeta2
                                       - 21.333333333333332 * H0 * x * zeta2
                                       + 21.333333333333332 * x2 * zeta2
                                       + 10.666666666666666 * zeta3
                                       + 10.666666666666666 * x * zeta3))
                                 / x))
           + CF * nf
                 * (LQm3
                        * (0.8888888888888888 + 1.1851851851851851 / x
                           - 0.8888888888888888 * x - 1.1851851851851851 * x2
                           + H0 * (1.7777777777777777 + 1.7777777777777777 * x))
                    + LQm2
                          * (-31.11111111111111
                             + H00
                                   * (-10.666666666666666
                                      - 10.666666666666666 * x)
                             + 1.7777777777777777 / x + 20.444444444444443 * x
                             + 8.88888888888889 * x2
                             + H0
                                   * (-14.222222222222221
                                      + x
                                            * (-19.555555555555554
                                               + 7.111111111111111 * x)))
                    + Lmmu2
                          * (LQm
                                 * (2.6666666666666665 + 3.5555555555555554 / x
                                    - 2.6666666666666665 * x
                                    - 3.5555555555555554 * x2
                                    + H0
                                          * (5.333333333333333
                                             + 5.333333333333333 * x))
                             + (1.7777777777777777
                                + H1
                                      * (-3.5555555555555554
                                         + x
                                               * (-2.6666666666666665
                                                  + x
                                                        * (2.6666666666666665
                                                           + 3.5555555555555554
                                                                 * x)))
                                + x
                                      * (-17.77777777777778
                                         - 5.333333333333333 * H00
                                         - 5.333333333333333 * H01
                                         + 1.7777777777777777 * x
                                         + H0
                                               * (-2.6666666666666665
                                                  + x
                                                        * (-13.333333333333332
                                                           + 3.5555555555555554
                                                                 * x))
                                         + 5.333333333333333 * zeta2
                                         + x
                                               * (-5.333333333333333 * H00
                                                  - 5.333333333333333 * H01
                                                  + 14.222222222222221 * x
                                                  + 5.333333333333333 * zeta2)))
                                   / x)
                    + (32. * LQm
                       * (0.41975308641975306 - 0.1111111111111111 * H11
                          - 0.2222222222222222 * Hm10
                          + (3.7962962962962963 + 3.518518518518518 * H0
                             + 1.611111111111111 * H00 + H000
                             + 0.4444444444444444 * H01
                             - 0.16666666666666666 * H011)
                                * x
                          + H1
                                * (0.18518518518518517
                                   + x
                                         * (0.25
                                            + (-0.08333333333333333
                                               - 0.35185185185185186 * x)
                                                  * x))
                          - 0.2222222222222222 * zeta2
                          + x
                                * (-0.6666666666666666 * Hm10
                                   + (-3.018518518518518
                                      + 1.6296296296296295 * H0
                                      + 3.2777777777777777 * H00 + H000
                                      + 0.611111111111111 * H01
                                      - 0.16666666666666666 * H011)
                                         * x
                                   + H11
                                         * (-0.08333333333333333
                                            + (0.08333333333333333
                                               + 0.1111111111111111 * x)
                                                  * x)
                                   - 0.4444444444444444 * zeta2
                                   + x
                                         * (Hm10
                                                * (-0.6666666666666666
                                                   - 0.2222222222222222 * x)
                                            + x
                                                  * (-1.1975308641975309
                                                     - 1.4629629629629628 * H0
                                                     - 0.6666666666666666 * H00
                                                     + 0.1111111111111111 * H01
                                                     - 0.1111111111111111
                                                           * zeta2)
                                            - 1.2777777777777777 * zeta2
                                            + 0.16666666666666666 * zeta3)
                                   + 0.16666666666666666 * zeta3)))
                          / x
                    + (-16.50480109739369 + 7.111111111111111 * H101
                       - 1.9753086419753085 * H11 + 7.111111111111111 * H110
                       - 1.1851851851851851 * H111
                       + (60.905349794238695 + 25.31687242798354 * H0
                          + 27.061728395061728 * H00 + 14.518518518518517 * H000
                          + 8.88888888888889 * H0000 - 7.111111111111111 * H0001
                          - 5.925925925925926 * H001
                          - 10.666666666666666 * H0010
                          + 3.5555555555555554 * H0011
                          - 19.753086419753085 * H01 - 23.11111111111111 * H010
                          - 10.666666666666666 * H0100
                          + 10.666666666666666 * H0101
                          + 2.962962962962963 * H011
                          + 10.666666666666666 * H0110
                          - 1.7777777777777777 * H0111
                          - 26.666666666666664 * H10)
                             * x
                       + H100
                             * (-7.111111111111111
                                + x
                                      * (-5.333333333333333
                                         + x
                                               * (5.333333333333333
                                                  + 7.111111111111111 * x)))
                       + H1
                             * (-11.851851851851851 - 7.111111111111111 * zeta2
                                + x
                                      * (-8.493827160493826
                                         - 5.333333333333333 * zeta2
                                         + x
                                               * (-9.28395061728395
                                                  + 5.333333333333333 * zeta2
                                                  + x
                                                        * (29.629629629629623
                                                           + 7.111111111111111
                                                                 * zeta2))))
                       + 9.481481481481481 * zeta3
                       + x
                             * (5.333333333333333 * H110
                                - 0.8888888888888888 * H111
                                + (-81.64609053497944 + 31.83539094650206 * H0
                                   + 5.728395061728395 * H00
                                   + 3.851851851851851 * H000
                                   + 8.88888888888889 * H0000
                                   - 7.111111111111111 * H0001
                                   + 4.7407407407407405 * H001
                                   - 10.666666666666666 * H0010
                                   + 3.5555555555555554 * H0011
                                   - 12.641975308641975 * H01
                                   - 12.444444444444443 * H010
                                   - 10.666666666666666 * H0100
                                   + 10.666666666666666 * H0101
                                   - 2.3703703703703702 * H011
                                   + 10.666666666666666 * H0110
                                   - 1.7777777777777777 * H0111
                                   + 37.33333333333333 * H10)
                                      * x
                                - 5.333333333333333 * H110 * x
                                + 0.8888888888888888 * H111 * x
                                + 37.24554183813443 * x2
                                - 9.876543209876543 * H0 * x2
                                + 10.666666666666666 * H00 * x2
                                - 7.111111111111111 * H000 * x2
                                + 7.111111111111111 * H001 * x2
                                - 21.925925925925924 * H01 * x2
                                + 7.111111111111111 * H010 * x2
                                - 3.5555555555555554 * H011 * x2
                                - 10.666666666666666 * H10 * x2
                                - 7.111111111111111 * H110 * x2
                                + 1.1851851851851851 * H111 * x2
                                + H101
                                      * (5.333333333333333
                                         + (-5.333333333333333
                                            - 7.111111111111111 * x)
                                               * x)
                                + H11
                                      * (6.222222222222221
                                         + x
                                               * (-11.555555555555555
                                                  + 7.308641975308641 * x))
                                + 19.753086419753085 * zeta2
                                + 5.925925925925926 * H0 * zeta2
                                + 7.111111111111111 * H00 * zeta2
                                - 10.666666666666666 * H01 * zeta2
                                + 12.641975308641975 * x * zeta2
                                - 4.7407407407407405 * H0 * x * zeta2
                                + 7.111111111111111 * H00 * x * zeta2
                                - 10.666666666666666 * H01 * x * zeta2
                                + 21.925925925925924 * x2 * zeta2
                                - 7.111111111111111 * H0 * x2 * zeta2
                                + 16. * zeta2_2 + 16. * x * zeta2_2
                                - 36.148148148148145 * zeta3
                                - 3.5555555555555554 * H0 * zeta3
                                - 34.37037037037037 * x * zeta3
                                - 3.5555555555555554 * H0 * x * zeta3
                                + 1.1851851851851851 * x2 * zeta3))
                          / x
                    + Lmmu
                          * (LQm2
                                 * (2.6666666666666665 + 3.5555555555555554 / x
                                    - 2.6666666666666665 * x
                                    - 3.5555555555555554 * x2
                                    + H0
                                          * (5.333333333333333
                                             + 5.333333333333333 * x))
                             + LQm
                                   * (-62.22222222222222
                                      + H00
                                            * (-21.333333333333332
                                               - 21.333333333333332 * x)
                                      + 3.5555555555555554 / x
                                      + 40.888888888888886 * x
                                      + 17.77777777777778 * x2
                                      + H0
                                            * (-28.444444444444443
                                               + x
                                                     * (-39.11111111111111
                                                        + 14.222222222222221 * x
                                                     )))
                             + (-7.111111111111111 * Hm10
                                + (97.77777777777777 + 78.22222222222221 * H0
                                   + 28.444444444444443 * H00
                                   + 21.333333333333336 * H000
                                   + 10.666666666666668 * H001
                                   + 23.111111111111114 * H01
                                   + 10.666666666666668 * H010
                                   - 10.666666666666668 * H011
                                   + 26.666666666666668 * H1)
                                      * x
                                + H10
                                      * (7.111111111111111
                                         + x
                                               * (5.333333333333333
                                                  + (-5.333333333333333
                                                     - 7.111111111111111 * x)
                                                        * x))
                                + H11
                                      * (-7.111111111111111
                                         + x
                                               * (-5.333333333333333
                                                  + x
                                                        * (5.333333333333333
                                                           + 7.111111111111111
                                                                 * x)))
                                - 7.111111111111111 * zeta2
                                + x
                                      * (Hm10
                                             * (-21.333333333333332
                                                + (-21.333333333333332
                                                   - 7.111111111111111 * x)
                                                      * x)
                                         - 23.111111111111114 * zeta2
                                         - 10.666666666666666 * H0 * zeta2
                                         + x2
                                               * (-28.444444444444443
                                                  - 46.22222222222223 * H0
                                                  - 14.222222222222221 * H00
                                                  - 7.111111111111111 * H01
                                                  + 10.666666666666666 * H1
                                                  + 7.111111111111111 * zeta2)
                                         + 21.333333333333332 * zeta3
                                         + x
                                               * (-69.33333333333333
                                                  + 53.33333333333333 * H0
                                                  + 92.44444444444446 * H00
                                                  + 21.333333333333332 * H000
                                                  + 10.666666666666666 * H001
                                                  + 12.444444444444443 * H01
                                                  + 10.666666666666666 * H010
                                                  - 10.666666666666666 * H011
                                                  - 37.33333333333333 * H1
                                                  - 33.77777777777778 * zeta2
                                                  - 10.666666666666666 * H0
                                                        * zeta2
                                                  + 21.333333333333332 * zeta3))
                               ) / x))
           + CF * CF
                 * (-92. * H000 - 100 * H000 * x + 24 * H0000 * x
                    - 8 * H00001 * x - 16 * H00010 * x - 16 * H00011 * x
                    + 103.33333333333333 * H001 * x + 88 * H0010 * x
                    - 16 * H00101 * x + 88 * H0011 * x + 16 * H00110 * x
                    + 16 * H00111 * x - 226.66666666666666 * H01 * x
                    - 202.66666666666666 * H010 * x + 96 * H0100 * x
                    + 48 * H01000 * x + 32 * H01001 * x + 80 * H0101 * x
                    + 32 * H01010 * x + 32 * H01011 * x
                    - 38.666666666666664 * H011 * x + 32 * H0110 * x
                    + 32 * H01100 * x - 8 * H0111 * x + 32 * H01110 * x
                    + 16 * H01111 * x + 53.925925925925924 * H1 * x
                    - 20 * H10 * x - 82.66666666666666 * H100 * x
                    - 24 * H1000 * x - 16 * H1001 * x - 60 * H101 * x
                    - 16 * H1010 * x - 16 * H1011 * x
                    - 33.33333333333333 * H11 * x
                    - 54.666666666666664 * H110 * x - 16 * H1100 * x
                    + 18.666666666666664 * H111 * x - 16 * H1110 * x
                    - 8 * H1111 * x + 97.18518518518518 * x2
                    - 94.22222222222221 * H000 * x2
                    + 10.666666666666666 * H0000 * x2
                    + 21.333333333333332 * H0001 * x2
                    - 67.55555555555556 * H001 * x2
                    + 21.333333333333332 * H0010 * x2
                    + 21.333333333333332 * H0011 * x2
                    - 89.48148148148148 * H01 * x2
                    - 131.55555555555554 * H010 * x2
                    - 21.333333333333332 * H0100 * x2
                    - 115.55555555555554 * H011 * x2
                    - 21.333333333333332 * H0110 * x2
                    - 10.666666666666666 * H0111 * x2
                    + 210.66666666666666 * H1 * x2
                    + 139.85185185185185 * H10 * x2
                    - 58.666666666666664 * H100 * x2 - 32 * H1000 * x2
                    - 21.333333333333332 * H1001 * x2 - 32 * H101 * x2
                    - 21.333333333333332 * H1010 * x2
                    - 21.333333333333332 * H1011 * x2
                    + 86.51851851851852 * H11 * x2 - 32 * H110 * x2
                    - 21.333333333333332 * H1100 * x2 - 16 * H111 * x2
                    - 21.333333333333332 * H1110 * x2
                    - 10.666666666666666 * H1111 * x2 + 20 * H000 * zeta2
                    + 32 * H001 * zeta2 - 52 * H01 * zeta2 - 16 * H010 * zeta2
                    + 8 * H011 * zeta2 - 35.33333333333333 * H1 * zeta2
                    - 8 * H10 * zeta2 + 4 * H11 * zeta2 - (40 * zeta2) / x
                    - (29.333333333333332 * H1 * zeta2) / x
                    - (10.666666666666666 * H10 * zeta2) / x
                    + (5.333333333333333 * H11 * zeta2) / x
                    + 224.66666666666666 * x * zeta2 + 20 * H000 * x * zeta2
                    + 32 * H001 * x * zeta2 - 60 * H01 * x * zeta2
                    - 16 * H010 * x * zeta2 + 8 * H011 * x * zeta2
                    + 59.33333333333333 * H1 * x * zeta2 + 8 * H10 * x * zeta2
                    - 4 * H11 * x * zeta2 + 113.48148148148148 * x2 * zeta2
                    - 5.333333333333333 * H01 * x2 * zeta2
                    + 5.333333333333333 * H1 * x2 * zeta2
                    + 10.666666666666666 * H10 * x2 * zeta2
                    - 5.333333333333333 * H11 * x2 * zeta2 - 32 * zeta2_2
                    + 36.800000000000004 * x * zeta2_2
                    + 39.46666666666667 * x2 * zeta2_2
                    + LQm3
                          * (-10.222222222222221 - 10.666666666666666 * H01
                             + H1 * (-5.333333333333333 - 7.111111111111111 / x)
                             + H00
                                   * (-5.333333333333333 - 5.333333333333333 * x
                                   )
                             + 10.222222222222221 * x
                             + 5.333333333333333 * H0 * x
                             + 10.666666666666666 * zeta2
                             + x
                                   * (-10.666666666666666 * H01
                                      + 7.111111111111111 * H0 * x
                                      + H1
                                            * (5.333333333333333
                                               + 7.111111111111111 * x)
                                      + 10.666666666666666 * zeta2))
                    + Lmmu2
                          * ((-8. * LQm
                              * ((1.9166666666666665 + H00 + 2. * H01) * x
                                 + H1
                                       * (1.3333333333333333 + x - 1. * x2
                                          - 1.3333333333333333 * x3)
                                 + x
                                       * (x
                                              * (-1.9166666666666665 + H00
                                                 + 2. * H01
                                                 + H0
                                                       * (-1.
                                                          - 1.3333333333333333
                                                                * x)
                                                 - 2. * zeta2)
                                          - 2. * zeta2)))
                                 / x
                             + (24.
                                * (0.8888888888888888 * H11
                                   + (1.8333333333333333
                                      + 1.1388888888888888 * H0
                                      + 0.3333333333333333 * H000 + 1. * H001
                                      + 0.3333333333333333 * H01
                                      + 0.6666666666666666 * H010
                                      + 1.3333333333333333 * H011)
                                         * x
                                   + H1
                                         * (-0.2222222222222222
                                            + x
                                                  * (2.8611111111111107
                                                     + (-0.861111111111111
                                                        - 1.7777777777777777 * x
                                                       ) * x))
                                   + H10
                                         * (0.4444444444444444
                                            + x
                                                  * (0.3333333333333333
                                                     + (-0.3333333333333333
                                                        - 0.4444444444444444 * x
                                                       ) * x))
                                   + x
                                         * (H11
                                                * (0.6666666666666666
                                                   + (-0.6666666666666666
                                                      - 0.8888888888888888 * x)
                                                         * x)
                                            + x2
                                                  * (-1.7777777777777777 * H0
                                                     - 0.4444444444444444 * H00
                                                     - 0.8888888888888888 * H01
                                                     + 0.8888888888888888
                                                           * zeta2)
                                            - 0.3333333333333333 * zeta2
                                            - 1. * H0 * zeta2
                                            + x
                                                  * (-1.8333333333333333
                                                     - 1.472222222222222 * H0
                                                     + 0.6666666666666666 * H00
                                                     + 0.3333333333333333 * H000
                                                     + H001
                                                     + 1.3333333333333333 * H01
                                                     + 0.6666666666666666 * H010
                                                     + 1.3333333333333333 * H011
                                                     - 1.3333333333333333
                                                           * zeta2
                                                     - 1. * H0 * zeta2
                                                     - 1. * zeta3)
                                            - 1. * zeta3)))
                                   / x)
                    + (LQm2
                       * (8. + 32. * H11
                          + (56.222222222222214 + 88. * H0 - 16. * H00
                             + 48. * H000 + 80. * H001 + 48. * H01 + 48. * H010
                             + 48. * H011)
                                * x
                          + H10 * (32. + x * (24. + (-24. - 32. * x) * x))
                          + H1
                                * (-14.222222222222221
                                   + x
                                         * (164.
                                            + (-148. - 1.7777777777777777 * x)
                                                  * x))
                          + x
                                * (H11 * (24. + (-24. - 32. * x) * x)
                                   - 48. * zeta2 - 80. * H0 * zeta2
                                   + x2
                                         * (26.666666666666664
                                            - 1.7777777777777777 * H0
                                            - 74.66666666666666 * H00
                                            - 74.66666666666666 * H01
                                            + 74.66666666666666 * zeta2)
                                   + x
                                         * (-90.88888888888889 - 140. * H0
                                            - 8. * H00 + 48. * H000 + 80. * H001
                                            + 16. * H01 + 48. * H010
                                            + 48. * H011 - 16. * zeta2
                                            - 80. * H0 * zeta2 - 32. * zeta3)
                                   - 32. * zeta3)))
                          / x
                    + H00
                          * (-114.66666666666666 + 24. * zeta2
                             + x
                                   * (-140.66666666666666
                                      + x * (6.518518518518518 - 32. * zeta2)
                                      + 8. * zeta2 - 5.333333333333333 * zeta3)
                             - 5.333333333333333 * zeta3)
                    + 5.333333333333333 * H01 * zeta3
                    + 2.6666666666666665 * H1 * zeta3
                    + (3.5555555555555554 * H1 * zeta3) / x
                    - 475.1111111111111 * x * zeta3
                    + 5.333333333333333 * H01 * x * zeta3
                    - 2.6666666666666665 * H1 * x * zeta3 - 80 * x2 * zeta3
                    - 3.5555555555555554 * H1 * x2 * zeta3
                    - 61.33333333333333 * zeta2 * zeta3
                    - 61.33333333333333 * x * zeta2 * zeta3
                    + H0
                          * (-170.14814814814815
                             + (-46.33333333333333 - 40. * zeta2) * zeta2
                             + x2
                                   * (151.11111111111111
                                      + 40.888888888888886 * zeta2
                                      - 3.5555555555555554 * zeta3)
                             - 24. * zeta3
                             + x
                                   * (149.25925925925924
                                      + (-155.66666666666666 - 40. * zeta2)
                                            * zeta2
                                      + 85.33333333333333 * zeta3))
                    + (LQm
                       * (Hm10
                              * (1.4222222222222223
                                 + x2
                                       * (378.6666666666666
                                          + 392.88888888888886 * x - 12.8 * x3))
                          + x
                                * (60.08888888888889 + 42.666666666666664 * H10
                                   - 64. * H100 - 106.66666666666664 * H101
                                   + 26.66666666666666 * H11
                                   - 106.66666666666664 * H110
                                   - 74.66666666666666 * H111
                                   + 21.333333333333332 * Hm101
                                   + 21.333333333333332 * Hm1m10
                                   + (-76.42962962962963 - 342. * H00
                                      + 80. * H000 - 168. * H0000 - 288. * H0001
                                      - 160. * H001 - 240. * H0010
                                      - 224. * H0011 + 128. * H00m10
                                      - 642.6666666666666 * H01 - 112. * H010
                                      - 64. * H0100 - 160. * H0101 - 152. * H011
                                      - 160. * H0110 - 112. * H0111
                                      + 160. * H0m10 + 64. * H0m100
                                      - 128. * H0m1m10)
                                         * x
                                   - 10.666666666666666 * Hm1 * zeta2
                                   + H1
                                         * (-43.25925925925926 + 96. * zeta2
                                            + x
                                                  * (-660. - 16. * zeta2
                                                     + x
                                                           * (849.3333333333333
                                                              + x
                                                                    * (-146.07407407407408
                                                                       - 96. * zeta2
                                                                    )
                                                              + 16. * zeta2)))
                                   - 42.666666666666664 * zeta3
                                   + x
                                         * (-433.3333333333333 * H11
                                            - 80. * H110 - 56. * H111
                                            + 128. * Hm100 + 64. * Hm101
                                            - 192. * Hm1m10
                                            + (305.762962962963
                                               + 67.77777777777777 * H00
                                               - 117.33333333333333 * H000
                                               - 168. * H0000 - 288. * H0001
                                               - 272. * H001 - 240. * H0010
                                               - 224. * H0011 + 296. * H01
                                               + 32. * H010 - 64. * H0100
                                               - 160. * H0101 - 40. * H011
                                               - 160. * H0110 - 112. * H0111
                                               + 53.33333333333332 * H0m10
                                               - 64. * H0m100 + 128. * H0m1m10)
                                                  * x
                                            + 417.3333333333333 * H11 * x
                                            + 80. * H110 * x + 56. * H111 * x
                                            + 128. * Hm100 * x + 64. * Hm101 * x
                                            - 192. * Hm1m10 * x
                                            - 289.4222222222222 * x2
                                            + 163.55555555555554 * H00 * x2
                                            + 277.3333333333333 * H000 * x2
                                            + 298.66666666666663 * H001 * x2
                                            + 195.55555555555554 * H01 * x2
                                            + 234.66666666666666 * H010 * x2
                                            + 202.66666666666669 * H011 * x2
                                            - 21.333333333333332 * H0m10 * x2
                                            + 64. * H100 * x2
                                            - 10.666666666666666 * H11 * x2
                                            + 106.66666666666664 * H110 * x2
                                            + 74.66666666666666 * H111 * x2
                                            + 21.333333333333332 * Hm101 * x2
                                            + 21.333333333333332 * Hm1m10 * x2
                                            + 12.8 * H00 * x3
                                            + H10
                                                  * (-441.3333333333333
                                                     + (441.3333333333333
                                                        - 42.666666666666664 * x
                                                       ) * x)
                                            + H101
                                                  * (-80.
                                                     + x
                                                           * (80.
                                                              + 106.66666666666664
                                                                    * x))
                                            + 642.6666666666666 * zeta2
                                            + 288. * H00 * zeta2
                                            + 96. * H01 * zeta2
                                            - 64. * H0m1 * zeta2
                                            - 160. * Hm1 * zeta2
                                            + 96.88888888888889 * x * zeta2
                                            + 288. * H00 * x * zeta2
                                            + 96. * H01 * x * zeta2
                                            + 64. * H0m1 * x * zeta2
                                            - 160. * Hm1 * x * zeta2
                                            - 195.55555555555554 * x2 * zeta2
                                            - 10.666666666666666 * Hm1 * x2
                                                  * zeta2
                                            - 12.8 * x3 * zeta2
                                            - 73.60000000000001 * zeta2_2
                                            - 105.60000000000001 * x * zeta2_2
                                            + 136. * zeta3 + 584. * x * zeta3)
                                   + H0
                                         * (-1.4222222222222223
                                            + x
                                                  * (-17.955555555555556
                                                     + 160. * zeta2
                                                     + x
                                                           * (973.6
                                                              + x
                                                                    * (-186.6074074074074
                                                                       - 298.66666666666663
                                                                             * zeta2
                                                                    )
                                                              + 325.3333333333333
                                                                    * zeta2
                                                              - 96. * zeta3)
                                                     + 32. * zeta3)))))
                          / x2
                    + Lmmu
                          * ((-16. * LQm2
                              * ((1.9166666666666665 + H00 + 2. * H01) * x
                                 + H1
                                       * (1.3333333333333333 + x - 1. * x2
                                          - 1.3333333333333333 * x3)
                                 + x
                                       * (x
                                              * (-1.9166666666666665 + H00
                                                 + 2. * H01
                                                 + H0
                                                       * (-1.
                                                          - 1.3333333333333333
                                                                * x)
                                                 - 2. * zeta2)
                                          - 2. * zeta2)))
                                 / x
                             + (80. * LQm
                                * (0.2 + 0.8 * H11
                                   + (1.661111111111111 + 1.75 * H0 - 0.2 * H00
                                      + 1. * H000 + 1.6 * H001 + 1.4 * H01
                                      + 1.2 * H010 + 1.2 * H011)
                                         * x
                                   + H1
                                         * (0.08888888888888888
                                            + x
                                                  * (3.7
                                                     + (-2.5
                                                        - 1.2888888888888888 * x
                                                       ) * x))
                                   + H10
                                         * (0.8
                                            + x * (0.6 + (-0.6 - 0.8 * x) * x))
                                   + x
                                         * (H11 * (0.6 + (-0.6 - 0.8 * x) * x)
                                            - 1.4 * zeta2
                                            - 1.5999999999999999 * H0 * zeta2
                                            + x2
                                                  * (0.3333333333333333
                                                     - 1.2888888888888888 * H0
                                                     - 1.3333333333333333 * H00
                                                     - 1.3333333333333333 * H01
                                                     + 1.3333333333333333
                                                           * zeta2)
                                            + x
                                                  * (-2.194444444444444
                                                     - 2.95 * H0 + 0.4 * H00
                                                     + H000 + 1.6 * H001
                                                     + 1.4 * H01 + 1.2 * H010
                                                     + 1.2 * H011 - 1.4 * zeta2
                                                     - 1.5999999999999999 * H0
                                                           * zeta2
                                                     - 0.4 * zeta3)
                                            - 0.4 * zeta3)))
                                   / x
                             + (Hm10
                                    * (1.4222222222222225
                                       + x2
                                             * (378.66666666666663
                                                + 392.8888888888889 * x
                                                - 12.8 * x3))
                                + x
                                      * (46.75555555555556
                                         - 39.11111111111111 * H10 - 64. * H100
                                         - 106.66666666666666 * H101
                                         - 14.222222222222221 * H11
                                         - 85.33333333333333 * H110 - 64. * H111
                                         + 21.333333333333332 * Hm101
                                         + 21.333333333333332 * Hm1m10
                                         + (-185.6888888888889 - 212. * H00
                                            + 16. * H000 - 80. * H0000
                                            - 176. * H0001 - 80. * H001
                                            - 160. * H0010 - 192. * H0011
                                            + 128. * H00m10
                                            - 366.66666666666663 * H01
                                            - 144. * H010 - 64. * H0100
                                            - 160. * H0101 - 192. * H011
                                            - 128. * H0110 - 96. * H0111
                                            + 160. * H0m10 + 64. * H0m100
                                            - 128. * H0m1m10)
                                               * x
                                         - 10.666666666666666 * Hm1 * zeta2
                                         + H1
                                               * (-25.777777777777775
                                                  + 96. * zeta2
                                                  + x
                                                        * (-355.1111111111111
                                                           - 16. * zeta2
                                                           + x
                                                                 * (283.1111111111111
                                                                    + x
                                                                          * (97.77777777777777
                                                                             - 96. * zeta2
                                                                          )
                                                                    + 16. * zeta2
                                                                 )))
                                         + 21.333333333333332 * zeta3
                                         + x
                                               * (-424. * H11 - 64. * H110
                                                  - 48. * H111 + 128. * Hm100
                                                  + 64. * Hm101 - 192. * Hm1m10
                                                  + (344.8
                                                     - 170.2222222222222 * H00
                                                     - 181.33333333333331 * H000
                                                     - 80. * H0000
                                                     - 176. * H0001
                                                     - 368. * H001
                                                     - 160. * H0010
                                                     - 192. * H0011 + 148. * H01
                                                     - 96. * H010 - 64. * H0100
                                                     - 160. * H0101
                                                     - 176. * H011
                                                     - 128. * H0110
                                                     - 96. * H0111
                                                     + 53.33333333333333 * H0m10
                                                     - 64. * H0m100
                                                     + 128. * H0m1m10)
                                                        * x
                                                  + 312. * H11 * x
                                                  + 64. * H110 * x
                                                  + 48. * H111 * x
                                                  + 128. * Hm100 * x
                                                  + 64. * Hm101 * x
                                                  - 192. * Hm1m10 * x
                                                  - 205.86666666666667 * x2
                                                  + 289.77777777777777 * H00
                                                        * x2
                                                  + 106.66666666666666 * H000
                                                        * x2
                                                  + 170.66666666666666 * H001
                                                        * x2
                                                  + 296.88888888888886 * H01
                                                        * x2
                                                  + 128. * H010 * x2
                                                  + 149.33333333333331 * H011
                                                        * x2
                                                  - 21.333333333333332 * H0m10
                                                        * x2
                                                  + 64. * H100 * x2
                                                  + 126.22222222222224 * H11
                                                        * x2
                                                  + 85.33333333333333 * H110
                                                        * x2
                                                  + 64. * H111 * x2
                                                  + 21.333333333333332 * Hm101
                                                        * x2
                                                  + 21.333333333333332 * Hm1m10
                                                        * x2
                                                  + 12.8 * H00 * x3
                                                  + H101
                                                        * (-80.
                                                           + x
                                                                 * (80.
                                                                    + 106.66666666666666
                                                                          * x))
                                                  + H10
                                                        * (-306.66666666666663
                                                           + x
                                                                 * (226.66666666666666
                                                                    + 119.1111111111111
                                                                          * x))
                                                  + 366.66666666666663 * zeta2
                                                  + 176. * H00 * zeta2
                                                  + 96. * H01 * zeta2
                                                  - 64. * H0m1 * zeta2
                                                  - 160. * Hm1 * zeta2
                                                  + 244.8888888888889 * x
                                                        * zeta2
                                                  + 176. * H00 * x * zeta2
                                                  + 96. * H01 * x * zeta2
                                                  + 64. * H0m1 * x * zeta2
                                                  - 160. * Hm1 * x * zeta2
                                                  - 296.88888888888886 * x2
                                                        * zeta2
                                                  - 10.666666666666666 * Hm1
                                                        * x2 * zeta2
                                                  - 12.8 * x3 * zeta2
                                                  - 16. * zeta2_2
                                                  - 48. * x * zeta2_2
                                                  + 80. * zeta3
                                                  + 512. * x * zeta3
                                                  - 96. * x2 * zeta3)
                                         + H0
                                               * (-1.4222222222222225
                                                  + x
                                                        * (-104.17777777777778
                                                           + 80. * zeta2
                                                           + 144. * zeta3
                                                           + x
                                                                 * (521.1555555555556
                                                                    + x
                                                                          * (110.57777777777778
                                                                             - 170.66666666666666
                                                                                   * zeta2
                                                                          )
                                                                    + 421.3333333333333
                                                                          * zeta2
                                                                    + 16. * zeta3
                                                                 )))))
                                   / x2)
                    + 120 * x * zeta5
                    + (-51.77777777777777 + 54.81481481481481 * H10
                       + 26.66666666666666 * H100 + 31.999999999999996 * H1000
                       + 21.333333333333332 * H1001 + 31.999999999999996 * H101
                       + 21.333333333333332 * H1010 + 21.333333333333332 * H1011
                       - 10.518518518518517 * H11 + 31.999999999999996 * H110
                       + 21.333333333333332 * H1100 - 15.999999999999998 * H111
                       + 21.333333333333332 * H1110 + 10.666666666666666 * H1111
                       + H1 * (-55.77777777777777 - 208.81481481481478 * x)
                       + (-210.07407407407408 + 8. * H0000 - 8. * H00001
                          - 16. * H0001 - 16. * H00010 - 16. * H00011
                          + 20.666666666666664 * H001
                          - 23.999999999999993 * H0010 - 16. * H00101
                          - 8. * H0011 + 16. * H00110 + 16. * H00111
                          - 78.66666666666666 * H01 - 45.33333333333333 * H010
                          + 64. * H0100 + 47.999999999999986 * H01000
                          + 32. * H01001 + 64. * H0101 + 32. * H01010
                          + 32. * H01011 + 2.6666666666666665 * H011
                          + 47.999999999999986 * H0110 + 32. * H01100
                          - 23.999999999999993 * H0111 + 32. * H01110
                          + 16. * H01111)
                             * x
                       + x
                             * (-174.66666666666666 * H10
                                + 114.66666666666666 * H100
                                + 23.999999999999996 * H1000 + 16. * H1001
                                + 59.999999999999986 * H101 + 16. * H1010
                                + 16. * H1011 - 42.666666666666664 * H11
                                + 54.66666666666666 * H110 + 16. * H1100
                                + 13.33333333333333 * H111 + 16. * H1110
                                + 8. * H1111 + 164.66666666666666 * x
                                + 96.66666666666666 * zeta2
                                - 108.88888888888889 * zeta3
                                + 119.99999999999997 * zeta5))
                          / x)
           + CA * CF
                 * (726.3950617283949 + 126.66666666666666 * H00
                    - 48.00000000000001 * H000 + 45.33333333333333 * H0000
                    - 15.999999999999998 * H00000 + 31.999999999999996 * H00010
                    - 90.66666666666666 * H0010 + 31.999999999999996 * H00100
                    - 33.777777777777786 * H01 + 133.33333333333331 * H010
                    - 50.66666666666667 * H0100 + 32. * H01010 + 16. * H01011
                    - 18.666666666666668 * H011 - 16. * H01100 + 32. * H01101
                    + 8. * H0111 + 32. * H01110 - 16. * H01111 + 16. * H0m1000
                    + 32. * H0m1010 + 32. * H0m1m100 - 64. * H0m1m1m10
                    - 24.296296296296294 * H1 - 236.44444444444446 * H10
                    + 117.33333333333333 * H100 + 10.666666666666666 * H101
                    + 16. * H1010 + 8. * H1011 - 36.44444444444445 * H11
                    + 10.666666666666666 * H110 - 8. * H1100 + 16. * H1101
                    - 13.333333333333332 * H111 + 16. * H1110 - 8. * H1111
                    - 44.44444444444445 * Hm10 - 10.666666666666666 * Hm100
                    + 8. * Hm1000 + 16. * Hm1010 + 21.333333333333332 * Hm1m10
                    + 16. * Hm1m100 - 32. * Hm1m1m10 - 537.4074074074075 / x
                    + 182.2222222222222 * zeta2
                    + (LQm3
                       * (-19.555555555555554
                          + H0
                                * (-3.5555555555555554
                                   + (-9.777777777777777 - 7.111111111111111 * x
                                     ) * x)
                          + H1
                                * (-3.5555555555555554
                                   + x
                                         * (-2.6666666666666665
                                            + x
                                                  * (2.6666666666666665
                                                     + 3.5555555555555554 * x)))
                          + x
                                * (20.
                                   + H00
                                         * (5.333333333333333
                                            - 10.666666666666666 * x)
                                   + H01
                                         * (-5.333333333333333
                                            - 5.333333333333333 * x)
                                   + 5.333333333333333 * zeta2
                                   + x
                                         * (-20. + 19.555555555555554 * x
                                            + 5.333333333333333 * zeta2))))
                          / x
                    + Lmmu2
                          * ((-16. * LQm
                              * (3.6666666666666665
                                 + H0
                                       * (0.6666666666666666
                                          + x
                                                * (1.8333333333333333
                                                   + 1.3333333333333333 * x))
                                 + H1
                                       * (0.6666666666666666
                                          + x
                                                * (0.5
                                                   + (-0.5
                                                      - 0.6666666666666666 * x)
                                                         * x))
                                 + x
                                       * (-3.75 + H01 + 3.75 * x
                                          + H00 * (-1. + 2. * x)
                                          + x
                                                * (H01 - 3.6666666666666665 * x
                                                   - 1. * zeta2)
                                          - 1. * zeta2)))
                                 / x
                             + (32.
                                * (-0.9583333333333333
                                   + 0.3333333333333333 * H01
                                   + 1.6666666666666665 * H1
                                   + 0.3333333333333333 * H10
                                   + 0.6666666666666666 * H11 + 0.625 * x
                                   + H0
                                         * (-0.16666666666666666
                                            + x
                                                  * (-2.625
                                                     + x
                                                           * (4.375
                                                              - 1.8333333333333333
                                                                    * x
                                                              - 1.5 * zeta2)))
                                   - 0.3333333333333333 * zeta2
                                   + x
                                         * (-0.5 * H000
                                            + 1.1666666666666665 * H01
                                            + 0.5 * H010 + H011
                                            - 0.20833333333333331 * H1
                                            + 0.25 * H10 + 0.5 * H11 + 8.125 * x
                                            + H00
                                                  * (0.9166666666666666
                                                     + 3.6666666666666665 * x)
                                            - 1.1666666666666665 * zeta2
                                            + x
                                                  * (H000 + 1.5 * H001
                                                     + 0.5 * H010 + H011
                                                     + 1.7083333333333333 * H1
                                                     - 0.25 * H10 - 0.5 * H11
                                                     + H01
                                                           * (1.9166666666666665
                                                              - 0.3333333333333333
                                                                    * x)
                                                     - 7.791666666666666 * x
                                                     + x
                                                           * (-3.1666666666666665
                                                                  * H1
                                                              - 0.3333333333333333
                                                                    * H10
                                                              - 0.6666666666666666
                                                                    * H11
                                                              + 0.3333333333333333
                                                                    * zeta2)
                                                     - 1.9166666666666665
                                                           * zeta2
                                                     - 1.5 * zeta3))))
                                   / x)
                    + (64. * LQm2
                       * (-2.085648148148148 + 0.3333333333333333 * H01
                          + 0.611111111111111 * H1 + 0.3333333333333333 * H10
                          + 0.3333333333333333 * H11 + 0.3333333333333333 * Hm10
                          + 3.0416666666666665 * x
                          + H0
                                * (-0.3611111111111111
                                   + x
                                         * (-2.2777777777777777
                                            + x
                                                  * (1.472222222222222
                                                     - 2.8055555555555554 * x
                                                     - 0.5 * zeta2)
                                            - 0.25 * zeta2))
                          + x
                                * (-0.75 * H000 + 0.25 * H001 + 0.125 * H01
                                   + 0.5 * H010 + 0.5 * H011 - 0.5 * H0m10
                                   - 0.7083333333333333 * H1 + 0.25 * H10
                                   + 0.25 * H11 - 0.25 * Hm10 - 1.375 * x
                                   + H00
                                         * (1.1666666666666665
                                            + 1.1666666666666665 * x)
                                   - 0.125 * zeta2
                                   + x
                                         * (H000 + H001 + 0.5 * H010
                                            + 0.5 * H011 + 0.5 * H0m10
                                            - 0.25 * H10 - 0.25 * H11
                                            - 0.25 * Hm10
                                            + H1
                                                  * (0.9583333333333333
                                                     - 0.861111111111111 * x)
                                            + 0.41898148148148145 * x
                                            - 0.3333333333333333 * H01 * x
                                            + x
                                                  * (-0.3333333333333333 * H10
                                                     - 0.3333333333333333 * H11
                                                     + 0.3333333333333333 * Hm10
                                                     + 0.3333333333333333
                                                           * zeta2)
                                            - 0.25 * zeta2 + 0.25 * zeta3)
                                   - 0.5 * zeta3)))
                          / x
                    + 305.3333333333333 * zeta3
                    + Lmmu
                          * (-251.4074074074074 + 408.8888888888889 * H00
                             - 149.33333333333331 * H000 + 96. * H0000
                             + 32. * H0001 - 106.66666666666663 * H001
                             - 32. * H0010 - 96. * H0011 + 32. * H00m10
                             + 190.2222222222222 * H01
                             - 82.66666666666666 * H010 - 80. * H0100
                             - 96. * H0101 + 18.666666666666664 * H011
                             - 128. * H0110 - 96. * H0111 - 80. * H0m10
                             + 80. * H0m100 + 64. * H0m101 - 32. * H0m1m10
                             - 300.4444444444444 * H1 + 103.99999999999999 * H10
                             - 55.99999999999999 * H100 - 48. * H101
                             + 61.33333333333333 * H11 - 64. * H110 - 48. * H111
                             + 202.66666666666663 * Hm10 + 24. * Hm100
                             + 32. * Hm101 + 16. * Hm1m10
                             - 105.48148148148148 / x
                             - (16. * LQm2
                                * (3.6666666666666665
                                   + H0
                                         * (0.6666666666666666
                                            + x
                                                  * (1.8333333333333333
                                                     + 1.3333333333333333 * x))
                                   + H1
                                         * (0.6666666666666666
                                            + x
                                                  * (0.5
                                                     + (-0.5
                                                        - 0.6666666666666666 * x
                                                       ) * x))
                                   + x
                                         * (-3.75 + H01 + 3.75 * x
                                            + H00 * (-1. + 2. * x)
                                            + x
                                                  * (H01
                                                     - 3.6666666666666665 * x
                                                     - 1. * zeta2)
                                            - 1. * zeta2)))
                                   / x
                             - 190.2222222222222 * zeta2
                             + (128. * LQm
                                * (-2.085648148148148 + 0.3333333333333333 * H01
                                   + 0.611111111111111 * H1
                                   + 0.3333333333333333 * H10
                                   + 0.3333333333333333 * H11
                                   + 0.3333333333333333 * Hm10
                                   + 3.0416666666666665 * x
                                   + H0
                                         * (-0.3611111111111111
                                            + x
                                                  * (-2.2777777777777777
                                                     + x
                                                           * (1.472222222222222
                                                              - 2.8055555555555554
                                                                    * x
                                                              - 0.5 * zeta2)
                                                     - 0.25 * zeta2))
                                   + x
                                         * (-0.75 * H000 + 0.25 * H001
                                            + 0.125 * H01 + 0.5 * H010
                                            + 0.5 * H011 - 0.5 * H0m10
                                            - 0.7083333333333333 * H1
                                            + 0.25 * H10 + 0.25 * H11
                                            - 0.25 * Hm10 - 1.375 * x
                                            + H00
                                                  * (1.1666666666666665
                                                     + 1.1666666666666665 * x)
                                            - 0.125 * zeta2
                                            + x
                                                  * (H000 + H001 + 0.5 * H010
                                                     + 0.5 * H011 + 0.5 * H0m10
                                                     - 0.25 * H10 - 0.25 * H11
                                                     - 0.25 * Hm10
                                                     + H1
                                                           * (0.9583333333333333
                                                              - 0.861111111111111
                                                                    * x)
                                                     + 0.41898148148148145 * x
                                                     - 0.3333333333333333 * H01
                                                           * x
                                                     + x
                                                           * (-0.3333333333333333
                                                                  * H10
                                                              - 0.3333333333333333
                                                                    * H11
                                                              + 0.3333333333333333
                                                                    * Hm10
                                                              + 0.3333333333333333
                                                                    * zeta2)
                                                     - 0.25 * zeta2
                                                     + 0.25 * zeta3)
                                            - 0.5 * zeta3)))
                                   / x
                             - 125.33333333333333 * zeta3
                             + (56.888888888888886 * H01 - 64. * H010
                                - 64. * H011 + 21.333333333333332 * H0m10
                                + 244.74074074074076 * H1
                                - 179.55555555555554 * H10
                                - 42.666666666666664 * H100 - 64. * H101
                                - 44.44444444444444 * H11
                                - 85.33333333333333 * H110 - 64. * H111
                                + 99.55555555555554 * Hm10 - 64. * Hm100
                                - 42.666666666666664 * Hm101
                                + 42.666666666666664 * Hm1m10
                                + 42.666666666666664 * zeta2
                                + 42.666666666666664 * H1 * zeta2
                                + 64. * Hm1 * zeta2
                                + x
                                      * (-32. * H00 + 80. * H01 - 80. * H0m1
                                         + 56. * H1 - 24. * Hm1 - 41.6 * zeta2)
                                      * zeta2
                                + x2
                                      * (1572.7407407407409
                                         - 421.3333333333333 * H000
                                         - 128. * H0000
                                         - 191.99999999999997 * H0001
                                         - 202.66666666666666 * H001
                                         - 224. * H0010
                                         - 191.99999999999997 * H0011
                                         - 32. * H00m10
                                         - 15.11111111111111 * H01
                                         - 122.66666666666666 * H010
                                         - 80. * H0100 - 96. * H0101
                                         + 42.666666666666664 * H011
                                         - 128. * H0110 - 96. * H0111
                                         - 96. * H0m10 - 80. * H0m100
                                         - 64. * H0m101 + 32. * H0m1m10
                                         + 292.4444444444444 * H1 - 240. * H10
                                         + 56. * H100 + 48. * H101
                                         - 117.33333333333333 * H11 + 64. * H110
                                         + 48. * H111
                                         + 10.666666666666666 * Hm10
                                         + 24. * Hm100 + 32. * Hm101
                                         + 16. * Hm1m10
                                         + 25.777777777777775 * zeta2
                                         + (80. * H01 + 80. * H0m1 - 56. * H1
                                            - 24. * Hm1 - 134.4 * zeta2)
                                               * zeta2
                                         + H00
                                               * (-308.4444444444444
                                                  + 160. * zeta2)
                                         - 301.3333333333333 * zeta3)
                                + x3
                                      * (-1215.8518518518517
                                         + 400.8888888888889 * H00
                                         + 42.666666666666664 * H001
                                         + 434.66666666666663 * H01 + 64. * H010
                                         + 85.33333333333333 * H011
                                         - 42.666666666666664 * H0m10
                                         - 236.74074074074076 * H1
                                         + 315.55555555555554 * H10
                                         + 42.666666666666664 * H100
                                         + 64. * H101 + 100.44444444444446 * H11
                                         + 85.33333333333333 * H110 + 64. * H111
                                         - 92.44444444444446 * Hm10
                                         - 64. * Hm100
                                         - 42.666666666666664 * Hm101
                                         + 42.666666666666664 * Hm1m10
                                         - 434.66666666666663 * zeta2
                                         - 42.666666666666664 * H1 * zeta2
                                         + 64. * Hm1 * zeta2 - 64. * zeta3)
                                + H0
                                      * (-17.77777777777778
                                         + 21.333333333333332 * zeta2
                                         + x
                                               * (-496.8888888888889
                                                  + 106.66666666666666 * zeta2
                                                  + x
                                                        * (39.11111111111111
                                                           + x
                                                                 * (966.8148148148148
                                                                    - 42.666666666666664
                                                                          * zeta2
                                                                 )
                                                           + 106.66666666666666
                                                                 * zeta2
                                                           - 32. * zeta3)
                                                  + 96. * zeta3)))
                                   / x)
                    + LQm
                          * (266.0740740740741 + 576.4444444444445 * H00
                             - 272. * H000 + 176. * H0000 + 24. * H0001
                             - 28. * H001 - 64. * H0010 - 80. * H0011
                             + 112. * H00m10 + 400.44444444444446 * H01
                             - 16. * H010 - 136. * H0100 - 80. * H0101
                             - 2.666666666666666 * H011 - 80. * H0110
                             - 80. * H0111 - 104. * H0m10 + 112. * H0m100
                             + 64. * H0m101 + 32. * H0m1m10
                             - 134.66666666666666 * H1 + 21.33333333333333 * H10
                             - 91.99999999999999 * H100 - 40. * H101 + 36. * H11
                             - 40. * H110 - 40. * H111
                             + 277.3333333333333 * Hm10 + 24. * Hm100
                             + 48. * Hm1m10 - 682.0740740740739 / x
                             - 400.44444444444446 * zeta2
                             - 37.33333333333333 * zeta3
                             + (92.44444444444446 * H01
                                - 42.66666666666666 * H010
                                - 53.33333333333333 * H011
                                + 21.33333333333333 * H0m10
                                + 327.4074074074074 * H1
                                - 55.11111111111112 * H10
                                - 74.66666666666666 * H100
                                - 53.33333333333333 * H101
                                - 47.111111111111114 * H11
                                - 53.33333333333333 * H110
                                - 53.33333333333333 * H111
                                + 209.77777777777774 * Hm10 - 96. * Hm100
                                - 64. * Hm101 + 117.33333333333333 * zeta2
                                + 53.33333333333333 * H1 * zeta2
                                + 64. * Hm1 * zeta2
                                + x * zeta2
                                      * (-24. * H00 + 96. * H01 - 48. * H0m1
                                         + 64. * H1 + 24. * Hm1 + 28. * zeta2)
                                + x3
                                      * (494.0740740740741
                                         + 670.2222222222222 * H00
                                         + 53.33333333333333 * H001 + 376. * H01
                                         + 74.66666666666666 * H010
                                         + 74.66666666666666 * H011
                                         - 106.66666666666666 * H0m10
                                         - 143.4074074074074 * H1
                                         + 135.11111111111111 * H10
                                         + 74.66666666666666 * H100
                                         + 53.33333333333333 * H101
                                         + 127.1111111111111 * H11
                                         + 53.33333333333333 * H110
                                         + 53.33333333333333 * H111
                                         + 129.77777777777777 * Hm10
                                         - 96. * Hm100 - 64. * Hm101
                                         - 376. * zeta2
                                         - 53.33333333333333 * H1 * zeta2
                                         + 64. * Hm1 * zeta2
                                         - 202.66666666666666 * zeta3)
                                + x2
                                      * (-78.07407407407408 - 328. * H000
                                         - 192. * H0000 - 200. * H0001
                                         - 116. * H001 - 160. * H0010
                                         - 176. * H0011 - 48. * H00m10
                                         + 32.44444444444444 * H01 - 32. * H010
                                         - 136. * H0100 - 80. * H0101
                                         - 10.666666666666666 * H011
                                         - 80. * H0110 - 80. * H0111
                                         + 8. * H0m10 - 112. * H0m100
                                         - 64. * H0m101 - 32. * H0m1m10
                                         - 49.33333333333333 * H1
                                         - 101.33333333333333 * H10 + 92. * H100
                                         + 40. * H101 - 116. * H11 + 40. * H110
                                         + 40. * H111
                                         + 197.33333333333331 * Hm10
                                         + 24. * Hm100 + 48. * Hm1m10
                                         + 164.88888888888889 * zeta2
                                         + zeta2
                                               * (96. * H01 + 48. * H0m1
                                                  - 64. * H1 + 24. * Hm1
                                                  + 32.800000000000004 * zeta2)
                                         + H00
                                               * (-462.2222222222222
                                                  + 152. * zeta2)
                                         - 93.33333333333333 * zeta3)
                                + 128. * zeta3
                                + H0
                                      * (-84.14814814814814
                                         + 21.333333333333332 * zeta2
                                         + x
                                               * (-1187.2592592592591
                                                  + 28. * zeta2 + 240. * zeta3
                                                  + x
                                                        * (-581.037037037037
                                                           + x
                                                                 * (-138.07407407407408
                                                                    - 53.33333333333333
                                                                          * zeta2
                                                                 )
                                                           + 124. * zeta2
                                                           + 224. * zeta3))))
                                   / x)
                    - 120.00000000000001 * zeta5
                    + (35.55555555555556 * H010 - 21.333333333333332 * H0100
                       - 43.72839506172839 * H1 + 199.7037037037037 * H10
                       - 100.44444444444444 * H100 + 17.77777777777778 * H101
                       + 21.333333333333332 * H1010 + 10.666666666666666 * H1011
                       - 9.925925925925924 * H11 + 17.77777777777778 * H110
                       - 10.666666666666666 * H1100 + 21.333333333333332 * H1101
                       - 0.8888888888888887 * H111 + 21.333333333333332 * H1110
                       - 10.666666666666666 * H1111 + 27.851851851851848 * Hm10
                       + 17.77777777777778 * Hm100 - 10.666666666666666 * Hm1000
                       - 21.333333333333332 * Hm1010
                       - 35.55555555555556 * Hm1m10
                       - 21.333333333333332 * Hm1m100
                       + 42.666666666666664 * Hm1m1m10
                       - 76.29629629629629 * zeta2
                       + 5.333333333333333 * H01 * zeta2
                       - 10.222222222222221 * H1 * zeta2
                       + 5.333333333333333 * H10 * zeta2
                       - 5.333333333333333 * H11 * zeta2
                       - 17.77777777777778 * Hm1 * zeta2 + 16. * Hm10 * zeta2
                       + 21.333333333333332 * Hm1m1 * zeta2
                       + 40.53333333333333 * zeta2_2 + 33.77777777777778 * zeta3
                       + 17.77777777777778 * H1 * zeta3 - 32. * Hm1 * zeta3
                       + x3
                             * (1864.5925925925926 + 315.85185185185185 * H00
                                - 48. * H000 + 156.44444444444443 * H010
                                - 8.88888888888889 * H011
                                - 455.7037037037037 * H10
                                + 108.44444444444444 * H100
                                - 33.77777777777778 * H101
                                - 21.333333333333332 * H1010
                                - 10.666666666666666 * H1011
                                + 48.59259259259259 * H11
                                - 33.77777777777778 * H110
                                + 10.666666666666666 * H1100
                                - 21.333333333333332 * H1101
                                + 24.888888888888886 * H111
                                - 21.333333333333332 * H1110
                                + 10.666666666666666 * H1111
                                + 107.85185185185186 * Hm10
                                + 33.77777777777778 * Hm100
                                - 10.666666666666666 * Hm1000
                                - 21.333333333333332 * Hm1010
                                - 67.55555555555556 * Hm1m10
                                - 21.333333333333332 * Hm1m100
                                + 42.666666666666664 * Hm1m1m10
                                + H01
                                      * (-16.59259259259259
                                         - 5.333333333333333 * zeta2)
                                + 272.7407407407407 * zeta2
                                - 5.333333333333333 * H10 * zeta2
                                + 5.333333333333333 * H11 * zeta2
                                - 33.77777777777778 * Hm1 * zeta2
                                + 16. * Hm10 * zeta2
                                + 21.333333333333332 * Hm1m1 * zeta2
                                + 5.333333333333333 * zeta2_2
                                + H1
                                      * (61.7283950617284
                                         + 34.22222222222222 * zeta2
                                         - 17.77777777777778 * zeta3)
                                + 341.3333333333333 * zeta3 - 32. * Hm1 * zeta3)
                       + x
                             * (45.33333333333333 * H00 * zeta2
                                - 24. * H000 * zeta2 + 8. * H001 * zeta2
                                - 14.666666666666666 * H01 * zeta2
                                + 8. * H010 * zeta2 - 8. * H011 * zeta2
                                - 24. * H0m10 * zeta2 - 32. * H0m1m1 * zeta2
                                - 26. * H1 * zeta2 + 4. * H10 * zeta2
                                - 4. * H11 * zeta2
                                + 10.666666666666666 * Hm1 * zeta2
                                - 12. * Hm10 * zeta2 - 16. * Hm1m1 * zeta2
                                - 36.53333333333333 * zeta2_2
                                + 69.33333333333333 * H00 * zeta3
                                + 26.66666666666666 * H01 * zeta3
                                + 48. * H0m1 * zeta3
                                + 13.33333333333333 * H1 * zeta3
                                + 24. * Hm1 * zeta3
                                + 49.33333333333333 * zeta2 * zeta3)
                       + H0
                             * (-64.79012345679013 - 17.77777777777778 * zeta2
                                - 3.555555555555555 * zeta3
                                + x
                                      * (-421.7777777777777
                                         + x2
                                               * (-1338.3703703703702
                                                  - 105.33333333333333 * zeta2)
                                         + (-65.55555555555556 - 8. * zeta2)
                                               * zeta2
                                         - 191.11111111111111 * zeta3
                                         + x
                                               * (-230.59259259259258
                                                  + zeta2
                                                        * (-26.888888888888886
                                                           + 60.800000000000004
                                                                 * zeta2)
                                                  + 11.555555555555555 * zeta3))
                             )
                       + x2
                             * (-2053.58024691358 - 10.666666666666666 * H0000
                                + 16. * H00000 - 32. * H00010
                                + 5.333333333333333 * H0010 - 64. * H00100
                                - 8. * H0011 - 10.222222222222221 * H01
                                + 213.3333333333333 * H010
                                - 50.666666666666664 * H0100 + 32. * H01010
                                + 16. * H01011 - 41.33333333333333 * H011
                                - 16. * H01100 + 32. * H01101 - 8. * H0111
                                + 32. * H01110 - 16. * H01111 + 48. * H0m10
                                - 16. * H0m1000 - 32. * H0m1010 - 32. * H0m1m100
                                + 64. * H0m1m1m10 + 6.296296296296296 * H1
                                + 492.4444444444444 * H10
                                - 125.33333333333331 * H100
                                + 5.333333333333333 * H101 - 16. * H1010
                                - 8. * H1011 - 2.2222222222222223 * H11
                                + 5.333333333333333 * H110 + 8. * H1100
                                - 16. * H1101 - 10.666666666666666 * H111
                                - 16. * H1110 + 8. * H1111
                                + 35.55555555555556 * Hm10
                                + 5.333333333333333 * Hm100 + 8. * Hm1000
                                + 16. * Hm1010 - 10.666666666666666 * Hm1m10
                                + 16. * Hm1m100 - 32. * Hm1m1m10
                                - 254.66666666666666 * zeta2 + 8. * H001 * zeta2
                                - 34.666666666666664 * H01 * zeta2
                                + 8. * H010 * zeta2 - 8. * H011 * zeta2
                                + 24. * H0m10 * zeta2 + 32. * H0m1m1 * zeta2
                                + 2. * H1 * zeta2 - 4. * H10 * zeta2
                                + 4. * H11 * zeta2
                                - 5.333333333333333 * Hm1 * zeta2
                                - 12. * Hm10 * zeta2 - 16. * Hm1m1 * zeta2
                                + 90.66666666666666 * zeta2_2
                                + H000 * (-62.66666666666666 + 24. * zeta2)
                                + H00
                                      * (-100.88888888888889
                                         - 26.66666666666666 * zeta2
                                         - 74.66666666666666 * zeta3)
                                + 514.6666666666666 * zeta3
                                + 26.66666666666666 * H01 * zeta3
                                - 48. * H0m1 * zeta3
                                - 13.33333333333333 * H1 * zeta3
                                + 24. * Hm1 * zeta3
                                + 57.333333333333336 * zeta2 * zeta3
                                + 280. * zeta5))
                          / x)
           + a_muindep_->MuIndependentNfIndependentTerm(x)
           + 1. / (1 + nf) * massless_as3_->MuIndependentTerms(x, 1 + nf);
}
