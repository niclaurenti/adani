#include "adani/MasslessCoefficientFunction.h"
#include "adani/Constants.h"
#include "adani/SpecialFunctions.h"

#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

//==========================================================================================//
//  MasslessCoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

MasslessCoefficientFunction::MasslessCoefficientFunction(
    const int &order, const char &kind, const char &channel
)
    : CoefficientFunction(order, kind, channel) {
    SetFunctions();
}

//==========================================================================================//
//  MasslessCoefficientFunction: mu-dependent + mu-independent terms
//------------------------------------------------------------------------------------------//

double MasslessCoefficientFunction::fx(
    double x, double m2Q2, double m2mu2, int nf
) const {
    return MuIndependentTerms(x, nf) + MuDependentTerms(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  MasslessCoefficientFunction: mu-dependent terms
//------------------------------------------------------------------------------------------//

double MasslessCoefficientFunction::MuDependentTerms(
    double /*x*/, double /*m2Q2*/, double /*m2mu2*/, int /*nf*/
) const {
    cout << "Error: mu dependent terms of the massless coefficient functions "
            "are not implemented!"
         << endl;
    exit(-1);
}

//==========================================================================================//
//  MasslessCoefficientFunction: done beacuse it is required by the abstract
//  base class CoefficientFunction
//------------------------------------------------------------------------------------------//

double MasslessCoefficientFunction::
    MuIndependentTerms(double /*x*/, double /*m2Q2*/, int /*nf*/) const {
    cout << "Error: massless coefficient functions do not depend on m^2/Q^2!"
         << endl;
    cout << "Call MasslessCoefficientFunction::MuIndependentTerms(double x, "
            "int nf)"
         << endl;
    exit(-1);
}

//==========================================================================================//
//  MasslessCoefficientFunction: mu-independent terms
//------------------------------------------------------------------------------------------//

double MasslessCoefficientFunction::MuIndependentTerms(double x, int nf) const {

    return (this->*mu_indep_)(x, nf);
}

//==========================================================================================//
//  MasslessCoefficientFunction: done beacuse it is required by the abstract
//  base class CoefficientFunction. It returns three identical values
//------------------------------------------------------------------------------------------//

Value MasslessCoefficientFunction::fxBand(
    double x, double m2Q2, double m2mu2, int nf
) const {
    return Value(fx(x, m2Q2, m2mu2, nf));
}

//==========================================================================================//
//  MasslessCoefficientFunction: function that sets the pointer mu_indep_ to the
//  correct function
//------------------------------------------------------------------------------------------//

void MasslessCoefficientFunction::SetFunctions() {
    if (GetOrder() == 1 && GetKind() == '2' && GetChannel() == 'g')
        mu_indep_ = &MasslessCoefficientFunction::C2_g1_massless;
    else if (GetOrder() == 1 && GetKind() == 'L' && GetChannel() == 'g')
        mu_indep_ = &MasslessCoefficientFunction::CL_g1_massless;

    else if (GetOrder() == 2 && GetKind() == '2' && GetChannel() == 'g')
        mu_indep_ = &MasslessCoefficientFunction::C2_g2_massless;
    else if (GetOrder() == 2 && GetKind() == '2' && GetChannel() == 'q')
        mu_indep_ = &MasslessCoefficientFunction::C2_ps2_massless;
    else if (GetOrder() == 2 && GetKind() == 'L' && GetChannel() == 'g')
        mu_indep_ = &MasslessCoefficientFunction::CL_g2_massless;
    else if (GetOrder() == 2 && GetKind() == 'L' && GetChannel() == 'q')
        mu_indep_ = &MasslessCoefficientFunction::CL_ps2_massless;

    else if (GetOrder() == 3 && GetKind() == '2' && GetChannel() == 'g')
        mu_indep_ = &MasslessCoefficientFunction::C2_g3_massless;
    else if (GetOrder() == 3 && GetKind() == '2' && GetChannel() == 'q')
        mu_indep_ = &MasslessCoefficientFunction::C2_ps3_massless;
    else if (GetOrder() == 3 && GetKind() == 'L' && GetChannel() == 'g')
        mu_indep_ = &MasslessCoefficientFunction::CL_g3_massless;
    else if (GetOrder() == 3 && GetKind() == 'L' && GetChannel() == 'q')
        mu_indep_ = &MasslessCoefficientFunction::CL_ps3_massless;
    else {
        cout << "Error: something has gone wrong in "
                "MasslessCoefficientFunction::SetFunctions!"
             << endl;
        exit(-1);
    }
}

//==========================================================================================//
//  Massless gluon coefficient functions for F2 at
//  O(as) for mu=Q.
//
//  Eq. (4.4) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double MasslessCoefficientFunction::C2_g1_massless(double x, int nf) const {

    return 4. * nf * TR
           * (-8. * x * x + 8. * x - 1.
              + log((1. - x) / x) * (2. * x * x - 2. * x + 1.));
}

//==========================================================================================//
//  Massless gluon coefficient functions for FL at
//  O(as) for mu=Q.
//
// Eq. (3) of Ref. [arXiv:hep-ph/0411112v2]
//------------------------------------------------------------------------------------------//

double MasslessCoefficientFunction::CL_g1_massless(double x, int nf) const {
    return 16. * nf * TR * x * (1. - x);
}

//==========================================================================================//
//  Massless gluon coefficient functions for F2 at
//  O(as^2) for mu=Q.
//
//  Eq. (B.6) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double MasslessCoefficientFunction::C2_g2_massless(double x, int nf) const {

    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;
    double x5 = x4 * x;

    double Hm1 = H_m1(x);
    double H0 = H_0(x);
    double H1 = H_1(x);

    double H00 = H_00(x);
    double H10 = H_10(x);
    double Hm10 = H_m10(x);
    double H01 = H_01(x);
    double H11 = H_11(x);

    double Hm1m10 = H_m1m10(x);
    double H0m10 = H_0m10(x);
    double Hm100 = H_m100(x);
    double H000 = H_000(x);
    double H100 = H_100(x);
    double H010 = H_010(x);
    double H110 = H_110(x);
    double Hm101 = H_m101(x);
    double H001 = H_001(x);
    double H101 = H_101(x);
    double H011 = H_011(x);
    double H111 = H_111(x);

    return (CA * nf
            * (12.740740740740742 + 5.333333333333333 * H11
               - 5.333333333333333 * Hm10
               + (26.555555555555554 + 58.00000000000001 * H0 - 2. * H00
                  + 20. * H000 + 8. * H001 - 8. * H01)
                     * x
               + H10
                     * (5.333333333333333
                        + x * (-4. + (80. - 89.33333333333333 * x) * x))
               - 5.333333333333333 * zeta2
               + H1
                     * (-11.555555555555555
                        + x
                              * (20.666666666666664 + 8. * zeta2
                                 + x
                                       * (151.33333333333331 - 16. * zeta2
                                          + x
                                                * (-174.44444444444443
                                                   + 8. * zeta2))))
               + x
                     * (-4. * H11 - 12. * H110 - 4. * H111 - 24. * Hm10
                        + 8. * Hm100 + 8. * Hm101 + 8. * Hm1m10
                        + (119.11111111111111 + 194.66666666666666 * H0
                           + 176. * H00 + 56. * H000 + 64. * H001 + 144. * H01
                           + 48. * H010 + 48. * H011)
                              * x
                        + 72. * H11 * x + 24. * H110 * x + 8. * H111 * x
                        + 16. * Hm100 * x + 16. * Hm101 * x + 16. * Hm1m10 * x
                        - 166.4074074074074 * x2 - 232.2222222222222 * H0 * x2
                        - 129.33333333333331 * H00 * x2 - 16. * H001 * x2
                        - 148. * H01 * x2 - 16. * H010 * x2 - 16. * H011 * x2
                        + 16. * H0m10 * x2 - 81.33333333333333 * H11 * x2
                        - 24. * H110 * x2 - 8. * H111 * x2
                        + 26.666666666666664 * Hm10 * x2 + 24. * Hm100 * x2
                        + 16. * Hm101 * x2 + H100 * (-12. + (24. - 16. * x) * x)
                        + H101 * (-4. + (8. - 8. * x) * x) + 8. * zeta2
                        - 8. * H0 * zeta2 - 4. * Hm1 * zeta2 - 144. * x * zeta2
                        - 64. * H0 * x * zeta2 - 8. * Hm1 * x * zeta2
                        + 148. * x2 * zeta2 + 16. * H0 * x2 * zeta2
                        - 16. * Hm1 * x2 * zeta2 + 4. * zeta3 - 48. * x * zeta3
                        + 24. * x2 * zeta3)))
               / x
           + (48. * CF * nf
              * (Hm10
                     * (0.011111111111111112 + x2 + 0.4444444444444444 * x3
                        + 0.4000000000000001 * x5)
                 + x
                       * (0.011111111111111112
                          + H0
                                * (-0.011111111111111112
                                   + x
                                         * (-0.3277777777777778
                                            + x
                                                  * (0.4708333333333334
                                                     + x * (-0.9 + zeta2)
                                                     - 0.6666666666666666
                                                           * zeta2)
                                            + 0.3333333333333333 * zeta2))
                          + x4
                                * (-0.4000000000000001 * H00
                                   + 0.4000000000000001 * zeta2)
                          + x2
                                * (0.9958333333333335 + 0.3055555555555555 * H00
                                   + 0.4166666666666667 * H000
                                   + 0.6666666666666666 * H001
                                   + 1.1666666666666667 * H01 + 0.5 * H010
                                   + 0.6666666666666666 * H011
                                   + 0.8333333333333334 * H1
                                   + 1.6666666666666667 * H10
                                   + 0.16666666666666666 * H100 + H101
                                   + 1.6666666666666667 * H11
                                   + 0.6666666666666666 * H110
                                   + 0.8333333333333334 * H111
                                   + 0.6666666666666666 * Hm100
                                   - 1.3333333333333333 * Hm1m10
                                   - 0.7222222222222222 * zeta2
                                   - 0.3333333333333333 * H1 * zeta2
                                   - 0.6666666666666666 * Hm1 * zeta2)
                          + x
                                * (-0.8986111111111111 - 0.0625 * H00
                                   - 0.20833333333333334 * H000
                                   - 0.3333333333333333 * H001
                                   - 0.3333333333333333 * H01 - 0.25 * H010
                                   - 0.3333333333333333 * H011
                                   + 0.6666666666666666 * H0m10
                                   - 0.2916666666666667 * H1
                                   - 0.5416666666666666 * H10
                                   - 0.08333333333333333 * H100 - 0.5 * H101
                                   - 0.5416666666666666 * H11
                                   - 0.3333333333333333 * H110
                                   - 0.4166666666666667 * H111
                                   + 0.3333333333333333 * Hm100
                                   - 0.6666666666666666 * Hm1m10
                                   + 0.3333333333333333 * zeta2
                                   + 0.16666666666666666 * H1 * zeta2
                                   - 0.3333333333333333 * Hm1 * zeta2
                                   + 0.6666666666666666 * zeta3)
                          + x3
                                * (-0.15 - 1.5 * H00 - 0.8333333333333334 * H000
                                   - 1. * H001 - 1.5 * H01
                                   - 0.6666666666666666 * H010
                                   - 0.8333333333333334 * H011
                                   + 0.6666666666666666 * H0m10 - 0.5 * H1
                                   - 1.5 * H10 - 0.5 * H100 - 1. * H101
                                   - 1.5 * H11 - 0.6666666666666666 * H110
                                   - 0.8333333333333334 * H111
                                   + 0.3333333333333333 * Hm100
                                   - 0.6666666666666666 * Hm1m10 + 1.5 * zeta2
                                   + 0.6666666666666666 * H1 * zeta2
                                   - 0.3333333333333333 * Hm1 * zeta2
                                   + 1.5 * zeta3))))
                 / x2;
}

//==========================================================================================//
//  Massless quark coefficient functions for F2 at
//  O(as^2) for mu=Q.
//
//  Eq. (B.7) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double MasslessCoefficientFunction::C2_ps2_massless(double x, int nf) const {

    double x2 = x * x;

    double H0 = H_0(x);
    double H1 = H_1(x);

    double H00 = H_00(x);
    double H10 = H_10(x);
    double Hm10 = H_m10(x);
    double H01 = H_01(x);
    double H11 = H_11(x);

    double H000 = H_000(x);
    double H010 = H_010(x);
    double H001 = H_001(x);
    double H011 = H_011(x);

    return (20. * CF * nf
            * (0.6370370370370371 + 0.26666666666666666 * H11
               - 0.26666666666666666 * Hm10
               + (0.8777777777777777 + 2.8 * H0 - 0.1 * H00 + H000 + 0.8 * H001
                  + 0.4 * H010 + 0.4 * H011)
                     * x
               + H10
                     * (0.26666666666666666
                        + x * (0.2 + (-0.2 - 0.26666666666666666 * x) * x))
               + H1
                     * (-0.5777777777777777
                        + x
                              * (1.7333333333333332
                                 + (-1.3333333333333333
                                    + 0.17777777777777776 * x)
                                       * x))
               - 0.26666666666666666 * zeta2
               + x
                     * ((-2.344444444444444 - 1.4666666666666666 * H0
                         + 1.4999999999999998 * H00 + 1. * H000 + 0.8 * H001
                         + 0.4 * H010 + 0.4 * H011)
                            * x
                        + 0.8296296296296296 * x2 - 0.711111111111111 * H0 * x2
                        - 1.0666666666666667 * H00 * x2 - 0.8 * H01 * x2
                        + Hm10 * (-0.8 + (-0.8 - 0.26666666666666666 * x) * x)
                        + H11 * (0.2 + (-0.2 - 0.26666666666666666 * x) * x)
                        - 0.8 * H0 * zeta2 - 0.8 * x * zeta2
                        - 0.8 * H0 * x * zeta2 + 0.8 * x2 * zeta2 - 0.4 * zeta3
                        - 0.4 * x * zeta3)))
           / x;
}

//==========================================================================================//
//  Massless gluon coefficient functions for FL at
//  O(as^2) for mu=Q.
//
//  Eq. (B.14) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double MasslessCoefficientFunction::CL_g2_massless(double x, int nf) const {

    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;
    double x5 = x4 * x;

    double H0 = H_0(x);
    double H1 = H_1(x);

    double H00 = H_00(x);
    double H10 = H_10(x);
    double Hm10 = H_m10(x);
    double H01 = H_01(x);
    double H11 = H_11(x);

    return CA * nf
               * (5.333333333333333 + 16. * H0 - 1.7777777777777777 / x
                  - (5.333333333333333 * H1) / x
                  + H1 * (16. + (144. - 154.66666666666666 * x) * x)
                  + x
                        * (90.66666666666666 + 128. * H0 + 96. * H00 + 96. * H01
                           + 32. * H10 + 32. * H11 + 32. * Hm10 - 64. * zeta2)
                  + x2
                        * (-94.22222222222221 - 208. * H0 - 32. * H01
                           - 32. * H10 - 32. * H11 + 32. * Hm10 + 32. * zeta2))
           + (CF * nf
              * (Hm10
                     * (2.1333333333333333 - 10.666666666666666 * x3 + 12.8 * x5
                     )
                 + x
                       * (2.1333333333333333
                          + H0
                                * (-2.1333333333333333
                                   + x
                                         * (-6.933333333333334
                                            + x
                                                  * (-41.6
                                                     + 19.200000000000003 * x)))
                          + x
                                * (-8.533333333333333
                                   + H1 * (-8. + x * (-24. + 32. * x))
                                   + x
                                         * (-60.800000000000004 - 16. * H01
                                            + H00
                                                  * (-21.333333333333332
                                                     - 12.8 * x2)
                                            + 5.333333333333333 * zeta2
                                            + x * (67.2 + 12.8 * x * zeta2))))))
                 / x2;
}

//==========================================================================================//
//  Massless gluon coefficient functions for FL at
//  O(as^2) for mu=Q.
//
//  Eq. (B.15) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double MasslessCoefficientFunction::CL_ps2_massless(double x, int nf) const {

    double x2 = x * x;
    double x3 = x2 * x;

    double H0 = H_0(x);
    double H1 = H_1(x);

    double H00 = H_00(x);
    double H01 = H_01(x);

    return (CF * nf
            * (-1.7777777777777777
               + H1 * (-5.333333333333333 + 16. * x - 10.666666666666666 * x3)
               + x
                     * (5.333333333333333 + H0 * (16. + (-16. - 32. * x) * x)
                        + x
                              * (-21.333333333333332 + 32. * H00 + 16. * H01
                                 + 17.77777777777778 * x - 16. * zeta2))))
           / x;
}

//==========================================================================================//
//  Massless gluon coefficient functions for F2 at
//  O(as^3) for mu=Q.
//  The term fl_g_11 is put to zero for the reason explained
//  in page 15 of arXiv:1205.5727
//
//  Eq. (B.9) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double MasslessCoefficientFunction::C2_g3_massless(double x, int nf) const {

    // the commented varibles are needed for the flav term

    double x2 = x * x;
    double x3 = x2 * x;

    // double xm = 1 - x;
    // double xm2 = xm * xm;
    // double xm3 = xm2 * xm;
    // double xm4 = xm3 * xm;
    // double xm5 = xm4 * xm;

    // double xp = 1 + x;
    // double xp2 = xp * xp;
    // double xp3 = xp2 * xp;
    // double xp4 = xp3 * xp;
    // double xp5 = xp4 * xp;

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
    const double H10m1 = Hr3[5];
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
    const double H10m1m10 = Hr5[86];
    const double Hm1m10m10 = Hr5[90];
    const double H0m10m10 = Hr5[91];
    const double Hm100m10 = Hr5[93];
    const double H000m10 = Hr5[94];
    const double H100m10 = Hr5[95];
    const double H010m10 = Hr5[97];
    const double H110m10 = Hr5[98];
    const double Hm1m1m100 = Hr5[108];
    const double H0m1m100 = Hr5[109];
    const double Hm10m100 = Hr5[111];
    const double H00m100 = Hr5[112];
    const double H10m100 = Hr5[113];
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
    const double H00m101 = Hr5[193];
    const double H10m101 = Hr5[194];
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
    const double H0m1011 = Hr5[226];
    const double Hm10011 = Hr5[228];
    const double H00011 = Hr5[229];
    const double H10011 = Hr5[230];
    const double H01011 = Hr5[232];
    const double H11011 = Hr5[233];
    const double Hm10111 = Hr5[237];
    const double H00111 = Hr5[238];
    const double H10111 = Hr5[239];
    const double H01111 = Hr5[241];
    const double H11111 = Hr5[242];

    delete[] Hr1;
    delete[] Hr2;
    delete[] Hr3;
    delete[] Hr4;
    delete[] Hr5;

    // double flav = 5. / 48 * fl11g(nf) * nf * nf;
    // this contribution is neglected
    // If you want to use it please check the 5/48 since
    // I'm not sure about it

    return /*flav
               * (-119.46666666666667 * H0001 + 494.93333333333334 * H0001 * x
                  + 1885.8666666666666 * H001 * x + 102.4 * H00m10 * x
                  - 1862.2577777777778 * H01 * x
                  - 1109.3333333333333 * H0100 * x + 2560 * H01001 * x
                  - 2560 * H01100 * x - 529.0666666666666 * H0m10 * x
                  + 614.4000000000001 * H0m100 * x
                  + 1809.0666666666666 * H0m101 * x
                  + 580.2666666666667 * H0m1m10 * x - 2666.951111111111 * H1 * x
                  - 2218.6666666666665 * H10 * x - 460.8 * H100 * x
                  + 853.3333333333333 * H1000 * x + 2560 * H10001 * x
                  + 4608 * H1001 * x - 2560 * H10100 * x - 512 * H10m10 * x
                  - 4437.333333333333 * H11 * x - 3413.333333333333 * H1100 * x
                  + 5120 * H11001 * x - 5120 * H11100 * x
                  - 2529.5644444444442 * Hm10 * x
                  - 1425.0666666666666 * Hm100 * x
                  - 853.3333333333333 * Hm1000 * x
                  + 1604.2666666666667 * Hm101 * x
                  + 1194.6666666666665 * Hm10m10 * x
                  + 4454.400000000001 * Hm1m10 * x
                  + 1194.6666666666665 * Hm1m100 * x
                  - 2389.333333333333 * Hm1m1m10 * x + 692.9066666666668 * x2
                  - 2926.9333333333334 * H0001 * x2
                  - 4017.7777777777774 * H001 * x2
                  + 1297.0666666666666 * H00m10 * x2
                  + 2234.0266666666666 * H01 * x2 + 3456 * H0100 * x2
                  - 2304 * H01001 * x2 + 2304 * H01100 * x2
                  - 3147.3777777777777 * H0m10 * x2
                  - 1109.3333333333333 * H0m100 * x2 - 1638.4 * H0m101 * x2
                  + 580.2666666666667 * H0m1m10 * x2
                  + 2199.8933333333334 * H1 * x2 + 1941.3333333333333 * H10 * x2
                  + 1510.4 * H100 * x2 - 938.6666666666666 * H1000 * x2
                  - 2304 * H10001 * x2 - 4565.333333333333 * H1001 * x2
                  + 2304 * H10100 * x2 + 768 * H10m10 * x2
                  + 3882.6666666666665 * H11 * x2 + 3456 * H1100 * x2
                  - 4608 * H11001 * x2 + 4608 * H11100 * x2
                  - 1120.4266666666667 * Hm10 * x2 - 640 * Hm100 * x2
                  - 938.6666666666666 * Hm1000 * x2
                  + 1867.3777777777777 * Hm101 * x2
                  + 1109.3333333333333 * Hm10m10 * x2
                  + 3147.3777777777777 * Hm1m10 * x2
                  + 1109.3333333333333 * Hm1m100 * x2
                  - 2218.6666666666665 * Hm1m1m10 * x2 + 501.76 * H000 * x3
                  + 2150.4 * H0001 * x3 + 501.76 * H001 * x3
                  - 2150.4 * H0100 * x3 - 501.76 * H0m10 * x3
                  + 2150.4 * H1001 * x3 - 2150.4 * H1100 * x3
                  - 230.4 * Hm10 * x3 - 501.76 * Hm100 * x3
                  - 501.76 * Hm101 * x3 + 501.76 * Hm1m10 * x3
                  + (34.13333333333333 * H0001) / xm5
                  - (34.13333333333333 * H0m100) / xm5
                  - (68.26666666666667 * H0m101) / xm5
                  + (17.066666666666666 * H000) / xm3
                  - (17.066666666666666 * H0000) / xm3
                  - (34.13333333333333 * H0001) / xm3
                  + (34.13333333333333 * H001) / xm3
                  + (34.13333333333333 * H0m100) / xm3
                  + (68.26666666666667 * H0m101) / xm3
                  - (34.13333333333333 * Hm100) / xm3
                  - (68.26666666666667 * Hm101) / xm3
                  - (8.533333333333333 * H000) / xm2
                  + (59.733333333333334 * H0000) / xm2
                  + (119.46666666666667 * H0001) / xm2
                  - (17.066666666666666 * H001) / xm2
                  - (119.46666666666667 * H0m100) / xm2
                  - (238.93333333333334 * H0m101) / xm2
                  + (17.066666666666666 * Hm100) / xm2
                  + (34.13333333333333 * Hm101) / xm2
                  + (56.888888888888886 * H000) / xm
                  + (113.77777777777777 * H001) / xm
                  + (5.688888888888889 * H01) / xm
                  - (113.77777777777777 * Hm100) / xm
                  - (227.55555555555554 * Hm101) / xm
                  - (17.066666666666666 * H0000) / xp5
                  + (34.13333333333333 * H00m10) / xp5
                  - (17.066666666666666 * H000) / xp3
                  + (17.066666666666666 * H0000) / xp3
                  - (34.13333333333333 * H00m10) / xp3
                  + (34.13333333333333 * H0m10) / xp3
                  + (8.533333333333333 * H000) / xp2
                  - (59.733333333333334 * H0000) / xp2
                  + (119.46666666666667 * H00m10) / xp2
                  - (17.066666666666666 * H0m10) / xp2
                  + (17.066666666666666 * Hm10) / xp2 - 2.8444444444444446 / xp
                  - (56.888888888888886 * H000) / xp
                  + (113.77777777777777 * H0m10) / xp + 230.4 * H01 * zeta2
                  + 512 * H010 * zeta2 + 452.26666666666665 * H0m1 * zeta2
                  + 769.4222222222222 * H1 * zeta2 + 1920 * H10 * zeta2
                  + 512 * H100 * zeta2 + 341.3333333333333 * H11 * zeta2
                  + 1024 * H110 * zeta2 + 1406.5777777777778 * Hm1 * zeta2
                  + 170.66666666666666 * Hm10 * zeta2
                  - 341.3333333333333 * Hm1m1 * zeta2
                  + (6.9688888888888885 * H1 * zeta2) / x2
                  + (59.733333333333334 * H10 * zeta2) / x2
                  + (20.906666666666666 * Hm1 * zeta2) / x2
                  + (31.857777777777777 * zeta2) / x
                  - 667.3066666666667 * x * zeta2
                  - 290.1333333333333 * H01 * x * zeta2
                  - 2560 * H010 * x * zeta2
                  - 1518.9333333333334 * H0m1 * x * zeta2
                  - 2227.2000000000003 * H1 * x * zeta2
                  - 5461.333333333333 * H10 * x * zeta2
                  - 2560 * H100 * x * zeta2
                  - 1194.6666666666665 * H11 * x * zeta2
                  - 5120 * H110 * x * zeta2
                  + 622.9333333333333 * Hm1 * x * zeta2
                  + 853.3333333333333 * Hm10 * x * zeta2
                  - 1194.6666666666665 * Hm1m1 * x * zeta2
                  - 2234.0266666666666 * x2 * zeta2
                  + 290.1333333333333 * H01 * x2 * zeta2
                  + 2304 * H010 * x2 * zeta2
                  + 1928.5333333333333 * H0m1 * x2 * zeta2
                  + 1573.6888888888889 * H1 * x2 * zeta2
                  + 5504 * H10 * x2 * zeta2 + 2304 * H100 * x2 * zeta2
                  + 1109.3333333333333 * H11 * x2 * zeta2
                  + 4608 * H110 * x2 * zeta2
                  - 293.6888888888889 * Hm1 * x2 * zeta2
                  + 938.6666666666666 * Hm10 * x2 * zeta2
                  - 1109.3333333333333 * Hm1m1 * x2 * zeta2 - 230.4 * x3 * zeta2
                  - 250.88 * H1 * x3 * zeta2 - 2150.4 * H10 * x3 * zeta2
                  + 752.64 * Hm1 * x3 * zeta2
                  + (68.26666666666667 * H0m1 * zeta2) / xm5
                  - (68.26666666666667 * H0m1 * zeta2) / xm3
                  + (68.26666666666667 * Hm1 * zeta2) / xm3
                  - (8.533333333333333 * zeta2) / xm2
                  + (238.93333333333334 * H0m1 * zeta2) / xm2
                  - (34.13333333333333 * Hm1 * zeta2) / xm2
                  - (5.688888888888889 * zeta2) / xm
                  + (227.55555555555554 * Hm1 * zeta2) / xm
                  + (8.533333333333333 * zeta2) / xp2
                  - 101.54666666666667 * zeta2_2 - 409.6 * H1 * zeta2_2
                  + 1419.9466666666667 * x * zeta2_2 + 2048 * H1 * x * zeta2_2
                  - 1508.6933333333334 * x2 * zeta2_2
                  - 1843.2 * H1 * x2 * zeta2_2 + 1720.32 * x3 * zeta2_2
                  - (6.826666666666667 * zeta2_2) / xm5
                  + (6.826666666666667 * zeta2_2) / xm3
                  - (23.893333333333334 * zeta2_2) / xm2
                  + (35.84 * zeta2_2) / xp5 - (35.84 * zeta2_2) / xp3
                  + (125.44 * zeta2_2) / xp2
                  + H00
                        * (28.58666666666667 + 13.937777777777777 / x
                           + 2.8444444444444446 / xm - 8.533333333333333 / xp2
                           + x3 * (230.4 - 2150.4 * zeta2)
                           + x
                                 * (2061.6533333333336
                                    - 392.5333333333333 * zeta2)
                           + 119.46666666666667 * zeta2
                           + ((-51.2 + (51.2 - 179.20000000000002 * xm) * xm2)
                              * zeta2)
                                 / xm5
                           + ((17.066666666666666
                               + xp2
                                     * (-17.066666666666666
                                        + 59.733333333333334 * xp))
                              * zeta2)
                                 / xp5
                           + x2 * (911.36 + 2926.9333333333334 * zeta2))
                  - 512 * H01 * zeta3 - 896 * H1 * zeta3 - 512 * H10 * zeta3
                  - 1024 * H11 * zeta3 + 341.3333333333333 * Hm1 * zeta3
                  - (59.733333333333334 * H1 * zeta3) / x2
                  + (59.733333333333334 * zeta3) / x
                  - 1373.8666666666666 * x * zeta3 + 2560 * H01 * x * zeta3
                  + 853.3333333333333 * H1 * x * zeta3 + 2560 * H10 * x * zeta3
                  + 5120 * H11 * x * zeta3
                  + 1194.6666666666665 * Hm1 * x * zeta3
                  - 1408.7111111111112 * x2 * zeta3 - 2304 * H01 * x2 * zeta3
                  - 640 * H1 * x2 * zeta3 - 2304 * H10 * x2 * zeta3
                  - 4608 * H11 * x2 * zeta3
                  + 1109.3333333333333 * Hm1 * x2 * zeta3 - 1254.4 * x3 * zeta3
                  + 2150.4 * H1 * x3 * zeta3 - (51.2 * zeta3) / xm3
                  + (25.6 * zeta3) / xm2 - (170.66666666666666 * zeta3) / xm
                  + (34.13333333333333 * zeta3) / xp3
                  - (17.066666666666666 * zeta3) / xp2
                  + (113.77777777777777 * zeta3) / xp
                  + (59.73333333333334 * H1100 - 13.937777777777777 * Hm100
                     - 13.937777777777777 * Hm101 + 13.937777777777777 * Hm1m10
                     + Hm10 * (-6.4 - 13.937777777777777 * x)
                     + x
                           * (7.537777777777777 + 59.73333333333334 * H001
                              - 45.79555555555555 * H01 + 45.79555555555555 * H1
                              - 59.73333333333334 * H100
                              + (-16.213333333333335 - 100.97777777777779 * H001
                                 - 119.46666666666668 * H00m10
                                 - 27.59111111111111 * H01 - 512. * H01001
                                 + 512. * H01100 - 770.8444444444444 * H0m10
                                 - 341.3333333333333 * H0m100
                                 - 221.86666666666667 * H0m101
                                 + 460.79999999999995 * H0m1m10
                                 + 358.68444444444447 * H1
                                 + 405.3333333333333 * H10
                                 - 349.8666666666667 * H100
                                 - 170.66666666666666 * H1000 - 512. * H10001)
                                    * x)
                     + H1001 * (-59.73333333333335 - 1749.3333333333333 * x2)
                     + x2
                           * (512. * H10100 + 810.6666666666666 * H11
                              + 1408. * H1100 - 1024. * H11001 + 1024. * H11100
                              - 2100.7644444444445 * Hm10 - 1088. * Hm100
                              - 170.66666666666666 * Hm1000
                              - 637.1555555555556 * Hm101
                              + 341.3333333333333 * Hm10m10
                              + 1538.8444444444444 * Hm1m10
                              + 341.3333333333333 * Hm1m100
                              - 682.6666666666666 * Hm1m1m10
                              + (-602.4533333333334 + 529.0666666666666 * H000
                                 - 102.4 * H0000)
                                    * x
                              + (17.066666666666666 * H0000) / xm5
                              + 27.59111111111111 * zeta2
                              + 95.28888888888889 * zeta3))
                        / x2
                  + H0
                        * ((-7.5377777777777775 - 59.733333333333334 * zeta2)
                               / x
                           + x2
                                 * (1928.5333333333333
                                    + 4017.7777777777774 * zeta2
                                    - 2568.5333333333333 * zeta3)
                           - (51.2 * zeta3) / xm5
                           + (34.13333333333333 * zeta3) / xp5
                           + x
                                 * (-383.43111111111114
                                    - 2414.9333333333334 * zeta2
                                    + 1877.3333333333333 * zeta3)
                           + x3 * (-1003.52 * zeta2 + 2150.4 * zeta3)
                           + (17.066666666666666 * zeta2
                              - 34.13333333333333 * zeta3
                              + xp
                                    * (-8.533333333333333
                                       - 8.533333333333333 * zeta2
                                       + 56.888888888888886 * xp * zeta2
                                       + 119.46666666666667 * zeta3))
                                 / xp3
                           + (-51.2 * zeta2 + 51.2 * zeta3
                              + xm
                                    * (25.6 * zeta2 - 179.20000000000002 * zeta3
                                       + xm
                                             * (-170.66666666666666 * zeta2
                                                + xm
                                                      * (0.8533333333333334
                                                         + 100.97777777777779
                                                               * zeta2
                                                         + 59.73333333333334
                                                               * zeta3))))
                                 / xm3))*/

           + CA * CF * nf
                 * (-1112.9304732510288 + 134.05185185185186 * H000
                    - 181.77777777777777 * H0000 + 120 * H00000
                    + 218.66666666666666 * H00001 - 94.66666666666666 * H0001
                    + 157.33333333333331 * H00010 + 168 * H00011
                    - 725.3333333333333 * H000m10 + 29.644444444444446 * H001
                    - 81.33333333333333 * H0010 - 101.33333333333333 * H00100
                    + 93.33333333333333 * H00101 - 118.66666666666666 * H0011
                    + 24 * H00110 + 82.66666666666666 * H00111
                    - 842.6666666666666 * H00m10 - 314.66666666666663 * H00m100
                    + 250.66666666666666 * H00m101
                    + 58.666666666666664 * H00m1m10 - 428.3930864197531 * H01
                    - 24.48888888888889 * H010 - 352 * H0100
                    - 114.66666666666666 * H01000 - 880 * H01001
                    - 150.66666666666666 * H0101 - 26.666666666666664 * H01010
                    + 37.33333333333333 * H01011 - 501.3333333333333 * H010m10
                    + 25.881481481481483 * H011 - 169.33333333333331 * H0110
                    + 808 * H01100 - 34.666666666666664 * H01101
                    - 145.77777777777777 * H0111 - 26.666666666666664 * H01110
                    + 40 * H01111 - 322.66666666666663 * H0m100
                    + 490.66666666666663 * H0m1000
                    + 1338.6666666666665 * H0m1001 + 293.3333333333333 * H0m101
                    + 256. * H0m1010 + 288. * H0m1011
                    + 69.33333333333333 * H0m1m10 - 896.0000000000001 * H0m1m100
                    - 1717.333333333333 * H0m1m101
                    + 149.33333333333331 * H0m1m1m10 - 900.7326748971193 * H1
                    - 358.4469135802469 * H10 - 550.6074074074074 * H100
                    - 354.22222222222223 * H1000 + 312. * H10000
                    + 301.3333333333333 * H10001 - 2369.333333333333 * H1001
                    - 26.666666666666664 * H10010 + 98.66666666666664 * H10011
                    - 293.3333333333333 * H100m10 - 44.22222222222222 * H101
                    - 302.66666666666663 * H1010 + 42.666666666666664 * H10100
                    - 114.66666666666667 * H10101 - 209.33333333333331 * H1011
                    - 128. * H10110 - 1013.3333333333333 * H10m10
                    - 272. * H10m100 - 426.66666666666663 * H10m101
                    - 181.33333333333331 * H10m1m10 - 170.61975308641976 * H11
                    - 34. * H110 + 1477.3333333333333 * H1100
                    + 285.3333333333333 * H11000 - 383.99999999999994 * H11001
                    - 221.33333333333331 * H1101 - 184. * H11010 - 32. * H11011
                    - 469.3333333333333 * H110m10 + 75.25925925925925 * H111
                    - 220. * H1110 + 277.3333333333333 * H11100 - 104. * H11101
                    - 116.44444444444444 * H1111 - 136. * H11110
                    - 4491.182222222223 * Hm10 - 2615.555555555555 * Hm100
                    + 795.9999999999999 * Hm1000 + 552. * Hm10000
                    + 309.3333333333333 * Hm10001 + 2250.6666666666665 * Hm1001
                    + 208. * Hm10010 + 240. * Hm10011
                    - 298.66666666666663 * Hm100m10 - 1102.6666666666665 * Hm101
                    + 402.66666666666663 * Hm1010 + 570.6666666666666 * Hm10100
                    + 90.66666666666666 * Hm10101 + 450.66666666666663 * Hm1011
                    + 80. * Hm10110 + 85.33333333333333 * Hm10111
                    + 170.66666666666666 * Hm10m10
                    - 1058.6666666666665 * Hm10m100 - 1184. * Hm10m101
                    + 432. * Hm10m1m10 + 1521.5111111111112 * Hm1m10
                    - 1464. * Hm1m100 - 968. * Hm1m1000
                    - 1469.3333333333333 * Hm1m1001
                    - 3053.333333333333 * Hm1m101
                    - 250.66666666666666 * Hm1m1010 - 288. * Hm1m1011
                    + 453.33333333333337 * Hm1m10m10
                    - 269.3333333333333 * Hm1m1m10
                    + 1781.3333333333333 * Hm1m1m100
                    + 2000.0000000000002 * Hm1m1m101
                    - 474.66666666666663 * Hm1m1m1m10
                    - (12.764444444444445 * H00) / x - (12.8 * H000) / x
                    + (1.0666666666666667 * H001) / x
                    - (16.675555555555555 * H01) / x
                    - (2.1333333333333333 * H010) / x
                    - (2.1333333333333333 * H011) / x
                    + (8.533333333333333 * H0m10) / x
                    + (11.024526748971194 * H1) / x
                    + (15.071604938271605 * H10) / x
                    + (19.081481481481482 * H100) / x
                    - (46.22222222222222 * H1000) / x
                    + (28.444444444444443 * H1001) / x
                    - (45.03703703703704 * H101) / x
                    + (78.22222222222221 * H1010) / x
                    + (85.33333333333333 * H1011) / x
                    + (5.738271604938271 * H11) / x
                    - (13.037037037037036 * H110) / x
                    + (71.11111111111111 * H1100) / x
                    + (78.22222222222221 * H1101) / x
                    - (29.037037037037035 * H111) / x
                    + (78.22222222222221 * H1110) / x
                    + (81.77777777777777 * H1111) / x
                    - (8.924444444444445 * Hm10) / x
                    - (62.22222222222222 * Hm100) / x
                    + (85.33333333333333 * Hm1000) / x
                    + (42.666666666666664 * Hm1001) / x
                    - (10.666666666666666 * Hm101) / x
                    - (85.33333333333333 * Hm10m10) / x
                    + (122.31111111111112 * Hm1m10) / x
                    - (170.66666666666666 * Hm1m100) / x
                    - (85.33333333333333 * Hm1m101) / x
                    + (85.33333333333333 * Hm1m1m10) / x
                    - 3590.3574897119342 * x + 1395.8074074074075 * H000 * x
                    - 207.1111111111111 * H0000 * x
                    + 389.3333333333333 * H00000 * x
                    + 693.3333333333333 * H00001 * x + 400 * H0001 * x
                    + 762.6666666666666 * H00010 * x
                    + 869.3333333333333 * H00011 * x
                    + 42.666666666666664 * H000m10 * x
                    + 1185.2444444444445 * H001 * x + 816 * H0010 * x
                    + 298.66666666666663 * H00100 * x
                    + 869.3333333333333 * H00101 * x
                    + 741.3333333333333 * H0011 * x
                    + 794.6666666666666 * H00110 * x + 816 * H00111 * x
                    + 656 * H00m10 * x + 522.6666666666666 * H00m100 * x
                    - 522.6666666666666 * H00m101 * x
                    - 1418.6666666666665 * H00m1m10 * x
                    - 2875.2375308641977 * H01 * x
                    - 262.0444444444445 * H010 * x - 1560 * H0100 * x
                    - 1040 * H01000 * x - 1877.3333333333333 * H01001 * x
                    + 1144 * H0101 * x + 789.3333333333333 * H01010 * x
                    + 693.3333333333333 * H01011 * x
                    + 42.666666666666664 * H010m10 * x
                    - 491.6740740740741 * H011 * x
                    + 1026.6666666666665 * H0110 * x
                    + 2138.6666666666665 * H01100 * x
                    + 741.3333333333333 * H01101 * x
                    + 971.5555555555555 * H0111 * x
                    + 789.3333333333333 * H01110 * x + 656 * H01111 * x
                    - 255.82222222222222 * H0m10 * x
                    + 2842.6666666666665 * H0m100 * x
                    + 1578.6666666666665 * H0m1000 * x + 736 * H0m1001 * x
                    + 3040 * H0m101 * x + 85.33333333333333 * H0m1010 * x
                    + 106.66666666666666 * H0m1011 * x
                    - 1749.3333333333333 * H0m10m10 * x
                    - 2218.6666666666665 * H0m1m10 * x
                    - 3050.6666666666665 * H0m1m100 * x
                    - 1301.3333333333333 * H0m1m101 * x
                    + 2133.333333333333 * H0m1m1m10 * x
                    - 3118.0245267489713 * H1 * x - 353.19506172839505 * H10 * x
                    + 907.7925925925927 * H100 * x
                    + 175.11111111111111 * H1000 * x - 624 * H10000 * x
                    - 602.6666666666666 * H10001 * x
                    + 981.3333333333333 * H1001 * x
                    + 53.33333333333333 * H10010 * x
                    - 197.33333333333331 * H10011 * x
                    + 586.6666666666666 * H100m10 * x
                    + 584.8888888888888 * H101 * x
                    + 1210.6666666666665 * H1010 * x
                    - 85.33333333333333 * H10100 * x
                    + 229.33333333333331 * H10101 * x + 880 * H1011 * x
                    + 256 * H10110 * x + 608 * H10m10 * x + 544 * H10m100 * x
                    + 853.3333333333333 * H10m101 * x
                    + 362.66666666666663 * H10m1m10 * x
                    - 1389.441975308642 * H11 * x + 424.4444444444444 * H110 * x
                    + 808 * H1100 * x - 570.6666666666666 * H11000 * x
                    + 768 * H11001 * x + 936 * H1101 * x + 368 * H11010 * x
                    + 64 * H11011 * x + 938.6666666666666 * H110m10 * x
                    + 93.92592592592592 * H111 * x
                    + 965.3333333333333 * H1110 * x
                    - 554.6666666666666 * H11100 * x + 208 * H11101 * x
                    + 608.8888888888888 * H1111 * x + 272 * H11110 * x
                    - 6122.9451851851845 * Hm10 * x - 3628 * Hm100 * x
                    + 1482.6666666666665 * Hm1000 * x + 1104 * Hm10000 * x
                    - 149.33333333333331 * Hm10001 * x
                    + 2493.333333333333 * Hm1001 * x + 416 * Hm10010 * x
                    + 480 * Hm10011 * x - 597.3333333333333 * Hm100m10 * x
                    - 2983.111111111111 * Hm101 * x
                    + 469.3333333333333 * Hm1010 * x
                    + 1909.3333333333333 * Hm10100 * x
                    + 181.33333333333331 * Hm10101 * x
                    + 549.3333333333333 * Hm1011 * x + 160 * Hm10110 * x
                    + 170.66666666666666 * Hm10111 * x
                    - 1045.3333333333333 * Hm10m10 * x
                    - 2117.333333333333 * Hm10m100 * x - 2368 * Hm10m101 * x
                    + 864 * Hm10m1m10 * x + 441.06666666666666 * Hm1m10 * x
                    - 3378.6666666666665 * Hm1m100 * x - 1936 * Hm1m1000 * x
                    - 2938.6666666666665 * Hm1m1001 * x - 3744 * Hm1m101 * x
                    - 501.3333333333333 * Hm1m1010 * x - 576 * Hm1m1011 * x
                    + 906.6666666666666 * Hm1m10m10 * x
                    + 1082.6666666666665 * Hm1m1m10 * x
                    + 3562.6666666666665 * Hm1m1m100 * x + 4000 * Hm1m1m101 * x
                    - 949.3333333333333 * Hm1m1m1m10 * x
                    + 4795.588148148148 * x2 - 1676.6222222222223 * H000 * x2
                    - 1712 * H0000 * x2 - 53.33333333333333 * H00001 * x2
                    - 1255.111111111111 * H0001 * x2 - 160 * H00010 * x2
                    - 181.33333333333331 * H00011 * x2
                    - 597.3333333333333 * H000m10 * x2
                    - 595.9555555555556 * H001 * x2
                    - 2135.111111111111 * H0010 * x2
                    - 394.66666666666663 * H00100 * x2
                    - 266.66666666666663 * H00101 * x2
                    - 2334.222222222222 * H0011 * x2
                    - 373.3333333333333 * H00110 * x2
                    - 245.33333333333331 * H00111 * x2
                    + 2074.6666666666665 * H00m10 * x2
                    + 309.3333333333333 * H00m100 * x2
                    + 810.6666666666666 * H00m101 * x2
                    - 533.3333333333333 * H00m1m10 * x2
                    + 516.5046913580247 * H01 * x2
                    - 1598.5037037037039 * H010 * x2
                    - 2389.333333333333 * H0100 * x2
                    - 90.66666666666666 * H01000 * x2 - 1680 * H01001 * x2
                    - 1882.6666666666665 * H0101 * x2
                    - 325.3333333333333 * H01010 * x2 - 176 * H01011 * x2
                    - 533.3333333333333 * H010m10 * x2
                    - 1190.4296296296297 * H011 * x2
                    - 1650.6666666666665 * H0110 * x2
                    + 1322.6666666666665 * H01100 * x2
                    - 309.3333333333333 * H01101 * x2
                    - 1644.4444444444443 * H0111 * x2
                    - 325.3333333333333 * H01110 * x2 - 160 * H01111 * x2
                    + 1129.5111111111112 * H0m10 * x2 + 3048 * H0m100 * x2
                    + 880 * H0m1000 * x2 + 1754.6666666666665 * H0m1001 * x2
                    + 1434.6666666666665 * H0m101 * x2
                    + 373.3333333333333 * H0m1010 * x2 + 448 * H0m1011 * x2
                    - 405.3333333333333 * H0m10m10 * x2
                    - 2229.333333333333 * H0m1m10 * x2
                    - 1722.6666666666665 * H0m1m100 * x2 - 2272 * H0m1m101 * x2
                    + 640 * H0m1m1m10 * x2 + 4143.4034567901235 * H1 * x2
                    + 447.61975308641973 * H10 * x2
                    - 761.0074074074074 * H100 * x2
                    + 164.44444444444443 * H1000 * x2 + 240 * H10000 * x2
                    + 346.66666666666663 * H10001 * x2
                    + 435.55555555555554 * H1001 * x2
                    - 117.33333333333333 * H10010 * x2
                    + 133.33333333333331 * H10011 * x2
                    - 458.66666666666663 * H100m10 * x2
                    - 726.2962962962963 * H101 * x2
                    - 1291.5555555555554 * H1010 * x2
                    - 170.66666666666666 * H10100 * x2
                    - 229.33333333333331 * H10101 * x2
                    - 1034.6666666666665 * H1011 * x2 - 256 * H10110 * x2
                    + 21.333333333333332 * H10m10 * x2 - 416 * H10m100 * x2
                    - 469.3333333333333 * H10m101 * x2
                    - 234.66666666666666 * H10m1m10 * x2
                    + 1405.6444444444444 * H11 * x2
                    - 789.8518518518518 * H110 * x2
                    - 1975.111111111111 * H1100 * x2
                    + 122.66666666666666 * H11000 * x2 - 448 * H11001 * x2
                    - 1062.2222222222222 * H1101 * x2 - 368 * H11010 * x2
                    - 64 * H11011 * x2 - 554.6666666666666 * H110m10 * x2
                    - 395.1111111111111 * H111 * x2
                    - 1070.2222222222222 * H1110 * x2
                    + 106.66666666666666 * H11100 * x2 - 208 * H11101 * x2
                    - 757.3333333333333 * H1111 * x2 - 272 * H11110 * x2
                    - 969.28 * Hm10 * x2 - 383.1111111111111 * Hm100 * x2
                    + 1384 * Hm1000 * x2 + 720 * Hm10000 * x2
                    - 277.3333333333333 * Hm10001 * x2
                    + 1069.3333333333333 * Hm1001 * x2 + 352 * Hm10010 * x2
                    + 416 * Hm10011 * x2 - 213.33333333333331 * Hm100m10 * x2
                    - 1521.3333333333333 * Hm101 * x2
                    + 245.33333333333331 * Hm1010 * x2
                    + 1525.3333333333333 * Hm10100 * x2
                    + 181.33333333333331 * Hm10101 * x2
                    + 277.3333333333333 * Hm1011 * x2 + 160 * Hm10110 * x2
                    + 170.66666666666666 * Hm10111 * x2 - 1520 * Hm10m10 * x2
                    - 1221.3333333333333 * Hm10m100 * x2 - 1472 * Hm10m101 * x2
                    + 480 * Hm10m1m10 * x2 - 1407.9111111111113 * Hm1m10 * x2
                    - 2909.333333333333 * Hm1m100 * x2 - 1104 * Hm1m1000 * x2
                    - 1850.6666666666665 * Hm1m1001 * x2
                    - 1690.6666666666665 * Hm1m101 * x2
                    - 373.3333333333333 * Hm1m1010 * x2 - 448 * Hm1m1011 * x2
                    + 522.6666666666666 * Hm1m10m10 * x2 + 1696 * Hm1m1m10 * x2
                    + 2090.6666666666665 * Hm1m1m100 * x2
                    + 2464 * Hm1m1m101 * x2
                    - 565.3333333333333 * Hm1m1m1m10 * x2 - 875.52 * H000 * x3
                    - 460.8 * H0000 * x3 + 38.400000000000006 * H0001 * x3
                    - 693.12 * H001 * x3 - 76.80000000000001 * H0010 * x3
                    - 76.80000000000001 * H0011 * x3
                    + 307.20000000000005 * H00m10 * x3 + 176 * H010 * x3
                    - 614.4000000000001 * H0100 * x3 + 517.12 * H0m10 * x3
                    + 576 * H0m100 * x3 + 768 * H0m101 * x3
                    - 76.80000000000001 * H0m1m10 * x3 - 115.2 * H1000 * x3
                    + 614.4000000000001 * H1001 * x3 - 176 * H101 * x3
                    + 230.4 * H10m10 * x3 + 176 * H110 * x3
                    - 614.4000000000001 * H1100 * x3
                    + 783.1466666666668 * Hm10 * x3 + 875.52 * Hm100 * x3
                    + 345.6 * Hm1000 * x3 + 576 * Hm1001 * x3
                    + 517.12 * Hm101 * x3 + 76.80000000000001 * Hm1010 * x3
                    + 76.80000000000001 * Hm1011 * x3
                    - 76.80000000000001 * Hm10m10 * x3 - 517.12 * Hm1m10 * x3
                    - 576 * Hm1m100 * x3 - 768 * Hm1m101 * x3
                    + 76.80000000000001 * Hm1m1m10 * x3
                    + (9.955555555555556 * H0m10 + 4.266666666666667 * H0m100
                       + 4.266666666666667 * H0m101
                       - 4.266666666666667 * H0m1m10 + 3.2 * H1000
                       - 17.066666666666666 * H1001 - 6.4 * H10m10
                       + 17.066666666666666 * H1100 + 20.628148148148142 * Hm10
                       + 21.297777777777778 * Hm100 + 9.600000000000001 * Hm1000
                       + 16. * Hm1001 + 11.342222222222222 * Hm101
                       + 2.1333333333333333 * Hm1010
                       + 2.1333333333333333 * Hm1011
                       - 2.1333333333333333 * Hm10m10
                       - 11.342222222222222 * Hm1m10 - 16. * Hm1m100
                       - 21.333333333333332 * Hm1m101
                       + 2.1333333333333333 * Hm1m1m10 + 59.611851851851846 * x
                       - 2298.311111111111 * H0m10 * x2)
                          / x2
                    + 428.3930864197531 * zeta2
                    - 218.66666666666666 * H000 * zeta2 - 64 * H001 * zeta2
                    - 221.33333333333331 * H00m1 * zeta2
                    + 185.33333333333331 * H01 * zeta2
                    + 629.3333333333333 * H010 * zeta2 - 40 * H011 * zeta2
                    - 258.66666666666663 * H0m1 * zeta2
                    - 1589.3333333333333 * H0m10 * zeta2 + 1792 * H0m1m1 * zeta2
                    + 804.9777777777778 * H1 * zeta2 + 1948 * H10 * zeta2
                    - 597.3333333333333 * H100 * zeta2
                    - 101.33333333333333 * H101 * zeta2 + 336 * H10m1 * zeta2
                    + 356 * H11 * zeta2 - 77.33333333333333 * H110 * zeta2
                    - 133.33333333333331 * H111 * zeta2
                    + 1863.4222222222222 * Hm1 * zeta2 - 2672 * Hm10 * zeta2
                    - 605.3333333333333 * Hm100 * zeta2 + 1400 * Hm10m1 * zeta2
                    + 2918.6666666666665 * Hm1m1 * zeta2
                    + 1930.6666666666665 * Hm1m10 * zeta2
                    - 2237.333333333333 * Hm1m1m1 * zeta2
                    - (2.1333333333333333 * H01 * zeta2) / x2
                    - (6.4 * H0m1 * zeta2) / x2
                    - (5.671111111111111 * H1 * zeta2) / x2
                    + (12.8 * H10 * zeta2) / x2
                    - (1.0666666666666667 * H11 * zeta2) / x2
                    - (17.013333333333335 * Hm1 * zeta2) / x2
                    - (20.266666666666666 * Hm10 * zeta2) / x2
                    + (22.400000000000002 * Hm1m1 * zeta2) / x2
                    + (7.751111111111111 * zeta2) / x
                    - (16.11851851851852 * H1 * zeta2) / x
                    + (14.222222222222221 * H10 * zeta2) / x
                    - (35.55555555555556 * H11 * zeta2) / x
                    + (71.82222222222222 * Hm1 * zeta2) / x
                    - (85.33333333333333 * Hm10 * zeta2) / x
                    + (128 * Hm1m1 * zeta2) / x - 3247.707654320988 * x * zeta2
                    - 650.6666666666666 * H000 * x * zeta2
                    - 160 * H001 * x * zeta2
                    - 186.66666666666666 * H00m1 * x * zeta2
                    - 34.666666666666664 * H01 * x * zeta2
                    + 2773.333333333333 * H010 * x * zeta2
                    + 325.3333333333333 * H011 * x * zeta2
                    - 4149.333333333333 * H0m1 * x * zeta2
                    - 1632 * H0m10 * x * zeta2 + 2368 * H0m1m1 * x * zeta2
                    - 805.4222222222222 * H1 * x * zeta2
                    - 154.66666666666666 * H10 * x * zeta2
                    + 1194.6666666666665 * H100 * x * zeta2
                    + 202.66666666666666 * H101 * x * zeta2
                    - 672 * H10m1 * x * zeta2
                    - 394.66666666666663 * H11 * x * zeta2
                    + 154.66666666666666 * H110 * x * zeta2
                    + 266.66666666666663 * H111 * x * zeta2
                    + 3203.6444444444446 * Hm1 * x * zeta2
                    - 3320 * Hm10 * x * zeta2
                    - 442.66666666666663 * Hm100 * x * zeta2
                    + 2800 * Hm10m1 * x * zeta2
                    + 4285.333333333333 * Hm1m1 * x * zeta2
                    + 3861.333333333333 * Hm1m10 * x * zeta2
                    - 4474.666666666666 * Hm1m1m1 * x * zeta2
                    - 516.5046913580247 * x2 * zeta2
                    + 53.33333333333333 * H000 * x2 * zeta2
                    - 1077.3333333333333 * H00m1 * x2 * zeta2
                    + 768 * H01 * x2 * zeta2
                    + 1210.6666666666665 * H010 * x2 * zeta2
                    - 10.666666666666666 * H011 * x2 * zeta2
                    - 2549.333333333333 * H0m1 * x2 * zeta2
                    - 2224 * H0m10 * x2 * zeta2 + 2592 * H0m1m1 * x2 * zeta2
                    + 22.340740740740742 * H1 * x2 * zeta2
                    - 1184.888888888889 * H10 * x2 * zeta2
                    - 682.6666666666666 * H100 * x2 * zeta2
                    - 10.666666666666666 * H101 * x2 * zeta2
                    + 352 * H10m1 * x2 * zeta2
                    + 214.2222222222222 * H11 * x2 * zeta2
                    - 90.66666666666666 * H110 * x2 * zeta2
                    - 74.66666666666666 * H111 * x2 * zeta2
                    + 817.3777777777779 * Hm1 * x2 * zeta2
                    - 1818.6666666666665 * Hm10 * x2 * zeta2
                    - 58.666666666666664 * Hm100 * x2 * zeta2
                    - 64 * Hm101 * x2 * zeta2 + 1712 * Hm10m1 * x2 * zeta2
                    + 2538.6666666666665 * Hm1m1 * x2 * zeta2
                    + 2389.333333333333 * Hm1m10 * x2 * zeta2
                    - 2746.6666666666665 * Hm1m1m1 * x2 * zeta2
                    + 783.1466666666668 * x3 * zeta2
                    + 38.400000000000006 * H01 * x3 * zeta2
                    - 806.4000000000001 * H0m1 * x3 * zeta2
                    + 434.56 * H1 * x3 * zeta2 - 460.8 * H10 * x3 * zeta2
                    + 38.400000000000006 * H11 * x3 * zeta2
                    - 775.6800000000001 * Hm1 * x3 * zeta2
                    - 729.6 * Hm10 * x3 * zeta2
                    + 806.4000000000001 * Hm1m1 * x3 * zeta2
                    - 221.95555555555558 * zeta2_2
                    - 150.93333333333334 * H1 * zeta2_2
                    - 191.46666666666667 * Hm1 * zeta2_2 + (25.6 * zeta2_2) / x
                    + 1112.9777777777779 * x * zeta2_2
                    + 301.8666666666667 * H1 * x * zeta2_2
                    - 997.3333333333333 * Hm1 * x * zeta2_2
                    + 536.1777777777778 * x2 * zeta2_2
                    - 113.06666666666666 * H1 * x2 * zeta2_2
                    - 706.1333333333333 * Hm1 * x2 * zeta2_2
                    + 480 * x3 * zeta2_2 + 94.16296296296296 * zeta3
                    - 1329.3333333333333 * H01 * zeta3
                    - 1534.6666666666665 * H0m1 * zeta3
                    - 2221.5555555555557 * H1 * zeta3 - 1328 * H10 * zeta3
                    - 1514.6666666666665 * H11 * zeta3
                    - 2416.6666666666665 * Hm1 * zeta3
                    - 1645.3333333333333 * Hm10 * zeta3
                    + 1909.3333333333333 * Hm1m1 * zeta3
                    - (25.066666666666666 * H1 * zeta3) / x2
                    - (18.666666666666664 * Hm1 * zeta3) / x2
                    - (45.15555555555556 * zeta3) / x
                    + (49.77777777777777 * H1 * zeta3) / x
                    - (106.66666666666666 * Hm1 * zeta3) / x
                    - 4454.014814814815 * x * zeta3 + 696 * H01 * x * zeta3
                    - 2088 * H0m1 * x * zeta3
                    + 1601.7777777777776 * H1 * x * zeta3
                    + 2656 * H10 * x * zeta3
                    + 3029.333333333333 * H11 * x * zeta3
                    - 3613.333333333333 * Hm1 * x * zeta3
                    - 4058.6666666666665 * Hm10 * x * zeta3
                    + 3818.6666666666665 * Hm1m1 * x * zeta3
                    + 1579.8222222222223 * x2 * zeta3 - 2672 * H01 * x2 * zeta3
                    - 2272 * H0m1 * x2 * zeta3
                    - 767.5555555555555 * H1 * x2 * zeta3
                    - 1568 * H10 * x2 * zeta3
                    - 1781.3333333333333 * H11 * x2 * zeta3
                    - 2222.6666666666665 * Hm1 * x2 * zeta3
                    - 2842.6666666666665 * Hm10 * x2 * zeta3
                    + 2378.6666666666665 * Hm1m1 * x2 * zeta3
                    + 764.8000000000001 * x3 * zeta3
                    + 902.4000000000001 * H1 * x3 * zeta3
                    - 672 * Hm1 * x3 * zeta3 - 176 * zeta2 * zeta3
                    - 389.3333333333333 * x * zeta2 * zeta3
                    - 533.3333333333333 * x2 * zeta2 * zeta3
                    + H00
                          * (-165.14123456790125 + 94.66666666666666 * zeta2
                             - 106.66666666666666 * zeta3
                             + x
                                   * (3360.9624691358026 + 256. * zeta2
                                      + 677.3333333333333 * zeta3
                                      + x
                                            * (457.68
                                               + 1255.111111111111 * zeta2
                                               + x
                                                     * (-783.1466666666668
                                                        + 268.8 * zeta2)
                                               + 74.66666666666666 * zeta3)))
                    + 191.99999999999997 * zeta4 - 24 * H1 * zeta4
                    + (32 * zeta4) / x - 48 * x * zeta4 + 48 * H1 * x * zeta4
                    - 248 * x2 * zeta4 - 48 * H1 * x2 * zeta4
                    + (H0
                       * (-13.552592592592593 + 7.466666666666667 * zeta2
                          + x
                                * (-380.9769547325103
                                   + (-29.644444444444446
                                      - 110.93333333333334 * zeta2)
                                         * zeta2
                                   + 146.2222222222222 * zeta3
                                   + x3 * (1210.24 * zeta2 + 1574.4 * zeta3)
                                   + 96. * zeta4
                                   + x2
                                         * (2736.8434567901236
                                            + (595.9555555555556
                                               - 228.26666666666668 * zeta2)
                                                  * zeta2
                                            + 2376.8888888888887 * zeta3
                                            + 48. * zeta4)
                                   + x
                                         * (-2610.383374485597
                                            + zeta2
                                                  * (-1441.0666666666666
                                                     + 61.86666666666667 * zeta2
                                                  )
                                            + 6987.555555555555 * zeta3
                                            + 240. * zeta4))))
                          / x
                    + 308.66666666666663 * zeta5 + 8201.333333333332 * x * zeta5
                    + 237.33333333333331 * x2 * zeta5)
           + CF * CF * nf
                 * (248.37166666666667 - 59.33333333333333 * H000
                    + 14.666666666666666 * H0000 - 60 * H00000
                    - 174.66666666666666 * H00001 - 94.66666666666666 * H0001
                    - 202.66666666666666 * H00010 - 269.3333333333333 * H00011
                    + 533.3333333333333 * H000m10 - 100.73333333333333 * H001
                    - 177.33333333333331 * H0010 - 74.66666666666666 * H00100
                    - 250.66666666666666 * H00101 - 248 * H0011
                    - 234.66666666666666 * H00110 - 290.66666666666663 * H00111
                    + 1162.6666666666665 * H00m10 + 629.3333333333333 * H00m100
                    + 224 * H00m101 - 394.66666666666663 * H00m1m10
                    + 102.17777777777778 * H01 - 154.66666666666666 * H010
                    + 46.666666666666664 * H0100 - 114.66666666666666 * H01000
                    + 440 * H01001 - 302.66666666666663 * H0101
                    - 234.66666666666666 * H01010 - 322.66666666666663 * H01011
                    + 384 * H010m10 - 209.33333333333331 * H011
                    - 321.3333333333333 * H0110 - 797.3333333333333 * H01100
                    - 269.3333333333333 * H01101 - 364 * H0111
                    - 269.3333333333333 * H01110 - 301.3333333333333 * H01111
                    + 1530.6666666666665 * H0m100 + 197.33333333333331 * H0m1000
                    - 96.00000000000001 * H0m1001 + 624. * H0m101
                    + 32. * H0m1010 + 32. * H0m1011 - 256. * H0m10m10
                    - 1072. * H0m1m10 - 224. * H0m1m100 + 352. * H0m1m101
                    + 288. * H0m1m1m10 + 731.2444444444445 * H1
                    - 79.66666666666666 * H10 - 71.06666666666666 * H100
                    + 30.666666666666668 * H1000 - 277.3333333333333 * H10000
                    - 336. * H10001 + 1262.6666666666665 * H1001
                    - 389.3333333333333 * H10010 - 496.00000000000006 * H10011
                    + 469.3333333333333 * H100m10 - 471.33333333333326 * H101
                    - 490.6666666666667 * H1010 - 437.33333333333337 * H10100
                    - 432. * H10101 - 550.6666666666666 * H1011
                    - 394.66666666666663 * H10110 - 469.3333333333333 * H10111
                    + 821.3333333333333 * H10m10 + 490.6666666666667 * H10m100
                    + 384.00000000000006 * H10m101 - 128. * H10m1m10
                    - 18.333333333333332 * H11 - 482. * H110 - 1604. * H1100
                    - 218.66666666666669 * H11000 + 21.333333333333332 * H11001
                    - 478.6666666666666 * H1101 - 341.3333333333333 * H11010
                    - 442.66666666666663 * H11011 + 341.3333333333333 * H110m10
                    - 522. * H111 - 505.3333333333333 * H1110 - 656. * H11100
                    - 378.66666666666663 * H11101 - 529.3333333333333 * H1111
                    - 357.3333333333333 * H11110 - 400. * H11111
                    + 4989.688888888889 * Hm10 + 3834.1333333333328 * Hm100
                    + 640. * Hm1000 - 37.33333333333333 * Hm10000
                    + 384.00000000000006 * Hm10001 - 96.00000000000001 * Hm1001
                    + 32. * Hm10010 + 32. * Hm10011
                    - 213.33333333333331 * Hm100m10 + 2092.2666666666664 * Hm101
                    + 64. * Hm1010 - 458.6666666666666 * Hm10100 + 64. * Hm1011
                    - 847.9999999999999 * Hm10m10 - 288. * Hm10m100
                    + 128. * Hm10m101 + 234.66666666666666 * Hm10m1m10
                    - 2429.333333333333 * Hm1m10 - 984. * Hm1m100
                    - 320. * Hm1m1000 + 32. * Hm1m1001 + 512. * Hm1m101
                    - 64. * Hm1m1010 - 64. * Hm1m1011
                    + 213.33333333333331 * Hm1m10m10 + 1040. * Hm1m1m10
                    + 192.00000000000003 * Hm1m1m100
                    - 384.00000000000006 * Hm1m1m101
                    - 213.33333333333331 * Hm1m1m1m10
                    + (3.9111111111111114 * H00) / x
                    + (4.266666666666667 * H000) / x
                    - (9.600000000000001 * H001) / x
                    + (13.511111111111111 * H01) / x
                    - (4.266666666666667 * H0m10) / x
                    - (13.511111111111111 * H1) / x
                    + (13.866666666666667 * H100) / x
                    + (0.35555555555555557 * Hm10) / x
                    - (4.266666666666667 * Hm100) / x
                    - (8.533333333333333 * Hm101) / x + 1331.3366666666668 * x
                    - 1655.0222222222224 * H000 * x
                    + 81.33333333333333 * H0000 * x
                    - 61.33333333333333 * H00000 * x
                    + 221.33333333333331 * H00001 * x + 156 * H0001 * x
                    + 341.3333333333333 * H00010 * x
                    + 474.66666666666663 * H00011 * x
                    - 960.7555555555556 * H001 * x
                    + 626.6666666666666 * H0010 * x
                    + 490.66666666666663 * H00100 * x
                    + 501.3333333333333 * H00101 * x
                    + 845.3333333333333 * H0011 * x
                    + 469.3333333333333 * H00110 * x
                    + 581.3333333333333 * H00111 * x - 384 * H00m10 * x
                    - 128 * H00m100 * x + 320 * H00m101 * x + 576 * H00m1m10 * x
                    + 2542.733333333333 * H01 * x + 1592 * H010 * x
                    + 2184 * H0100 * x + 826.6666666666666 * H01000 * x
                    + 2341.333333333333 * H01001 * x + 1064 * H0101 * x
                    + 448 * H01010 * x + 645.3333333333333 * H01011 * x
                    - 128 * H010m10 * x + 1757.3333333333333 * H011 * x
                    + 1032 * H0110 * x - 1050.6666666666665 * H01100 * x
                    + 560 * H01101 * x + 1162.6666666666665 * H0111 * x
                    + 517.3333333333333 * H01110 * x
                    + 602.6666666666666 * H01111 * x
                    + 661.6888888888889 * H0m10 * x
                    - 1013.3333333333333 * H0m100 * x
                    - 202.66666666666666 * H0m1000 * x + 64 * H0m1001 * x
                    - 1216 * H0m101 * x + 64 * H0m1010 * x + 64 * H0m1011 * x
                    + 426.66666666666663 * H0m10m10 * x + 1024 * H0m1m10 * x
                    + 661.3333333333333 * H0m1m100 * x + 192 * H0m1m101 * x
                    - 448 * H0m1m1m10 * x + 2181.2000000000003 * H1 * x
                    + 1216.6666666666665 * H10 * x
                    + 621.0666666666666 * H100 * x
                    + 650.6666666666666 * H1000 * x
                    + 554.6666666666666 * H10000 * x + 672 * H10001 * x
                    + 1461.3333333333333 * H1001 * x
                    + 778.6666666666666 * H10010 * x + 992 * H10011 * x
                    - 938.6666666666666 * H100m10 * x + 1476 * H101 * x
                    + 1418.6666666666665 * H1010 * x
                    + 874.6666666666666 * H10100 * x + 864 * H10101 * x
                    + 1642.6666666666665 * H1011 * x
                    + 789.3333333333333 * H10110 * x
                    + 938.6666666666666 * H10111 * x
                    - 853.3333333333333 * H10m10 * x
                    - 981.3333333333333 * H10m100 * x - 768 * H10m101 * x
                    + 256 * H10m1m10 * x + 1399.3333333333333 * H11 * x
                    + 1521.3333333333333 * H110 * x
                    + 597.3333333333333 * H1100 * x
                    + 437.3333333333333 * H11000 * x
                    - 42.666666666666664 * H11001 * x + 1440 * H1101 * x
                    + 682.6666666666666 * H11010 * x
                    + 885.3333333333333 * H11011 * x
                    - 682.6666666666666 * H110m10 * x
                    + 1649.3333333333333 * H111 * x + 1488 * H1110 * x
                    + 1312 * H11100 * x + 757.3333333333333 * H11101 * x
                    + 1562.6666666666665 * H1111 * x
                    + 714.6666666666666 * H11110 * x + 800 * H11111 * x
                    + 6253.688888888889 * Hm10 * x
                    + 4819.022222222223 * Hm100 * x
                    + 1029.3333333333333 * Hm1000 * x
                    - 74.66666666666666 * Hm10000 * x + 1536 * Hm10001 * x
                    + 149.33333333333331 * Hm1001 * x + 64 * Hm10010 * x
                    + 64 * Hm10011 * x - 426.66666666666663 * Hm100m10 * x
                    + 4018.488888888889 * Hm101 * x + 160 * Hm1010 * x
                    - 1685.3333333333333 * Hm10100 * x + 160 * Hm1011 * x
                    - 885.3333333333333 * Hm10m10 * x - 576 * Hm10m100 * x
                    + 256 * Hm10m101 * x + 469.3333333333333 * Hm10m1m10 * x
                    - 1768.8888888888887 * Hm1m10 * x - 1072 * Hm1m100 * x
                    - 640 * Hm1m1000 * x + 64 * Hm1m1001 * x
                    + 501.3333333333333 * Hm1m101 * x - 128 * Hm1m1010 * x
                    - 128 * Hm1m1011 * x + 426.66666666666663 * Hm1m10m10 * x
                    + 1216 * Hm1m1m10 * x + 384 * Hm1m1m100 * x
                    - 768 * Hm1m1m101 * x - 426.66666666666663 * Hm1m1m1m10 * x
                    - 1541.1466666666668 * x2 - 1128 * H000 * x2
                    - 1226.6666666666665 * H0000 * x2 - 480 * H00000 * x2
                    - 874.6666666666666 * H00001 * x2
                    - 2778.6666666666665 * H0001 * x2
                    - 842.6666666666666 * H00010 * x2
                    - 1045.3333333333333 * H00011 * x2
                    + 1365.3333333333333 * H000m10 * x2
                    - 2703.2000000000003 * H001 * x2 - 1280 * H0010 * x2
                    - 522.6666666666666 * H00100 * x2 - 864 * H00101 * x2
                    - 1504 * H0011 * x2 - 800 * H00110 * x2
                    - 938.6666666666666 * H00111 * x2
                    + 186.66666666666666 * H00m10 * x2
                    + 1557.3333333333333 * H00m100 * x2 + 512 * H00m101 * x2
                    - 725.3333333333333 * H00m1m10 * x2
                    - 962.9333333333333 * H01 * x2
                    - 1089.3333333333333 * H010 * x2
                    + 253.33333333333331 * H0100 * x2
                    - 309.3333333333333 * H01000 * x2 + 800 * H01001 * x2
                    - 1261.3333333333333 * H0101 * x2
                    - 682.6666666666666 * H01010 * x2
                    - 874.6666666666666 * H01011 * x2
                    + 426.66666666666663 * H010m10 * x2 - 1212 * H011 * x2
                    - 1266.6666666666665 * H0110 * x2
                    - 2069.333333333333 * H01100 * x2
                    - 757.3333333333333 * H01101 * x2
                    - 1365.3333333333333 * H0111 * x2
                    - 725.3333333333333 * H01110 * x2 - 800 * H01111 * x2
                    - 574.4 * H0m10 * x2 - 93.33333333333333 * H0m100 * x2
                    + 768 * H0m1000 * x2 + 192 * H0m1001 * x2 - 96 * H0m101 * x2
                    + 128 * H0m1010 * x2 + 128 * H0m1011 * x2
                    - 426.66666666666663 * H0m10m10 * x2 + 80 * H0m1m10 * x2
                    - 640 * H0m1m100 * x2 + 256 * H0m1m101 * x2
                    + 426.66666666666663 * H0m1m1m10 * x2
                    - 2983.6000000000004 * H1 * x2
                    - 1186.6666666666665 * H10 * x2
                    - 520.5333333333333 * H100 * x2 - 1024 * H1000 * x2
                    - 426.66666666666663 * H10000 * x2 - 736 * H10001 * x2
                    - 2458.6666666666665 * H1001 * x2
                    - 778.6666666666666 * H10010 * x2 - 992 * H10011 * x2
                    + 810.6666666666666 * H100m10 * x2 - 1012 * H101 * x2
                    - 1184 * H1010 * x2 - 682.6666666666666 * H10100 * x2
                    - 864 * H10101 * x2 - 1408 * H1011 * x2
                    - 789.3333333333333 * H10110 * x2
                    - 938.6666666666666 * H10111 * x2
                    + 405.3333333333333 * H10m10 * x2
                    + 853.3333333333333 * H10m100 * x2 + 512 * H10m101 * x2
                    - 256 * H10m1m10 * x2 - 1449.3333333333333 * H11 * x2
                    - 1089.3333333333333 * H110 * x2
                    + 253.33333333333331 * H1100 * x2
                    - 309.3333333333333 * H11000 * x2
                    - 341.3333333333333 * H11001 * x2
                    - 1261.3333333333333 * H1101 * x2
                    - 682.6666666666666 * H11010 * x2
                    - 885.3333333333333 * H11011 * x2
                    + 426.66666666666663 * H110m10 * x2 - 1212 * H111 * x2
                    - 1266.6666666666665 * H1110 * x2 - 928 * H11100 * x2
                    - 757.3333333333333 * H11101 * x2
                    - 1365.3333333333333 * H1111 * x2
                    - 714.6666666666666 * H11110 * x2 - 800 * H11111 * x2
                    + 1135.4666666666667 * Hm10 * x2
                    + 953.0666666666666 * Hm100 * x2
                    + 202.66666666666666 * Hm1000 * x2
                    + 53.33333333333333 * Hm10000 * x2 + 1280 * Hm10001 * x2
                    + 64 * Hm1001 * x2 + 64 * Hm10010 * x2 + 64 * Hm10011 * x2
                    - 554.6666666666666 * Hm100m10 * x2
                    + 2036.8000000000002 * Hm101 * x2 + 96 * Hm1010 * x2
                    - 1301.3333333333333 * Hm10100 * x2 + 96 * Hm1011 * x2
                    - 37.33333333333333 * Hm10m10 * x2 - 704 * Hm10m100 * x2
                    + 469.3333333333333 * Hm10m1m10 * x2 + 728 * Hm1m10 * x2
                    + 93.33333333333333 * Hm1m100 * x2 - 768 * Hm1m1000 * x2
                    - 192 * Hm1m1001 * x2 + 352 * Hm1m101 * x2
                    - 128 * Hm1m1010 * x2 - 128 * Hm1m1011 * x2
                    + 426.66666666666663 * Hm1m10m10 * x2 + 176 * Hm1m1m10 * x2
                    + 640 * Hm1m1m100 * x2 - 256 * Hm1m1m101 * x2
                    - 426.66666666666663 * Hm1m1m1m10 * x2
                    + 294.40000000000003 * H000 * x3
                    + 153.60000000000002 * H0000 * x3 - 345.6 * H0001 * x3
                    + 140.8 * H001 * x3 - 153.60000000000002 * H00m10 * x3
                    + 499.20000000000005 * H0100 * x3 - 140.8 * H0m10 * x3
                    - 153.60000000000002 * H0m100 * x3
                    - 307.20000000000005 * H0m101 * x3
                    + 76.80000000000001 * H1000 * x3
                    - 499.20000000000005 * H1001 * x3
                    - 153.60000000000002 * H10m10 * x3
                    + 499.20000000000005 * H1100 * x3 - 341.12 * Hm10 * x3
                    - 294.40000000000003 * Hm100 * x3
                    - 76.80000000000001 * Hm1000 * x3
                    - 153.60000000000002 * Hm1001 * x3 - 140.8 * Hm101 * x3
                    + 140.8 * Hm1m10 * x3 + 153.60000000000002 * Hm1m100 * x3
                    + 307.20000000000005 * Hm1m101 * x3
                    + (-2.1333333333333333 * H1000 + 13.866666666666667 * H1001
                       + 4.266666666666667 * H10m10 - 13.866666666666667 * H1100
                       - 8.231111111111112 * Hm10 - 8.177777777777779 * Hm100
                       - 2.1333333333333333 * Hm1000
                       - 4.266666666666667 * Hm1001 - 3.9111111111111114 * Hm101
                       + 3.9111111111111114 * Hm1m10
                       + 4.266666666666667 * Hm1m100
                       + 8.533333333333333 * Hm1m101 - 12.853333333333333 * x
                       + H0m10 * (-4.266666666666667 + 2605.866666666667 * x2))
                          / x2
                    - 102.17777777777776 * zeta2
                    + 174.66666666666666 * H000 * zeta2
                    + 53.33333333333333 * H001 * zeta2
                    - 421.3333333333333 * H00m1 * zeta2
                    - 233.33333333333331 * H01 * zeta2 - 376 * H010 * zeta2
                    + 125.33333333333333 * H011 * zeta2 - 1160 * H0m1 * zeta2
                    + 160 * H0m10 * zeta2 - 208 * H0m1m1 * zeta2
                    - 743.3333333333333 * H1 * zeta2 - 1276 * H10 * zeta2
                    + 464 * H100 * zeta2 + 314.66666666666663 * H101 * zeta2
                    - 448 * H10m1 * zeta2 - 41.33333333333333 * H11 * zeta2
                    + 42.666666666666664 * H110 * zeta2 + 272 * H111 * zeta2
                    - 3306.9333333333334 * Hm1 * zeta2
                    + 82.66666666666666 * Hm10 * zeta2 - 256 * Hm100 * zeta2
                    + 64 * Hm101 * zeta2 - 10.666666666666666 * Hm10m1 * zeta2
                    + 8 * Hm1m1 * zeta2 - 96 * Hm1m10 * zeta2
                    + 277.3333333333333 * Hm1m1m1 * zeta2
                    + (1.9555555555555557 * H1 * zeta2) / x2
                    - (11.733333333333333 * H10 * zeta2) / x2
                    + (5.866666666666666 * Hm1 * zeta2) / x2
                    + (6.4 * Hm10 * zeta2) / x2
                    - (8.533333333333333 * Hm1m1 * zeta2) / x2
                    - (13.155555555555557 * zeta2) / x
                    + (8.533333333333333 * Hm1 * zeta2) / x
                    + 3710.9555555555557 * x * zeta2
                    - 221.33333333333331 * H000 * x * zeta2
                    - 789.3333333333333 * H001 * x * zeta2
                    - 32 * H00m1 * x * zeta2 - 1576 * H01 * x * zeta2
                    - 2618.6666666666665 * H010 * x * zeta2
                    - 784 * H011 * x * zeta2 + 1728 * H0m1 * x * zeta2
                    + 213.33333333333331 * H0m10 * x * zeta2
                    - 416 * H0m1m1 * x * zeta2
                    - 591.5555555555555 * H1 * x * zeta2
                    - 1445.3333333333333 * H10 * x * zeta2
                    - 928 * H100 * x * zeta2
                    - 629.3333333333333 * H101 * x * zeta2
                    + 896 * H10m1 * x * zeta2 - 832 * H11 * x * zeta2
                    - 85.33333333333333 * H110 * x * zeta2
                    - 544 * H111 * x * zeta2
                    - 4902.933333333333 * Hm1 * x * zeta2
                    - 165.33333333333331 * Hm10 * x * zeta2
                    - 1280 * Hm100 * x * zeta2 + 128 * Hm101 * x * zeta2
                    - 21.333333333333332 * Hm10m1 * x * zeta2
                    + 106.66666666666666 * Hm1m1 * x * zeta2
                    - 192 * Hm1m10 * x * zeta2
                    + 554.6666666666666 * Hm1m1m1 * x * zeta2
                    + 962.9333333333333 * x2 * zeta2
                    + 874.6666666666666 * H000 * x2 * zeta2
                    + 501.3333333333333 * H001 * x2 * zeta2
                    - 874.6666666666666 * H00m1 * x2 * zeta2
                    + 1301.3333333333333 * H01 * x2 * zeta2
                    - 800 * H010 * x2 * zeta2 + 544 * H011 * x2 * zeta2
                    + 136 * H0m1 * x2 * zeta2 - 192 * H0m10 * x2 * zeta2
                    - 42.666666666666664 * H0m1m1 * x2 * zeta2
                    + 1376 * H1 * x2 * zeta2
                    + 2642.6666666666665 * H10 * x2 * zeta2
                    + 864 * H100 * x2 * zeta2
                    + 629.3333333333333 * H101 * x2 * zeta2
                    - 640 * H10m1 * x2 * zeta2
                    + 1173.3333333333333 * H11 * x2 * zeta2
                    + 341.3333333333333 * H110 * x2 * zeta2
                    + 544 * H111 * x2 * zeta2
                    - 1672.8000000000002 * Hm1 * x2 * zeta2
                    + 120 * Hm10 * x2 * zeta2 - 1152 * Hm100 * x2 * zeta2
                    + 128 * Hm101 * x2 * zeta2
                    + 234.66666666666666 * Hm10m1 * x2 * zeta2
                    - 264 * Hm1m1 * x2 * zeta2 + 192 * Hm1m10 * x2 * zeta2
                    + 42.666666666666664 * Hm1m1m1 * x2 * zeta2
                    - 341.12 * x3 * zeta2
                    + 307.20000000000005 * H0m1 * x3 * zeta2
                    - 70.4 * H1 * x3 * zeta2
                    + 422.40000000000003 * H10 * x3 * zeta2
                    + 211.20000000000002 * Hm1 * x3 * zeta2
                    + 230.4 * Hm10 * x3 * zeta2
                    - 307.20000000000005 * Hm1m1 * x3 * zeta2
                    + 115.73333333333333 * zeta2_2
                    + 43.733333333333334 * H1 * zeta2_2
                    + 311.46666666666664 * Hm1 * zeta2_2
                    + 122.93333333333334 * x * zeta2_2
                    - 87.46666666666667 * H1 * x * zeta2_2
                    + 1237.3333333333333 * Hm1 * x * zeta2_2
                    - 1345.0666666666666 * x2 * zeta2_2
                    - 155.73333333333332 * H1 * x2 * zeta2_2
                    + 840.5333333333333 * Hm1 * x2 * zeta2_2
                    - 506.88 * x3 * zeta2_2 + 740.4666666666667 * zeta3
                    + 1050.6666666666665 * H01 * zeta3 + 96 * H0m1 * zeta3
                    + 1637.3333333333333 * H1 * zeta3
                    + 965.3333333333333 * H10 * zeta3
                    + 906.6666666666666 * H11 * zeta3 - 184 * Hm1 * zeta3
                    + 634.6666666666666 * Hm10 * zeta3
                    - 133.33333333333331 * Hm1m1 * zeta3
                    + (18.133333333333333 * H1 * zeta3) / x2
                    + (6.4 * Hm1 * zeta3) / x2 - (24.53333333333333 * zeta3) / x
                    + 5948.488888888889 * x * zeta3
                    - 53.33333333333333 * H01 * x * zeta3
                    + 320 * H0m1 * x * zeta3
                    - 1029.3333333333333 * H1 * x * zeta3
                    - 1930.6666666666665 * H10 * x * zeta3
                    - 1813.3333333333333 * H11 * x * zeta3
                    - 352 * Hm1 * x * zeta3
                    + 2037.3333333333333 * Hm10 * x * zeta3
                    - 266.66666666666663 * Hm1m1 * x * zeta3
                    + 263.73333333333335 * x2 * zeta3
                    + 2325.333333333333 * H01 * x2 * zeta3
                    - 117.33333333333333 * H0m1 * x2 * zeta3
                    + 794.6666666666666 * H1 * x2 * zeta3
                    + 1418.6666666666665 * H10 * x2 * zeta3
                    + 1173.3333333333333 * H11 * x2 * zeta3
                    + 104 * Hm1 * x2 * zeta3
                    + 1333.3333333333333 * Hm10 * x2 * zeta3
                    + 117.33333333333333 * Hm1m1 * x2 * zeta3 - 352 * x3 * zeta3
                    - 652.8000000000001 * H1 * x3 * zeta3
                    + 230.4 * Hm1 * x3 * zeta3 - 160 * zeta2 * zeta3
                    + 85.33333333333333 * x * zeta2 * zeta3
                    - 832 * x2 * zeta2 * zeta3
                    + H00
                          * (-96.88888888888889 + 94.66666666666666 * zeta2
                             + 284. * zeta3
                             + x
                                   * (-4094.2222222222217 - 540. * zeta2
                                      - 760. * zeta3
                                      + x
                                            * (-1669.3333333333335
                                               + 2778.6666666666665 * zeta2
                                               + x * (341.12 + 192. * zeta2)
                                               + 1237.3333333333333 * zeta3)))
                    + (H0
                       * (8.586666666666668 + 5.333333333333333 * zeta2
                          + x
                                * (60.41777777777778
                                   + zeta2
                                         * (100.73333333333333
                                            + 7.733333333333333 * zeta2)
                                   + x3 * (-281.6 * zeta2 - 883.2 * zeta3)
                                   + 300. * zeta3
                                   + x2
                                         * (-2476.08
                                            + (2703.2000000000003
                                               - 89.60000000000001 * zeta2)
                                                  * zeta2
                                            + 818.6666666666666 * zeta3)
                                   + x
                                         * (502.1066666666667
                                            + zeta2
                                                  * (1622.4444444444443
                                                     + 460.26666666666665
                                                           * zeta2)
                                            - 4094.6666666666665 * zeta3
                                            - 48. * zeta4)
                                   - 24. * zeta4)))
                          / x
                    - 108. * zeta4 + 144 * x * zeta4
                    - 123.99999999999997 * zeta5 - 7677.333333333331 * x * zeta5
                    + 367.99999999999994 * x2 * zeta5)
           + CA * CA * nf
                 * (557.9650205761317 - 1208.074074074074 * H000
                    + 159.55555555555554 * H0000 - 240 * H00000 - 56 * H00001
                    + 32 * H0001 + 36 * H00010 + 12 * H00011
                    - 109.33333333333333 * H000m10 - 752.1777777777778 * H001
                    - 4 * H0010 + 153.33333333333331 * H00100
                    + 29.333333333333332 * H00101 - 33.33333333333333 * H0011
                    + 98.66666666666666 * H00110 + 32 * H00111 + 384 * H00m10
                    - 90.66666666666666 * H00m100 - 82.66666666666666 * H00m101
                    - 104 * H00m1m10 + 542.6666666666666 * H01
                    - 321.3333333333333 * H010 + 72 * H0100
                    + 154.66666666666666 * H01000 + 297.3333333333333 * H01001
                    - 13.333333333333332 * H0101 + 41.33333333333333 * H01010
                    - 18.666666666666664 * H01011 + 181.33333333333331 * H010m10
                    - 492.88888888888886 * H011 - 24 * H0110
                    - 217.33333333333331 * H01100 - 8 * H01110
                    - 26.666666666666664 * H01111 - 162.4 * H0m10
                    + 274.66666666666663 * H0m100 - 146.66666666666666 * H0m1000
                    - 214.66666666666666 * H0m1001 - 5.333333333333333 * H0m101
                    - 24 * H0m1010 - 10.666666666666666 * H0m1011
                    - 61.33333333333333 * H0m10m10
                    - 117.33333333333333 * H0m1m10
                    + 57.33333333333333 * H0m1m100 + 296 * H0m1m101
                    + 16 * H0m1m1m10 - 1283.1522633744858 * H1
                    + 22.296296296296294 * H10 + 139.1851851851852 * H100
                    - 229.33333333333331 * H10000 - 290.66666666666663 * H10001
                    + 645.7777777777777 * H1001 - 149.33333333333331 * H10010
                    - 194.66666666666666 * H10011 + 117.33333333333333 * H100m10
                    - 114.29629629629629 * H101 - 45.77777777777778 * H1010
                    - 162.66666666666666 * H10100 - 93.33333333333333 * H10101
                    - 57.77777777777777 * H1011 - 117.33333333333333 * H10110
                    - 106.66666666666666 * H10111 + 336. * H10m10
                    + 96. * H10m100 + 101.33333333333333 * H10m101
                    - 5.333333333333333 * H10m1m10 - 277.2098765432099 * H11
                    - 85.7037037037037 * H110 - 653.7777777777777 * H1100
                    - 298.66666666666663 * H11000 - 85.33333333333333 * H11001
                    - 92.44444444444444 * H1101 - 114.66666666666666 * H11010
                    - 101.33333333333333 * H11011 + 144. * H110m10
                    - 159.25925925925924 * H111 - 114.22222222222221 * H1110
                    - 261.3333333333333 * H11100 - 93.33333333333333 * H11101
                    - 95.55555555555556 * H1111 - 82.66666666666666 * H11110
                    - 80. * H11111 + 1203.3604938271606 * Hm10
                    - 40.85925925925926 * Hm100 - 218.2222222222222 * Hm1000
                    + 10.666666666666666 * Hm10000
                    + 146.66666666666666 * Hm10001 - 598.6666666666666 * Hm1001
                    + 112. * Hm10010 + 133.33333333333331 * Hm10011
                    + 48. * Hm100m10 + 3.9111111111111114 * Hm101
                    - 15.999999999999998 * Hm1010 - 80. * Hm10100
                    + 26.666666666666664 * Hm10101 + 7.999999999999999 * Hm1011
                    + 37.33333333333333 * Hm10110 + 10.666666666666666 * Hm10111
                    - 37.33333333333333 * Hm10m10 + 269.3333333333333 * Hm10m100
                    + 272. * Hm10m101 - 170.66666666666666 * Hm10m1m10
                    - 384.88888888888886 * Hm1m10 + 528. * Hm1m100
                    + 80. * Hm1m1000 + 178.66666666666666 * Hm1m1001
                    + 976. * Hm1m101 - 80. * Hm1m1010
                    - 85.33333333333333 * Hm1m1011
                    - 170.66666666666666 * Hm1m10m10 + 80. * Hm1m1m10
                    - 533.3333333333333 * Hm1m1m100 - 560. * Hm1m1m101
                    + 208. * Hm1m1m1m10 - (2.8000000000000003 * H00) / x
                    + (2.1333333333333333 * H000) / x - (3.2 * H001) / x
                    + (148.70617283950617 * H01) / x
                    - (92.44444444444444 * H010) / x
                    + (42.666666666666664 * H0100) / x
                    + (53.33333333333333 * H0101) / x
                    - (121.48148148148148 * H011) / x
                    + (53.33333333333333 * H0110) / x
                    + (56.888888888888886 * H0111) / x
                    + (44.08888888888889 * H0m10) / x
                    - (85.33333333333333 * H0m100) / x
                    - (42.666666666666664 * H0m101) / x
                    + (42.666666666666664 * H0m1m10) / x
                    + (756.1827160493827 * H1) / x
                    - (309.679012345679 * H10) / x
                    + (45.33333333333333 * H100) / x
                    + (120.88888888888889 * H1000) / x
                    + (99.55555555555554 * H1001) / x
                    + (70.81481481481481 * H101) / x
                    + (85.33333333333333 * H1010) / x + (64 * H1011) / x
                    + (14.222222222222221 * H10m10) / x - (332 * H11) / x
                    + (35.85185185185185 * H110) / x
                    + (99.55555555555554 * H1100) / x
                    + (71.11111111111111 * H1101) / x
                    + (80.29629629629629 * H111) / x
                    + (71.11111111111111 * H1110) / x
                    + (46.22222222222222 * H1111) / x
                    + (347.6888888888889 * Hm10) / x
                    - (496.94814814814816 * Hm100) / x
                    + (74.66666666666666 * Hm1000) / x
                    + (42.666666666666664 * Hm1001) / x
                    - (251.97037037037038 * Hm101) / x
                    + (42.666666666666664 * Hm1010) / x
                    + (56.888888888888886 * Hm1011) / x
                    + (49.77777777777777 * Hm10m10) / x
                    + (132.74074074074073 * Hm1m10) / x
                    + (85.33333333333333 * Hm1m100) / x
                    + (71.11111111111111 * Hm1m101) / x
                    - (56.888888888888886 * Hm1m1m10) / x
                    + 11575.300411522634 * x + 5504.118518518519 * H000 * x
                    + 3741.333333333333 * H0000 * x + 1056 * H00000 * x
                    + 1445.3333333333333 * H00001 * x
                    + 4310.666666666666 * H0001 * x + 1256 * H00010 * x
                    + 1432 * H00011 * x + 133.33333333333331 * H000m10 * x
                    + 5990.044444444445 * H001 * x + 3316 * H0010 * x
                    + 1181.3333333333333 * H00100 * x
                    + 1114.6666666666665 * H00101 * x
                    + 3713.333333333333 * H0011 * x
                    + 1189.3333333333333 * H00110 * x
                    + 1066.6666666666665 * H00111 * x
                    + 69.33333333333333 * H00m100 * x
                    + 346.66666666666663 * H00m101 * x
                    + 506.66666666666663 * H00m1m10 * x
                    + 9308.22222222222 * H01 * x + 4810.222222222222 * H010 * x
                    + 3209.333333333333 * H0100 * x
                    + 1205.3333333333333 * H01000 * x + 1512 * H01001 * x
                    + 2480 * H0101 * x + 690.6666666666666 * H01010 * x
                    + 613.3333333333333 * H01011 * x
                    + 21.333333333333332 * H010m10 * x
                    + 5154.222222222222 * H011 * x + 2672 * H0110 * x
                    + 440 * H01100 * x + 650.6666666666666 * H01101 * x
                    + 2421.333333333333 * H0111 * x
                    + 645.3333333333333 * H01110 * x
                    + 469.3333333333333 * H01111 * x + 598.4 * H0m10 * x
                    - 214.66666666666666 * H0m100 * x + 432 * H0m1000 * x
                    + 546.6666666666666 * H0m1001 * x - 488 * H0m101 * x
                    + 304 * H0m1010 * x + 362.66666666666663 * H0m1011 * x
                    + 570.6666666666666 * H0m10m10 * x
                    + 725.3333333333333 * H0m1m10 * x
                    + 621.3333333333333 * H0m1m100 * x + 240 * H0m1m101 * x
                    - 800 * H0m1m1m10 * x + 10918.576131687243 * H1 * x
                    + 5420.666666666666 * H10 * x + 3401.185185185185 * H100 * x
                    + 1194.6666666666665 * H1000 * x
                    + 458.66666666666663 * H10000 * x
                    + 581.3333333333333 * H10001 * x
                    + 1037.7777777777778 * H1001 * x
                    + 298.66666666666663 * H10010 * x
                    + 389.3333333333333 * H10011 * x
                    - 234.66666666666666 * H100m10 * x
                    + 3555.7037037037035 * H101 * x
                    + 1092.888888888889 * H1010 * x
                    + 325.3333333333333 * H10100 * x
                    + 186.66666666666666 * H10101 * x
                    + 1048.888888888889 * H1011 * x
                    + 234.66666666666666 * H10110 * x
                    + 213.33333333333331 * H10111 * x
                    - 250.66666666666666 * H10m10 * x - 192 * H10m100 * x
                    - 202.66666666666666 * H10m101 * x
                    + 10.666666666666666 * H10m1m10 * x
                    + 6270.197530864197 * H11 * x
                    + 3711.4074074074074 * H110 * x
                    + 1572.888888888889 * H1100 * x
                    + 597.3333333333333 * H11000 * x
                    + 170.66666666666666 * H11001 * x
                    + 1144.888888888889 * H1101 * x
                    + 229.33333333333331 * H11010 * x
                    + 202.66666666666666 * H11011 * x - 288 * H110m10 * x
                    + 3508.7407407407404 * H111 * x
                    + 1161.7777777777776 * H1110 * x
                    + 522.6666666666666 * H11100 * x
                    + 186.66666666666666 * H11101 * x
                    + 943.1111111111111 * H1111 * x
                    + 165.33333333333331 * H11110 * x + 160 * H11111 * x
                    + 1827.520987654321 * Hm10 * x
                    + 1103.437037037037 * Hm100 * x
                    + 1027.5555555555554 * Hm1000 * x
                    + 21.333333333333332 * Hm10000 * x
                    + 485.3333333333333 * Hm10001 * x + 628 * Hm1001 * x
                    + 224 * Hm10010 * x + 266.66666666666663 * Hm10011 * x
                    + 96 * Hm100m10 * x + 991.4666666666667 * Hm101 * x
                    + 461.3333333333333 * Hm1010 * x - 352 * Hm10100 * x
                    + 53.33333333333333 * Hm10101 * x + 528 * Hm1011 * x
                    + 74.66666666666666 * Hm10110 * x
                    + 21.333333333333332 * Hm10111 * x
                    + 157.33333333333331 * Hm10m10 * x
                    + 538.6666666666666 * Hm10m100 * x + 544 * Hm10m101 * x
                    - 341.3333333333333 * Hm10m1m10 * x
                    - 85.33333333333333 * Hm1m10 * x
                    + 193.33333333333331 * Hm1m100 * x + 160 * Hm1m1000 * x
                    + 357.3333333333333 * Hm1m1001 * x
                    + 578.6666666666666 * Hm1m101 * x - 160 * Hm1m1010 * x
                    - 170.66666666666666 * Hm1m1011 * x
                    - 341.3333333333333 * Hm1m10m10 * x - 304 * Hm1m1m10 * x
                    - 1066.6666666666665 * Hm1m1m100 * x - 1120 * Hm1m1m101 * x
                    + 416 * Hm1m1m1m10 * x - 10949.21083676269 * x2
                    - 5658.962962962963 * H000 * x2 - 360 * H0000 * x2
                    - 192 * H00001 * x2 - 1633.7777777777776 * H0001 * x2
                    - 160 * H00010 * x2 - 288 * H00011 * x2
                    + 85.33333333333333 * H000m10 * x2
                    - 6529.125925925926 * H001 * x2 - 1296 * H0010 * x2
                    - 32 * H00100 * x2 - 213.33333333333331 * H00101 * x2
                    - 1504.4444444444443 * H0011 * x2
                    - 170.66666666666666 * H00110 * x2 - 192 * H00111 * x2
                    + 939.5555555555555 * H00m10 * x2
                    + 149.33333333333331 * H00m100 * x2
                    + 85.33333333333333 * H00m101 * x2
                    - 213.33333333333331 * H00m1m10 * x2
                    - 6122.898765432099 * H01 * x2 - 4906 * H010 * x2
                    - 1023.5555555555555 * H0100 * x2
                    - 117.33333333333333 * H01000 * x2
                    + 154.66666666666666 * H01001 * x2
                    - 1294.6666666666665 * H0101 * x2 - 176 * H01010 * x2
                    - 229.33333333333331 * H01011 * x2
                    + 202.66666666666666 * H010m10 * x2
                    - 5191.7037037037035 * H011 * x2
                    - 1675.5555555555554 * H0110 * x2
                    - 618.6666666666666 * H01100 * x2
                    - 213.33333333333331 * H01101 * x2
                    - 1200.888888888889 * H0111 * x2
                    - 229.33333333333331 * H01110 * x2 - 192 * H01111 * x2
                    + 141.61481481481482 * H0m10 * x2
                    + 1565.7777777777776 * H0m100 * x2
                    + 170.66666666666666 * H0m1000 * x2
                    + 261.3333333333333 * H0m1001 * x2 + 1024 * H0m101 * x2
                    + 138.66666666666666 * H0m1010 * x2
                    + 170.66666666666666 * H0m1011 * x2
                    - 10.666666666666666 * H0m10m10 * x2
                    - 402.66666666666663 * H0m1m10 * x2
                    - 122.66666666666666 * H0m1m100 * x2
                    + 53.33333333333333 * H0m1m101 * x2 - 64 * H0m1m1m10 * x2
                    - 10677.721810699588 * H1 * x2
                    - 5370.913580246914 * H10 * x2
                    - 3737.4074074074074 * H100 * x2
                    - 1472.888888888889 * H1000 * x2
                    - 298.66666666666663 * H10000 * x2
                    - 437.3333333333333 * H10001 * x2
                    - 1690.6666666666665 * H1001 * x2
                    - 266.66666666666663 * H10010 * x2
                    - 357.3333333333333 * H10011 * x2
                    + 202.66666666666666 * H100m10 * x2
                    - 3678.7407407407404 * H101 * x2
                    - 1270.2222222222222 * H1010 * x2
                    - 245.33333333333331 * H10100 * x2
                    - 186.66666666666666 * H10101 * x2
                    - 1162.2222222222222 * H1011 * x2
                    - 234.66666666666666 * H10110 * x2
                    - 213.33333333333331 * H10111 * x2
                    + 47.11111111111111 * H10m10 * x2 + 160 * H10m100 * x2
                    + 74.66666666666666 * H10m101 * x2
                    - 74.66666666666666 * H10m1m10 * x2
                    - 6001.456790123457 * H11 * x2
                    - 3986.148148148148 * H110 * x2
                    - 1376.4444444444443 * H1100 * x2
                    - 405.3333333333333 * H11000 * x2
                    - 234.66666666666666 * H11001 * x2 - 1260 * H1101 * x2
                    - 229.33333333333331 * H11010 * x2
                    - 202.66666666666666 * H11011 * x2 + 160 * H110m10 * x2
                    - 3656.148148148148 * H111 * x2
                    - 1294.2222222222222 * H1110 * x2
                    - 394.66666666666663 * H11100 * x2
                    - 186.66666666666666 * H11101 * x2 - 1024 * H1111 * x2
                    - 165.33333333333331 * H11110 * x2 - 160 * H11111 * x2
                    + 1204.083950617284 * Hm10 * x2
                    + 998.9777777777778 * Hm100 * x2
                    + 1531.5555555555554 * Hm1000 * x2
                    + 181.33333333333331 * Hm10000 * x2
                    + 613.3333333333333 * Hm10001 * x2
                    + 1458.6666666666665 * Hm1001 * x2 + 256 * Hm10010 * x2
                    + 298.66666666666663 * Hm10011 * x2 - 64 * Hm100m10 * x2
                    + 885.1407407407407 * Hm101 * x2
                    + 581.3333333333333 * Hm1010 * x2 - 256 * Hm10100 * x2
                    + 53.33333333333333 * Hm10101 * x2
                    + 638.2222222222222 * Hm1011 * x2
                    + 74.66666666666666 * Hm10110 * x2
                    + 21.333333333333332 * Hm10111 * x2
                    + 167.11111111111111 * Hm10m10 * x2
                    + 122.66666666666666 * Hm10m100 * x2 + 160 * Hm10m101 * x2
                    - 149.33333333333331 * Hm10m1m10 * x2
                    + 208.07407407407408 * Hm1m10 * x2
                    - 470.66666666666663 * Hm1m100 * x2 - 224 * Hm1m1000 * x2
                    - 122.66666666666666 * Hm1m1001 * x2
                    - 398.2222222222222 * Hm1m101 * x2 - 224 * Hm1m1010 * x2
                    - 234.66666666666666 * Hm1m1011 * x2
                    - 149.33333333333331 * Hm1m10m10 * x2
                    - 379.55555555555554 * Hm1m1m10 * x2
                    - 394.66666666666663 * Hm1m1m100 * x2 - 480 * Hm1m1m101 * x2
                    + 224 * Hm1m1m1m10 * x2 - 24 * H000 * x3
                    + 76.80000000000001 * H0000 * x3 - 115.2 * H0001 * x3
                    - 24 * H001 * x3 - 76.80000000000001 * H00m10 * x3
                    + 192 * H0100 * x3 + 24 * H0m10 * x3
                    - 76.80000000000001 * H0m100 * x3
                    - 153.60000000000002 * H0m101 * x3
                    + 38.400000000000006 * H1000 * x3 - 192 * H1001 * x3
                    - 76.80000000000001 * H10m10 * x3 + 192 * H1100 * x3
                    + 45.6 * Hm10 * x3 + 24 * Hm100 * x3
                    - 38.400000000000006 * Hm1000 * x3
                    - 76.80000000000001 * Hm1001 * x3 + 24 * Hm101 * x3
                    - 24 * Hm1m10 * x3 + 76.80000000000001 * Hm1m100 * x3
                    + 153.60000000000002 * Hm1m101 * x3
                    + (5.333333333333333 * H1001 + 2.1333333333333333 * H10m10
                       - 5.333333333333333 * H1100 + 1.2666666666666666 * Hm10
                       + 0.6666666666666666 * Hm100
                       - 1.0666666666666667 * Hm1000
                       - 2.1333333333333333 * Hm1001
                       + 0.6666666666666666 * Hm101
                       - 0.6666666666666666 * Hm1m10
                       + 2.1333333333333333 * Hm1m100
                       + 4.266666666666667 * Hm1m101 - 1376.4743484224966 * x
                       + H1000 * (-1.0666666666666667 + 21.333333333333332 * x2)
                      ) / x2
                    - 542.6666666666666 * zeta2 + 56 * H000 * zeta2
                    - 81.33333333333333 * H001 * zeta2
                    + 30.666666666666664 * H00m1 * zeta2
                    - 45.33333333333333 * H01 * zeta2
                    - 237.33333333333331 * H010 * zeta2 - 8 * H011 * zeta2
                    - 53.33333333333333 * H0m1 * zeta2
                    + 274.66666666666663 * H0m10 * zeta2 - 288 * H0m1m1 * zeta2
                    - 78.14814814814814 * H1 * zeta2
                    - 496.4444444444444 * H10 * zeta2
                    + 373.3333333333333 * H100 * zeta2
                    + 178.66666666666666 * H101 * zeta2 - 104 * H10m1 * zeta2
                    + 52.44444444444444 * H11 * zeta2
                    + 242.66666666666666 * H110 * zeta2
                    + 197.33333333333331 * H111 * zeta2
                    - 196.35555555555555 * Hm1 * zeta2 + 748 * Hm10 * zeta2
                    - 64 * Hm100 * zeta2 - 24 * Hm101 * zeta2
                    - 357.3333333333333 * Hm10m1 * zeta2 - 936 * Hm1m1 * zeta2
                    - 336 * Hm1m10 * zeta2 + 664 * Hm1m1m1 * zeta2
                    - (0.3333333333333333 * H1 * zeta2) / x2
                    - (4.266666666666667 * H10 * zeta2) / x2
                    - (Hm1 * zeta2) / x2 + (3.2 * Hm10 * zeta2) / x2
                    - (4.266666666666667 * Hm1m1 * zeta2) / x2
                    + (198.9827160493827 * zeta2) / x
                    - (74.66666666666666 * H01 * zeta2) / x
                    + (64 * H0m1 * zeta2) / x
                    - (137.1851851851852 * H1 * zeta2) / x
                    - (117.33333333333333 * H10 * zeta2) / x
                    - (99.55555555555554 * H11 * zeta2) / x
                    + (318.3407407407407 * Hm1 * zeta2) / x
                    - (24.888888888888886 * Hm10 * zeta2) / x
                    - (99.55555555555554 * Hm1m1 * zeta2) / x
                    - 7480.701234567901 * x * zeta2 - 1312 * H000 * x * zeta2
                    - 1368 * H001 * x * zeta2
                    - 93.33333333333333 * H00m1 * x * zeta2
                    - 2842.6666666666665 * H01 * x * zeta2
                    - 1786.6666666666665 * H010 * x * zeta2
                    - 1050.6666666666665 * H011 * x * zeta2
                    + 850.6666666666666 * H0m1 * x * zeta2
                    - 272 * H0m10 * x * zeta2 - 640 * H0m1m1 * x * zeta2
                    - 3513.037037037037 * H1 * x * zeta2
                    - 1241.7777777777776 * H10 * x * zeta2
                    - 746.6666666666666 * H100 * x * zeta2
                    - 357.3333333333333 * H101 * x * zeta2
                    + 208 * H10m1 * x * zeta2
                    - 1296.888888888889 * H11 * x * zeta2
                    - 485.3333333333333 * H110 * x * zeta2
                    - 394.66666666666663 * H111 * x * zeta2
                    - 1034.1333333333332 * Hm1 * x * zeta2
                    - 424 * Hm10 * x * zeta2 - 320 * Hm100 * x * zeta2
                    - 48 * Hm101 * x * zeta2
                    - 714.6666666666666 * Hm10m1 * x * zeta2
                    - 730.6666666666666 * Hm1m1 * x * zeta2
                    - 672 * Hm1m10 * x * zeta2 + 1328 * Hm1m1m1 * x * zeta2
                    + 6122.898765432099 * x2 * zeta2 + 192 * H000 * x2 * zeta2
                    + 106.66666666666666 * H001 * x2 * zeta2
                    - 192 * H00m1 * x2 * zeta2
                    + 1093.3333333333333 * H01 * x2 * zeta2
                    - 58.666666666666664 * H010 * x2 * zeta2
                    + 245.33333333333331 * H011 * x2 * zeta2
                    - 1225.3333333333333 * H0m1 * x2 * zeta2
                    - 165.33333333333331 * H0m10 * x2 * zeta2
                    - 85.33333333333333 * H0m1m1 * x2 * zeta2
                    + 3782.7777777777774 * H1 * x2 * zeta2
                    + 1797.7777777777776 * H10 * x2 * zeta2
                    + 506.66666666666663 * H100 * x2 * zeta2
                    + 261.3333333333333 * H101 * x2 * zeta2
                    - 112 * H10m1 * x2 * zeta2
                    + 1449.7777777777776 * H11 * x2 * zeta2
                    + 389.3333333333333 * H110 * x2 * zeta2
                    + 298.66666666666663 * H111 * x2 * zeta2
                    - 781.1037037037038 * Hm1 * x2 * zeta2
                    - 1351.5555555555554 * Hm10 * x2 * zeta2
                    - 544 * Hm100 * x2 * zeta2 - 16 * Hm101 * x2 * zeta2
                    - 234.66666666666666 * Hm10m1 * x2 * zeta2
                    + 208.44444444444443 * Hm1m1 * x2 * zeta2
                    - 32 * Hm1m10 * x2 * zeta2 + 592 * Hm1m1m1 * x2 * zeta2
                    + 45.6 * x3 * zeta2 + 153.60000000000002 * H0m1 * x3 * zeta2
                    + 12 * H1 * x3 * zeta2
                    + 153.60000000000002 * H10 * x3 * zeta2
                    - 36 * Hm1 * x3 * zeta2 + 115.2 * Hm10 * x3 * zeta2
                    - 153.60000000000002 * Hm1m1 * x3 * zeta2
                    + 63.86666666666667 * zeta2_2 + 3.2 * H1 * zeta2_2
                    + 45.33333333333333 * Hm1 * zeta2_2
                    + (100.97777777777779 * zeta2_2) / x
                    + 824.6222222222223 * x * zeta2_2 - 6.4 * H1 * x * zeta2_2
                    + 244.26666666666665 * Hm1 * x * zeta2_2
                    - 757.3333333333333 * x2 * zeta2_2
                    - 27.200000000000003 * H1 * x2 * zeta2_2
                    + 197.86666666666667 * Hm1 * x2 * zeta2_2
                    - 207.36 * x3 * zeta2_2
                    + H00
                          * (1044.6098765432098 - 32. * zeta2
                             + x
                                   * (7912.325925925927
                                      - 4310.666666666666 * zeta2
                                      + x
                                            * (-5109.614814814815
                                               + 1633.7777777777776 * zeta2
                                               + x
                                                     * (-45.6
                                                        + 38.400000000000006
                                                              * zeta2))
                                      - 1477.3333333333333 * zeta3)
                             - 378.66666666666663 * zeta3)
                    - 353.46666666666664 * zeta3
                    + 413.3333333333333 * H01 * zeta3
                    + 225.33333333333331 * H0m1 * zeta3
                    + 814.6666666666666 * H1 * zeta3
                    + 621.3333333333333 * H10 * zeta3 + 656 * H11 * zeta3
                    + 706.6666666666666 * Hm1 * zeta3
                    + 221.33333333333331 * Hm10 * zeta3
                    - 485.3333333333333 * Hm1m1 * zeta3
                    + (7.466666666666667 * H1 * zeta3) / x2
                    + (3.2 * Hm1 * zeta3) / x2 + (70.51851851851852 * zeta3) / x
                    - (104.88888888888889 * H1 * zeta3) / x
                    + (48 * Hm1 * zeta3) / x - 3365.7777777777774 * x * zeta3
                    - 1146.6666666666665 * H01 * x * zeta3
                    + 338.66666666666663 * H0m1 * x * zeta3
                    - 2050.6666666666665 * H1 * x * zeta3
                    - 1242.6666666666665 * H10 * x * zeta3
                    - 1312 * H11 * x * zeta3
                    + 257.3333333333333 * Hm1 * x * zeta3
                    + 634.6666666666666 * Hm10 * x * zeta3
                    - 970.6666666666666 * Hm1m1 * x * zeta3
                    + 5519.933333333333 * x2 * zeta3 + 992 * H01 * x2 * zeta3
                    - 37.33333333333333 * H0m1 * x2 * zeta3
                    + 1815.5555555555554 * H1 * x2 * zeta3
                    + 826.6666666666666 * H10 * x2 * zeta3
                    + 848 * H11 * x2 * zeta3 - 532 * Hm1 * x2 * zeta3
                    + 202.66666666666666 * Hm10 * x2 * zeta3
                    - 346.66666666666663 * Hm1m1 * x2 * zeta3 + 60 * x3 * zeta3
                    - 268.8 * H1 * x3 * zeta3 + 115.2 * Hm1 * x3 * zeta3
                    + 214.66666666666666 * zeta2 * zeta3
                    + 1448 * x * zeta2 * zeta3
                    - 229.33333333333331 * x2 * zeta2 * zeta3
                    + (H0
                       * (-160.96872427983539 + 47.28888888888889 * zeta2
                          - 14.222222222222221 * zeta3
                          + x
                                * (-2650.351440329218
                                   + (752.1777777777778
                                      - 52.266666666666666 * zeta2)
                                         * zeta2
                                   + x3 * (48. * zeta2 - 384. * zeta3)
                                   - 46.22222222222222 * zeta3
                                   + x
                                         * (10963.58683127572
                                            + zeta2
                                                  * (-5391.644444444445
                                                     + 741.3333333333333 * zeta2
                                                  )
                                            - 5904.888888888889 * zeta3
                                            - 192. * zeta4)
                                   + x2
                                         * (-10913.509465020576
                                            + (6529.125925925926
                                               - 103.46666666666667 * zeta2)
                                                  * zeta2
                                            + 2166.222222222222 * zeta3
                                            - 48. * zeta4)
                                   - 72. * zeta4)))
                          / x
                    - 84. * zeta4 + 24 * H1 * zeta4 - (32 * zeta4) / x
                    - 96 * x * zeta4 - 48 * H1 * x * zeta4 + 248 * x2 * zeta4
                    + 48 * H1 * x2 * zeta4 - 870.6666666666666 * zeta5
                    - 1624. * x * zeta5 - 581.3333333333331 * x2 * zeta5)
           + CF * nf * nf
                 * (2262.7246502057615 + 658.5925925925926 * H000
                    + 181.77777777777777 * H0000 + 120 * H00000 + 88 * H00001
                    + 145.33333333333331 * H0001 + 48 * H00010 + 48 * H00011
                    + 489.77777777777777 * H001 + 88 * H0010 + 16 * H00100
                    + 16 * H00101 + 88 * H0011 + 16 * H00110 + 16 * H00111
                    - 160 * H00m10 + 697.7975308641975 * H01
                    + 270.66666666666663 * H010 + 24 * H0100
                    + 34.666666666666664 * H0101 + 273.6296296296296 * H011
                    + 34.666666666666664 * H0110 + 33.77777777777778 * H0111
                    - 170.66666666666666 * H0m100 - 53.33333333333333 * H0m101
                    + 74.66666666666666 * H0m1m10 + 1110.1711934156378 * H1
                    + 416.7654320987654 * H10 + 114.07407407407408 * H100
                    + 8.88888888888889 * H1000 + 21.333333333333332 * H1001
                    + 141.33333333333331 * H101 + 32. * H1010
                    + 26.666666666666664 * H1011 - 21.333333333333332 * H10m10
                    + 416.8641975308642 * H11 + 116.44444444444446 * H110
                    + 5.333333333333333 * H1100 + 10.666666666666666 * H1101
                    + 119.4074074074074 * H111 + 21.333333333333332 * H1110
                    + 15.11111111111111 * H1111 - 518.0444444444445 * Hm10
                    - 273.77777777777777 * Hm100 - 32. * Hm1000
                    - 10.666666666666666 * Hm1001 - 80. * Hm101
                    + 42.666666666666664 * Hm10m10 + 147.55555555555554 * Hm1m10
                    + 74.66666666666666 * Hm1m100 + 21.333333333333332 * Hm1m101
                    - 42.666666666666664 * Hm1m1m10
                    + (1.0666666666666667 * H00) / x
                    + (0.3555555555555556 * H01 - 34.33086419753087 * H1
                       + 18.962962962962962 * H10 - 7.111111111111109 * H100
                       - 7.111111111111109 * H101 + 18.962962962962962 * H11
                       - 7.111111111111109 * H110 - 7.111111111111109 * H111
                       - 19.31851851851852 * Hm10 + 21.333333333333332 * Hm100
                       + 7.111111111111109 * Hm101 - 7.111111111111109 * Hm1m10)
                          / x
                    - 3215.919094650206 * x - 1007.4074074074074 * H000 * x
                    - 803.5555555555555 * H0000 * x - 240 * H00000 * x
                    - 176 * H00001 * x - 442.66666666666663 * H0001 * x
                    - 96 * H00010 * x - 96 * H00011 * x - 372 * H001 * x
                    - 160 * H0010 * x - 32 * H00100 * x - 32 * H00101 * x
                    - 160 * H0011 * x - 32 * H00110 * x - 32 * H00111 * x
                    + 64 * H00m10 * x - 506.49876543209876 * H01 * x
                    - 133.33333333333331 * H010 * x
                    - 21.333333333333332 * H0101 * x
                    - 142.8148148148148 * H011 * x
                    - 21.333333333333332 * H0110 * x
                    - 19.555555555555554 * H0111 * x
                    + 241.77777777777777 * H0m10 * x
                    + 170.66666666666666 * H0m100 * x + 64 * H0m101 * x
                    - 21.333333333333332 * H0m1m10 * x
                    - 1506.3522633744856 * H1 * x - 722.2716049382716 * H10 * x
                    - 239.7037037037037 * H100 * x
                    - 17.77777777777778 * H1000 * x
                    - 42.666666666666664 * H1001 * x
                    - 290.66666666666663 * H101 * x - 64 * H1010 * x
                    - 53.33333333333333 * H1011 * x
                    + 42.666666666666664 * H10m10 * x
                    - 685.8765432098766 * H11 * x
                    - 219.55555555555554 * H110 * x
                    - 10.666666666666666 * H1100 * x
                    - 21.333333333333332 * H1101 * x
                    - 229.03703703703704 * H111 * x
                    - 42.666666666666664 * H1110 * x
                    - 30.22222222222222 * H1111 * x
                    - 282.7851851851852 * Hm10 * x
                    - 110.22222222222221 * Hm100 * x - 64 * Hm1000 * x
                    - 21.333333333333332 * Hm1001 * x
                    - 14.222222222222221 * Hm101 * x
                    + 85.33333333333333 * Hm10m10 * x
                    + 149.33333333333331 * Hm1m10 * x
                    + 149.33333333333331 * Hm1m100 * x
                    + 42.666666666666664 * Hm1m101 * x
                    - 85.33333333333333 * Hm1m1m10 * x + 881.8127297668038 * x2
                    + 7.111111111111111 * H000 * x2
                    - 103.1111111111111 * H0000 * x2
                    - 53.33333333333333 * H0001 * x2
                    + 101.33333333333333 * H001 * x2
                    - 10.666666666666666 * H0010 * x2
                    - 10.666666666666666 * H0011 * x2 - 128 * H00m10 * x2
                    + 204.89876543209877 * H01 * x2
                    + 174.2222222222222 * H010 * x2
                    + 10.666666666666666 * H0101 * x2
                    + 146.37037037037035 * H011 * x2
                    + 10.666666666666666 * H0110 * x2
                    + 8.88888888888889 * H0111 * x2
                    - 19.555555555555554 * H0m10 * x2
                    - 138.66666666666666 * H0m100 * x2
                    - 42.666666666666664 * H0m101 * x2 + 64 * H0m1m10 * x2
                    + 417.8534979423868 * H1 * x2 + 329.9753086419753 * H10 * x2
                    + 197.92592592592592 * H100 * x2
                    + 49.77777777777777 * H1000 * x2 + 64 * H1001 * x2
                    + 177.77777777777777 * H101 * x2
                    + 53.33333333333333 * H1010 * x2
                    + 53.33333333333333 * H1011 * x2
                    - 21.333333333333332 * H10m10 * x2
                    + 288.69135802469134 * H11 * x2
                    + 181.33333333333331 * H110 * x2
                    + 21.333333333333332 * H1100 * x2 + 32 * H1101 * x2
                    + 153.48148148148147 * H111 * x2 + 32 * H1110 * x2
                    + 30.22222222222222 * H1111 * x2
                    + 154.9037037037037 * Hm10 * x2
                    + 131.55555555555554 * Hm100 * x2 - 32 * Hm1000 * x2
                    - 10.666666666666666 * Hm1001 * x2
                    + 55.11111111111111 * Hm101 * x2
                    + 42.666666666666664 * Hm10m10 * x2
                    + 12.444444444444443 * Hm1m10 * x2
                    + 74.66666666666666 * Hm1m100 * x2
                    + 21.333333333333332 * Hm1m101 * x2
                    - 42.666666666666664 * Hm1m1m10 * x2
                    + 38.400000000000006 * H000 * x3
                    + 44.800000000000004 * H001 * x3 - 32 * H010 * x3
                    - 12.8 * H0m10 * x3 + 32 * H101 * x3 - 32 * H110 * x3
                    - 45.653333333333336 * Hm10 * x3
                    - 38.400000000000006 * Hm100 * x3 - 12.8 * Hm101 * x3
                    + 12.8 * Hm1m10 * x3
                    + (-1.5051851851851852 * Hm10 - 1.0666666666666667 * Hm100
                       - 0.3555555555555556 * Hm101
                       + 0.3555555555555556 * Hm1m10 + 60.37245541838134 * x
                       + H0m10 * (-0.7111111111111112 - 275.55555555555554 * x2)
                      ) / x2
                    - 697.7975308641976 * zeta2 - 88 * H000 * zeta2
                    - 16 * H001 * zeta2 + 2.6666666666666665 * H01 * zeta2
                    + 90.66666666666666 * H0m1 * zeta2
                    - 67.55555555555556 * H1 * zeta2
                    - 10.666666666666666 * H10 * zeta2
                    + 10.666666666666666 * H11 * zeta2
                    + 153.77777777777777 * Hm1 * zeta2
                    + 21.333333333333332 * Hm10 * zeta2
                    - 42.666666666666664 * Hm1m1 * zeta2
                    + (0.17777777777777778 * H1 * zeta2) / x2
                    + (0.5333333333333333 * Hm1 * zeta2) / x2
                    - (19.674074074074074 * zeta2) / x
                    + (10.666666666666666 * H1 * zeta2) / x
                    - (10.666666666666666 * Hm1 * zeta2) / x
                    + 223.71358024691358 * x * zeta2 + 176 * H000 * x * zeta2
                    + 32 * H001 * x * zeta2 + 32 * H01 * x * zeta2
                    - 74.66666666666666 * H0m1 * x * zeta2
                    + 216 * H1 * x * zeta2
                    + 21.333333333333332 * H10 * x * zeta2
                    - 21.333333333333332 * H11 * x * zeta2
                    + 88.88888888888889 * Hm1 * x * zeta2
                    + 42.666666666666664 * Hm10 * x * zeta2
                    - 85.33333333333333 * Hm1m1 * x * zeta2
                    - 204.89876543209877 * x2 * zeta2
                    + 21.333333333333332 * H01 * x2 * zeta2
                    + 74.66666666666666 * H0m1 * x2 * zeta2
                    - 171.55555555555554 * H1 * x2 * zeta2
                    - 53.33333333333333 * H10 * x2 * zeta2
                    - 10.666666666666666 * H11 * x2 * zeta2
                    - 48.888888888888886 * Hm1 * x2 * zeta2
                    + 21.333333333333332 * Hm10 * x2 * zeta2
                    - 42.666666666666664 * Hm1m1 * x2 * zeta2
                    - 45.653333333333336 * x3 * zeta2
                    - 38.400000000000006 * H1 * x3 * zeta2
                    + 19.200000000000003 * Hm1 * x3 * zeta2
                    - 20.177777777777777 * zeta2_2
                    - 150.04444444444445 * x * zeta2_2
                    - 28.08888888888889 * x2 * zeta2_2
                    - 386.5185185185185 * zeta3
                    - 12.444444444444443 * H1 * zeta3
                    + 37.33333333333333 * Hm1 * zeta3
                    + (4.7407407407407405 * zeta3) / x
                    + 680.1481481481482 * x * zeta3
                    + 24.888888888888886 * H1 * x * zeta3
                    + 74.66666666666666 * Hm1 * x * zeta3
                    - 370.66666666666663 * x2 * zeta3
                    - 30.22222222222222 * H1 * x2 * zeta3
                    + 37.33333333333333 * Hm1 * x2 * zeta3 + 64 * x3 * zeta3
                    + 32 * zeta2 * zeta3 - 64 * x * zeta2 * zeta3
                    + H00 * (951.2938271604938 - 145.33333333333331 * zeta2 - 50.666666666666664 * zeta3 + x * (-686.1876543209877 + 506.66666666666663 * zeta2 + x * (-69.35308641975308 + 45.653333333333336 * x + 53.33333333333333 * zeta2) + 101.33333333333333 * zeta3))
                    + (H0
                       * (1.1496296296296296
                          + x
                                * (1746.1478189300412 - 57.6 * x3 * zeta2
                                   + zeta2
                                         * (-489.7777777777777
                                            + 11.200000000000001 * zeta2)
                                   + x2
                                         * (518.1175308641975
                                            - 101.33333333333333 * zeta2
                                            - 21.333333333333332 * zeta3)
                                   - 140.88888888888889 * zeta3
                                   + x
                                         * (-1422.6373662551441
                                            + (613.7777777777777
                                               - 22.400000000000002 * zeta2)
                                                  * zeta2
                                            + 375.1111111111111 * zeta3))))
                          / x
                    - 24. * zeta5 + 48 * x * zeta5)
           + (CA * nf * nf
              * (Hm10
                     * (-0.26666666666666666
                        + x
                              * (5.925925925925926
                                 + x
                                       * (72.69135802469135 + 8. * zeta2
                                          + x
                                                * (32.79012345679012
                                                   + 16. * zeta2
                                                   + x
                                                         * (-62.17283950617284
                                                            - 9.600000000000001
                                                                  * x
                                                            + 26.66666666666666
                                                                  * zeta2)))))
                 + x
                       * (-1.0513031550068588 - 5.925925925925926 * H10
                          - 3.5555555555555554 * H100
                          - 3.5555555555555554 * H101 - 5.925925925925926 * H11
                          - 3.5555555555555554 * H110
                          - 3.5555555555555554 * H111
                          + 10.666666666666666 * Hm100
                          + 3.5555555555555554 * Hm101
                          - 3.5555555555555554 * Hm1m10
                          + (-64.23786008230452 - 13.432098765432098 * H00
                             - 19.703703703703702 * H000
                             - 0.8888888888888888 * H0000
                             + 2.6666666666666665 * H0001
                             - 3.1111111111111107 * H001
                             + 10.666666666666666 * H00m10
                             + 11.11111111111111 * H01
                             + 2.6666666666666665 * H010
                             + 5.333333333333333 * H0100
                             - 5.333333333333333 * H0101
                             + 2.6666666666666665 * H011
                             + 5.333333333333333 * H0110
                             + 37.33333333333333 * H0m10
                             + 5.333333333333333 * H0m100
                             - 10.666666666666666 * H0m1m10)
                                * x
                          + 5.925925925925926 * zeta2
                          - 5.333333333333333 * Hm1 * zeta2
                          + 2.3703703703703702 * zeta3
                          + x
                                * (18.666666666666664 * H1000
                                   + 0.8888888888888888 * H1001
                                   + 0.2962962962962963 * H101
                                   + 4.444444444444445 * H1010
                                   - 0.8888888888888888 * H1011
                                   + 10.666666666666666 * H10m10
                                   + 13.728395061728394 * H11
                                   + 22.814814814814813 * H110
                                   + 20.444444444444443 * H1100
                                   + 9.777777777777777 * H1101
                                   + 8.592592592592592 * H111
                                   + 6.222222222222221 * H1110
                                   + 6.222222222222221 * H1111
                                   + 44.148148148148145 * Hm100
                                   - 12.444444444444443 * Hm1000
                                   - 10.666666666666666 * Hm1001
                                   + 7.111111111111111 * Hm101
                                   - 5.333333333333333 * Hm1010
                                   - 5.333333333333333 * Hm1011 - 16. * Hm10m10
                                   - 40.888888888888886 * Hm1m10
                                   - 21.333333333333332 * Hm1m100
                                   - 5.333333333333333 * Hm1m101
                                   + 16. * Hm1m1m10
                                   + (-1186.1201646090535
                                      - 804.8888888888888 * H00
                                      - 380.7407407407407 * H000
                                      - 90.66666666666666 * H0000 - 80. * H0001
                                      - 295.55555555555554 * H001
                                      - 53.33333333333332 * H0010
                                      - 53.33333333333332 * H0011
                                      - 565.4814814814813 * H01
                                      - 172.88888888888889 * H010
                                      - 42.666666666666664 * H0100
                                      - 21.333333333333332 * H0101
                                      - 175.11111111111111 * H011
                                      - 42.666666666666664 * H0110 - 32. * H0111
                                      + 15.11111111111111 * H0m10
                                      + 10.666666666666666 * H0m100
                                      - 21.333333333333332 * H0m1m10)
                                         * x
                                   - 37.33333333333333 * H1000 * x
                                   - 1.7777777777777777 * H1001 * x
                                   - 25.48148148148148 * H101 * x
                                   - 8.88888888888889 * H1010 * x
                                   + 1.7777777777777777 * H1011 * x
                                   - 21.333333333333332 * H10m10 * x
                                   - 269.53086419753083 * H11 * x
                                   - 116.74074074074073 * H110 * x
                                   - 40.888888888888886 * H1100 * x
                                   - 19.555555555555554 * H1101 * x
                                   - 66.96296296296296 * H111 * x
                                   - 12.444444444444443 * H1110 * x
                                   - 12.444444444444443 * H1111 * x
                                   + 4.7407407407407405 * Hm100 * x
                                   - 24.888888888888886 * Hm1000 * x
                                   - 21.333333333333332 * Hm1001 * x
                                   - 12.444444444444443 * Hm101 * x
                                   - 10.666666666666666 * Hm1010 * x
                                   - 10.666666666666666 * Hm1011 * x
                                   - 32. * Hm10m10 * x
                                   - 55.11111111111111 * Hm1m10 * x
                                   - 42.666666666666664 * Hm1m100 * x
                                   - 10.666666666666666 * Hm1m101 * x
                                   + 32. * Hm1m1m10 * x
                                   + 1291.9031550068587 * x2
                                   + 251.20987654320987 * H00 * x2
                                   + 69.62962962962962 * H000 * x2
                                   + 10.666666666666666 * H0001 * x2
                                   + 105.77777777777777 * H001 * x2
                                   + 10.666666666666666 * H0010 * x2
                                   + 10.666666666666666 * H0011 * x2
                                   - 10.666666666666666 * H00m10 * x2
                                   + 360. * H01 * x2
                                   + 93.33333333333333 * H010 * x2
                                   + 16. * H0100 * x2
                                   + 93.33333333333333 * H011 * x2
                                   + 21.333333333333332 * H0110 * x2
                                   + 10.666666666666666 * H0111 * x2
                                   - 30.22222222222222 * H0m10 * x2
                                   - 26.66666666666666 * H0m100 * x2
                                   - 10.666666666666666 * H0m101 * x2
                                   + 21.333333333333332 * H1000 * x2
                                   - 8.88888888888889 * H1001 * x2
                                   + 31.703703703703702 * H101 * x2
                                   + 14.222222222222221 * H1010 * x2
                                   - 1.7777777777777777 * H1011 * x2
                                   + 10.666666666666666 * H10m10 * x2
                                   + 288.7901234567902 * H11 * x2
                                   + 122.96296296296296 * H110 * x2
                                   + 35.55555555555556 * H1100 * x2
                                   + 14.222222222222221 * H1101 * x2
                                   + 73.18518518518518 * H111 * x2
                                   + 17.77777777777778 * H1110 * x2
                                   + 12.444444444444443 * H1111 * x2
                                   - 65.48148148148148 * Hm100 * x2
                                   - 40.888888888888886 * Hm1000 * x2
                                   - 26.66666666666666 * Hm1001 * x2
                                   - 30.22222222222222 * Hm101 * x2
                                   - 10.666666666666666 * Hm1010 * x2
                                   - 10.666666666666666 * Hm1011 * x2
                                   - 10.666666666666666 * Hm10m10 * x2
                                   - 3.5555555555555554 * Hm1m10 * x2
                                   - 5.333333333333333 * Hm1m100 * x2
                                   + 10.666666666666666 * Hm1m1m10 * x2
                                   + 9.600000000000001 * H00 * x3
                                   + H100
                                         * (25.48148148148148
                                            + x
                                                  * (-100.74074074074073
                                                     + 90.07407407407408 * x))
                                   - 11.11111111111111 * zeta2
                                   - 2.6666666666666665 * H00 * zeta2
                                   - 5.333333333333333 * H0m1 * zeta2
                                   - 17.77777777777778 * H11 * zeta2
                                   - 27.555555555555554 * Hm1 * zeta2
                                   + 13.33333333333333 * Hm1m1 * zeta2
                                   + 598.2716049382716 * x * zeta2
                                   + 80. * H00 * x * zeta2
                                   + 32. * H01 * x * zeta2
                                   - 10.666666666666666 * H0m1 * x * zeta2
                                   + 35.55555555555556 * H11 * x * zeta2
                                   - 15.11111111111111 * Hm1 * x * zeta2
                                   + 26.66666666666666 * Hm1m1 * x * zeta2
                                   - 360. * x2 * zeta2
                                   - 10.666666666666666 * H00 * x2 * zeta2
                                   + 10.666666666666666 * H0m1 * x2 * zeta2
                                   - 19.555555555555554 * H11 * x2 * zeta2
                                   + 28.444444444444443 * Hm1 * x2 * zeta2
                                   + 5.333333333333333 * Hm1m1 * x2 * zeta2
                                   - 9.600000000000001 * x3 * zeta2
                                   + 3.7333333333333334 * zeta2_2 - 15.28888888888889 * x * zeta2_2 + 6.933333333333334 * x2 * zeta2_2 + H10 * (11.851851851851851 - 3.5555555555555554 * zeta2 + x * (-265.48148148148147 + 7.111111111111111 * zeta2 + x * (282.96296296296293 + 8.88888888888889 * zeta2)))
                                   - 29.777777777777775 * zeta3
                                   - 8. * Hm1 * zeta3
                                   + 307.55555555555554 * x * zeta3
                                   - 16. * Hm1 * x * zeta3
                                   - 110.81481481481481 * x2 * zeta3
                                   + 2.6666666666666665 * Hm1 * x2 * zeta3)
                          + H0
                                * (0.26666666666666666
                                   + x
                                         * (-48.07572016460905
                                            + 3.1111111111111107 * zeta2
                                            - 7.111111111111111 * zeta3
                                            + x
                                                  * (-1122.4460905349795
                                                     + 310.66666666666663
                                                           * zeta2
                                                     + x
                                                           * (704.9086419753087
                                                              - 105.77777777777777
                                                                    * zeta2
                                                              - 16. * zeta3)
                                                     + 78.22222222222221 * zeta3
                                                  )))
                          + H1
                                * (2.962962962962963 + 5.333333333333333 * zeta2
                                   + x
                                         * (3.547325102880659
                                            - 20.74074074074074 * zeta2
                                            + 8. * zeta3
                                            + x
                                                  * (-745.5144032921811
                                                     + 53.03703703703703 * zeta2
                                                     - 16. * zeta3
                                                     + x
                                                           * (783.4403292181071
                                                              - 33.48148148148148
                                                                    * zeta2
                                                              + 18.666666666666664
                                                                    * zeta3)))))
              )) / x2;
}

//==========================================================================================//
//  Massless quark coefficient functions for F2 at
//  O(as^3) for mu=Q.
//  The term fl_ps_11 is put to zero for the reason
//  explained in page 15 of arXiv:1205.5727
//
//  Eq. (B.10) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double MasslessCoefficientFunction::C2_ps3_massless(double x, int nf) const {

    // the commented varibles are needed for the flav term

    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;

    // double xm = 1 - x;

    // double xp = 1 + x;

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
    const double H0m1m1m10 = Hr5[82];
    const double H00m1m10 = Hr5[85];
    const double H0m10m10 = Hr5[91];
    const double H000m10 = Hr5[94];
    const double H010m10 = Hr5[97];
    const double H0m1m100 = Hr5[109];
    const double H00m100 = Hr5[112];
    const double H0m1000 = Hr5[118];
    const double H00000 = Hr5[121];
    const double H01000 = Hr5[124];
    // const double Hm10100 = Hr5[129];
    const double H00100 = Hr5[130];
    const double H01100 = Hr5[133];
    const double H0m1010 = Hr5[145];
    const double H00010 = Hr5[148];
    const double H01010 = Hr5[151];
    const double H00110 = Hr5[157];
    const double H01110 = Hr5[160];
    const double H0m1m101 = Hr5[190];
    const double H00m101 = Hr5[193];
    const double H0m1001 = Hr5[199];
    // const double Hm10001 = Hr5[201];
    const double H00001 = Hr5[202];
    const double H01001 = Hr5[205];
    const double H00101 = Hr5[211];
    const double H01101 = Hr5[214];
    const double H0m1011 = Hr5[226];
    const double H00011 = Hr5[229];
    const double H01011 = Hr5[232];
    const double H00111 = Hr5[238];
    const double H01111 = Hr5[241];

    delete[] Hr1;
    delete[] Hr2;
    delete[] Hr3;
    delete[] Hr4;
    delete[] Hr5;

    // double flav = 5. / 18 * fl11ps(nf) * nf;
    // this contribution is neglected
    // If you want to use it please check the 5/48 since
    // I'm not sure about it

    return /*flav
               * (-38.400000000000006 + 34.13333333333333 * H001 - 512. * H01001
                  + 512. * H01100 - 384.00000000000006 * H0m10
                  - 699.7333333333333 * H1 + 238.93333333333334 * H100
                  + 1024. * H1100 - 426.66666666666663 * Hm10
                  - 401.06666666666666 * Hm100 - 512. * Hm10001
                  - 674.1333333333333 * Hm101 + 512. * Hm10100 + 128. * Hm1m10
                  + 938.6666666666666 * H001 * x - 1843.2 * H01 * x
                  - 597.3333333333333 * H0100 * x - 1024 * H01001 * x
                  + 1024 * H01100 * x - 42.666666666666664 * H0m10 * x
                  + 85.33333333333333 * H0m100 * x
                  + 426.66666666666663 * H0m101 * x + 256 * H0m1m10 * x
                  - 153.60000000000002 * H1 * x - 25.6 * H100 * x
                  - 128 * H1000 * x + 597.3333333333333 * H1001 * x
                  - 597.3333333333333 * H1100 * x
                  - 1066.6666666666665 * Hm10 * x
                  - 913.0666666666666 * Hm100 * x + 128 * Hm1000 * x
                  - 1024 * Hm10001 * x + 341.3333333333333 * Hm1001 * x
                  - 1399.4666666666667 * Hm101 * x + 1024 * Hm10100 * x
                  + 256 * Hm10m10 * x + 426.66666666666663 * Hm1m10 * x
                  - 341.3333333333333 * Hm1m100 * x
                  - 682.6666666666666 * Hm1m101 * x + 230.4 * x2
                  + 768 * H0001 * x2 - 537.6 * H001 * x2 + 768 * H01 * x2
                  - 768 * H0100 * x2 - 1152 * H0m10 * x2 + 768 * H1 * x2
                  - 307.20000000000005 * H100 * x2 + 768 * H1001 * x2
                  - 768 * H1100 * x2 - 768 * Hm10 * x2
                  - 307.20000000000005 * Hm100 * x2 + 537.6 * Hm101 * x2
                  + 1152 * Hm1m10 * x2 + 768 * H000 * x3 + 768 * H001 * x3
                  + 307.20000000000005 * H0100 * x3 - 768 * H0m10 * x3
                  - 307.20000000000005 * H0m100 * x3
                  - 614.4000000000001 * H0m101 * x3
                  - 307.20000000000005 * H1001 * x3
                  + 307.20000000000005 * H1100 * x3
                  - 998.4000000000001 * Hm10 * x3 - 768 * Hm100 * x3
                  - 307.20000000000005 * Hm1001 * x3 - 768 * Hm101 * x3
                  + 768 * Hm1m10 * x3 + 307.20000000000005 * Hm1m100 * x3
                  + 614.4000000000001 * Hm1m101 * x3
                  + (-85.33333333333333 * H01 + 85.33333333333333 * H1
                     - 34.13333333333333 * H100 + 85.33333333333333 * Hm10
                     + 34.13333333333333 * Hm100 + 68.26666666666667 * Hm101
                     + 503.4666666666667 * H01 * x
                     + (-345.6 + 42.666666666666664 * H000 + 512. * H0001) * x2)
                        / x
                  + (-34.13333333333333 * H1001 + 34.13333333333333 * H1100
                     + 110.93333333333335 * Hm10 + 85.33333333333333 * Hm100
                     + 34.13333333333333 * Hm1001 + 85.33333333333333 * Hm101
                     - 85.33333333333333 * Hm1m10 - 34.13333333333333 * Hm1m100
                     - 68.26666666666667 * Hm1m101 + 25.6 * x
                     - 1024. * H1001 * x2)
                        / x2
                  - (64 * H000) / xm + (128 * H0m10) / xm + (64 * H000) / xp
                  - (512 * H01) / xp
                  + H00
                        * (170.66666666666666 - 85.33333333333333 / x + 128 / xm
                           - 256 / xp
                           + x
                                 * (-392.5333333333333
                                    + x
                                          * (768. + 998.4000000000001 * x
                                             - 768. * zeta2)
                                    - 512. * zeta2))
                  - 503.4666666666667 * zeta2 + 512. * H010 * zeta2
                  + 64 * H1 * zeta2 + 1024 * H10 * zeta2
                  + 738.1333333333333 * Hm1 * zeta2 + 512 * Hm100 * zeta2
                  - (42.666666666666664 * H1 * zeta2) / x2
                  + (34.13333333333333 * H10 * zeta2) / x2
                  - (128 * Hm1 * zeta2) / x2
                  - (34.13333333333333 * Hm10 * zeta2) / x2
                  + (68.26666666666667 * Hm1m1 * zeta2) / x2
                  + (170.66666666666666 * zeta2) / x
                  - (68.26666666666667 * Hm1 * zeta2) / x
                  + 776.5333333333333 * x * zeta2 - 128 * H01 * x * zeta2
                  + 1024 * H010 * x * zeta2
                  - 298.66666666666663 * H0m1 * x * zeta2
                  - 213.33333333333331 * H1 * x * zeta2
                  - 469.3333333333333 * H10 * x * zeta2
                  + 1612.8000000000002 * Hm1 * x * zeta2
                  - 469.3333333333333 * Hm10 * x * zeta2
                  + 1024 * Hm100 * x * zeta2
                  + 682.6666666666666 * Hm1m1 * x * zeta2 - 768 * x2 * zeta2
                  + 576 * H1 * x2 * zeta2 - 768 * H10 * x2 * zeta2
                  + 38.400000000000006 * Hm1 * x2 * zeta2
                  - 998.4000000000001 * x3 * zeta2
                  + 614.4000000000001 * H0m1 * x3 * zeta2
                  - 384 * H1 * x3 * zeta2
                  + 307.20000000000005 * H10 * x3 * zeta2
                  + 1152 * Hm1 * x3 * zeta2
                  + 307.20000000000005 * Hm10 * x3 * zeta2
                  - 614.4000000000001 * Hm1m1 * x3 * zeta2 + (512 * zeta2) / xp
                  - 409.6 * Hm1 * zeta2_2 + 524.8000000000001 * x * zeta2_2
                  - 819.2 * Hm1 * x * zeta2_2 + 614.4000000000001 * x2 * zeta2_2
                  - 184.32 * x3 * zeta2_2 - 209.06666666666666 * zeta3
                  - 512 * H01 * zeta3 - 1024 * H1 * zeta3 - 512 * Hm10 * zeta3
                  - (17.066666666666666 * H1 * zeta3) / x2
                  - (51.2 * Hm1 * zeta3) / x2 + (68.26666666666667 * zeta3) / x
                  - 3447.4666666666667 * x * zeta3 - 1024 * H01 * x * zeta3
                  + 810.6666666666666 * H1 * x * zeta3 - 512 * Hm1 * x * zeta3
                  - 1024 * Hm10 * x * zeta3 + 268.8 * x2 * zeta3
                  + 768 * H1 * x2 * zeta3 - 1920 * x3 * zeta3
                  - 153.60000000000002 * H1 * x3 * zeta3
                  + 460.8 * Hm1 * x3 * zeta3 + (192 * zeta3) / xm
                  + 512 * x * zeta2 * zeta3
                  + H0
                        * (311.4666666666666 - 25.6 / x - 320. / xp
                           - 34.13333333333333 * zeta2 + (64. * zeta2) / xm
                           - (64. * zeta2) / xp
                           + x
                                 * (-76.80000000000001
                                    - 981.3333333333333 * zeta2
                                    + 682.6666666666666 * zeta3
                                    + x
                                          * (998.4000000000001
                                             + (537.6 - 1536. * x) * zeta2
                                             + 767.9999999999999 * zeta3
                                             - 614.4000000000001 * x * zeta3)))
                  + 5120 * x * zeta5)*/
           + CA * CF * nf
                 * (2009.7366255144034 - 867.2592592592592 * H000
                    + 379.55555555555554 * H0000 - 240 * H00000
                    - 106.66666666666666 * H00001 + 102.22222222222221 * H0001
                    + 8 * H00010 + 29.333333333333332 * H00011
                    - 165.33333333333331 * H000m10 - 754.4592592592593 * H001
                    + 28 * H0010 + 104 * H00100 + 69.33333333333333 * H00101
                    - 37.77777777777778 * H0011 + 122.66666666666666 * H00110
                    + 101.33333333333333 * H00111 + 301.3333333333333 * H00m10
                    - 106.66666666666666 * H00m100
                    - 58.666666666666664 * H00m101
                    - 58.666666666666664 * H00m1m10 + 812.6913580246913 * H01
                    - 309.77777777777777 * H010 + 22.666666666666664 * H0100
                    + 234.66666666666666 * H01000 + 216 * H01001
                    - 13.333333333333332 * H0101 + 98.66666666666666 * H01010
                    + 80 * H01011 + 10.666666666666666 * H010m10
                    - 437.037037037037 * H011 + 13.333333333333332 * H0110
                    + 56 * H01100 + 90.66666666666666 * H01101
                    + 0.8888888888888888 * H0111 + 90.66666666666666 * H01110
                    + 74.66666666666666 * H01111 - 568 * H0m10
                    + 378.66666666666663 * H0m100 - 106.66666666666666 * H0m1000
                    - 114.66666666666666 * H0m1001 + 168 * H0m101
                    - 58.666666666666664 * H0m1010 - 64 * H0m1011
                    - 101.33333333333333 * H0m10m10 - 184 * H0m1m10
                    - 114.66666666666666 * H0m1m100
                    - 37.33333333333333 * H0m1m101
                    + 117.33333333333333 * H0m1m1m10 - 1103.3333333333333 * H1
                    + 251.25925925925924 * H10 + 54.31111111111111 * H100
                    + 202.66666666666666 * H1000 - 71.11111111111111 * H101
                    + 52. * H1010 + 40. * H1011 - 5.333333333333333 * H10m10
                    + 24.74074074074074 * H11 + 54.22222222222222 * H110
                    - 73.33333333333331 * H1100 + 45.333333333333336 * H1101
                    - 13.333333333333332 * H111 + 45.333333333333336 * H1110
                    + 37.33333333333333 * H1111 + 163.52592592592592 * Hm10
                    - 328. * Hm100 + 48. * Hm1000 - 76. * Hm1001
                    - 84.44444444444444 * Hm101 + 8. * Hm1010
                    + 21.333333333333332 * Hm1011 - 109.33333333333333 * Hm10m10
                    - 143.1111111111111 * Hm1m10 - 76. * Hm1m100
                    + 146.66666666666663 * Hm1m101
                    + 122.66666666666666 * Hm1m1m10
                    + (1.4222222222222223 * H00) / x
                    + (4.2666666666666675 * H001 + 143.3283950617284 * H01
                       - 92.44444444444444 * H010 + 42.666666666666664 * H0100
                       + 53.33333333333332 * H0101 - 121.48148148148148 * H011
                       + 53.33333333333332 * H0110 + 56.88888888888889 * H0111
                       + 46.22222222222222 * H0m10 - 85.33333333333333 * H0m100
                       - 42.666666666666664 * H0m101
                       + 42.666666666666664 * H0m1m10 + 826.5810699588478 * H1
                       - 341.0864197530863 * H10 + 39.288888888888884 * H100
                       + 99.55555555555551 * H1000 + 78.22222222222221 * H1001
                       + 54.81481481481481 * H101 + 63.99999999999999 * H1010
                       + 53.33333333333332 * H1011 + 14.222222222222223 * H10m10
                       - 362.32098765432096 * H11 + 19.85185185185185 * H110
                       + 78.22222222222221 * H1100 + 60.44444444444443 * H1101
                       + 28.74074074074074 * H111 + 60.44444444444443 * H1110
                       + 49.77777777777776 * H1111 + 330.42962962962963 * Hm10
                       - 491.25925925925924 * Hm100
                       + 138.66666666666669 * Hm1000
                       + 106.66666666666664 * Hm1001
                       - 251.25925925925927 * Hm101 + 63.99999999999999 * Hm1010
                       + 78.22222222222221 * Hm1011
                       + 28.444444444444446 * Hm10m10
                       + 136.29629629629628 * Hm1m10
                       + 21.333333333333332 * Hm1m100
                       + 49.77777777777776 * Hm1m101
                       - 35.55555555555556 * Hm1m1m10)
                          / x
                    - 538.40329218107 * x + 579.4074074074074 * H000 * x
                    + 595.5555555555555 * H0000 * x + 240 * H00000 * x
                    + 266.66666666666663 * H00001 * x
                    + 286.22222222222223 * H0001 * x + 232 * H00010 * x
                    + 253.33333333333331 * H00011 * x
                    + 37.33333333333333 * H000m10 * x
                    - 100.05925925925926 * H001 * x
                    + 158.66666666666666 * H0010 * x
                    + 221.33333333333331 * H00100 * x + 176 * H00101 * x
                    + 119.55555555555554 * H0011 * x + 208 * H00110 * x
                    + 197.33333333333331 * H00111 * x
                    + 2.6666666666666665 * H00m10 * x
                    + 37.33333333333333 * H00m100 * x
                    + 90.66666666666666 * H00m101 * x
                    + 122.66666666666666 * H00m1m10 * x
                    + 357.4320987654321 * H01 * x
                    - 25.777777777777775 * H010 * x
                    + 241.33333333333331 * H0100 * x
                    + 234.66666666666666 * H01000 * x + 216 * H01001 * x
                    + 13.333333333333332 * H0101 * x
                    + 98.66666666666666 * H01010 * x + 80 * H01011 * x
                    + 10.666666666666666 * H010m10 * x
                    - 83.7037037037037 * H011 * x
                    + 18.666666666666664 * H0110 * x + 56 * H01100 * x
                    + 90.66666666666666 * H01101 * x
                    + 19.555555555555554 * H0111 * x
                    + 90.66666666666666 * H01110 * x
                    + 74.66666666666666 * H01111 * x
                    + 139.55555555555554 * H0m10 * x
                    - 165.33333333333331 * H0m100 * x
                    + 106.66666666666666 * H0m1000 * x
                    + 114.66666666666666 * H0m1001 * x
                    - 242.66666666666666 * H0m101 * x
                    + 58.666666666666664 * H0m1010 * x + 64 * H0m1011 * x
                    + 101.33333333333333 * H0m10m10 * x + 88 * H0m1m10 * x
                    + 114.66666666666666 * H0m1m100 * x
                    + 37.33333333333333 * H0m1m101 * x
                    - 117.33333333333333 * H0m1m1m10 * x
                    + 601.8518518518518 * H1 * x - 77.48148148148148 * H10 * x
                    + 132.35555555555555 * H100 * x
                    - 202.66666666666666 * H1000 * x - 268 * H1001 * x
                    + 183.11111111111111 * H101 * x - 52 * H1010 * x
                    - 40 * H1011 * x + 5.333333333333333 * H10m10 * x
                    + 185.48148148148147 * H11 * x
                    + 25.777777777777775 * H110 * x
                    + 94.66666666666666 * H1100 * x
                    - 45.33333333333333 * H1101 * x
                    + 109.33333333333333 * H111 * x
                    - 45.33333333333333 * H1110 * x
                    - 37.33333333333333 * H1111 * x
                    - 237.36296296296297 * Hm10 * x
                    - 86.22222222222221 * Hm100 * x + 48 * Hm1000 * x
                    - 76 * Hm1001 * x + 77.33333333333333 * Hm101 * x
                    + 8 * Hm1010 * x + 21.333333333333332 * Hm1011 * x
                    - 109.33333333333333 * Hm10m10 * x
                    - 144.88888888888889 * Hm1m10 * x - 76 * Hm1m100 * x
                    + 146.66666666666666 * Hm1m101 * x
                    + 122.66666666666666 * Hm1m1m10 * x - 140.0488340192044 * x2
                    - 930.6666666666666 * H000 * x2
                    + 14.222222222222221 * H0001 * x2
                    - 679.525925925926 * H001 * x2
                    - 85.33333333333333 * H0010 * x2
                    - 92.44444444444444 * H0011 * x2
                    + 184.88888888888889 * H00m10 * x2
                    + 305.50123456790124 * H01 * x2
                    - 431.1111111111111 * H010 * x2
                    - 174.2222222222222 * H0100 * x2
                    - 74.66666666666666 * H0101 * x2
                    - 390.8148148148148 * H011 * x2
                    - 110.22222222222221 * H0110 * x2
                    - 88.88888888888889 * H0111 * x2
                    - 198.5185185185185 * H0m10 * x2
                    + 199.1111111111111 * H0m100 * x2
                    + 106.66666666666666 * H0m101 * x2
                    - 21.333333333333332 * H0m1m10 * x2
                    - 325.09958847736624 * H1 * x2
                    + 167.30864197530863 * H10 * x2
                    - 225.95555555555558 * H100 * x2
                    - 99.55555555555554 * H1000 * x2
                    - 14.222222222222221 * H1001 * x2
                    - 166.8148148148148 * H101 * x2 - 64 * H1010 * x2
                    - 53.33333333333333 * H1011 * x2
                    - 14.222222222222221 * H10m10 * x2
                    + 152.09876543209876 * H11 * x2
                    - 99.85185185185185 * H110 * x2
                    - 142.22222222222223 * H1100 * x2
                    - 60.44444444444444 * H1101 * x2
                    - 124.74074074074073 * H111 * x2
                    - 60.44444444444444 * H1110 * x2
                    - 49.77777777777777 * H1111 * x2
                    - 66.90370370370371 * Hm10 * x2
                    - 235.25925925925924 * Hm100 * x2
                    + 138.66666666666666 * Hm1000 * x2
                    + 106.66666666666666 * Hm1001 * x2
                    - 75.25925925925925 * Hm101 * x2 + 64 * Hm1010 * x2
                    + 78.22222222222221 * Hm1011 * x2
                    + 28.444444444444443 * Hm10m10 * x2
                    + 120.29629629629629 * Hm1m10 * x2
                    + 21.333333333333332 * Hm1m100 * x2
                    + 49.77777777777777 * Hm1m101 * x2
                    - 35.55555555555556 * Hm1m1m10 * x2 - 12.8 * H000 * x3
                    - 38.400000000000006 * H0001 * x3 - 12.8 * H001 * x3
                    + 38.400000000000006 * H0100 * x3 + 12.8 * H0m10 * x3
                    - 38.400000000000006 * H1001 * x3
                    + 38.400000000000006 * H1100 * x3 + 3.2 * Hm10 * x3
                    + 12.8 * Hm100 * x3 + 12.8 * Hm101 * x3 - 12.8 * Hm1m10 * x3
                    + (4.266666666666667 * H1100 - 0.35555555555555557 * Hm10
                       - 1.4222222222222223 * Hm100 - 1.4222222222222223 * Hm101
                       + 1.4222222222222223 * Hm1m10 - 1331.284499314129 * x
                       + H1001 * (-4.266666666666667 + 246.66666666666666 * x2))
                          / x2
                    - 812.6913580246913 * zeta2
                    + 106.66666666666666 * H000 * zeta2
                    - 98.66666666666666 * H001 * zeta2
                    + 29.333333333333332 * H00m1 * zeta2
                    - 78.66666666666666 * H01 * zeta2
                    - 261.3333333333333 * H010 * zeta2
                    - 149.33333333333331 * H011 * zeta2 - 260 * H0m1 * zeta2
                    + 69.33333333333333 * H0m10 * zeta2 + 96 * H0m1m1 * zeta2
                    - 0.4444444444444444 * H1 * zeta2 - 304 * H10 * zeta2
                    - 106.66666666666666 * H11 * zeta2
                    + 12.888888888888888 * Hm1 * zeta2
                    + 18.666666666666664 * Hm10 * zeta2
                    - 85.33333333333333 * Hm1m1 * zeta2
                    + (0.7111111111111111 * H1 * zeta2) / x2
                    + (4.266666666666667 * H10 * zeta2) / x2
                    + (2.1333333333333333 * Hm1 * zeta2) / x2
                    + (187.10123456790123 * zeta2) / x
                    - (74.66666666666666 * H01 * zeta2) / x
                    + (64 * H0m1 * zeta2) / x
                    - (122.96296296296296 * H1 * zeta2) / x
                    - (85.33333333333333 * H10 * zeta2) / x
                    - (78.22222222222221 * H11 * zeta2) / x
                    + (319.4074074074074 * Hm1 * zeta2) / x
                    - (99.55555555555554 * Hm10 * zeta2) / x
                    - (67.55555555555556 * Hm1m1 * zeta2) / x
                    - 594.795061728395 * x * zeta2
                    - 229.33333333333331 * H000 * x * zeta2
                    - 237.33333333333331 * H001 * x * zeta2
                    - 29.333333333333332 * H00m1 * x * zeta2
                    - 57.33333333333333 * H01 * x * zeta2
                    - 261.3333333333333 * H010 * x * zeta2
                    - 149.33333333333331 * H011 * x * zeta2
                    + 286.66666666666663 * H0m1 * x * zeta2
                    - 69.33333333333333 * H0m10 * x * zeta2
                    - 96 * H0m1m1 * x * zeta2
                    - 110.66666666666666 * H1 * x * zeta2
                    + 325.3333333333333 * H10 * x * zeta2
                    + 106.66666666666666 * H11 * x * zeta2
                    - 149.77777777777777 * Hm1 * x * zeta2
                    + 18.666666666666664 * Hm10 * x * zeta2
                    - 85.33333333333333 * Hm1m1 * x * zeta2
                    - 305.50123456790124 * x2 * zeta2 + 64 * H01 * x2 * zeta2
                    - 117.33333333333333 * H0m1 * x2 * zeta2
                    + 226.96296296296296 * H1 * x2 * zeta2
                    + 21.333333333333332 * H10 * x2 * zeta2
                    + 78.22222222222221 * H11 * x2 * zeta2
                    + 135.4074074074074 * Hm1 * x2 * zeta2
                    - 99.55555555555554 * Hm10 * x2 * zeta2
                    - 67.55555555555556 * Hm1m1 * x2 * zeta2 + 3.2 * x3 * zeta2
                    + 6.4 * H1 * x3 * zeta2
                    + 38.400000000000006 * H10 * x3 * zeta2
                    - 19.200000000000003 * Hm1 * x3 * zeta2
                    + 71.46666666666667 * zeta2_2
                    + (100.97777777777779 * zeta2_2) / x
                    - 60.53333333333333 * x * zeta2_2
                    + 82.13333333333333 * x2 * zeta2_2 - 30.72 * x3 * zeta2_2
                    + H00
                          * (1979.2987654320987 - 102.22222222222221 * zeta2
                             + x
                                   * (1157.9654320987654
                                      - 283.55555555555554 * zeta2
                                      + x
                                            * (55.32839506172839
                                               - 14.222222222222221 * zeta2
                                               + x
                                                     * (-3.2
                                                        + 38.400000000000006
                                                              * zeta2))
                                      - 352.00000000000006 * zeta3)
                             - 341.3333333333333 * zeta3)
                    - 542.4592592592593 * zeta3 - 168 * H01 * zeta3
                    - 42.666666666666664 * H0m1 * zeta3 - 140 * H1 * zeta3
                    + 40 * Hm1 * zeta3 - (4.266666666666667 * H1 * zeta3) / x2
                    + (87.82222222222222 * zeta3) / x
                    - (69.33333333333333 * H1 * zeta3) / x
                    + (5.333333333333333 * Hm1 * zeta3) / x
                    + 652.8296296296296 * x * zeta3 - 168 * H01 * x * zeta3
                    + 42.666666666666664 * H0m1 * x * zeta3
                    + 118.66666666666666 * H1 * x * zeta3 + 40 * Hm1 * x * zeta3
                    + 426.2518518518519 * x2 * zeta3
                    + 133.33333333333331 * H1 * x2 * zeta3
                    + 5.333333333333333 * Hm1 * x2 * zeta3 + 32 * x3 * zeta3
                    - 38.400000000000006 * H1 * x3 * zeta3 + 296 * zeta2 * zeta3
                    + 301.3333333333333 * x * zeta2 * zeta3
                    + (H0
                       * (-163.56872427983538 + 41.955555555555556 * zeta2
                          - 14.222222222222221 * zeta3
                          + x
                                * (-1342.1201646090535
                                   + (754.4592592592593 - 56. * zeta2) * zeta2
                                   + x3
                                         * (25.6 * zeta2
                                            - 38.400000000000006 * zeta3)
                                   - 35.55555555555556 * zeta3
                                   + x2
                                         * (-986.9069958847737
                                            + 679.525925925926 * zeta2
                                            + 387.55555555555554 * zeta3)
                                   + x
                                         * (697.4106995884773
                                            + zeta2
                                                  * (239.61481481481485
                                                     + 122.13333333333333
                                                           * zeta2)
                                            - 768.8888888888888 * zeta3
                                            - 48. * zeta4)
                                   - 48. * zeta4)))
                          / x
                    - 24. * zeta4 - (32 * zeta4) / x + 24 * x * zeta4
                    + 32 * x2 * zeta4 - 620. * zeta5 - 132 * x * zeta5)
           + CF * CF * nf
                 * (92.3446913580247 + 575.5555555555555 * H000
                    - 146.66666666666666 * H0000 + 240. * H00000
                    + 442.66666666666663 * H00001 + 170.66666666666666 * H0001
                    + 357.3333333333333 * H00010 + 400. * H00011
                    - 384. * H000m10 + 1067.288888888889 * H001
                    + 194.66666666666666 * H0010 + 181.3333333333334 * H00100
                    + 287.99999999999994 * H00101 + 293.3333333333333 * H0011
                    + 245.33333333333331 * H00110 + 234.66666666666663 * H00111
                    - 352.00000000000006 * H00m10 - 575.9999999999999 * H00m100
                    - 256. * H00m101 + 384. * H00m1m10 + 914.6370370370369 * H01
                    + 695.1111111111111 * H010 + 149.33333333333331 * H0100
                    - 165.33333333333331 * H01000 - 85.33333333333333 * H01001
                    + 194.66666666666666 * H0101 + 149.33333333333331 * H01010
                    + 143.99999999999997 * H01011 + 846.6666666666667 * H011
                    + 162.66666666666666 * H0110 + 234.66666666666663 * H01100
                    + 133.33333333333331 * H01101 + 192. * H0111
                    + 133.33333333333331 * H01110 + 117.33333333333331 * H01111
                    - 736.0000000000002 * H0m100 - 256. * H0m1000
                    - 64. * H0m1001 - 416. * H0m101 + 256. * H0m10m10
                    + 416. * H0m1m10 + 448. * H0m1m100 + 128. * H0m1m101
                    - 256. * H0m1m1m10 + 1466.5086419753088 * H1
                    + 698.8148148148148 * H10 + 373.15555555555557 * H100
                    - 210.66666666666663 * H1000 - 266.66666666666663 * H1001
                    + 571.1111111111111 * H101 + 74.66666666666666 * H1010
                    + 71.99999999999999 * H1011 + 964.8888888888887 * H11
                    + 432.44444444444434 * H110 + 277.3333333333333 * H1100
                    + 66.66666666666666 * H1101 + 433.77777777777766 * H111
                    + 66.66666666666666 * H1110 + 58.66666666666666 * H1111
                    - 1049.2444444444445 * Hm10 - 1797.333333333333 * Hm100
                    - 448. * Hm1000 - 192. * Hm1001 - 976. * Hm101
                    - 64. * Hm1010 - 64. * Hm1011 + 320. * Hm10m10
                    + 1104. * Hm1m10 + 575.9999999999999 * Hm1m100
                    + 64. * Hm1m101 - 320. * Hm1m1m10
                    + (2.1333333333333333 * H00) / x
                    + (-8.533333333333333 * H001 + 7.822222222222221 * H01
                       - 127.50123456790124 * H1 + 70.41975308641975 * H10
                       + 30.459259259259266 * H100 - 24.888888888888882 * H1000
                       + 49.777777777777764 * H1001 - 43.25925925925927 * H101
                       + 99.55555555555553 * H1010 + 96. * H1011
                       + 68.09876543209876 * H11 - 11.25925925925926 * H110
                       + 92.44444444444446 * H1100 + 88.88888888888886 * H1101
                       + 7.111111111111106 * H111 + 88.88888888888886 * H1110
                       + 78.22222222222221 * H1111 - 15.288888888888879 * Hm10
                       - 64. * Hm100 + 21.333333333333332 * Hm1000
                       - 21.333333333333332 * Hm1001 - 17.77777777777778 * Hm101
                       - 21.333333333333332 * Hm1010
                       - 21.333333333333332 * Hm1011 - 64. * Hm10m10
                       + 110.22222222222229 * Hm1m10
                       - 106.66666666666659 * Hm1m100 - 64. * Hm1m101
                       + 64. * Hm1m1m10)
                          / x
                    - 537.9002469135803 * x + 1070.6666666666665 * H000 * x
                    + 293.3333333333333 * H0000 * x + 240 * H00000 * x
                    + 442.66666666666663 * H00001 * x + 720 * H0001 * x
                    + 357.3333333333333 * H00010 * x + 400 * H00011 * x
                    + 912.2666666666667 * H001 * x + 216 * H0010 * x
                    + 181.33333333333331 * H00100 * x + 288 * H00101 * x
                    + 320 * H0011 * x + 245.33333333333331 * H00110 * x
                    + 234.66666666666666 * H00111 * x
                    - 53.33333333333333 * H00m10 * x + 64 * H00m100 * x
                    - 128 * H00m101 * x - 256 * H00m1m10 * x
                    - 690.1037037037038 * H01 * x - 348.4444444444444 * H010 * x
                    - 304 * H0100 * x - 165.33333333333331 * H01000 * x
                    - 85.33333333333333 * H01001 * x
                    + 45.33333333333333 * H0101 * x
                    + 149.33333333333331 * H01010 * x + 144 * H01011 * x
                    - 245.33333333333331 * H011 * x
                    + 13.333333333333332 * H0110 * x
                    + 234.66666666666666 * H01100 * x
                    + 133.33333333333331 * H01101 * x
                    + 53.33333333333333 * H0111 * x
                    + 133.33333333333331 * H01110 * x
                    + 117.33333333333333 * H01111 * x
                    - 748.4444444444445 * H0m10 * x - 32 * H0m100 * x
                    + 256 * H0m1000 * x + 64 * H0m1001 * x
                    + 202.66666666666666 * H0m101 * x - 256 * H0m10m10 * x
                    + 53.33333333333333 * H0m1m10 * x - 448 * H0m1m100 * x
                    - 128 * H0m1m101 * x + 256 * H0m1m1m10 * x
                    - 1308.5086419753086 * H1 * x - 905.037037037037 * H10 * x
                    - 506.4888888888889 * H100 * x
                    + 210.66666666666666 * H1000 * x
                    + 309.3333333333333 * H1001 * x
                    - 587.1111111111111 * H101 * x
                    - 74.66666666666666 * H1010 * x - 72 * H1011 * x
                    - 1084.4444444444443 * H11 * x
                    - 416.4444444444444 * H110 * x - 320 * H1100 * x
                    - 66.66666666666666 * H1101 * x
                    - 433.77777777777777 * H111 * x
                    - 66.66666666666666 * H1110 * x
                    - 58.666666666666664 * H1111 * x
                    - 995.3185185185185 * Hm10 * x
                    - 1850.6666666666665 * Hm100 * x - 448 * Hm1000 * x
                    - 192 * Hm1001 * x - 1032.888888888889 * Hm101 * x
                    - 64 * Hm1010 * x - 64 * Hm1011 * x + 320 * Hm10m10 * x
                    + 1096.888888888889 * Hm1m10 * x + 576 * Hm1m100 * x
                    + 64 * Hm1m101 * x - 320 * Hm1m1m10 * x
                    + 409.7837037037037 * x2 - 216.88888888888889 * H000 * x2
                    - 423.1111111111111 * H0000 * x2
                    - 647.1111111111111 * H0001 * x2
                    - 293.6888888888889 * H001 * x2
                    - 369.77777777777777 * H0010 * x2
                    - 387.55555555555554 * H0011 * x2
                    + 106.66666666666666 * H00m10 * x2
                    - 143.98024691358023 * H01 * x2
                    - 137.48148148148147 * H010 * x2 - 64 * H0100 * x2
                    - 266.66666666666663 * H0101 * x2
                    - 158.8148148148148 * H011 * x2 - 224 * H0110 * x2
                    - 206.2222222222222 * H0111 * x2
                    - 110.22222222222221 * H0m10 * x2
                    + 106.66666666666666 * H0m100 * x2
                    + 21.333333333333332 * H0m101 * x2
                    - 106.66666666666666 * H0m1m10 * x2
                    - 30.498765432098764 * H1 * x2
                    + 135.80246913580245 * H10 * x2
                    + 102.87407407407407 * H100 * x2
                    + 24.888888888888886 * H1000 * x2
                    - 177.77777777777777 * H1001 * x2
                    + 59.25925925925925 * H101 * x2
                    - 99.55555555555554 * H1010 * x2 - 96 * H1011 * x2
                    + 51.456790123456784 * H11 * x2
                    - 4.7407407407407405 * H110 * x2
                    + 35.55555555555556 * H1100 * x2
                    - 88.88888888888889 * H1101 * x2
                    - 7.111111111111111 * H111 * x2
                    - 88.88888888888889 * H1110 * x2
                    - 78.22222222222221 * H1111 * x2
                    + 46.93333333333333 * Hm10 * x2 - 96 * Hm100 * x2
                    + 21.333333333333332 * Hm1000 * x2
                    - 21.333333333333332 * Hm1001 * x2
                    - 81.77777777777777 * Hm101 * x2
                    - 21.333333333333332 * Hm1010 * x2
                    - 21.333333333333332 * Hm1011 * x2 - 64 * Hm10m10 * x2
                    + 110.22222222222221 * Hm1m10 * x2
                    - 106.66666666666666 * Hm1m100 * x2 - 64 * Hm1m101 * x2
                    + 64 * Hm1m1m10 * x2 - 19.200000000000003 * H000 * x3
                    + 76.80000000000001 * H0001 * x3 + 6.4 * H001 * x3
                    - 76.80000000000001 * H0100 * x3 - 6.4 * H0m10 * x3
                    + 76.80000000000001 * H1001 * x3
                    - 76.80000000000001 * H1100 * x3
                    + 7.253333333333334 * Hm10 * x3
                    + 19.200000000000003 * Hm100 * x3 - 6.4 * Hm101 * x3
                    + 6.4 * Hm1m10 * x3
                    + (8.533333333333333 * H1001 - 8.533333333333333 * H1100
                       - 1.0429629629629629 * Hm10 - 2.1333333333333333 * Hm100
                       + 0.7111111111111111 * Hm101
                       - 0.7111111111111111 * Hm1m10 + 35.77185185185185 * x
                       + H0m10 * (-2.8444444444444446 - 912.0000000000002 * x2))
                          / x2
                    - 914.6370370370369 * zeta2
                    - 442.66666666666663 * H000 * zeta2 - 96 * H001 * zeta2
                    + 448 * H00m1 * zeta2 + 13.333333333333332 * H01 * zeta2
                    + 213.33333333333331 * H010 * zeta2
                    - 5.333333333333333 * H011 * zeta2 + 624 * H0m1 * zeta2
                    + 192 * H0m10 * zeta2 - 256 * H0m1m1 * zeta2
                    - 19.11111111111111 * H1 * zeta2
                    + 426.66666666666663 * H10 * zeta2
                    + 93.33333333333333 * H11 * zeta2 + 1528 * Hm1 * zeta2
                    + 352 * Hm10 * zeta2 - 224 * Hm1m1 * zeta2
                    - (0.35555555555555557 * H1 * zeta2) / x2
                    - (8.533333333333333 * H10 * zeta2) / x2
                    - (1.0666666666666667 * Hm1 * zeta2) / x2
                    - (23.11111111111111 * zeta2) / x
                    - (11.851851851851851 * H1 * zeta2) / x
                    - (17.77777777777778 * H10 * zeta2) / x
                    - (56.888888888888886 * H11 * zeta2) / x
                    + (72.88888888888889 * Hm1 * zeta2) / x
                    - (10.666666666666666 * Hm10 * zeta2) / x
                    + (96 * Hm1m1 * zeta2) / x - 305.2148148148148 * x * zeta2
                    - 442.66666666666663 * H000 * x * zeta2
                    - 160 * H001 * x * zeta2 - 72 * H01 * x * zeta2
                    + 213.33333333333331 * H010 * x * zeta2
                    - 5.333333333333333 * H011 * x * zeta2
                    - 176 * H0m1 * x * zeta2 - 192 * H0m10 * x * zeta2
                    + 256 * H0m1m1 * x * zeta2
                    + 38.666666666666664 * H1 * x * zeta2
                    - 469.3333333333333 * H10 * x * zeta2
                    - 93.33333333333333 * H11 * x * zeta2
                    + 1581.3333333333333 * Hm1 * x * zeta2
                    + 352 * Hm10 * x * zeta2 - 224 * Hm1m1 * x * zeta2
                    + 143.98024691358023 * x2 * zeta2
                    + 213.33333333333331 * H01 * x2 * zeta2
                    - 74.66666666666666 * H0m1 * x2 * zeta2
                    - 4.148148148148148 * H1 * x2 * zeta2
                    + 145.77777777777777 * H10 * x2 * zeta2
                    + 56.888888888888886 * H11 * x2 * zeta2
                    + 136.88888888888889 * Hm1 * x2 * zeta2
                    - 10.666666666666666 * Hm10 * x2 * zeta2
                    + 96 * Hm1m1 * x2 * zeta2 + 7.253333333333334 * x3 * zeta2
                    - 3.2 * H1 * x3 * zeta2
                    - 76.80000000000001 * H10 * x3 * zeta2
                    + 9.600000000000001 * Hm1 * x3 * zeta2 + 64 * zeta2_2
                    + (25.6 * zeta2_2) / x + 276.8 * x * zeta2_2
                    - 209.77777777777777 * x2 * zeta2_2 + 61.44 * x3 * zeta2_2
                    + H00
                          * (-39.48148148148148 - 170.66666666666666 * zeta2
                             + x
                                   * (63.33333333333333
                                      - 773.3333333333331 * zeta2
                                      + x
                                            * (10.76543209876543
                                               + x
                                                     * (-7.253333333333334
                                                        - 76.80000000000001
                                                              * zeta2)
                                               + 647.1111111111111 * zeta2)
                                      - 85.33333333333333 * zeta3)
                             - 213.33333333333331 * zeta3)
                    - 394.4888888888889 * zeta3
                    + 165.33333333333331 * H01 * zeta3 + 224 * H0m1 * zeta3
                    + 226.66666666666666 * H1 * zeta3 + 256 * Hm1 * zeta3
                    + (8.533333333333333 * H1 * zeta3) / x2
                    - (97.42222222222223 * zeta3) / x
                    + (14.222222222222221 * H1 * zeta3) / x
                    - (64 * Hm1 * zeta3) / x - 2453.0666666666666 * x * zeta3
                    + 165.33333333333331 * H01 * x * zeta3
                    - 224 * H0m1 * x * zeta3 - 184 * H1 * x * zeta3
                    + 256 * Hm1 * x * zeta3 - 198.28148148148148 * x2 * zeta3
                    - 142.22222222222223 * H1 * x2 * zeta3
                    - 64 * Hm1 * x2 * zeta3 - 16 * x3 * zeta3
                    + 76.80000000000001 * H1 * x3 * zeta3 + 144 * zeta2 * zeta3
                    + 16 * x * zeta2 * zeta3 + 24. * zeta4 + (32 * zeta4) / x
                    - 24 * x * zeta4 - 32 * x2 * zeta4
                    + (H0
                       * (1.754074074074074 + 8.533333333333333 * zeta2
                          + x
                                * (254.14024691358023
                                   + zeta2
                                         * (-1067.288888888889
                                            + 77.33333333333333 * zeta2)
                                   - 26.666666666666664 * zeta3
                                   + x3
                                         * (-12.8 * zeta2
                                            + 76.80000000000001 * zeta3)
                                   + x2
                                         * (149.80345679012345
                                            + 293.6888888888889 * zeta2
                                            + 120.88888888888887 * zeta3)
                                   + 48. * zeta4
                                   + x
                                         * (-1610.9313580246915
                                            + zeta2
                                                  * (-1660.711111111111
                                                     + 109.33333333333333
                                                           * zeta2)
                                            + 282.66666666666663 * zeta3
                                            + 48. * zeta4))))
                          / x
                    + 176.00000000000003 * zeta5 + 176 * x * zeta5)
           + (CF * nf * nf
              * ((29.620850480109738 + 0.7111111111111111 * H0
                  + 1.382716049382716 * H1 - 5.925925925925926 * H10
                  + 3.5555555555555554 * H101 - 11.65432098765432 * H11
                  + 3.5555555555555554 * H110 + 4.7407407407407405 * H111)
                     * x
                 + Hm10
                       * (-0.7111111111111111
                          + x
                                * (-3.5555555555555554
                                   + x
                                         * (56.888888888888886
                                            + x
                                                  * (71.11111111111111
                                                     + x
                                                           * (17.77777777777778
                                                              + 6.4 * x)))))
                 + x
                       * (14.22222222222222 * Hm100 * (1. + x) * (1. + x)
                              * (1. + x)
                          - 3.5555555555555554 * zeta2
                          - 3.5555555555555554 * H1 * zeta2
                          + x4 * (-6.4 * H00 + 6.4 * zeta2)
                          + 4.7407407407407405 * zeta3
                          + x2
                                * (185.07325102880657 - 223.20987654320987 * H00
                                   - 162.07407407407408 * H000
                                   - 40.888888888888886 * H0000
                                   - 3.5555555555555554 * H0001
                                   - 36.74074074074074 * H001
                                   + 12.444444444444443 * H0011
                                   - 42.172839506172835 * H01
                                   - 19.555555555555554 * H010
                                   + 5.333333333333333 * H0101
                                   - 10.962962962962962 * H011
                                   + 5.333333333333333 * H0110
                                   + 7.111111111111111 * H0111
                                   + 26.17283950617284 * H1
                                   + 2.6666666666666665 * H10
                                   - 2.6666666666666665 * H101
                                   - 22.22222222222222 * H11
                                   - 2.6666666666666665 * H110
                                   - 3.5555555555555554 * H111
                                   + 113.28395061728395 * zeta2
                                   + zeta2
                                         * (3.5555555555555554 * H00
                                            - 5.333333333333333 * H01
                                            + 2.6666666666666665 * H1
                                            + 3.7333333333333334 * zeta2)
                                   + H0
                                         * (-42.26502057613169
                                            + 36.74074074074074 * zeta2
                                            - 1.7777777777777777 * zeta3)
                                   + 5.037037037037036 * zeta3)
                          + x3
                                * (42.082853223593965 + 64.5925925925926 * H00
                                   + 28.444444444444443 * H000
                                   - 3.5555555555555554 * H001
                                   + 19.555555555555554 * H01
                                   - 3.5555555555555554 * H010
                                   - 14.222222222222221 * H011
                                   + 14.222222222222221 * H0m10
                                   + 25.28395061728395 * H1
                                   + 11.25925925925926 * H10
                                   - 3.5555555555555554 * H101
                                   + 11.65432098765432 * H11
                                   - 3.5555555555555554 * H110
                                   - 4.7407407407407405 * H111
                                   - 19.555555555555554 * zeta2
                                   + 3.5555555555555554 * H1 * zeta2
                                   + H0
                                         * (115.279012345679
                                            + 3.5555555555555554 * zeta2)
                                   + 5.925925925925925 * zeta3)
                          + x
                                * (-256.77695473251026 - 155.6543209876543 * H00
                                   - 66.07407407407408 * H000
                                   - 40.888888888888886 * H0000
                                   - 3.5555555555555554 * H0001
                                   - 31.40740740740741 * H001
                                   + 12.444444444444443 * H0011
                                   - 22.61728395061728 * H01
                                   - 14.222222222222221 * H010
                                   + 5.333333333333333 * H0101
                                   - 5.62962962962963 * H011
                                   + 5.333333333333333 * H0110
                                   + 7.111111111111111 * H0111
                                   + 42.666666666666664 * H0m10
                                   - 52.839506172839506 * H1 - 8. * H10
                                   + 2.6666666666666665 * H101
                                   + 22.22222222222222 * H11
                                   + 2.6666666666666665 * H110
                                   + 3.5555555555555554 * H111
                                   + 22.61728395061728 * zeta2
                                   + zeta2
                                         * (3.5555555555555554 * H00
                                            - 5.333333333333333 * H01
                                            - 2.6666666666666665 * H1
                                            + 3.7333333333333334 * zeta2)
                                   + H0
                                         * (-180.93168724279835
                                            + 31.40740740740741 * zeta2
                                            - 1.7777777777777777 * zeta3)
                                   + 12.148148148148145 * zeta3))))
                 / x2;
}

//==========================================================================================//
//  Massless gluon coefficient functions for FL at
//  O(as^3) for mu=Q.
//  The term fl_ps_11 is put to zero for the reason
//  explained in page 15 of arXiv:1205.5727
//
//  Eq. (B.17) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double MasslessCoefficientFunction::CL_g3_massless(double x, int nf) const {

    // the commented varibles are needed for the flav term

    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;

    // double xm = 1 - x;
    // double xm2 = xm * xm;
    // double xm3 = xm2 * xm;

    // double xp = 1 + x;
    // double xp2 = xp * xp;
    // double xp3 = xp2 * xp;

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
    const double H0m1m1m10 = Hr5[82];
    const double H00m1m10 = Hr5[85];
    const double H0m10m10 = Hr5[91];
    const double H0m1m100 = Hr5[109];
    const double H00m100 = Hr5[112];
    const double H0m1000 = Hr5[118];
    const double H01000 = Hr5[124];
    const double Hm10100 = Hr5[129];
    const double H00100 = Hr5[130];
    // const double H10100 = Hr5[131];
    const double H01100 = Hr5[133];
    // const double H11100 = Hr5[134];
    const double H0m1m101 = Hr5[190];
    const double H00m101 = Hr5[193];
    const double H0m1001 = Hr5[199];
    const double Hm10001 = Hr5[201];
    // const double H10001 = Hr5[203];
    const double H01001 = Hr5[205];
    // const double H11001 = Hr5[206];

    delete[] Hr1;
    delete[] Hr2;
    delete[] Hr3;
    delete[] Hr4;
    delete[] Hr5;

    // double flav = 5. / 48 * fl11g(nf) * nf * nf;
    // this contribution is neglected
    // If you want to use it please check the 5/48 since
    // I'm not sure about it
    return /*flav*(-76.2311111111111 + 358.40000000000003*H001 -
     161.56444444444443*H01 - 17.066666666666666*H0m10 + 117.76*H1 -
     375.46666666666664*H100 + 768*H1100 - 959.1466666666668*Hm10 -
     426.66666666666663*Hm100 - 409.6*Hm101 + 443.73333333333335*Hm1m10 +
     (238.93333333333334*H001)/x - (183.1822222222222*H01)/x +
     (183.1822222222222*H1)/x - (238.93333333333334*H100)/x - (256*H1001)/x +
     (256*H1100)/x - (55.75111111111111*Hm10)/x - (182.04444444444445*Hm101)/x -
     (182.04444444444445*Hm1m10)/x - 814.08*x + 147.91111111111113*H000*x +
     136.53333333333333*H0000*x + 221.86666666666667*H0001*x +
     967.1111111111111*H001*x - 136.53333333333333*H00m10*x -
     1162.808888888889*H01*x - 597.3333333333333*H0100*x + 1536*H01001*x -
     1536*H01100*x - 147.91111111111113*H0m10*x + 375.46666666666664*H0m100*x +
     887.4666666666667*H0m101*x + 136.53333333333333*H0m1m10*x -
     2179.9822222222224*H1*x - 1408*H10*x - 307.20000000000005*H100*x +
     512*H1000*x + 1536*H10001*x + 2645.333333333333*H1001*x - 1536*H10100*x -
     512*H10m10*x - 2816*H11*x - 2133.333333333333*H1100*x + 3072*H11001*x -
     3072*H11100*x - 1408.568888888889*Hm10*x - 659.9111111111112*Hm100*x -
     512*Hm1000*x + 1132.088888888889*Hm101*x + 512*Hm10m10*x +
     2451.911111111111*Hm1m10*x + 512*Hm1m100*x - 1024*Hm1m1m10*x + 860.16*x2 -
     1723.7333333333333*H0001*x2 - 2514.488888888889*H001*x2 +
     1092.2666666666667*H00m10*x2 + 1716.9066666666668*H01*x2 + 2304*H0100*x2 -
     1536*H01001*x2 + 1536*H01100*x2 - 2104.8888888888887*H0m10*x2 -
     512*H0m100*x2 - 1092.2666666666667*H0m101*x2 - 68.26666666666667*H0m1m10*x2
     + 1879.04*H1*x2 + 1408*H10*x2 + 921.6*H100*x2 - 512*H1000*x2 -
     1536*H10001*x2 - 2816*H1001*x2 + 1536*H10100*x2 + 512*H10m10*x2 +
     2816*H11*x2 + 2304*H1100*x2 - 3072*H11001*x2 + 3072*H11100*x2 -
     633.1733333333334*Hm10*x2 - 512*Hm100*x2 - 512*Hm1000*x2 +
     1080.888888888889*Hm101*x2 + 512*Hm10m10*x2 + 2104.8888888888887*Hm1m10*x2
     + 512*Hm1m100*x2 - 1024*Hm1m1m10*x2 + 334.50666666666666*H000*x3 +
     1433.6000000000001*H0001*x3 + 334.50666666666666*H001*x3 -
     1433.6000000000001*H0100*x3 - 334.50666666666666*H0m10*x3 +
     1433.6000000000001*H1001*x3 - 1433.6000000000001*H1100*x3 -
     153.60000000000002*Hm10*x3 - 334.50666666666666*Hm100*x3 -
     334.50666666666666*Hm101*x3 + 334.50666666666666*Hm1m10*x3 -
     (768.*(-0.3111111111111111*H1100 + 0.03333333333333333*Hm10 +
     0.07259259259259258*Hm100 + 0.07259259259259258*Hm101 -
     0.07259259259259258*Hm1m10 - 0.03925925925925926*x +
     H1001*(0.3111111111111111 + x2)))/x2 + (17.066666666666666*H0000)/xm3 +
     (34.13333333333333*H0001)/xm3 - (34.13333333333333*H0m100)/xm3 -
     (68.26666666666667*H0m101)/xm3 + (17.066666666666666*H000)/xm2 -
     (17.066666666666666*H0000)/xm2 - (34.13333333333333*H0001)/xm2 +
     (34.13333333333333*H001)/xm2 + (34.13333333333333*H0m100)/xm2 +
     (68.26666666666667*H0m101)/xm2 - (34.13333333333333*Hm100)/xm2 -
     (68.26666666666667*Hm101)/xm2 - (8.533333333333333*H000)/xm -
     (17.066666666666666*H001)/xm + (17.066666666666666*Hm100)/xm +
     (34.13333333333333*Hm101)/xm - (17.066666666666666*H0000)/xp3 +
     (34.13333333333333*H00m10)/xp3 - (17.066666666666666*H000)/xp2 +
     (17.066666666666666*H0000)/xp2 - (34.13333333333333*H00m10)/xp2 +
     (34.13333333333333*H0m10)/xp2 + (8.533333333333333*H000)/xp -
     (17.066666666666666*H0m10)/xp + (17.066666666666666*Hm10)/xp +
     161.56444444444443*zeta2 + 221.86666666666667*H1*zeta2 + 768*H10*zeta2 +
     631.4666666666667*Hm1*zeta2 + (27.875555555555554*H1*zeta2)/x2 +
     (238.93333333333334*H10*zeta2)/x2 + (83.62666666666667*Hm1*zeta2)/x2 +
     (127.43111111111111*zeta2)/x + (91.02222222222223*H1*zeta2)/x +
     (256*H10*zeta2)/x + (91.02222222222223*Hm1*zeta2)/x - 245.76*x*zeta2
     - 68.26666666666667*H01*x*zeta2 - 1536*H010*x*zeta2 - 819.2*H0m1*x*zeta2 -
     1225.9555555555555*H1*x*zeta2 - 3157.333333333333*H10*x*zeta2 -
     1536*H100*x*zeta2 - 512*H11*x*zeta2 - 3072*H110*x*zeta2
     + 93.86666666666666*Hm1*x*zeta2 + 512*Hm10*x*zeta2 - 512*Hm1m1*x*zeta2 -
     1716.9066666666668*x2*zeta2 - 34.13333333333333*H01*x2*zeta2 +
     1536*H010*x2*zeta2 + 1058.1333333333332*H0m1*x2*zeta2 +
     1052.4444444444443*H1*x2*zeta2 + 3328*H10*x2*zeta2 + 1536*H100*x2*zeta2 +
     512*H11*x2*zeta2 + 3072*H110*x2*zeta2 - 28.444444444444443*Hm1*x2*zeta2 +
     512*Hm10*x2*zeta2 - 512*Hm1m1*x2*zeta2 - 153.60000000000002*x3*zeta2 -
     167.25333333333333*H1*x3*zeta2 - 1433.6000000000001*H10*x3*zeta2 +
     501.76*Hm1*x3*zeta2 + (68.26666666666667*H0m1*zeta2)/xm3 -
     (68.26666666666667*H0m1*zeta2)/xm2 + (68.26666666666667*Hm1*zeta2)/xm2 -
     (8.533333333333333*zeta2)/xm - (34.13333333333333*Hm1*zeta2)/xm +
     (8.533333333333333*zeta2)/xp + 493.2266666666667*x*zeta2_2 +
     1228.8000000000002*H1*x*zeta2_2 - 812.3733333333333*x2*zeta2_2 -
     1228.8000000000002*H1*x2*zeta2_2 + 1146.88*x3*zeta2_2 -
     (6.826666666666667*zeta2_2)/xm3 + (6.826666666666667*zeta2_2)/xm2 +
     (35.84*zeta2_2)/xp3 - (35.84*zeta2_2)/xp2 + H00*(100.12444444444445
     + 55.75111111111111/x - 8.533333333333333/xp + x3*(153.60000000000002 -
     1433.6000000000001*zeta2) + x*(1104.2133333333334 -
     358.40000000000003*zeta2) + ((-51.2 + 51.2*xm)*zeta2)/xm3 +
     ((17.066666666666666 - 17.066666666666666*xp)*zeta2)/xp3 +
        x2*(607.5733333333334 + 1723.7333333333333*zeta2)) + 384*zeta3 -
     768*H1*zeta3 - (238.93333333333334*H1*zeta3)/x2 +
     (238.93333333333334*zeta3)/x - (256*H1*zeta3)/x -
     190.57777777777778*x*zeta3 + 1536*H01*x*zeta3 +
     597.3333333333333*H1*x*zeta3 + 1536*H10*x*zeta3 + 3072*H11*x*zeta3 +
     512*Hm1*x*zeta3 - 1277.1555555555556*x2*zeta3 - 1536*H01*x2*zeta3 -
     768*H1*x2*zeta3 - 1536*H10*x2*zeta3 - 3072*H11*x2*zeta3 + 512*Hm1*x2*zeta3
     - 836.2666666666667*x3*zeta3 + 1433.6000000000001*H1*x3*zeta3 -
     (51.2*zeta3)/xm2 + (25.6*zeta3)/xm + (34.13333333333333*zeta3)/xp2 -
     (17.066666666666666*zeta3)/xp + H0*((-30.15111111111111 -
     238.93333333333334*zeta2)/x + x2*(1698.1333333333332 +
     2514.488888888889*zeta2 - 1826.1333333333332*zeta3) +
        (xp*(xp*(-8.533333333333333 - 8.533333333333333*zeta2)
     + 17.066666666666666*zeta2 - 34.13333333333333*zeta3)
     + 34.13333333333333*zeta3)/xp3 + x*(-657.6355555555556 -
     1115.0222222222224*zeta2 + 768*zeta3) + x3*(-669.0133333333333*zeta2 +
     1433.6000000000001*zeta3) +
        (-51.2*zeta3 + xm*(xm2*(-22.186666666666667 - 358.40000000000003*zeta2)
     - 51.2*zeta2 + 25.6*xm*zeta2 + 51.2*zeta3))/xm3))*/
        +CA * CF * nf
            * (-204.1659259259259 + 92.80000000000001 * H000
               + 82.13333333333333 * H001 - 128. * H00m10
               + 205.79555555555552 * H01 + 4.266666666666667 * H010
               + 4.266666666666667 * H011 - 256. * H0m100 - 128. * H0m101
               + 128. * H0m1m10 - 128.4088888888889 * H1 - 172.8 * H10
               - 226.1333333333333 * H100 - 128. * H1000 - 480. * H1001
               - 16. * H101 - 172.8 * H11 - 47.99999999999999 * H110
               + 415.99999999999994 * H1100 - 16. * H111
               - 1725.0844444444444 * Hm10 - 992.0000000000001 * Hm100
               - 128. * Hm1000 - 32. * Hm1001 - 565.3333333333331 * Hm101
               + 128. * Hm10m10 + 910.9333333333333 * Hm1m10
               + 223.99999999999997 * Hm1m100 + 64. * Hm1m101 - 128. * Hm1m1m10
               - (55.324444444444445 * H00) / x - (51.20000000000001 * H000) / x
               + (4.266666666666667 * H001) / x - (70.9688888888889 * H01) / x
               - (8.533333333333333 * H010) / x - (8.533333333333333 * H011) / x
               + (34.13333333333333 * H0m10) / x + (105.81333333333333 * H1) / x
               + (40.53333333333333 * H10) / x - (46.93333333333333 * H100) / x
               + (128. * H1001) / x + (40.53333333333333 * H11) / x
               - (128. * H1100) / x + (21.90222222222222 * Hm10) / x
               + (42.666666666666664 * Hm100) / x
               + (46.22222222222222 * Hm101) / x
               - (4.977777777777778 * Hm1m10) / x - 3370.56 * x
               - 470.40000000000003 * H000 * x - 437.3333333333333 * H0000 * x
               - 351.99999999999994 * H0001 * x - 139.02222222222224 * H001 * x
               - 213.3333333333333 * H0010 * x - 64. * H00100 * x
               - 213.3333333333333 * H0011 * x - 128. * H00m10 * x
               + 64. * H00m100 * x - 128. * H00m101 * x - 256. * H00m1m10 * x
               - 2575.36 * H01 * x - 838.4000000000001 * H010 * x
               - 1162.6666666666665 * H0100 * x - 256. * H01000 * x
               - 960. * H01001 * x - 32. * H0101 * x
               - 843.7333333333332 * H011 * x - 95.99999999999999 * H0110 * x
               + 831.9999999999999 * H01100 * x - 32. * H0111 * x
               - 393.9555555555556 * H0m10 * x + 885.3333333333331 * H0m100 * x
               + 256. * H0m1000 * x + 64. * H0m1001 * x + 128. * H0m101 * x
               - 256. * H0m10m10 * x - 1856.0000000000002 * H0m1m10 * x
               - 447.99999999999994 * H0m1m100 * x - 128. * H0m1m101 * x
               + 256. * H0m1m1m10 * x - 3276.0355555555557 * H1 * x
               - 1197.8666666666666 * H10 * x - 435.2000000000001 * H100 * x
               + 95.99999999999999 * H1000 * x - 629.3333333333333 * H1001 * x
               - 149.3333333333333 * H101 * x - 576. * H10m10 * x
               - 1480.5333333333333 * H11 * x - 298.6666666666666 * H110 * x
               + 437.3333333333333 * H1100 * x - 207.99999999999997 * H111 * x
               - 2758.8029629629623 * Hm10 * x - 1390.2222222222224 * Hm100 * x
               - 223.99999999999997 * Hm1000 * x
               - 767.9999999999999 * Hm10001 * x
               - 223.99999999999997 * Hm1001 * x
               - 1809.7777777777776 * Hm101 * x
               - 42.666666666666664 * Hm1010 * x
               + 767.9999999999999 * Hm10100 * x
               - 42.666666666666664 * Hm1011 * x
               - 213.3333333333333 * Hm10m10 * x
               - 148.62222222222223 * Hm1m10 * x
               - 95.99999999999999 * Hm1m100 * x
               + 234.66666666666666 * Hm1m101 * x
               + 341.3333333333333 * Hm1m1m10 * x + 3425.70074074074 * x2
               + 38.400000000000006 * H000 * x2 + 746.6666666666666 * H0001 * x2
               + 744.1777777777778 * H001 * x2 + 810.6666666666666 * H00m10 * x2
               + 1528.32 * H01 * x2 + 82.13333333333333 * H010 * x2
               - 512. * H0100 * x2 - 767.9999999999999 * H01001 * x2
               + 204.80000000000004 * H011 * x2
               + 767.9999999999999 * H01100 * x2
               + 1142.0444444444445 * H0m10 * x2 + 512. * H0m100 * x2
               + 21.333333333333332 * H0m101 * x2
               - 490.6666666666666 * H0m1m10 * x2 + 3298.6311111111113 * H1 * x2
               + 1330.1333333333332 * H10 * x2 + 708.2666666666667 * H100 * x2
               + 95.99999999999999 * H1000 * x2 + 640. * H1001 * x2
               + 282.6666666666666 * H101 * x2
               + 447.99999999999994 * H10m10 * x2
               + 1612.8000000000002 * H11 * x2 + 229.3333333333333 * H110 * x2
               - 383.99999999999994 * H1100 * x2
               + 223.99999999999997 * H111 * x2 - 582.1866666666667 * Hm10 * x2
               + 138.66666666666666 * Hm100 * x2
               + 95.99999999999999 * Hm1000 * x2
               - 767.9999999999999 * Hm10001 * x2 + 128. * Hm1001 * x2
               - 903.111111111111 * Hm101 * x2
               + 767.9999999999999 * Hm10100 * x2
               - 383.99999999999994 * Hm10m10 * x2
               - 1359.6444444444444 * Hm1m10 * x2 - 640. * Hm1m100 * x2
               - 256. * Hm1m101 * x2 + 512. * Hm1m1m10 * x2
               - 583.6800000000001 * H000 * x3 - 307.20000000000005 * H0000 * x3
               + 25.600000000000005 * H0001 * x3 - 462.08 * H001 * x3
               - 51.20000000000001 * H0010 * x3 - 51.20000000000001 * H0011 * x3
               + 204.80000000000004 * H00m10 * x3
               + 117.33333333333333 * H010 * x3 - 409.6000000000001 * H0100 * x3
               + 344.74666666666667 * H0m10 * x3
               + 383.99999999999994 * H0m100 * x3 + 512. * H0m101 * x3
               - 51.20000000000001 * H0m1m10 * x3
               - 76.80000000000001 * H1000 * x3 + 409.6000000000001 * H1001 * x3
               - 117.33333333333333 * H101 * x3
               + 153.60000000000002 * H10m10 * x3
               + 117.33333333333333 * H110 * x3 - 409.6000000000001 * H1100 * x3
               + 522.0977777777778 * Hm10 * x3 + 583.6800000000001 * Hm100 * x3
               + 230.4 * Hm1000 * x3 + 383.99999999999994 * Hm1001 * x3
               + 344.74666666666667 * Hm101 * x3
               + 51.20000000000001 * Hm1010 * x3
               + 51.20000000000001 * Hm1011 * x3
               - 51.20000000000001 * Hm10m10 * x3
               - 344.74666666666667 * Hm1m10 * x3
               - 383.99999999999994 * Hm1m100 * x3 - 512. * Hm1m101 * x3
               + 51.20000000000001 * Hm1m1m10 * x3
               + (39.82222222222222 * H0m10 + 17.066666666666666 * H0m100
                  + 17.066666666666666 * H0m101 - 17.066666666666666 * H0m1m10
                  + 12.800000000000004 * H1000 - 68.26666666666667 * H1001
                  - 25.60000000000001 * H10m10 + 68.26666666666667 * H1100
                  + 92.46814814814813 * Hm10 + 89.45777777777778 * Hm100
                  + 38.400000000000006 * Hm1000 + 64. * Hm1001
                  + 49.635555555555555 * Hm101 + 8.533333333333333 * Hm1010
                  + 8.533333333333333 * Hm1011 - 8.533333333333333 * Hm10m10
                  - 49.635555555555555 * Hm1m10 - 64. * Hm1m100
                  - 85.33333333333333 * Hm1m101 + 8.533333333333333 * Hm1m1m10
                  + 149.02518518518517 * x - 652.8000000000001 * H0m10 * x2)
                     / x2
               - 205.79555555555552 * zeta2 + 64. * H01 * zeta2
               + 191.99999999999997 * H0m1 * zeta2
               + 471.4666666666666 * H1 * zeta2 + 544. * H10 * zeta2
               + 64. * H11 * zeta2 + 1020.8000000000001 * Hm1 * zeta2
               + 95.99999999999999 * Hm10 * zeta2 - 128. * Hm1m1 * zeta2
               - (8.533333333333333 * H01 * zeta2) / x2
               - (25.600000000000005 * H0m1 * zeta2) / x2
               - (24.817777777777778 * H1 * zeta2) / x2
               + (51.20000000000001 * H10 * zeta2) / x2
               - (4.266666666666667 * H11 * zeta2) / x2
               - (74.45333333333333 * Hm1 * zeta2) / x2
               - (81.06666666666666 * Hm10 * zeta2) / x2
               + (89.60000000000001 * Hm1m1 * zeta2) / x2
               + (92.8711111111111 * zeta2) / x
               + (2.488888888888889 * H1 * zeta2) / x - (128. * H10 * zeta2) / x
               - (48.71111111111111 * Hm1 * zeta2) / x
               - 183.44296296296295 * x * zeta2 + 128. * H001 * x * zeta2
               + 960. * H01 * x * zeta2 + 1088. * H010 * x * zeta2
               + 128. * H011 * x * zeta2 - 1056. * H0m1 * x * zeta2
               - 191.99999999999997 * H0m10 * x * zeta2
               + 256. * H0m1m1 * x * zeta2 + 223.64444444444445 * H1 * x * zeta2
               + 447.99999999999994 * H10 * x * zeta2
               + 170.66666666666666 * H11 * x * zeta2
               + 1735.4666666666665 * Hm1 * x * zeta2
               + 405.3333333333333 * Hm10 * x * zeta2
               + 767.9999999999999 * Hm100 * x * zeta2 - 64. * Hm1m1 * x * zeta2
               - 1528.32 * x2 * zeta2 - 245.3333333333333 * H01 * x2 * zeta2
               + 767.9999999999999 * H010 * x2 * zeta2
               - 266.6666666666666 * H0m1 * x2 * zeta2
               - 962.4888888888889 * H1 * x2 * zeta2 - 608. * H10 * x2 * zeta2
               - 256. * H11 * x2 * zeta2 + 223.28888888888892 * Hm1 * x2 * zeta2
               - 95.99999999999999 * Hm10 * x2 * zeta2
               + 767.9999999999999 * Hm100 * x2 * zeta2
               + 512. * Hm1m1 * x2 * zeta2 + 522.0977777777778 * x3 * zeta2
               + 25.600000000000005 * H01 * x3 * zeta2
               - 537.6 * H0m1 * x3 * zeta2 + 289.7066666666667 * H1 * x3 * zeta2
               - 307.20000000000005 * H10 * x3 * zeta2
               + 25.600000000000005 * H11 * x3 * zeta2
               - 517.12 * Hm1 * x3 * zeta2 - 486.4 * Hm10 * x3 * zeta2
               + 537.6 * Hm1m1 * x3 * zeta2 - 247.46666666666667 * x * zeta2_2
               - 614.4000000000001 * Hm1 * x * zeta2_2
               + 945.0666666666665 * x2 * zeta2_2
               - 614.4000000000001 * Hm1 * x2 * zeta2_2 + 320. * x3 * zeta2_2
               - 83.2 * zeta3 - 240. * H1 * zeta3
               + 111.99999999999999 * Hm1 * zeta3
               - (100.26666666666667 * H1 * zeta3) / x2
               - (74.66666666666664 * Hm1 * zeta3) / x2
               + (132.26666666666662 * zeta3) / x + (128. * H1 * zeta3) / x
               - 2702.577777777778 * x * zeta3 - 480. * H01 * x * zeta3
               - 223.99999999999997 * H0m1 * x * zeta3
               - 933.3333333333331 * H1 * x * zeta3
               + 37.33333333333332 * Hm1 * x * zeta3
               - 767.9999999999999 * Hm10 * x * zeta3
               + 1538.8444444444444 * x2 * zeta3
               - 767.9999999999999 * H01 * x2 * zeta3 + 544. * H1 * x2 * zeta3
               - 447.99999999999994 * Hm1 * x2 * zeta3
               - 767.9999999999999 * Hm10 * x2 * zeta3
               + 509.86666666666673 * x3 * zeta3 + 601.6 * H1 * x3 * zeta3
               - 447.99999999999994 * Hm1 * x3 * zeta3
               + 256. * x * zeta2 * zeta3
               + H00
                     * (-5.937777777777779
                        + x
                              * (-182.7081481481481 + 223.99999999999997 * zeta2
                                 + x
                                       * (1669.12 - 746.6666666666666 * zeta2
                                          + x
                                                * (-522.0977777777778
                                                   + 179.20000000000002 * zeta2)
                                       )
                                 + 128. * zeta3))
               + (H0
                  * (-59.899259259259225 + 29.86666666666667 * zeta2
                     + x
                           * (-48.68148148148148 - 82.13333333333333 * zeta2
                              + x
                                    * (-3091.413333333333
                                       + (-254.93333333333337
                                          - 38.400000000000006 * zeta2)
                                             * zeta2
                                       + 2208. * zeta3
                                       + x
                                             * (3126.257777777778
                                                - 744.1777777777778 * zeta2
                                                + 874.6666666666666 * zeta3)
                                       + x2
                                             * (806.8266666666667 * zeta2
                                                + 1049.6000000000001 * zeta3))))
                 ) / x
               + 3920. * x * zeta5)
        + CF * CF * nf
              * (20.053333333333335 - 40. * H000 - 4.266666666666667 * H001
                 + 256. * H00m10 - 257.06666666666666 * H01
                 - 47.99999999999999 * H010 - 64. * H011 + 512. * H0m100
                 + 256. * H0m101 - 256. * H0m1m10 + 146.75555555555556 * H1
                 + 104.00000000000001 * H10 - 84.26666666666667 * H100
                 + 256. * H1000 + 576. * H1001 - 95.99999999999999 * H101
                 + 88. * H11 - 64. * H110 - 448.00000000000006 * H1100
                 - 80. * H111 + 1754.3111111111111 * Hm10
                 + 1501.8666666666666 * Hm100 + 256. * Hm1000 + 64. * Hm1001
                 + 806.4 * Hm101 - 256. * Hm10m10 - 917.3333333333334 * Hm1m10
                 - 448.00000000000006 * Hm1m100 - 128. * Hm1m101
                 + 256. * Hm1m1m10 + (15.644444444444446 * H00) / x
                 + (17.066666666666666 * H000) / x
                 - (38.400000000000006 * H001) / x
                 + (54.04444444444446 * H01) / x
                 - (17.066666666666666 * H0m10) / x
                 - (54.04444444444446 * H1) / x + (55.46666666666667 * H100) / x
                 - (128. * H1001) / x + (128. * H1100) / x
                 + (1.4222222222222225 * Hm10) / x
                 - (17.066666666666666 * Hm100) / x
                 - (5.68888888888889 * Hm101) / x
                 + (28.444444444444443 * Hm1m10) / x + 1787.9466666666667 * x
                 - 374.7555555555556 * H000 * x + 165.33333333333331 * H0000 * x
                 - 53.33333333333332 * H0001 * x - 848.3555555555556 * H001 * x
                 - 224.00000000000003 * H0010 * x + 128. * H00100 * x
                 - 256. * H0011 * x - 85.33333333333333 * H00m10 * x
                 - 128. * H00m100 * x + 256. * H00m101 * x + 512. * H00m1m10 * x
                 + 1114.4888888888888 * H01 * x + 64. * H010 * x
                 + 650.6666666666666 * H0100 * x + 512. * H01000 * x
                 + 1152. * H01001 * x - 191.99999999999997 * H0101 * x
                 + 80. * H011 * x - 128. * H0110 * x
                 - 896.0000000000001 * H01100 * x - 160. * H0111 * x
                 + 65.42222222222223 * H0m10 * x
                 - 1301.3333333333333 * H0m100 * x - 512. * H0m1000 * x
                 - 128. * H0m1001 * x - 512. * H0m101 * x + 512. * H0m10m10 * x
                 + 1578.6666666666665 * H0m1m10 * x
                 + 896.0000000000001 * H0m1m100 * x + 256. * H0m1m101 * x
                 - 512. * H0m1m1m10 * x + 1556.3555555555556 * H1 * x
                 - 56.00000000000001 * H10 * x + 105.60000000000001 * H100 * x
                 + 85.33333333333333 * H1000 * x + 661.3333333333333 * H1001 * x
                 - 95.99999999999999 * H101 * x + 341.3333333333333 * H10m10 * x
                 + 40. * H11 * x - 64. * H110 * x
                 - 533.3333333333333 * H1100 * x - 80. * H111 * x
                 + 2378.311111111111 * Hm10 * x + 1484.088888888889 * Hm100 * x
                 - 85.33333333333333 * Hm1000 * x
                 + 767.9999999999999 * Hm10001 * x
                 - 106.66666666666664 * Hm1001 * x
                 + 1993.9555555555553 * Hm101 * x
                 - 767.9999999999999 * Hm10100 * x + 256. * Hm10m10 * x
                 + 220.44444444444443 * Hm1m10 * x
                 + 490.66666666666663 * Hm1m100 * x
                 + 213.3333333333333 * Hm1m101 * x - 256. * Hm1m1m10 * x
                 - 1752.32 * x2 + 352. * H000 * x2
                 - 746.6666666666666 * H0001 * x2
                 - 681.2444444444445 * H001 * x2
                 - 938.6666666666666 * H00m10 * x2
                 + 196.26666666666665 * H01 * x2 + 128. * H010 * x2
                 + 576. * H0100 * x2 + 767.9999999999999 * H01001 * x2
                 + 160. * H011 * x2 - 767.9999999999999 * H01100 * x2
                 - 1142.0444444444445 * H0m10 * x2 - 1024. * H0m100 * x2
                 - 341.3333333333333 * H0m101 * x2
                 + 682.6666666666666 * H0m1m10 * x2
                 - 1649.0666666666668 * H1 * x2 - 47.99999999999999 * H10 * x2
                 - 76.80000000000001 * H100 * x2
                 - 383.99999999999994 * H1000 * x2
                 - 832.0000000000001 * H1001 * x2
                 + 191.99999999999997 * H101 * x2 - 256. * H10m10 * x2
                 - 128. * H11 * x2 + 128. * H110 * x2 + 576. * H1100 * x2
                 + 160. * H111 * x2 + 435.2000000000001 * Hm10 * x2
                 - 198.4 * Hm100 * x2 - 383.99999999999994 * Hm1000 * x2
                 + 767.9999999999999 * Hm10001 * x2 - 256. * Hm1001 * x2
                 + 1103.6444444444444 * Hm101 * x2
                 - 767.9999999999999 * Hm10100 * x2 + 512. * Hm10m10 * x2
                 + 1244.4444444444443 * Hm1m10 * x2 + 1024. * Hm1m100 * x2
                 + 512. * Hm1m101 * x2 - 512. * Hm1m1m10 * x2
                 + 196.26666666666665 * H000 * x3 + 102.4 * H0000 * x3
                 - 230.39999999999998 * H0001 * x3
                 + 93.86666666666666 * H001 * x3 - 102.4 * H00m10 * x3
                 + 332.8 * H0100 * x3 - 93.86666666666666 * H0m10 * x3
                 - 102.4 * H0m100 * x3 - 204.8 * H0m101 * x3 + 51.2 * H1000 * x3
                 - 332.8 * H1001 * x3 - 102.4 * H10m10 * x3 + 332.8 * H1100 * x3
                 - 227.41333333333338 * Hm10 * x3
                 - 196.26666666666665 * Hm100 * x3 - 51.2 * Hm1000 * x3
                 - 102.4 * Hm1001 * x3 - 93.86666666666666 * Hm101 * x3
                 + 93.86666666666666 * Hm1m10 * x3 + 102.4 * Hm1m100 * x3
                 + 204.8 * Hm1m101 * x3
                 + (-8.533333333333333 * H1000 + 55.46666666666667 * H1001
                    + 17.066666666666666 * H10m10 - 55.46666666666667 * H1100
                    - 37.19111111111111 * Hm10 - 32.711111111111116 * Hm100
                    - 8.533333333333333 * Hm1000 - 17.066666666666666 * Hm1001
                    - 15.644444444444446 * Hm101 + 15.644444444444446 * Hm1m10
                    + 17.066666666666666 * Hm1m100 + 34.13333333333333 * Hm1m101
                    - 55.68000000000001 * x
                    + H0m10 * (-17.066666666666666 + 823.4666666666667 * x2))
                       / x2
                 + 257.06666666666666 * zeta2 - 128. * H01 * zeta2
                 - 383.99999999999994 * H0m1 * zeta2
                 - 362.66666666666663 * H1 * zeta2 - 704. * H10 * zeta2
                 - 128. * H11 * zeta2 - 1265.0666666666666 * Hm1 * zeta2
                 - 191.99999999999997 * Hm10 * zeta2 + 256. * Hm1m1 * zeta2
                 + (7.822222222222223 * H1 * zeta2) / x2
                 - (46.93333333333333 * H10 * zeta2) / x2
                 + (23.466666666666665 * Hm1 * zeta2) / x2
                 + (25.6 * Hm10 * zeta2) / x2
                 - (34.13333333333333 * Hm1m1 * zeta2) / x2
                 - (52.622222222222234 * zeta2) / x
                 - (14.222222222222221 * H1 * zeta2) / x
                 + (128. * H10 * zeta2) / x
                 + (19.91111111111111 * Hm1 * zeta2) / x
                 + 1263.8222222222223 * x * zeta2 - 256. * H001 * x * zeta2
                 - 597.3333333333333 * H01 * x * zeta2
                 - 1408. * H010 * x * zeta2 - 256. * H011 * x * zeta2
                 + 1301.3333333333333 * H0m1 * x * zeta2
                 + 383.99999999999994 * H0m10 * x * zeta2
                 - 512. * H0m1m1 * x * zeta2
                 - 14.222222222222221 * H1 * x * zeta2
                 - 618.6666666666666 * H10 * x * zeta2 - 128. * H11 * x * zeta2
                 - 1883.7333333333333 * Hm1 * x * zeta2 + 64. * Hm10 * x * zeta2
                 - 767.9999999999999 * Hm100 * x * zeta2
                 - 341.3333333333333 * Hm1m1 * x * zeta2
                 - 196.26666666666665 * x2 * zeta2
                 + 341.3333333333333 * H01 * x2 * zeta2
                 - 767.9999999999999 * H010 * x2 * zeta2
                 + 682.6666666666666 * H0m1 * x2 * zeta2
                 + 430.2222222222222 * H1 * x2 * zeta2 + 960. * H10 * x2 * zeta2
                 + 256. * H11 * x2 * zeta2
                 - 481.4222222222222 * Hm1 * x2 * zeta2
                 + 383.99999999999994 * Hm10 * x2 * zeta2
                 - 767.9999999999999 * Hm100 * x2 * zeta2
                 - 767.9999999999999 * Hm1m1 * x2 * zeta2
                 - 227.41333333333338 * x3 * zeta2 + 204.8 * H0m1 * x3 * zeta2
                 - 46.93333333333333 * H1 * x3 * zeta2
                 + 281.6 * H10 * x3 * zeta2 + 140.8 * Hm1 * x3 * zeta2
                 + 153.60000000000002 * Hm10 * x3 * zeta2
                 - 204.8 * Hm1m1 * x3 * zeta2 + 377.6 * x * zeta2_2
                 + 614.4000000000001 * Hm1 * x * zeta2_2
                 - 874.6666666666666 * x2 * zeta2_2
                 + 614.4000000000001 * Hm1 * x2 * zeta2_2
                 - 337.92 * x3 * zeta2_2
                 + (H0
                    * (38.61333333333334 + 21.333333333333332 * zeta2
                       + x
                             * (10.204444444444444 + 4.266666666666667 * zeta2
                                + x
                                      * (1307.1822222222222
                                         + zeta2
                                               * (913.7777777777777
                                                  + 76.80000000000001 * zeta2)
                                         + x
                                               * (-1310.72
                                                  + 681.2444444444445 * zeta2
                                                  - 874.6666666666666 * zeta3)
                                         + x2
                                               * (-187.73333333333332 * zeta2
                                                  - 588.8000000000001 * zeta3)
                                         - 1696. * zeta3))))
                       / x
                 + H00
                       * (-55.55555555555555
                          + x
                                * (-988.8888888888888 - 32. * zeta2
                                   + x
                                         * (-47.99999999999999
                                            + 746.6666666666666 * zeta2
                                            + x
                                                  * (227.41333333333338
                                                     + 128. * zeta2))
                                   - 256. * zeta3))
                 + 296.5333333333333 * zeta3 + 95.99999999999999 * H1 * zeta3
                 - 224.00000000000003 * Hm1 * zeta3
                 + (72.53333333333333 * H1 * zeta3) / x2
                 + (25.6 * Hm1 * zeta3) / x2 - (98.13333333333333 * zeta3) / x
                 - (128. * H1 * zeta3) / x + 2831.288888888889 * x * zeta3
                 + 191.99999999999997 * H01 * x * zeta3
                 + 448.00000000000006 * H0m1 * x * zeta3
                 + 714.6666666666666 * H1 * x * zeta3 + 288. * Hm1 * x * zeta3
                 + 767.9999999999999 * Hm10 * x * zeta3
                 - 1400.177777777778 * x2 * zeta3
                 + 767.9999999999999 * H01 * x2 * zeta3 - 320. * H1 * x2 * zeta3
                 + 640. * Hm1 * x2 * zeta3
                 + 767.9999999999999 * Hm10 * x2 * zeta3
                 - 234.66666666666666 * x3 * zeta3
                 - 435.2000000000001 * H1 * x3 * zeta3
                 + 153.60000000000002 * Hm1 * x3 * zeta3
                 - 128. * x * zeta2 * zeta3 - 4000. * x * zeta5)
        + CA * CA * nf
              * (239.7037037037037 - 160. * H000 - 102.93333333333334 * H001
                 + 6.222222222222221 * H01 + 64. * H010 + 96. * H011
                 - 84.26666666666667 * H0m10 - 356. * H1
                 + 229.33333333333331 * H10 + 370.66666666666663 * H100
                 + 96. * H1001 + 224. * H101 + 298.66666666666663 * H11
                 + 224. * H110 - 96. * H1100 + 192. * H111
                 + 441.68888888888887 * Hm10 + 4.266666666666667 * Hm100
                 + 120.53333333333333 * Hm101 - 272. * Hm1m10
                 - (4.266666666666667 * H1000) / x2
                 + (21.333333333333332 * H1001) / x2
                 + (8.533333333333333 * H10m10) / x2
                 - (21.333333333333332 * H1100) / x2
                 + (5.066666666666666 * Hm10) / x2
                 + (2.6666666666666665 * Hm100) / x2
                 - (4.266666666666667 * Hm1000) / x2
                 - (8.533333333333333 * Hm1001) / x2
                 + (2.6666666666666665 * Hm101) / x2
                 - (2.6666666666666665 * Hm1m10) / x2
                 + (8.533333333333333 * Hm1m100) / x2
                 + (17.066666666666666 * Hm1m101) / x2
                 - (11.200000000000001 * H00) / x
                 + (-550.5283950617284 + 8.533333333333333 * H000 - 12.8 * H001
                    - 4.0888888888888895 * H01 - 42.666666666666664 * H010
                    - 53.33333333333332 * H011 + 12.8 * H0m10
                    - 95.46666666666667 * H1 - 190.2222222222222 * H10
                    - 64. * H100 - 32. * H1001 - 74.66666666666666 * H101
                    - 199.1111111111111 * H11 - 74.66666666666666 * H110
                    + 32. * H1100 - 64. * H111 + 112.5333333333333 * Hm10
                    - 93.86666666666666 * Hm100 - 63.28888888888889 * Hm101
                    - 3.5555555555555554 * Hm1m10)
                       / x
                 + (11159.85185185185 + 3019.733333333333 * H000
                    + 1589.333333333333 * H0000 + 2128. * H0001
                    + 3686.3999999999996 * H001 + 1792.0000000000002 * H0010
                    + 2048. * H0011 + 330.6666666666666 * H00m10
                    + 4763.555555555556 * H01 + 3264. * H010
                    + 1477.3333333333333 * H0100 + 191.99999999999997 * H01001
                    + 1280. * H0101 + 3680. * H011 + 1280. * H0110
                    - 191.99999999999997 * H01100 + 1152. * H0111
                    + 366.93333333333334 * H0m10 + 1002.6666666666666 * H0m100
                    + 629.3333333333333 * H0m101 + 160. * H0m1m10
                    + 4531.851851851851 * H1 + 4451.555555555556 * H10
                    + 1669.3333333333333 * H100)
                       * x
                 + 154.66666666666666 * H1000 * x
                 + 586.6666666666666 * H1001 * x + 1888. * H101 * x
                 + 384. * H1010 * x + 448. * H1011 * x
                 + 74.66666666666666 * H10m10 * x + 4782.222222222222 * H11 * x
                 + 2122.6666666666665 * H110 * x + 309.3333333333333 * H1100 * x
                 + 448. * H1101 * x + 1845.3333333333333 * H111 * x
                 + 448. * H1110 * x + 384. * H1111 * x
                 + 78.57777777777778 * Hm10 * x + 2324.2666666666664 * Hm100 * x
                 + 549.3333333333333 * Hm1000 * x + 192. * Hm10001 * x
                 + 682.6666666666666 * Hm1001 * x + 1507.2 * Hm101 * x
                 + 128. * Hm1010 * x - 192. * Hm10100 * x + 128. * Hm1011 * x
                 - 192. * Hm10m10 * x - 752. * Hm1m10 * x
                 - 810.6666666666666 * Hm1m100 * x
                 - 725.3333333333333 * Hm1m101 * x + 128. * Hm1m1m10 * x
                 - 10849.027160493826 * x2 - 704. * H000 * x2
                 - 266.66666666666663 * H0001 * x2
                 - 2209.2444444444445 * H001 * x2 - 128. * H0010 * x2
                 - 256. * H0011 * x2 - 202.66666666666666 * H00m10 * x2
                 - 7092.977777777778 * H01 * x2 - 2048. * H010 * x2
                 - 32. * H0100 * x2 + 192. * H01001 * x2 - 384. * H0101 * x2
                 - 2250.6666666666665 * H011 * x2 - 384. * H0110 * x2
                 - 192. * H01100 * x2 - 384. * H0111 * x2
                 + 1689.4222222222224 * H0m10 * x2 + 448. * H0m100 * x2
                 + 490.66666666666663 * H0m101 * x2
                 - 405.3333333333333 * H0m1m10 * x2
                 - 4080.3851851851855 * H1 * x2 - 4490.666666666666 * H10 * x2
                 - 1976. * H100 * x2 - 176. * H1000 * x2 - 544. * H1001 * x2
                 - 2037.3333333333333 * H101 * x2 - 384. * H1010 * x2
                 - 448. * H1011 * x2 - 32. * H10m10 * x2
                 - 4881.777777777777 * H11 * x2 - 2272. * H110 * x2
                 - 352. * H1100 * x2 - 448. * H1101 * x2
                 - 1973.3333333333333 * H111 * x2 - 448. * H1110 * x2
                 - 384. * H1111 * x2 - 225.24444444444444 * Hm10 * x2
                 + 2239.4666666666667 * Hm100 * x2 + 528. * Hm1000 * x2
                 + 192. * Hm10001 * x2 + 640. * Hm1001 * x2
                 + 1336.7111111111112 * Hm101 * x2 + 128. * Hm1010 * x2
                 - 192. * Hm10100 * x2 + 128. * Hm1011 * x2
                 - 192. * Hm10m10 * x2 - 496.8888888888889 * Hm1m10 * x2
                 - 768. * Hm1m100 * x2 - 640. * Hm1m101 * x2
                 + 128. * Hm1m1m10 * x2 - 16. * H000 * x3 + 51.2 * H0000 * x3
                 - 76.80000000000001 * H0001 * x3 - 16. * H001 * x3
                 - 51.2 * H00m10 * x3 + 128. * H0100 * x3 + 16. * H0m10 * x3
                 - 51.2 * H0m100 * x3 - 102.4 * H0m101 * x3 + 25.6 * H1000 * x3
                 - 128. * H1001 * x3 - 51.2 * H10m10 * x3 + 128. * H1100 * x3
                 + 30.400000000000002 * Hm10 * x3 + 16. * Hm100 * x3
                 - 25.6 * Hm1000 * x3 - 51.2 * Hm1001 * x3 + 16. * Hm101 * x3
                 - 16. * Hm1m10 * x3 + 51.2 * Hm1m100 * x3
                 + 102.4 * Hm1m101 * x3 - 6.222222222222221 * zeta2
                 - 360. * H1 * zeta2 - 96. * H10 * zeta2
                 - 256.5333333333333 * Hm1 * zeta2
                 - (1.3333333333333333 * H1 * zeta2) / x2
                 - (17.066666666666666 * H10 * zeta2) / x2
                 - (4. * Hm1 * zeta2) / x2 + (12.8 * Hm10 * zeta2) / x2
                 - (17.066666666666666 * Hm1m1 * zeta2) / x2
                 + (116.62222222222222 * zeta2) / x
                 + (76.44444444444444 * H1 * zeta2) / x
                 + (32. * H10 * zeta2) / x
                 + (61.51111111111112 * Hm1 * zeta2) / x
                 - 4684.977777777778 * x * zeta2 - 1360. * H01 * x * zeta2
                 - 192. * H010 * x * zeta2
                 - 549.3333333333333 * H0m1 * x * zeta2 - 1512. * H1 * x * zeta2
                 - 453.3333333333333 * H10 * x * zeta2 - 384. * H11 * x * zeta2
                 - 1883.2 * Hm1 * x * zeta2 - 816. * Hm10 * x * zeta2
                 - 192. * Hm100 * x * zeta2
                 + 789.3333333333333 * Hm1m1 * x * zeta2
                 + 7092.977777777778 * x2 * zeta2
                 + 181.33333333333331 * H01 * x2 * zeta2
                 - 192. * H010 * x2 * zeta2
                 - 693.3333333333333 * H0m1 * x2 * zeta2
                 + 1788.8888888888887 * H1 * x2 * zeta2
                 + 432. * H10 * x2 * zeta2 + 384. * H11 * x2 * zeta2
                 - 1585.1555555555556 * Hm1 * x2 * zeta2
                 - 752. * Hm10 * x2 * zeta2 - 192. * Hm100 * x2 * zeta2
                 + 704. * Hm1m1 * x2 * zeta2 + 30.400000000000002 * x3 * zeta2
                 + 102.4 * H0m1 * x3 * zeta2 + 8. * H1 * x3 * zeta2
                 + 102.4 * H10 * x3 * zeta2 - 24. * Hm1 * x3 * zeta2
                 + 76.80000000000001 * Hm10 * x3 * zeta2
                 - 102.4 * Hm1m1 * x3 * zeta2 - 38.400000000000006 * x * zeta2_2
                 + 153.60000000000002 * Hm1 * x * zeta2_2
                 - 511.46666666666664 * x2 * zeta2_2
                 + 153.60000000000002 * Hm1 * x2 * zeta2_2
                 - 138.24 * x3 * zeta2_2
                 + H00
                       * (184.26666666666665
                          + x
                                * (4392.711111111112
                                   - 1797.3333333333333 * zeta2
                                   + x
                                         * (-8556.977777777778
                                            + 266.66666666666663 * zeta2
                                            + x
                                                  * (-30.400000000000002
                                                     + 25.6 * zeta2))))
                 - 328.5333333333333 * zeta3 + 96. * H1 * zeta3
                 + (29.866666666666664 * H1 * zeta3) / x2
                 + (12.8 * Hm1 * zeta3) / x2 + (64. * zeta3) / x
                 - (32. * H1 * zeta3) / x - 1720. * x * zeta3
                 + 192. * H01 * x * zeta3 - 250.66666666666666 * H1 * x * zeta3
                 - 704. * Hm1 * x * zeta3 + 192. * Hm10 * x * zeta3
                 + 3251.2888888888892 * x2 * zeta3 + 192. * H01 * x2 * zeta3
                 + 336. * H1 * x2 * zeta3 - 640. * Hm1 * x2 * zeta3
                 + 192. * Hm10 * x2 * zeta3 + 40. * x3 * zeta3
                 - 179.20000000000002 * H1 * x3 * zeta3
                 + 76.80000000000001 * Hm1 * x3 * zeta3
                 - 96. * x * zeta2 * zeta3
                 + (H0
                    * (-74.45925925925926 + 25.6 * zeta2
                       + x
                             * (-465.4222222222222 + 102.93333333333334 * zeta2
                                + x
                                      * (7496.948148148148
                                         - 3319.4666666666667 * zeta2
                                         - 1962.6666666666663 * zeta3
                                         + x
                                               * (-547.2888888888889
                                                  + (2209.2444444444445
                                                     + 32. * x)
                                                        * zeta2
                                                  - 250.66666666666666 * zeta3
                                                  - 256. * x * zeta3)))))
                       / x
                 - 960. * x * zeta5)
        + (CF * nf * nf
           * (-4.266666666666667 * Hm100 - 1.4222222222222223 * Hm101
              + 1.4222222222222223 * Hm1m10
              + Hm10
                    * (-6.73185185185185
                       + (-8.533333333333333 - 63.28888888888889 * x) * x)
              + x
                    * (15.747160493827161 + 5.3096296296296295 * H0
                       + 4.266666666666667 * H00 + 1.422222222222222 * H01
                       - 2.6074074074074076 * H1 + 7.111111111111111 * H10
                       + 7.111111111111111 * H11
                       + (582.957037037037 + 377.2562962962963 * H0
                          + 152.53333333333333 * H00 + 80. * H000 + 64. * H001
                          + 137.95555555555555 * H01 + 32. * H010 + 32. * H011
                          + 292.8 * H1 + 74.66666666666666 * H10
                          + 74.66666666666666 * H11)
                             * x
                       + (-1142.6607407407407 - 470.32888888888897 * H0
                          - 422.2814814814815 * H00 - 549.3333333333333 * H000
                          - 352. * H0000 - 192. * H0001
                          - 199.1111111111111 * H001 - 64. * H0010 - 64. * H0011
                          - 119.82222222222221 * H01 - 10.666666666666666 * H011
                         ) * x2)
              + H0m10
                    * (-2.8444444444444446
                       + x3
                             * (135.11111111111111
                                + (-42.666666666666664 - 8.533333333333333 * x)
                                      * x))
              + 0.7111111111111111 * H1 * zeta2
              + 2.1333333333333333 * Hm1 * zeta2 - 9.955555555555556 * x * zeta2
              + H1 * x3
                    * (-371.02222222222224 + 14.222222222222221 * zeta2
                       + x
                             * (80.82962962962964
                                + (10.666666666666666 - 25.6 * x) * zeta2))
              + x2
                    * (-170.66666666666666 * H11 * x
                       + 10.666666666666666 * H110 * x
                       + 81.3037037037037 * Hm10 * x
                       + 21.333333333333332 * Hm100 * x
                       + 7.111111111111111 * Hm101 * x
                       - 7.111111111111111 * Hm1m10 * x + 543.9565432098766 * x2
                       - 37.83111111111111 * H0 * x2
                       + 68.26666666666667 * H00 * x2 - 256. * H000 * x2
                       - 138.66666666666666 * H001 * x2
                       + 104.53333333333333 * H01 * x2 - 32. * H010 * x2
                       - 42.666666666666664 * H011 * x2
                       + 88.88888888888889 * H11 * x2
                       + 10.666666666666666 * H110 * x2
                       + 112.35555555555555 * Hm10 * x2
                       + 30.435555555555556 * H00 * x3 + 25.6 * H000 * x3
                       + 29.866666666666664 * H001 * x3
                       - 21.333333333333332 * H010 * x3
                       - 21.333333333333332 * H110 * x3
                       - 30.435555555555556 * Hm10 * x3 - 25.6 * Hm100 * x3
                       - 8.533333333333333 * Hm101 * x3
                       + 8.533333333333333 * Hm1m10 * x3
                       + H10 * x * (-181.33333333333331 + 99.55555555555554 * x)
                       + H101 * x
                             * (-10.666666666666666
                                + x
                                      * (-10.666666666666666
                                         + 21.333333333333332 * x))
                       - 137.95555555555555 * zeta2 - 64. * H0 * zeta2
                       + 201.12592592592594 * x * zeta2
                       + 334.22222222222223 * H0 * x * zeta2
                       + 192. * H00 * x * zeta2
                       - 10.666666666666666 * Hm1 * x * zeta2
                       - 104.53333333333333 * x2 * zeta2
                       + 138.66666666666666 * H0 * x2 * zeta2
                       - 30.435555555555556 * x3 * zeta2
                       - 38.400000000000006 * H0 * x3 * zeta2
                       + 12.8 * Hm1 * x3 * zeta2 + 6.4 * x * zeta2_2
                       + (-32.
                          + x
                                * (412.44444444444434 + 128. * H0
                                   + x
                                         * (-74.66666666666666
                                            + 42.666666666666664 * x)))
                             * zeta3)))
              / x2
        + (CA * nf * nf
           * ((28.85925925925926 + 1.0666666666666667 * H0
               + 14.814814814814813 * H1 + 3.5555555555555554 * H10
               + 3.5555555555555554 * H11)
                  * x
              + Hm10
                    * (-1.0666666666666667
                       + x
                             * (-3.5555555555555554
                                + x
                                      * (-5.333333333333333
                                         + x
                                               * (-56.888888888888886
                                                  + (-60.44444444444444
                                                     - 6.4 * x)
                                                        * x))))
              + x
                    * ((-47.14074074074074 - 9.422222222222222 * H0
                        + 5.333333333333333 * H00 - 21.333333333333332 * H1
                        - 10.666666666666666 * H10 - 10.666666666666666 * H11)
                           * x
                       + x4 * (6.4 * H00 - 6.4 * zeta2)
                       - 3.5555555555555554 * zeta2
                       + x3
                             * (739.5851851851853 + 138.66666666666666 * H00
                                + 21.333333333333332 * H001
                                + 170.66666666666666 * H01
                                + 21.333333333333332 * H010
                                + 21.333333333333332 * H011
                                - 42.666666666666664 * H0m10
                                + 452.1481481481481 * H1
                                + 138.66666666666666 * H10
                                + 21.333333333333332 * H100
                                + 138.66666666666666 * H11
                                + 42.666666666666664 * H110
                                + 21.333333333333332 * H111 - 64. * Hm100
                                - 21.333333333333332 * Hm101
                                + 21.333333333333332 * Hm1m10
                                + H0
                                      * (382.1037037037038
                                         - 21.333333333333332 * zeta2)
                                - 170.66666666666666 * zeta2
                                + 10.666666666666666 * H1 * zeta2
                                + 32. * Hm1 * zeta2 - 48. * zeta3)
                       + x2
                             * (-721.3037037037038 - 337.77777777777777 * H00
                                - 128. * H000 - 106.66666666666666 * H001
                                - 293.3333333333333 * H01 - 64. * H010
                                - 64. * H011 - 10.666666666666666 * H0m10
                                - 445.6296296296296 * H1
                                - 131.55555555555554 * H10
                                - 21.333333333333332 * H100
                                - 131.55555555555554 * H11
                                - 42.666666666666664 * H110
                                - 21.333333333333332 * H111 - 64. * Hm100
                                - 21.333333333333332 * Hm101
                                + 21.333333333333332 * Hm1m10
                                + 236.44444444444446 * zeta2
                                - 10.666666666666666 * H1 * zeta2
                                + 32. * Hm1 * zeta2
                                + H0 * (-643.7925925925927 + 96. * zeta2)
                                + 74.66666666666666 * zeta3))))
              / x2;
}

//==========================================================================================//
//  Massless gluon coefficient functions for FL at
//  O(as^3) for mu=Q.
//  The term fl_ps_11 is put to zero for the reason
//  explained in page 15 of arXiv:1205.5727
//
//  Eq. (B.18) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double MasslessCoefficientFunction::CL_ps3_massless(double x, int nf) const {

    // the commented varibles are needed for the flav term

    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;
    double x5 = x4 * x;

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
    // const double Hm1m1 = Hr2[0];
    const double H0m1 = Hr2[1];
    const double Hm10 = Hr2[3];
    const double H00 = Hr2[4];
    const double H10 = Hr2[5];
    const double H01 = Hr2[7];
    const double H11 = Hr2[8];

    // weight 3
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
    const double H0m1m10 = Hr4[28];
    const double H00m10 = Hr4[31];
    // const double H10m10 = Hr4[32];
    // const double Hm1m100 = Hr4[36];
    const double H0m100 = Hr4[37];
    // const double Hm1000 = Hr4[39];
    const double H0000 = Hr4[40];
    // const double H1000 = Hr4[41];
    const double H0100 = Hr4[43];
    const double H1100 = Hr4[44];
    const double H0010 = Hr4[49];
    const double H0110 = Hr4[52];
    // const double Hm1m101 = Hr4[63];
    const double H0m101 = Hr4[64];
    // const double Hm1001 = Hr4[66];
    const double H0001 = Hr4[67];
    const double H1001 = Hr4[68];
    const double H0101 = Hr4[70];
    const double H0011 = Hr4[76];
    const double H0111 = Hr4[79];

    //  weight 5
    // const double Hm10100 = Hr5[129];
    // const double H01100 = Hr5[133];
    // const double Hm10001 = Hr5[201];
    // const double H01001 = Hr5[205];

    delete[] Hr1;
    delete[] Hr2;
    delete[] Hr3;
    delete[] Hr4;
    delete[] Hr5;

    // double flav = 5. / 18 * fl11ps(nf) * nf;
    // this contribution is neglected
    // If you want to use it please check the 5/48 since
    // I'm not sure about it

    return /*flav*(-119.46666666666665*H001 + 298.66666666666663*H001*x -
     1228.8000000000002*H01*x - 85.33333333333333*H0100*x - 512.*H01001*x +
     512.*H01100*x + 85.33333333333333*H0m10*x + 85.33333333333333*H0m100*x +
     426.66666666666663*H0m101*x + 256.*H0m1m10*x - 230.4*H1*x + 25.6*H100*x -
     128.*H1000*x + 85.33333333333333*H1001*x + 256.*H10m10*x
     - 85.33333333333333*H1100*x - 554.6666666666666*Hm10*x -
     324.26666666666665*Hm100*x + 128.*Hm1000*x - 512.*Hm10001*x +
     341.3333333333333*Hm1001*x - 605.8666666666667*Hm101*x + 512.*Hm10100*x
     + 42.666666666666664*Hm1m10*x - 341.3333333333333*Hm1m100*x -
     682.6666666666666*Hm1m101*x + 153.60000000000002*x2 + 512.*H0001*x2 -
     358.40000000000003*H001*x2 + 512.*H01*x2 - 512.*H0100*x2 - 768.*H0m10*x2 +
     512.*H1*x2 - 204.8*H100*x2 + 512.*H1001*x2 - 512.*H1100*x2 - 512.*Hm10*x2 -
     204.8*Hm100*x2 + 358.40000000000003*Hm101*x2 + 768.*Hm1m10*x2 +
     512.*H000*x3 + 512.*H001*x3 + 204.8*H0100*x3 - 512.*H0m10*x3 -
     204.8*H0m100*x3 - 409.6*H0m101*x3 - 204.8*H1001*x3 + 204.8*H1100*x3 -
     665.6*Hm10*x3 - 512.*Hm100*x3 - 204.8*Hm1001*x3 - 512.*Hm101*x3 +
     512.*Hm1m10*x3 + 204.8*Hm1m100*x3 + 409.6*Hm1m101*x3 +
     (H00*(-341.3333333333333 + x*(170.66666666666666 + x*(-418.1333333333333 +
     x*(512. + 665.6*x - 512.*zeta2)))))/x - 128.*H1*zeta2 + 512.*H10*zeta2 -
     247.46666666666667*Hm1*zeta2 - (170.66666666666666*H1*zeta2)/x2 +
     (136.53333333333333*H10*zeta2)/x2 - (512.*Hm1*zeta2)/x2 -
     (136.53333333333333*Hm10*zeta2)/x2 + (273.06666666666666*Hm1m1*zeta2)/x2 +
     (682.6666666666666*zeta2)/x + (192.*H1*zeta2)/x - (256.*H10*zeta2)/x -
     (81.06666666666666*Hm1*zeta2)/x + 674.1333333333333*x*zeta2 -
     128.*H01*x*zeta2 + 512.*H010*x*zeta2 - 298.66666666666663*H0m1*x*zeta2 -
     21.333333333333332*H1*x*zeta2 + 42.666666666666664*H10*x*zeta2 +
     627.2*Hm1*x*zeta2 - 469.3333333333333*Hm10*x*zeta2 + 512.*Hm100*x*zeta2 +
     682.6666666666666*Hm1m1*x*zeta2 - 512.*x2*zeta2 + 384.*H1*x2*zeta2 -
     512.*H10*x2*zeta2 + 25.6*Hm1*x2*zeta2 - 665.6*x3*zeta2 +
     409.6*H0m1*x3*zeta2 - 256.*H1*x3*zeta2 + 204.8*H10*x3*zeta2 +
     768.*Hm1*x3*zeta2 + 204.8*Hm10*x3*zeta2 - 409.6*Hm1m1*x3*zeta2 +
     115.2*x*zeta2_2 - 409.6*Hm1*x*zeta2_2 + 409.6*x2*zeta2_2 -
     122.88*x3*zeta2_2 + (136.53333333333336*H1100 + 443.7333333333334*Hm10 +
     341.33333333333326*Hm100 + 136.53333333333336*Hm1001 +
        341.33333333333326*Hm101 - 341.33333333333326*Hm1m10 -
     136.53333333333336*Hm1m100 - 273.0666666666667*Hm1m101 +
        (102.39999999999999 - 341.33333333333326*H01 + 341.33333333333326*H1 -
     136.53333333333336*H100)*x + H1001*(-136.53333333333336 +
     (255.99999999999997 - 511.99999999999994*x)*x) + x*(341.33333333333326*Hm10
     + 136.53333333333336*Hm100 - 110.93333333333335*Hm101 - 384.*Hm1m10 +
           (-153.60000000000005 + 605.8666666666667*H01 - 622.9333333333334*H1 +
     187.73333333333335*H100)*x + H1100*(-255.99999999999997 +
     511.99999999999994*x) + x*(85.33333333333331*Hm10 +
     187.73333333333335*Hm100 + 119.46666666666665*Hm101 -
     255.99999999999997*Hm1m10 - 230.4*x - 85.33333333333331*H000*x -
              605.8666666666667*zeta2 - 324.2666666666666*zeta3)))/x2 -
     512.*H1*zeta3 - (68.26666666666667*H1*zeta3)/x2 - (204.8*Hm1*zeta3)/x2 +
     (273.06666666666666*zeta3)/x + (256.*H1*zeta3)/x -
     1565.8666666666666*x*zeta3 - 512.*H01*x*zeta3 +
     298.66666666666663*H1*x*zeta3 - 512.*Hm1*x*zeta3 - 512.*Hm10*x*zeta3 +
     179.20000000000002*x2*zeta3 + 512.*H1*x2*zeta3 - 1280.*x3*zeta3 -
     102.4*H1*x3*zeta3 + 307.20000000000005*Hm1*x3*zeta3 + 256.*x*zeta2*zeta3 +
     (H0*
        (-102.4 + x*(-34.13333333333333 + 119.46666666666665*zeta2 +
             x*(-179.20000000000002 - 213.33333333333331*zeta2 +
     170.66666666666666*zeta3 + x*(665.6 + (358.40000000000003 - 1024.*x)*zeta2
     + 512.*zeta3 - 409.6*x*zeta3))
             )))/x + 2560.*x*zeta5)*/
        +CA * CF * nf
            * (696.8888888888888 - 160. * H000 - 124.80000000000001 * H001
               + 110.2222222222222 * H01 + 31.99999999999999 * H011
               - 95.99999999999997 * H0m10 - 280.8888888888888 * H1
               + 85.33333333333333 * H10 + 284.79999999999995 * H100
               + 160. * H101 + 74.66666666666664 * H11 + 160. * H110
               + 160. * H111 + 216.17777777777775 * Hm10
               - 202.66666666666666 * Hm100 - 106.66666666666663 * Hm101
               - 85.33333333333333 * Hm1m10 - (17.06666666666666 * H1001) / x2
               + (17.06666666666666 * H1100) / x2
               - (1.422222222222222 * Hm10) / x2
               - (5.688888888888888 * Hm100) / x2
               - (5.688888888888888 * Hm101) / x2
               + (5.688888888888888 * Hm1m10) / x2 - 528.4740740740741 / x
               + (5.688888888888889 * H00) / x
               + (17.06666666666666 * H001 - 25.6 * H01
                  - 42.666666666666664 * H010 - 53.33333333333329 * H011
                  + 21.333333333333332 * H0m10 + 60.562962962962956 * H1
                  - 140.4444444444444 * H10 - 81.06666666666666 * H100
                  + 21.333333333333332 * H1001 - 53.33333333333329 * H101
                  - 129.77777777777766 * H11 - 53.33333333333329 * H110
                  - 21.333333333333332 * H1100 - 53.33333333333329 * H111
                  + 92.08888888888889 * Hm10 - 106.66666666666659 * Hm100
                  - 74.66666666666664 * Hm101 + 10.666666666666666 * Hm1m10)
                     / x
               + 1146.074074074074 * x + 392.88888888888886 * H000 * x
               + 405.3333333333333 * H0000 * x + 400. * H0001 * x
               - 20.977777777777778 * H001 * x + 320. * H0010 * x
               + 352. * H0011 * x + 74.66666666666666 * H00m10 * x
               - 158.2222222222222 * H01 * x + 32. * H010 * x
               + 293.3333333333333 * H0100 * x + 160. * H0101 * x
               + 5.333333333333333 * H011 * x + 160. * H0110 * x
               + 160. * H0111 * x - 184.88888888888889 * H0m10 * x
               + 202.66666666666666 * H0m100 * x
               + 106.66666666666666 * H0m101 * x
               + 85.33333333333333 * H0m1m10 * x - 380.4444444444444 * H1 * x
               + 250.66666666666666 * H10 * x - 50.13333333333333 * H100 * x
               - 21.333333333333332 * H1001 * x + 234.66666666666666 * H11 * x
               + 21.333333333333332 * H1100 * x - 226.48888888888894 * Hm10 * x
               + 103.11111111111109 * Hm100 * x + 103.11111111111109 * Hm101 * x
               - 103.11111111111109 * Hm1m10 * x - 1314.488888888889 * x2
               + 42.666666666666664 * H0001 * x2 - 102.4 * H001 * x2
               - 681.6 * H01 * x2 - 149.33333333333331 * H010 * x2
               - 42.666666666666664 * H0100 * x2
               - 149.33333333333331 * H011 * x2
               + 234.66666666666666 * H0m10 * x2 + 600.7703703703704 * H1 * x2
               - 195.55555555555554 * H10 * x2 - 153.60000000000002 * H100 * x2
               + 42.666666666666664 * H1001 * x2
               - 106.66666666666666 * H101 * x2 - 179.55555555555554 * H11 * x2
               - 106.66666666666666 * H110 * x2
               - 42.666666666666664 * H1100 * x2
               - 106.66666666666666 * H111 * x2 - 347.02222222222224 * Hm10 * x2
               + 213.33333333333331 * Hm100 * x2
               + 149.33333333333331 * Hm101 * x2
               - 21.333333333333332 * Hm1m10 * x2
               - 8.533333333333333 * H000 * x3 - 25.6 * H0001 * x3
               - 8.533333333333333 * H001 * x3 + 25.6 * H0100 * x3
               + 8.533333333333333 * H0m10 * x3 - 25.6 * H1001 * x3
               + 25.6 * H1100 * x3 + 2.1333333333333333 * Hm10 * x3
               + 8.533333333333333 * Hm100 * x3 + 8.533333333333333 * Hm101 * x3
               - 8.533333333333333 * Hm1m10 * x3 - 110.2222222222222 * zeta2
               - 202.66666666666666 * H1 * zeta2 + 64. * Hm1 * zeta2
               + (2.8444444444444446 * H1 * zeta2) / x2
               + (17.066666666666666 * H10 * zeta2) / x2
               + (8.533333333333333 * Hm1 * zeta2) / x2
               + (117.6888888888889 * zeta2) / x + (48. * H1 * zeta2) / x
               - (21.333333333333332 * H10 * zeta2) / x
               + (80. * Hm1 * zeta2) / x - 68.26666666666667 * x * zeta2
               - 202.66666666666666 * H01 * x * zeta2 - 64. * H0m1 * x * zeta2
               + 51.55555555555554 * H1 * x * zeta2
               + 21.333333333333332 * H10 * x * zeta2
               - 154.66666666666666 * Hm1 * x * zeta2 + 681.6 * x2 * zeta2
               + 96. * H1 * x2 * zeta2 - 42.666666666666664 * H10 * x2 * zeta2
               - 160. * Hm1 * x2 * zeta2 + 2.1333333333333333 * x3 * zeta2
               + 4.266666666666667 * H1 * x3 * zeta2 + 25.6 * H10 * x3 * zeta2
               - 12.8 * Hm1 * x3 * zeta2 - 92.26666666666667 * x * zeta2_2
               + 34.13333333333333 * x2 * zeta2_2 - 20.48 * x3 * zeta2_2
               + H00
                     * (341.6888888888889
                        + x
                              * (546.1333333333333 - 325.3333333333333 * zeta2
                                 + x
                                       * (-1366.0444444444445
                                          - 42.666666666666664 * zeta2
                                          + x
                                                * (-2.1333333333333333
                                                   + 25.6 * zeta2))))
               - 364.7999999999999 * zeta3
               + ((123.73333333333335 * x
                   + H1
                         * (-17.066666666666666 + 21.333333333333332 * x
                            - 21.333333333333332 * x3 + 42.666666666666664 * x4
                            - 25.6 * x5)
                   + x3
                         * (-97.42222222222223
                            + x * (441.6 + 21.333333333333332 * x)))
                  * zeta3)
                     / x2
               + (H0
                  * (-84.85925925925926 + 4.266666666666667 * zeta2
                     + x
                           * (-151.82222222222222 + 124.80000000000001 * zeta2
                              + x
                                    * (754.2518518518518
                                       - 163.91111111111113 * zeta2
                                       - 522.6666666666666 * zeta3
                                       + x
                                             * (1197.3925925925928
                                                + (102.4
                                                   + 17.066666666666666 * x)
                                                      * zeta2
                                                + 42.666666666666664 * zeta3
                                                - 25.6 * x * zeta3)))))
                     / x)
        + CF * CF * nf
              * (113.99111111111114 + 160. * H000 + 281.6 * H001
                 + 255.28888888888886 * H01 + 128. * H010 + 128. * H011
                 + 668.7999999999998 * H1 + 181.33333333333326 * H10
                 - 153.60000000000002 * H100 + 64. * H101
                 + 250.6666666666666 * H11 + 64. * H110 + 32. * H111
                 - 435.2 * Hm10 - 64. * Hm100 + 21.333333333333325 * Hm101
                 + 234.6666666666666 * Hm1m10 + (8.533333333333333 * H00) / x
                 + (-34.133333333333326 * H001 + 31.288888888888888 * H01
                    - 147.43703703703704 * H1 - 7.111111111111111 * H10
                    + 34.133333333333326 * H100 - 42.66666666666665 * H1001
                    - 21.333333333333325 * H101 - 23.1111111111111 * H11
                    - 21.333333333333325 * H110 + 42.66666666666665 * H1100
                    - 10.666666666666663 * H111 - 7.822222222222222 * Hm10
                    + 21.333333333333325 * Hm101 + 21.333333333333325 * Hm1m10)
                       / x
                 - 160.80592592592592 * x + 222.22222222222217 * H000 * x
                 + 64. * H0000 * x + 320. * H0001 * x
                 + 415.2888888888889 * H001 * x + 128. * H0010 * x
                 + 128. * H0011 * x - 21.333333333333332 * H00m10 * x
                 + 171.73333333333332 * H01 * x + 96. * H010 * x
                 - 170.66666666666666 * H0100 * x + 64. * H0101 * x
                 + 208. * H011 * x + 64. * H0110 * x + 32. * H0111 * x
                 + 35.55555555555556 * H0m10 * x + 64. * H0m100 * x
                 - 21.333333333333332 * H0m101 * x
                 - 234.66666666666666 * H0m1m10 * x - 307.9111111111111 * H1 * x
                 - 181.33333333333331 * H10 * x + 68.26666666666667 * H100 * x
                 + 42.666666666666664 * H1001 * x - 170.66666666666666 * H11 * x
                 - 42.666666666666664 * H1100 * x - 422.162962962963 * Hm10 * x
                 - 85.33333333333333 * Hm100 * x - 35.55555555555556 * Hm101 * x
                 + 163.55555555555554 * Hm1m10 * x + 1.6829629629629625 * x2
                 - 256. * H000 * x2 - 85.33333333333333 * H0001 * x2
                 - 307.20000000000005 * H001 * x2
                 - 132.26666666666665 * H01 * x2 - 128. * H010 * x2
                 + 85.33333333333333 * H0100 * x2
                 - 106.66666666666666 * H011 * x2
                 + 42.666666666666664 * H0m10 * x2
                 - 213.45185185185187 * H1 * x2 + 7.111111111111111 * H10 * x2
                 + 51.2 * H100 * x2 - 85.33333333333333 * H1001 * x2
                 - 42.666666666666664 * H101 * x2
                 - 56.888888888888886 * H11 * x2
                 - 42.666666666666664 * H110 * x2
                 + 85.33333333333333 * H1100 * x2
                 - 21.333333333333332 * H111 * x2
                 + 17.066666666666666 * Hm10 * x2
                 - 42.666666666666664 * Hm101 * x2
                 - 42.666666666666664 * Hm1m10 * x2 - 12.8 * H000 * x3
                 + 51.2 * H0001 * x3 + 4.266666666666667 * H001 * x3
                 - 51.2 * H0100 * x3 - 4.266666666666667 * H0m10 * x3
                 + 51.2 * H1001 * x3 - 51.2 * H1100 * x3
                 + 4.835555555555556 * Hm10 * x3 + 12.8 * Hm100 * x3
                 - 4.266666666666667 * Hm101 * x3
                 + 4.266666666666667 * Hm1m10 * x3
                 + (34.133333333333326 * H1001 - 34.133333333333326 * H1100
                    - 7.0162962962962965 * Hm10 - 8.533333333333331 * Hm100
                    + 2.8444444444444437 * Hm101 - 2.8444444444444437 * Hm1m10
                    + 45.13185185185186 * x
                    + H0m10 * (-11.377777777777775 - 213.33333333333331 * x2))
                       / x2
                 - 255.28888888888886 * zeta2 + 53.33333333333333 * H1 * zeta2
                 + 96. * Hm1 * zeta2 - (1.4222222222222223 * H1 * zeta2) / x2
                 - (34.13333333333333 * H10 * zeta2) / x2
                 - (4.266666666666667 * Hm1 * zeta2) / x2
                 - (39.11111111111111 * zeta2) / x
                 + (10.666666666666666 * H1 * zeta2) / x
                 + (42.666666666666664 * H10 * zeta2) / x
                 - (10.666666666666666 * Hm1 * zeta2) / x
                 - 593.8962962962963 * x * zeta2
                 + 53.33333333333333 * H01 * x * zeta2 - 96. * H0m1 * x * zeta2
                 - 81.77777777777777 * H1 * x * zeta2
                 - 42.666666666666664 * H10 * x * zeta2
                 + 117.33333333333333 * Hm1 * x * zeta2
                 + 132.26666666666665 * x2 * zeta2
                 + 21.333333333333332 * H1 * x2 * zeta2
                 + 85.33333333333333 * H10 * x2 * zeta2
                 + 21.333333333333332 * Hm1 * x2 * zeta2
                 + 4.835555555555556 * x3 * zeta2
                 - 2.1333333333333333 * H1 * x3 * zeta2
                 - 51.2 * H10 * x3 * zeta2 + 6.4 * Hm1 * x3 * zeta2
                 + 94.93333333333334 * x * zeta2_2
                 - 68.26666666666667 * x2 * zeta2_2 + 40.96 * x3 * zeta2_2
                 + H00
                       * (-37.33333333333333
                          + x
                                * (168.29629629629628
                                   - 341.3333333333333 * zeta2
                                   + x
                                         * (-64.
                                            + x
                                                  * (-4.835555555555556
                                                     - 51.2 * zeta2)
                                            + 85.33333333333333 * zeta2)))
                 + 25.6 * zeta3
                 + ((-76.80000000000001 * x
                     + x3
                           * (-514.4888888888889
                              + (-8.533333333333333 - 10.666666666666666 * x)
                                    * x)
                     + H1
                           * (34.13333333333333 - 42.666666666666664 * x
                              + 42.666666666666664 * x3 - 85.33333333333333 * x4
                              + 51.2 * x5))
                    * zeta3)
                       / x2
                 + (H0
                    * (9.86074074074074 + 34.13333333333333 * zeta2
                       + x
                             * (100.66962962962963 - 281.6 * zeta2
                                + x
                                      * (-476.1600000000001
                                         - 379.73333333333335 * zeta2
                                         + 245.33333333333326 * zeta3
                                         + x
                                               * (-19.176296296296297
                                                  + (307.20000000000005
                                                     - 8.533333333333333 * x)
                                                        * zeta2
                                                  - 85.33333333333333 * zeta3
                                                  + 51.2 * x * zeta3)))))
                       / x)
        + (CF * nf * nf
           * ((39.032098765432096 + 2.8444444444444446 * H0
               + 5.925925925925926 * H1 - 3.555555555555555 * H11)
                  * x
              + Hm10
                    * (-2.8444444444444446 - 7.11111111111111 * x
                       + 14.222222222222221 * x3 + 14.222222222222221 * x4
                       + 4.266666666666667 * x5)
              + x
                    * ((-116.62222222222222 - 51.2 * H0
                        - 21.333333333333332 * H00 - 24.888888888888886 * H1
                        + 10.666666666666666 * H11)
                           * x
                       - 7.11111111111111 * zeta2
                       + x4
                             * (-4.266666666666667 * H00
                                + 4.266666666666667 * zeta2)
                       + x3
                             * (51.04197530864197 + 60.91851851851852 * H0
                                + 42.666666666666664 * H00
                                - 7.111111111111111 * H01
                                + 22.51851851851852 * H1
                                - 7.111111111111111 * H11
                                + 7.111111111111111 * zeta2)
                       + x2
                             * (26.54814814814815 - 66.60740740740741 * H0
                                - 128. * H00 - 64. * H000
                                - 39.11111111111111 * H01
                                + 10.666666666666666 * H011
                                - 3.5555555555555554 * H1
                                + 53.33333333333332 * zeta2
                                - 10.666666666666666 * zeta3))))
              / x2;
}

//==========================================================================================//
//  Charge factors for the flavour topologies entering up to
//  three loops.
//
//  Table 2 of Ref. [arXiv:hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

// u, d, s, c, b, t
double charges[] = { 2. / 3., -1. / 3., -1. / 3., 2. / 3., -1. / 3., 2. / 3. };

double MasslessCoefficientFunction::fl11g(int nf) const {

    double eavg = 0., e2avg = 0.;

    for (int i = 0; i < nf; i++) {
        eavg += charges[i];
        e2avg += charges[i] * charges[i];
    }
    eavg /= nf;
    e2avg /= nf;

    return eavg * eavg / e2avg;
}

double MasslessCoefficientFunction::fl11ps(int nf) const {

    double eavg = 0., e2avg = 0.;

    for (int i = 0; i < nf; i++) {
        eavg += charges[i];
        e2avg += charges[i] * charges[i];
    }
    eavg /= nf;
    e2avg /= nf;

    return eavg * eavg / e2avg - 3 * eavg;
}

//==========================================================================================//
//  Massless gluon coefficient functions for F2 at O(as^2) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long)
//  result
//
//  Eq. (4.10) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

/*
double MasslessCoefficientFunction::C2_g2_massless_param(double x, int nf) const
{

    double x1 = 1. - x;
    double L0 = log(x);
    double L1 = log(x1);

    return nf
           * (58. / 9. * L1 * L1 * L1 - 24 * L1 * L1 - 34.88 * L1 + 30.586
              - (25.08 + 760.3 * x + 29.65 * L1 * L1 * L1) * x1
              + 1204 * x * L0 * L0 + L0 * L1 * (293.8 + 711.2 * x + 1043 * L0)
              + 115.6 * L0 - 7.109 * L0 * L0 + 70. / 9. * L0 * L0 * L0
              + 11.9033 * x1 / x);
}
*/
//==========================================================================================//
//  Massless quark coefficient functions for F2 at O(as^2) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long)
//  result
//
//  Eq. (4.9) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//
/*
double MasslessCoefficientFunction::C2_ps2_massless_param(double x, int nf)
const {

    double x1 = 1 - x;
    double L0 = log(x);
    double L02 = L0 * L0;
    double L03 = L02 * L0;

    double L1 = log(x1);
    double L12 = L1 * L1;
    double L13 = L12 * L1;

    return nf
           * ((8. / 3 * L12 - 32. / 3 * L1 + 9.8937) * x1
              + (9.57 - 13.41 * x + 0.08 * L13) * x1 * x1 + 5.667 * x * L03
              - L02 * L1 * (20.26 - 33.93 * x) + 43.36 * x1 * L0 - 1.053 * L02
              + 40. / 9 * L03 + 5.2903 / x * x1 * x1);
}
*/
//==========================================================================================//
//  Massless gluon coefficient functions for FL at O(as^2) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long)
//  result
//
//  Eq. (6) of Ref. [arXiv:hep-ph/0411112v2]
//------------------------------------------------------------------------------------------//
/*
double MasslessCoefficientFunction::CL_g2_massless_param(double x, int nf) const
{

    double L0 = log(x);
    double L02 = L0 * L0;

    double x1 = 1. - x;
    double L1 = log(x1);
    double L12 = L1 * L1;

    return nf
           * ((94.74 - 49.2 * x) * x1 * L12 + 864.8 * x1 * L1
              + 1161 * x * L1 * L0 + 60.06 * x * L02 + 39.66 * x1 * L0
              - 5.333 * (1. / x - 1));
}
*/
//==========================================================================================//
//  Massless gluon coefficient functions for FL at O(as^2) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long)
//  result
//
//  Eq. (5) of Ref. [arXiv:hep-ph/0411112v2]
//------------------------------------------------------------------------------------------//
/*
double MasslessCoefficientFunction::CL_ps2_massless_param(double x, int nf)
const {

    double L0 = log(x);
    double L02 = L0 * L0;

    double x1 = 1. - x;
    double L1 = log(x1);

    return nf
           * ((15.94 - 5.212 * x) * x1 * x1 * L1 + (0.421 + 1.520 * x) * L02
              + 28.09 * x1 * L0 - (2.370 / x - 19.27) * x1 * x1 * x1);
}
*/
//==========================================================================================//
//  Massless gluon coefficient functions for F2 at O(as^3) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long)
//  result. The term fl_g_11 is put to zero for the reason explained in page 15
//  of arXiv:1205.5727
//
//  Eq. (4.13) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//
/*
double MasslessCoefficientFunction::C2_g3_massless_param(
    double x, int nf
) const { // remember that there is a delta(x1) that has been omitted

    // double fl_g_11 = fl11g(nf) ;

    double x2 = x * x;
    // double x3 = x2 * x;

    double x1 = 1 - x;
    double L0 = log(x);
    double L1 = log(x1);

    double L02 = L0 * L0;
    double L03 = L02 * L0;
    double L04 = L03 * L0;
    double L05 = L04 * L0;

    double L12 = L1 * L1;
    double L13 = L12 * L1;
    double L14 = L13 * L1;
    double L15 = L14 * L1;

    double c_nf =
        (966. / 81. * L15 - 1871. / 18. * L14 + 89.31 * L13 + 979.2 * L12
         - 2405 * L1 + 1372 * x1 * L14 - 15729 - 310510 * x + 331570 * x2
         - 244150 * x * L02 - 253.3 * x * L05 + L0 * L1 * (138230 - 237010 * L0)
         - 11860 * L0 - 700.8 * L02 - 1440 * L03 + 4961. / 162. * L04
         - 134. / 9. * L05
         - (6362.54 + 932.089 * L0)
               / x // there's a typo in the paper: -932.089 -> +932.089
        );         // there is + 0.625*delta(x1)

    double c_nf2 =
        (131. / 81. * L14 - 14.72 * L13 + 3.607 * L12 - 226.1 * L1 + 4.762
         - 190 * x - 818.4 * x2 - 4019 * x * L02 - L0 * L1 * (791.5 + 4646 * L0)
         + 739.0 * L0 + 418.0 * L02 + 104.3 * L03 + 809. / 81. * L04
         + 12. / 9. * L05 + 84.423 / x);

    // double c_nf_fl = (
    //     3.211 * L12 + 19.04 * x * L1 + 0.623 * x1 * L13 - 64.47
    //     * x + 121.6 * x2 - 45.82 * x3 - x * L0 * L1 * (31.68 +
    //     37.24 * L0) + 11.27 * x2 * L03 - 82.40 * x * L0 - 16.08
    //     * x * L02 + 520./81. * x * L03 + 20./27. * x * L04
    // );

    return (
        c_nf * nf + c_nf2 * nf * nf
        //+ c_nf_fl * nf * nf * fl_g_11
    );
}
*/
//==========================================================================================//
//  Massless quark coefficient functions for F2 at O(as^3) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long)
//  result. The term fl_ps_11 is put to zero for the reason explained in page 15
//  of arXiv:1205.5727
//
//  Eq. (4.12) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//
/*
double MasslessCoefficientFunction::C2_ps3_massless_param(
    double x, int nf
) const { // remember that there is a delta(x1) that has been omitted

    // double fl_ps_11 = fl11ps(nf) ;

    double x2 = x * x;

    double x1 = 1. - x;
    double L0 = log(x);
    double L1 = log(x1);

    double L02 = L0 * L0;
    double L03 = L02 * L0;
    double L04 = L03 * L0;
    double L05 = L04 * L0;

    double L12 = L1 * L1;
    double L13 = L12 * L1;
    double L14 = L13 * L1;

    double c_nf =
        ((856. / 81 * L14 - 6032. / 81 * L13 + 130.57 * L12 - 542 * L1 + 8501
          - 4714 * x + 61.5 * x2)
             * x1
         + L0 * L1 * (8831 * L0 + 4162 * x1) - 15.44 * x * L05 + 3333 * x * L02
         + 1615 * L0 + 1208 * L02 - 333.73 * L03 + 4244. / 81 * L04
         - 40. / 9 * L05 - 1. / x * (2731.82 * x1 + 414.262 * L0));

    double c_nf2 =
        ((-64. / 81 * L13 + 208. / 81 * L12 + 23.09 * L1 - 220.27 + 59.80 * x
          - 177.6 * x2)
             * x1
         - L0 * L1 * (160.3 * L0 + 135.4 * x1) - 24.14 * x * L03
         - 215.4 * x * L02 - 209.8 * L0 - 90.38 * L02 - 3568. / 243 * L03
         - 184. / 81 * L04 + 40.2426 * x1 / x);

    // double c_fl_nf = (
    //     (126.42 - 50.29 * x - 50.15 * x2) * x1 - 26.717
    //     - 9.075 * x * x1 * L1 - x * L02 * (101.8 + 34.79 * L0 + 3.070 * L02)
    //     + 59.59 * L0 - 320./81 * L02 * (5 + L0)
    // ) * x ;

    return (
        c_nf * nf + c_nf2 * nf * nf
        //+ c_fl_nf * fl_ps_11 * nf
    );
}
*/
//==========================================================================================//
//  Massless gluon coefficient functions for FL at O(as^3) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long)
//  result. The term fl_ps_11 is put to zero for the reason explained in page 15
//  of arXiv:1205.5727
//
//  Eq. (10) of Ref. [arXiv:hep-ph/0411112v2]
//------------------------------------------------------------------------------------------//
/*
double MasslessCoefficientFunction::CL_g3_massless_param(
    double x, int nf
) const { // remember that there is a delta(x1) that has been omitted

    // double fl_g_11 = fl11g(nf) ;

    double x2 = x * x;

    double x1 = 1. - x;
    double L0 = log(x);
    double L1 = log(x1);

    double L02 = L0 * L0;
    double L03 = L02 * L0;
    // double L04 = L03 * L0;

    double L12 = L1 * L1;
    double L13 = L12 * L1;
    double L14 = L13 * L1;

    double c_nf =
        ((144. * L14 - 47024. / 27. * L13 + 6319. * L12 + 53160. * L1) * x1
         + 72549. * L0 * L1 + 88238. * L02 * L1
         + (3709. - 33514. * x - 9533. * x2) * x1 + 66773. * x * L02
         - 1117. * L0 + 45.37 * L02 - 5360. / 27. * L03
         - (2044.70 * x1 + 409.506 * L0) / x);

    double c_nf2 =
        ((32. / 3. * L13 - 1216. / 9. * L12 - 592.3 * L1 + 1511. * x * L1) * x1
         + 311.3 * L0 * L1 + 14.24 * L02 * L1 + (577.3 - 729. * x) * x1
         + 30.78 * x * L03 + 366. * L0 + 1000. / 9. * L02 + 160. / 9. * L03
         + 88.5037 / x * x1);

    // double c_nf_fl = (
    //     (
    //         -0.0105 * L13 + 1.55 * L12 + 19.72 * x * L1 - 66.745 * x
    //         + 0.615 * x2
    //     ) * x1
    //     + 20./27. * x * L04 + (280./81. + 2.26 * x) * x * L03
    //     - (15.4 - 2.201 * x) * x * L02 - (71.66 - 0.121 * x) * x * L0
    // ) ;

    return (
        c_nf * nf + c_nf2 * nf * nf
        //+ c_nf_fl * nf * nf * fl_g_11
    );
}
*/
//==========================================================================================//
//  Massless gluon coefficient functions for FL at O(as^3) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long)
//  result. The term fl_ps_11 is put to zero for the reason explained in page 15
//  of arXiv:1205.5727
//
//  Eq. (9) of Ref. [arXiv:hep-ph/0411112v2]
//------------------------------------------------------------------------------------------//
/*
double MasslessCoefficientFunction::CL_ps3_massless_param(
    double x, int nf
) const { // remember that there is a delta(x1) that has been omitted

    // double fl_ps_11 = fl11ps(nf) ;

    // double x2=x*x;

    double x1 = 1. - x;
    double L0 = log(x);
    double L1 = log(x1);

    double L02 = L0 * L0;
    double L03 = L02 * L0;

    double L12 = L1 * L1;
    double L13 = L12 * L1;

    double c_nf =
        ((1568. / 27 * L13 - 3968. / 9 * L12 + 5124 * L1) * x1 * x1
         + (2184 * L0 + 6059 * x1) * L0 * L1 - (795.6 + 1036 * x) * x1 * x1
         - 143.6 * x1 * L0 + 2848. / 9 * L02 - 1600. / 27 * L03
         - (885.53 * x1 + 182.00 * L0) / x * x1);

    double c_nf2 =
        ((-32. / 9 * L12 + 29.52 * L1) * x1 * x1
         + (35.18 * L0 + 73.06 * x1) * L0 * L1 - 35.24 * x * L02
         - (14.16 - 69.84 * x) * x1 * x1 - 69.41 * x1 * L0 - 128. / 9 * L02
         + 40.239 / x * x1 * x1);

    // double c_fl_nf = (
    //     ((107.0 + 321.05 * x - 54.62 * x2) * x1 - 26.717 + 9.773
    //     * L0 + (363.8 + 68.32 * L0) * x * L0 - 320./81 * L02 * (2
    //     + L0)) * x
    // );

    return (
        c_nf * nf + c_nf2 * nf * nf
        //+ c_fl_nf * fl_ps_11 * nf
    );
}
*/
