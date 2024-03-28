#include "adani/MatchingCondition.h"
#include "adani/Constants.h"
#include "adani/SpecialFunctions.h"

#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

//==========================================================================================//
//  MatchingCondition: constructor
//------------------------------------------------------------------------------------------//

MatchingCondition::MatchingCondition(
    const int &order, const char &entry1, const char &entry2,
    const string &version
) {
    // check order
    if (order != 3) {
        cout << "Error: only order = 3 is implemented. Got " << order << endl;
        exit(-1);
    }
    order_ = order;

    // check entry1
    if (entry1 != 'Q') {
        cout << "Error: only entry1 = 'Q' is implemented. Got " << entry1
             << endl;
        exit(-1);
    }
    entry1_ = entry1;

    // check entry2
    if (entry2 != 'g' && entry2 != 'q') {
        cout << "Error: entry2 must be g or q. Got " << entry2 << endl;
        exit(-1);
    }
    entry2_ = entry2;

    // check version
    if (version != "exact" && version != "abmp" && version != "klmv"
        && version != "gm") {
        cout << "Error: version must be 'exact', 'abmp', 'gm' or "
                "'klmv'! Got "
             << version << endl;
        exit(-1);
    }

    if (entry2 == 'g' && version == "exact") {
        cout << "Error: aQg channel doesn't have 'exact' version!" << endl;
        exit(-1);
    }

    if (entry2 == 'q' && (version == "abmp" || version == "gm")) {
        cout << "Error: aQq channel doesn't have 'abmp' or 'gm' "
                "version!"
             << endl;
        exit(-1);
    }

    version_ = version;
}

//==========================================================================================//
//  MatchingCondition: nf independent part of the a_Qi (i=q,q) term.
//  a_Qi is the mu independent part of the unrenormalized OMA
//------------------------------------------------------------------------------------------//

Value MatchingCondition::MuIndependentNfIndependentTerm(double x) const {

    double central, higher, lower;

    vector<double> res = NotOrdered(x);

    central = res[0];

    if (res[1] > res[2]) {
        higher = res[1];
        lower = res[2];
    } else {
        higher = res[2];
        lower = res[1];
    }

    return Value(central, higher, lower);
}

//==========================================================================================//
//  MatchingCondition: nf independent part of the a_Qi (i=q,q) term without
//  ordering the upper and lower bands
//------------------------------------------------------------------------------------------//

vector<double> MatchingCondition::NotOrdered(double x) const {

    double central, higher, lower;
    if (entry2_ == 'q') {
        if (version_ == "exact") {
            central = a_Qq_PS_30(x, 0);
            return { central, central, central };
        } else if (version_ == "klmv") {
            higher = a_Qq_PS_30(x, 1);
            lower = a_Qq_PS_30(x, -1);
            central = 0.5 * (higher + lower);

            return { central, higher, lower };
        } else {
            cout << "Error: something has gone wrong in "
                    "MatchingCondition::MuIndependentNfIndependentTerm"
                 << endl;
            exit(-1);
        }
    } else {
        int low_id;
        if (version_ == "exact") {
            central = a_Qg_30(x, 0);
            return { central, central, central };
        } else if (version_ == "abmp")
            low_id = -1;
        else if (version_ == "klmv")
            low_id = -12;
        else {
            cout << "Error: something has gone wrong in "
                    "MatchingCondition::NotOrdered"
                 << endl;
            exit(-1);
        }

        higher = a_Qg_30(x, 1);
        lower = a_Qg_30(x, low_id);
        central = 0.5 * (higher + lower);

        return { central, higher, lower };
    }
}

//==========================================================================================//
//  Matching condition Qg O(as)
//
//  Eq. (B.2) from Ref. [arXiv:hep-ph/9612398v1]
//------------------------------------------------------------------------------------------//

// double MatchingCondition::K_Qg1(double x, double m2mu2) const {

//     return 2 * TR * (x * x + (x - 1) * (x - 1)) * log(1. / m2mu2);
// }

//==========================================================================================//
//  Local part of the matching condition gg O(as)
//
//  Eq. (B.6) from Ref. [arXiv:hep-ph/9612398v1]
//------------------------------------------------------------------------------------------//

// double MatchingCondition::K_gg1_local(double m2mu2) const { return -4. / 3. *
// TR * log(1. / m2mu2); }

//==========================================================================================//
//  Matching condition Qg O(as^2)
//
//  Eq. (B.3) from Ref. [arXiv:hep-ph/9612398v1]
//------------------------------------------------------------------------------------------//

// double MatchingCondition::K_Qg2(double x, double m2mu2) const {

//     double x2 = x * x;

//     double L = log(x);
//     double L2 = L * L;
//     double L3 = L2 * L;

//     double Lm = log(1. - x);
//     double Lm2 = Lm * Lm;
//     double Lm3 = Lm2 * Lm;

//     double Lp = log(1. + x);
//     double Lp2 = Lp * Lp;

//     double Li2xm = Li2(1. - x);
//     double Li3xm = Li3(1. - x);
//     double Li2minus = Li2(-x);
//     double Li3minus = Li3(-x);
//     double S12xm = S12(1. - x);
//     double S12minus = S12(-x);

//     double Lmu = log(m2mu2);
//     double Lmu2 = Lmu * Lmu;

//     double logmu2_CFTR =
//         (Lm * (8. - 16. * x + 16 * x2) - L * (4. - 8. * x + 16. * x2)
//          - (2. - 8. * x));

//     double logmu2_CATR =
//         (-Lm * (8. - 16. * x + 16. * x2) - (8. + 32. * x) * L - 16. / 3. / x
//          - 4. - 32. * x + 124. / 3. * x2);

//     double logmu2_TR2 = -16. / 3. * (2. * x2 - 2. * x + 1.);

//     double logmu2 =
//         CF * TR * logmu2_CFTR + CA * TR * logmu2_CATR + TR * TR * logmu2_TR2;

//     double logmu_CFTR =
//         ((8. - 16. * x + 16. * x2) * (2. * L * Lm - Lm2 + 2. * zeta2)
//          - (4. - 8. * x + 16. * x2) * L2 - 32. * x * (1. - x) * Lm
//          - (12. - 16. * x + 32. * x2) * L - 56. + 116. * x - 80. * x2);

//     double logmu_CATR =
//         ((16. + 32. * x + 32. * x2) * (Li2minus + L * Lp)
//          + (8. - 16. * x + 16. * x2) * Lm2 + (8. + 16. * x) * L2
//          + 32. * x * zeta2 + 32. * x * (1. - x) * Lm
//          - (8. + 64. * x + 352. / 3. * x2) * L - 160. / 9. / x + 16. - 200. *
//          x
//          + 1744. / 9. * x2);

//     double logmu = CF * TR * logmu_CFTR + CA * TR * logmu_CATR;

//     double const_CFTR =
//         ((1. - 2. * x + 2. * x2)
//              * (8. * zeta3 + 4. / 3. * Lm3 - 8. * Lm * Li2xm + 8. * zeta2 * L
//                 - 4. * L * Lm2 + 2. / 3. * L3 - 8. * L * Li2xm + 8. * Li3xm
//                 - 24. * S12xm)
//          + x2 * (-16. * zeta2 * L + 4. / 3. * L3 + 16 * L * Li2xm + 32 *
//          S12xm)
//          - (4. + 96. * x - 64. * x2) * Li2xm - (4. - 48. * x + 40. * x2) *
//          zeta2
//          - (8. + 48. * x - 24. * x2) * L * Lm + (4. + 8. * x - 12. * x2) *
//          Lm2
//          - (1. + 12. * x - 20. * x2) * L2 - (52. * x - 48. * x2) * Lm
//          - (16. + 18. * x + 48. * x2) * L + 26. - 82. * x + 80. * x2);

//     double const_CATR =
//         ((1. - 2. * x + 2. * x2) * (-4. / 3. * Lm3 + 8 * Lm * Li2xm - 8 *
//         Li3xm)
//          + (1. + 2. * x + 2. * x2)
//                * (-8. * zeta2 * Lp - 16 * Lp * Li2minus - 8 * L * Lp2
//                   + 4 * L2 * Lp + 8 * L * Li2minus - 8 * Li3minus
//                   - 16 * S12minus)
//          + (16. + 64. * x) * (2. * S12xm + L * Li2xm)
//          - (4. / 3. + 8. / 3. * x) * L3 + (8. - 32. * x + 16. * x2) * zeta3
//          - (16. + 64. * x) * zeta2 * L
//          + (16. * x + 16. * x2) * (Li2minus + L * Lp)
//          // there is a typo here in the e-Print (16+16*x2 -> 16*x+16*x2)
//          + (32. / 3. / x + 12. + 64. * x - 272. / 3. * x2) * Li2xm
//          - (12 + 48 * x - 260. / 3. * x2 + 32. / 3. / x) * zeta2
//          - 4 * x2 * L * Lm - (2 + 8 * x - 10 * x2) * Lm2
//          + (2 + 8 * x + 46. / 3. * x2) * L2 + (4 + 16 * x - 16 * x2) * Lm
//          - (56. / 3. + 172. / 3. * x + 1600. / 9. * x2) * L - 448. / 27. / x
//          - 4. / 3. - 628. / 3. * x + 6352. / 27. * x2);

//     double const_tot = CF * TR * const_CFTR + CA * TR * const_CATR;

//     double tmp = const_tot + logmu * Lmu + logmu2 * Lmu2;

//     return 0.5 * tmp;
// }

//==========================================================================================//
//  Approximation of the nf-independent part of the mu-independent part of the
//  unrenormalized matching condition Qg at O(as^3).
//
//  v = 0 : exact result
//  v = 1 : Eq. (3.49) of Ref. [arXiv:1205.5727]
//  v = -1 : Eq. (16) Ref. of [arXiv:1701.05838]
//  v = -12 : Eq. (3.50) of Ref. [arXiv:1205.5727]
//  v = 2 : approximation from Giacomo Magni, based on the results of [arXiv:2403.00513]
//------------------------------------------------------------------------------------------//

double MatchingCondition::a_Qg_30(double x, int v) const {

    double L = log(x);
    double L2 = L * L;

    double x1 = 1. - x;

    double L1 = log(x1);
    double L12 = L1 * L1;
    double L13 = L12 * L1;

    if (v == 0) {
        cout << "a_Qg_30 exact is not known/implemented yet!" << endl;
        exit(-1);
    } else if (v == 1) {
        return (
            354.1002 * L13 + 479.3838 * L12 - 7856.784 * (2. - x)
            - 6233.530 * L2 + 9416.621 / x + 1548.891 / x * L
        );
    } else if (v == -1) { // Updated version w.r.t v==-12
        return (
            226.3840 * L13 - 652.2045 * L12 - 2686.387 * L1
            - 7714.786 * (2. - x) - 2841.851 * L2 + 7721.120 / x
            + 1548.891 / x * L
        );
    } else if (v == -12) { // Version of the paper (used only for benchamrk)
        return (
            -2658.323 * L12 - 7449.948 * L1 - 7460.002 * (2. - x)
            + 3178.819 * L2 + 4710.725 / x + 1548.891 / x * L
        );
    } else if (v == 2) {
        double L14 = L13 * L1;
        double L15 = L14 * L1;
        double L3 = L2 * L;
        double L4 = L3 * L;
        double L5 = L4 * L;
        return -5882.68 + 8956.65 / x + 10318.5 * x - 8363.19 * x * x
               + 737.165 * L1 - 332.537 * L12 + 4.3802 * L13 - 8.20988 * L14
               + 3.7037 * L15 + 11013.4 * L + (1548.89 * L) / x
               + 6558.74 * x * L - 720.048 * L2 + 514.091 * L3 - 21.7593 * L4
               + 4.84444 * L5 - 274.207 * (-L1 + L)
               - 274.207 * (-1 + x - x * L1 + x * L);
    } else {
        cout << "Error in MatchingCondition::a_Qg_30: Choose either v=0, v=1, "
                "v=-1, v=-12 or v=2"
             << endl;
        exit(-1);
    }
}

//==========================================================================================//
//  nf-independent part of the mu-independent part of the unrenormalized
//  matching condition Qq at O(as^3). Both the exact resul and the
//  approximate one are implemented. The latter is used just as benchmark for
//  the plots of the paper.
//
//  v = 0 : exact result from Eq. (5.41, 5.42, 5.45) of Ref. [arXiv:1409.1135]
//  v = 1 : approximation from Eq. (3.53) of [arXiv:1205.5727]
//  v = 2 : approximation from Eq. (3.53) of [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double MatchingCondition::a_Qq_PS_30(double x, int v) const {

    double x2 = x * x;

    if (v == 0) {

        double x3 = x2 * x;
        double x4 = x3 * x;
        double x5 = x4 * x;
        double x6 = x5 * x;

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
        const double H0m1m1m1 = Hr4[1];
        const double H00m1m1 = Hr4[4];
        const double H01m1m1 = Hr4[7];
        const double H000m1 = Hr4[13];
        const double H0m11m1 = Hr4[19];
        const double H001m1 = Hr4[22];
        const double H011m1 = Hr4[25];
        const double H0m1m11 = Hr4[55];
        const double H00m11 = Hr4[58];
        const double H01m11 = Hr4[61];
        const double H0m101 = Hr4[64];
        const double H0001 = Hr4[67];
        const double H0m111 = Hr4[73];
        const double H0011 = Hr4[76];
        const double H0111 = Hr4[79];

        //  weight 5
        const double H00m1m1m1 = Hr5[4];
        const double H0m10m1m1 = Hr5[10];
        const double H000m1m1 = Hr5[13];
        const double H00m10m1 = Hr5[31];
        const double H0000m1 = Hr5[40];
        const double H0010m1 = Hr5[49];
        const double H0001m1 = Hr5[67];
        const double H000m11 = Hr5[175];
        const double H0m1m101 = Hr5[190];
        const double H00m101 = Hr5[193];
        const double H00001 = Hr5[202];
        const double H00101 = Hr5[211];
        const double H00011 = Hr5[229];
        const double H01011 = Hr5[232];
        const double H00111 = Hr5[238];
        const double H01111 = Hr5[241];

        delete[] Hr1;
        delete[] Hr2;
        delete[] Hr3;
        delete[] Hr4;
        delete[] Hr5;

        wx = 1 - 2 * x;
        double *tildeHr1 = new double[sz];
        double *tildeHr2 = new double[sz * sz];
        double *tildeHr3 = new double[sz * sz * sz];
        double *tildeHr4 = new double[sz * sz * sz * sz];
        double *tildeHr5 = new double[sz * sz * sz * sz * sz];

        apf_hplog_(
            &wx, &nw, tildeHr1, tildeHr2, tildeHr3, tildeHr4, tildeHr5, &n1, &n2
        );

        const double tildeH0m1 = tildeHr2[1];
        const double tildeH01 = tildeHr2[7];

        const double tildeH0m1m1 = tildeHr3[1];
        const double tildeH01m1 = tildeHr3[7];
        const double tildeH0m11 = tildeHr3[19];
        const double tildeH011 = tildeHr3[25];

        const double tildeH0m1m1m1 = tildeHr4[1];
        const double tildeH01m1m1 = tildeHr4[7];
        const double tildeH0m11m1 = tildeHr4[19];
        const double tildeH011m1 = tildeHr4[25];
        const double tildeH0m1m11 = tildeHr4[55];
        const double tildeH01m11 = tildeHr4[61];
        const double tildeH0m111 = tildeHr4[73];
        const double tildeH0111 = tildeHr4[79];

        const double tildeH0m11m1m1 = tildeHr5[19];
        const double tildeH011m1m1 = tildeHr5[25];
        const double tildeH0m1m11m1 = tildeHr5[55];
        const double tildeH01m11m1 = tildeHr5[61];
        const double tildeH0m111m1 = tildeHr5[73];
        const double tildeH0111m1 = tildeHr5[79];
        const double tildeH0m1m1m11 = tildeHr5[163];
        const double tildeH01m1m11 = tildeHr5[169];
        const double tildeH0m1m111 = tildeHr5[217];
        const double tildeH01m111 = tildeHr5[223];
        const double tildeH0m1111 = tildeHr5[235];
        const double tildeH01111 = tildeHr5[241];

        double Li_412 = 0.5174790617; //=Li_4(1/2)=H_{0,0,0,1}(1/2)
        double Li_512 = 0.5084005792; //=Li_5(1/2)=H_{0,0,0,0,1}(1/2)

        delete[] tildeHr1;
        delete[] tildeHr2;
        delete[] tildeHr3;
        delete[] tildeHr4;
        delete[] tildeHr5;

        double H0_2 = H0 * H0;
        double H0_3 = H0_2 * H0;
        double H0_4 = H0_3 * H0;
        double H0_5 = H0_4 * H0;

        double H1_2 = H1 * H1;
        double H1_3 = H1_2 * H1;
        double H1_4 = H1_3 * H1;

        double ln2_2 = ln2 * ln2;
        double ln2_3 = ln2_2 * ln2;
        double ln2_4 = ln2_3 * ln2;
        double ln2_5 = ln2_4 * ln2;

        double B_4 = -4. * zeta2 * ln2_2 + 2. / 3 * ln2_4 - 13. / 2 * zeta4
                     + 16 * Li_412;

        double aQqPS30 = (
            // This contribution has been already introduced in
            // D2_ps3_highscale
            /*CF * nf * TR * TR * (
                16. / 81 * (x - 1.) / x * (4. * x2 + 7 * x + 4) * H1_3
                + 64./729 * (x - 1.) / x * (5165 * x2 - 2632 * x + 1736)
                - 64./243 * (1350 * x2 - 569 * x - 218) * H0
                + 64./81 * (111 * x2 + 64 * x + 100) * H0_2
                - 128./81 * (6 * x2 + 4 * x - 5) * H0_3
                + 64./27 * (x + 1.) * H0_4
                + (
                    64./9 * (x - 1.) / x * (4. * x2 + 7. * x + 4.) * H0_2
                    - 64./27 * (x - 1.) / x * (74. * x2 - 43. * x + 20) * H0
                    + 64./81 * (x - 1.) / x * (150. * x2 + 103. * x + 60.)
                ) * H1
                + (
                    -32./9 * (x - 1.) / x * (4. * x2 + 7. * x + 4) * H0
                    + 16./81 * (x - 1.) / x * (74. * x2 - 43. * x + 20)
                ) * H1_2
                + (
                    128./9 * 1./x * (2. * x3 + x2 -2 * x + 4.) * H0
                    + 64./81 * 1./x * (111. * x3 - 415. * x2 + 89. * x - 60)
                    - 128./3 * (x + 1.) * H0_2
                ) * H01
                + (
                    -128. / 27 * (2. * x + 3.) / x * (9. * x2 - 8. * x + 4.)
                    + 256./3 * (x + 1.) * H0
                ) * H001
                + (
                    64./27 * 1./x * (6. * x3 + 5. * x2 - 4. * x - 12.)
                    + 128./3 * (x + 1.) * H0
                ) * H011
                - 256./9 * (x + 1.) * H0001 - 640./9 * (x + 1.) * H0011
                - 64./9 * (x + 1.) * H0111
                + (
                    80./9 * (x - 1.) / x * (4. * x2 + 7. * x + 4.) * H1
                    + 16./81 * 1./x * (666. * x3 - 95. * x2 + 589. * x - 60.)
                    - 224./27 * (6. * x2 + 4. * x - 5.) * H0
                    + 224./9 * (x + 1.) * H0_2 - 160./3 * (x + 1) * H01
                ) * zeta2
                + (
                    32./27 * 1./x * (104 * x3 + 67 * x2 - 89 * x + 28)
                    - 320./3 * (x + 1) * H0
                ) * zeta3 + 560./3 * (x + 1.) * zeta4
            )*/

            +CF * TR * TR
                * (-32. / 81 * (x - 1.) / x * (4. * x2 + 7. * x + 4.) * H1_3
                   - 128. / 1215
                         * (18. * x4 - 171. * x3 + 3006. * x2 + 3502. * x + 775.
                         )
                         * H0
                   - 128. / 3645 * (x - 1.) / x
                         * (108. * x4 - 918. * x3 - 13889. * x2 + 145. * x
                            - 3035.)
                   + 64. / 405
                         * (6 * x5 - 60 * x4 + 30 * x3 + 630 * x2 - 985 * x
                            + 320)
                         * H0_2
                   - 128. / 81 * (6 * x2 + 4 * x - 5) * H0_3
                   + 64. / 27 * (x + 1) * H0_4
                   + (64. / 405 * (x - 1.) / x
                          * (12 * x4 - 102 * x3 - 698 * x2 - 255 * x - 290)
                      + 64. / 135 * (x - 1) / x
                            * (4 * x5 - 36 * x4 - 16 * x3 - 156 * x2 - 431 * x
                               - 100)
                            * H0)
                         * H1
                   + (-32. / 9 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H0
                      - 32. / 405 * (x - 1) / x
                            * (12 * x5 - 108 * x4 - 48 * x3 + 172 * x2
                               - 2183 * x - 200))
                         * H1_2
                   + (128. / 9 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H1
                      - 256. / 405 * 1. / x
                            * (6 * x6 - 60 * x5 + 30 * x4 - 1030 * x2 + 122 * x
                               + 75)
                      + 128. / 9 * (2 * x2 + 11 * x + 8) * H0)
                         * H01
                   - 128. / 3 * (x + 1.) * H01 * H01
                   - 128. / 27 * (6 * x2 + 62 * x + 53) * H001
                   + (-128. / 27 * (3 * x - 1) / x * (4 * x2 + 19 * x + 18)
                      + 128. / 3 * (x + 1) * H0)
                         * H011
                   - 256. / 9 * (x + 1) * H0001 + 128. / 9 * (x + 1) * H0011
                   + 128. / 9 * (x + 1) * H0111
                   + (-32. / 3 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H1
                      + 32. / 405 * 1. / x
                            * (24 * x6 - 240 * x5 + 120 * x4 + 2250 * x3
                               - 7265 * x2 - 1145 * x - 600)
                      - 32. / 27 * (60 * x2 + 127 * x + 37) * H0
                      + 224. / 9 * (x + 1) * H0_2 + 64 * (x + 1) * H01)
                         * zeta2
                   + (64. / 27 * 1. / x * (64 * x3 + 251 * x2 + 155 * x - 64)
                      - 128 * (x + 1) * H0)
                         * zeta3
                   - 128. / 3 * (x + 1) * zeta4)

            + CF * CF * TR
                  * (2. / 27 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H1_4
                     - 8. / 81 * (x - 1) / x * (5400 * x2 + 3811 * x + 4614)
                     + 4. / 81 * (4376 * x3 + 5311 * x2 - 9879 * x + 840) * H0_2
                           / (1 - x)
                     + 4. / 81 * (11500 * x2 - 1187 * x + 3989) * H0
                     - 2. / 27 * (704 * x2 - 313 * x + 151) * H0_3
                     + 4. / 27 * (76 * x2 + 15 * x + 33) * H0_4
                     - 6. / 5 * (x + 1) * H0_5
                     + (-80. / 27 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H0_3
                        + 16. / 81 * (x - 1) / x * (67 * x2 - 1320 * x - 1490)
                        + 4. / 27 * (x - 1) / x * (904 * x2 + 2683 * x - 392)
                              * H0_2
                        - 16. / 81 * (x - 1) / x * (1312 * x2 + 4357 * x - 722)
                              * H0)
                           * H1
                     + (-16. / 9 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H0_2
                        - 8. / 27 * (x - 1) / x * (122 * x2 - 55 * x - 40) * H0
                        + 8. / 81 * (x - 1) / x * (913 * x2 + 3928 * x - 320))
                           * H1_2
                     + (-(x - 1) / x * 16. / 27 * (4 * x2 + 7 * x + 4) * H0
                        + (x - 1) / x * 4. / 27 * (40 * x2 + 41 * x + 4))
                           * H1_3
                     + (-32. / 3 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H1_2
                        - 16. / 9 * 1. / x * (4 * x3 + 225 * x2 + 54 * x + 20)
                              * H0_2
                        - 16. / 27 * 1. / x
                              * (882 * x3 + 1527 * x2 - 2364 * x + 196) * H0
                        + 16. / 81 * 1. / x
                              * (3293 * x3 + 8316 * x2 - 5229 * x + 722)
                        + 160. / 9 * (x + 1) * H0_3
                        + (-448. / 9 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H0
                           + 64. / 9 * 1. / x
                                 * (2 * x3 - 102 * x2 + 102 * x - 11))
                              * H1)
                           * H01
                     + (-16. / 9 * 1. / x * (104 * x3 - 39 * x2 - 105 * x - 56)
                        + 448. / 3 * (x + 1) * H0)
                           * H01 * H01
                     + (896. / 9 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H1
                        - 32. / 9 * 1. / x * (28 * x3 - 537 * x2 - 123 * x - 20)
                              * H0
                        + 8. / 27 * 1. / x
                              * (1652 * x3 + 3273 * x2 - 7665 * x + 392)
                        - 48 * (x + 1) * H0_2 - 1792. / 3 * (x + 1) * H01)
                           * H001
                     + (608. / 9 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H1
                        + 32. / 3 * 1. / x * (36 * x3 + 23 * x2 - 32 * x - 40)
                              * H0
                        + 32. / 27 * 1. / x
                              * (209 * x3 + 1743 * x2 - 1572 * x + 152)
                        + 64. / 3 * (x + 1) * H0_2 + 128 * (x + 1) * H01)
                           * H011
                     + (32. / 9 * 1. / x * (156 * x3 - 876 * x2 - 255 * x - 20)
                        + 640. / 3 * (x + 1) * H0)
                           * H0001
                     + (32. / 9 * 1. / x * (8 * x3 - 372 * x2 - 81 * x + 120)
                        - 1792. / 3 * (x + 1) * H0)
                           * H0011
                     + (-16. / 9 * 1. / x
                            * (300 * x3 + 243 * x2 - 219 * x - 304)
                        + 64. / 3 * (x + 1) * H0)
                           * H0111
                     - 2560. / 3 * (x + 1) * H00001 + 3392 * (x + 1) * H00011
                     + 4288. / 3 * (x + 1) * H00101 - 832 * (x + 1) * H00111
                     - 1600. / 3 * (x + 1) * H01011 - 32. / 3 * (x + 1) * H01111
                     + (124. / 9 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H1_2
                        - 4. / 81 * 1. / x
                              * (8896 * x3 + 21003 * x2 + 129 * x - 1620)
                        + 2. / 9 * (1536 * x2 + 1879 * x + 1943) * H0
                        - 8. / 9 * (68 * x2 + 99 * x - 57) * H0_2
                        + 188. / 9 * (x + 1) * H0_3
                        + (112. / 9 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H0
                           + 4. / 27 * 1. / x
                                 * (752 * x3 + 4197 * x2 - 5169 * x + 652))
                              * H1
                        + (8. / 9 * 1. / x
                               * (196 * x3 - 327 * x2 - 213 * x + 56)
                           - 224. / 3 * (x + 1) * H0)
                              * H01
                        - 608. / 3 * (x + 1) * H001 - 496. / 3 * (x + 1) * H011)
                           * zeta2
                     + 1024. / 3 * zeta2 * ln2 * ln2 * ln2 * (x + 1)
                     + (-592. / 9 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H1
                        - 4. / 27 * 1. / x
                              * (2824 * x3 - 4125 * x2 - 17331 * x + 1938)
                        - 16. / 3 * (76 * x2 - 275 * x - 112) * H0
                        + 920. / 3 * (x + 1) * H0_2 + 1184. / 3 * (x + 1) * H01)
                           * zeta3
                     + 1792 * (x + 1) * zeta3 * ln2 * ln2
                     + (-896. / 3 * 1. / x * (4 * x3 - 9 * x2 - 6 * x - 2)
                        + 1792 * (x + 1) * H0)
                           * zeta3 * ln2
                     + (-4. / 9 * 1. / x
                            * (628 * x3 - 3669 * x2 - 351 * x + 540)
                        + 928. / 3 * (x + 1) * H0)
                           * zeta4
                     + 1664 * (x + 1) * zeta4 * ln2
                     + (-32. / 3 * 1. / x * (8 * x3 - 30 * x2 - 15 * x - 2)
                        + 128 * (x + 1) * H0)
                           * B_4
                     + 256 * (x + 1) * B_4 * ln2
                     + 944. / 3 * (x + 1) * zeta2 * zeta3
                     - 3544 * (x + 1) * zeta5 + 4096 * (x + 1) * Li_512
                     - 512. / 15 * ln2_5 * (x + 1))

            + CA * CF * TR
                  * (-2. / 27 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H1_4
                     - 16. / 9 * (x + 1) / x * (32 * x2 + 61 * x + 32)
                           * (H0m1m11 + H0m11m1 + H01m1m1)
                     - 8. / 729 * (x - 1) / x
                           * (1024228 * x2 - 83309 * x + 274870)
                     + (-112. / 27 * (x + 1) / x * (4 * x2 - 7 * x + 4) * Hm1
                            * Hm1 * Hm1
                        + 32. / 27 * (x + 1) / x * (188 * x2 - 83 * x + 80)
                              * Hm1 * Hm1
                        - 8. / 27 * (x + 1) / x * (4200 * x2 - 1577 * x + 1236)
                              * Hm1
                        + 4. / 243 * 1. / x
                              * (503464 * x3 + 110993 * x2 + 171290 * x + 20992)
                       ) * H0
                     + (-2. / 3 * (x + 1) / x * (40 * x2 - 37 * x + 40) * Hm1
                            * Hm1
                        + 4. / 27 * 1. / x * (202 * x3 - 18 * x2 + 9 * x + 256)
                              * Hm1)
                           * H0_2
                     + 4. / 81 * 1. / (x + 1)
                           * (22712 * x4 - 5914 * x3 - 9627 * x2 + 5671 * x
                              - 13652)
                           * H0_2 / (1 - x)
                     + (4. / 81 * (1290 * x2 + 1213 * x + 1045)
                        + 28. / 9 * (x + 1) / x * (8 * x2 - 11 * x + 8) * Hm1)
                           * H0_3
                     + 2. / 27 * (95 * x - 172) * H0_4
                     - 4. / 15 * (4 * x - 5) * H0_5
                     + (4. / 27 * (x - 1) / x * (152 * x2 + 203 * x + 152)
                            * H0_3
                        - 4. / 243 * (x - 1) / x
                              * (6776 * x2 + 15425 * x - 11926)
                        + 8. / 81 * (x - 1) / x * (21512 * x2 - 1057 * x + 7436)
                              * H0
                        - 4. / 9 * 1. / x
                              * (936 * x3 - 640 * x2 + 379 * x - 684) * H0_2)
                           * H1
                     + (-2. / 9 * (x - 1) / x * (136 * x2 + 283 * x + 136)
                            * H0_2
                        + 8. / 27 * (x - 1) / x * (481 * x2 + 340 * x + 184)
                              * H0
                        - 4. / 81 * (x - 1) / x * (1754 * x2 + 4169 * x - 658))
                           * H1_2
                     + (64. / 27 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H0
                        - 4. / 81 * (x - 1) / x * (154 * x2 + 163 * x + 46))
                           * H1_3
                     + (-32. / 9 * (x - 1) / x * (4 * x2 - 11 * x + 4) * H0 * H1
                        + 112. / 9 * (x + 1) / x * (4 * x2 - 7 * x + 4) * Hm1
                              * Hm1
                        - 64. / 27 * (x + 1) / x * (188 * x2 - 83 * x + 80)
                              * Hm1
                        + 8. / 27 * (x + 1) / x * (4200 * x2 - 1577 * x + 1236)
                        - 4. / 9 * 1. / x * (112 * x3 - 27 * x2 - 57 * x + 168)
                              * H0_2
                        + (-8. / 9 * (x + 1) / x * (40 * x2 - 151 * x + 40)
                               * Hm1
                           + 8. / 27 * 1. / x
                                 * (1102 * x3 - 2154 * x2 + 765 * x - 256))
                              * H0
                        + 280. / 9 * (x - 1) * H0_3)
                           * H0m1
                     + (-8. / 9 * 1. / x * (8 * x3 + 81 * x2 - 87 * x + 80)
                        - 152. / 3 * (x - 1) * H0)
                           * H0m1 * H0m1
                     + (-8. / 9 * (x + 1) / x * (32 * x2 + 61 * x + 32) * Hm1
                            * Hm1
                        - 16. / 27 * 1. / x * (82 * x3 - 78 * x2 - 267 * x - 80)
                              * Hm1
                        - 4. / 9 * 1. / x
                              * (144 * x3 - 555 * x2 - 471 * x - 344) * H0_2
                        - 8. / 81 * 1. / (x + 1) / x
                              * (20970 * x4 + 2819 * x3 - 10430 * x2 + 857 * x
                                 - 6540)
                        + (8. / 9 * 1. / x
                               * (154 * x3 - 1068 * x2 - 217 * x - 844)
                           + 64 * (x + 1) / x * (2 * x2 - 3 * x + 2) * Hm1)
                              * H0
                        - 248. / 9 * (x + 1) * H0_3
                        + (8. / 9 * (x - 1) / x * (184 * x2 + 349 * x + 184)
                               * H0
                           - 16. / 27 * 1. / x
                                 * (167 * x3 - 711 * x2 + 657 * x - 140))
                              * H1
                        + (-32. / 9 * (x - 1) / x * (4 * x2 - 11 * x + 4)
                           - 64. / 3 * (x + 1) * H0)
                              * H0m1)
                           * H01
                     + (4. / 9 * 1. / x * (224 * x3 - 201 * x2 - 105 * x - 112)
                        - 392. / 3 * (x + 1) * H0)
                           * H01 * H01
                     + (-224. / 9 * (x + 1) / x * (4 * x2 - 7 * x + 4) * Hm1
                        + 64. / 27 * (x + 1) / x * (188 * x2 - 83 * x + 80)
                        + 8. / 9 * 1. / x * (56 * x3 + 51 * x2 - 285 * x + 200)
                              * H0
                        - 152. / 3 * (x - 1) * H0_2 + 448. / 3 * (x - 1) * H0m1)
                           * H0m1m1
                     + (16. / 9 * (x + 1) / x * (32 * x2 + 61 * x + 32) * Hm1
                        - 32. / 9 * 1. / x * (32 * x3 - 3 * x2 - 33 * x + 40)
                              * H0
                        + 16. / 27 * 1. / x * (82 * x3 - 78 * x2 - 267 * x - 80)
                       ) * H0m11
                     + (64. / 9 * (x - 1) / x * (4 * x2 - 11 * x + 4) * H1
                        + 8. / 9 * (x + 1) / x * (200 * x2 - 413 * x + 200)
                              * Hm1
                        - 8. / 3 * 1. / x * (8 * x3 - 9 * x2 + 27 * x - 56) * H0
                        - 8. / 27 * 1. / x
                              * (2406 * x3 - 4326 * x2 + 1539 * x - 256)
                        - 32. / 3 * (15 * x - 17) * H0_2 + 304 * (x - 1) * H0m1
                        + 128. / 3 * (x + 1) * H01)
                           * H00m1
                     + (-8. / 3 * (x + 1) / x * (88 * x2 - 169 * x + 88) * Hm1
                        - 8. / 9 * (x - 1) / x * (296 * x2 + 491 * x + 296) * H1
                        + 8. / 9 * 1. / x
                              * (136 * x3 - 1605 * x2 - 591 * x - 536) * H0
                        + 8. / 27 * 1. / x
                              * (1936 * x3 + 3913 * x2 + 3019 * x + 3012)
                        + 16. / 3 * (69 * x + 25) * H0_2
                        - 1136. / 3 * (x - 1) * H0m1 + 1136. / 3 * (x + 1) * H01
                       ) * H001
                     + (16. / 9 * (x + 1) / x * (32 * x2 + 61 * x + 32) * Hm1
                        - 32. / 9 * 1. / x * (32 * x3 - 3 * x2 - 33 * x + 40)
                              * H0
                        + 16. / 27 * 1. / x * (82 * x3 - 78 * x2 - 267 * x - 80)
                       ) * H01m1
                     + (-176. / 9 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H1
                        - 16. / 9 * 1. / x
                              * (148 * x3 + 189 * x2 - 45 * x - 140) * H0
                        - 8. / 27 * 1. / x
                              * (356 * x3 + 3725 * x2 - 3469 * x + 272)
                        + 64. / 9 * (x + 1) / x * (2 * x2 + x + 2) * Hm1
                        + 104 * (x + 1) * H0_2)
                           * H011
                     + (224. / 9 * (x + 1) / x * (4 * x2 - 7 * x + 4)
                        - 448. / 3 * (x - 1) * H0)
                           * H0m1m1m1
                     + (8. / 9 * 1. / x * (136 * x3 - 183 * x2 - 75 * x + 104)
                        + 64. / 3 * (9 * x - 7) * H0)
                           * H0m101
                     - 64 * (x + 1) * (2 * x2 + x + 2) / (9 * x) * H0m111
                     + (-8. / 9 * (x + 1) / x * (200 * x2 - 413 * x + 200)
                        + 320 * (x - 1) * H0)
                           * H00m1m1
                     + (8. / 3 * (x + 1) / x * (88 * x2 - 169 * x + 88)
                        + 128. / 3 * (x + 1) * H0)
                           * H00m11
                     + (8. / 3 * 1. / x * (80 * x3 - 33 * x2 + 45 * x - 56)
                        + 1168. / 3 * (x - 1) * H0)
                           * H000m1
                     + (-16. / 9 * 1. / x
                            * (68 * x3 - 1577 * x2 - 260 * x - 364)
                        - 496 * (3 * x + 1) * H0)
                           * H0001
                     + (8. / 3 * (x + 1) / x * (88 * x2 - 169 * x + 88)
                        + 128. / 3 * (x + 1) * H0)
                           * H001m1
                     + (16. / 9 * 1. / x * (40 * x3 + 617 * x2 + 119 * x - 164)
                        + 16. / 3 * (35 * x + 37) * H0)
                           * H0011
                     - 64. / 9 * (x + 1) / x * (2 * x2 + x + 2)
                           * (H01m11 + H011m1)
                     + (16. / 9 * 1. / x * (104 * x3 + 101 * x2 - 64 * x - 104)
                        - 256. / 3 * (x + 1) * H0)
                           * H0111
                     + 160. / 3 * (x - 1) * H0m1m101
                     - 896. / 3 * (x - 1) * H0m10m1m1
                     - 1792. / 3 * (x - 1) * H00m1m1m1
                     - 624 * (x - 1) * H00m10m1
                     + 32. / 3 * (33 * x - 41) * H00m101
                     - 1872 * (x - 1) * H000m1m1 + 16 * (63 * x - 79) * H000m11
                     - 128 * (3 * x - 1) * H0000m1
                     + 16. / 3 * (411 * x + 197) * H00001
                     + 16 * (63 * x - 79) * H0001m1
                     - 16. / 3 * (341 * x + 347) * H00011
                     + 16. / 3 * (63 * x - 79) * H0010m1
                     - 32. / 3 * (68 * x + 67) * H00101 + 160 * (x + 1) * H00111
                     + 352. / 3 * (x + 1) * H01011 + 32. / 3 * (x + 1) * H01111
                     + (16. / 9 * (x + 1) / x * (2 * x2 + 55 * x + 2) * Hm1
                            * Hm1
                        - 20. / 9 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H1_2
                        + 16. / 27 * 1. / x
                              * (458 * x3 + 132 * x2 - 273 * x + 80) * Hm1
                        - 8. / 81 * 1. / (x + 1) / x
                              * (5729 * x4 + 3096 * x3 + 1145 * x2 + 4805 * x
                                 + 703)
                        + (-8. / 9 * (x + 1) / x * (148 * x2 - 97 * x + 148)
                               * Hm1
                           + 4. / 27 * 1. / x
                                 * (1318 * x3 - 2183 * x2 - 275 * x + 240))
                              * H0
                        + 4. / 9 * (8 * x2 + 23 * x - 211) * H0_2
                        - 16. / 9 * (5 * x - 4) * H0_3
                        + (-4. / 27 * 1. / x
                               * (214 * x3 + 3147 * x2 - 3579 * x + 326)
                           - (32 * (x - 1) * (x2 + x + 1) * H0) / (3 * x))
                              * H1
                        + (8. / 9 * 1. / x * (76 * x3 + 135 * x2 - 21 * x + 148)
                           - 304. / 3 * (x - 1) * H0)
                              * H0m1
                        + (-8. / 3 * 1. / x * (32 * x3 - 124 * x2 - 55 * x + 24)
                           + 32. / 3 * (x + 1) * H0)
                              * H01
                        - 128 * (x - 1) * H0m1m1 + 544. / 3 * (x - 1) * H00m1
                        + 16. / 3 * (7 * x + 3) * H001
                        + 80. / 3 * (x + 1) * H011)
                           * zeta2
                     - 512. / 3 * (x + 1) * zeta2 * ln2_3
                     + (8. / 9 * (x - 1) / x * (88 * x2 + 235 * x + 88) * H1
                        - 16. / 9 * 1. / x * (52 * x3 + 819 * x2 - 144 * x + 28)
                              * H0
                        - 4. / 27 * 1. / x
                              * (4612 * x3 + 15262 * x2 + 8524 * x + 2559)
                        + (16 * (x - 4) * (x + 1) * (4 * x - 1) * Hm1) / x
                        + 64. / 3 * (5 * x - 6) * H0_2
                        + 608. / 3 * (x - 1) * H0m1 - 496. / 3 * (x + 1) * H01)
                           * zeta3
                     - 896 * (x + 1) * zeta3 * ln2_2
                     + (448. / 3 * 1. / x * (4 * x3 - 9 * x2 - 6 * x - 2)
                        - 896 * (x + 1) * H0)
                           * zeta3 * ln2
                     + (-2. / 9 * 1. / x
                            * (1752 * x3 + 11325 * x2 + 1401 * x + 1828)
                        - 76. / 3 * (17 * x - 15) * H0)
                           * zeta4
                     - 832 * (x + 1) * zeta4 * ln2
                     + (16. / 3 * 1. / x * (8 * x3 - 30 * x2 - 15 * x - 2)
                        - 64 * (x + 1) * H0)
                           * B_4
                     - 128 * (x + 1) * B_4 * ln2
                     + 8. / 3 * (31 * x - 127) * zeta2 * zeta3
                     - 12 * (47 * x - 145) * zeta5 - 2048 * Li_512 * (x + 1)
                     + 256. / 15 * (x + 1) * ln2_5)
        );

        double tildeaQqPS30 =
            (CF * TR * (CA / 2 - CF)
             * (-64 * (x + 1) * H0 * H0 * tildeH0m1m1
                + 64. / 3 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H1
                      * tildeH0m1m1
                + 16. / 3 * 1. / x * (32 * x3 + 200 * x2 - 104 * x + 1)
                      * tildeH0m1m1
                + 64. / 3 * (4 * x2 - 21 * x - 9) * H0 * tildeH0m1m1
                - 128 * (x + 1) * H01 * tildeH0m1m1
                + (-64. / 3 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H1
                   - 16. / 3 * 1. / x * (32 * x3 + 200 * x2 - 104 * x + 1)
                   - 64. / 3 * (4 * x2 - 21 * x - 9) * H0 + 64 * (x + 1) * H0_2
                   + 128 * (x + 1) * H01)
                      * tildeH0m11
                + (64. / 3 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H1
                   + 16. / 3 * 1. / x * (32 * x3 + 200 * x2 - 104 * x + 1)
                   + 64. / 3 * (4 * x2 - 21 * x - 9) * H0 - 64 * (x + 1) * H0_2
                   - 128 * (x + 1) * H01)
                      * tildeH01m1
                + (-64. / 3 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H1
                   - 16. / 3 * 1. / x * (32 * x3 + 200 * x2 - 104 * x + 1)
                   - 64. / 3 * (4 * x2 - 21 * x - 9) * H0 + 64 * (x + 1) * H0_2
                   + 128 * (x + 1) * H01)
                      * tildeH011
                + (64 * (x - 1) * (4 * x2 + 7 * x + 4) / x - 384 * (x + 1) * H0)
                      * tildeH0m1m1m1
                + (-64. / 3 * 1. / x * (4 * x3 + 27 * x2 + 3 * x - 8)
                   + 128 * (x + 1) * H0)
                      * tildeH0m1m11
                + (64. / 3 * (4 * x2 - 21 * x - 9) - 128 * (x + 1) * H0)
                      * tildeH0m11m1
                + (-64. / 3 * 1. / x * (12 * x3 - 39 * x2 - 21 * x - 4)
                   + 384 * (x + 1) * H0)
                      * tildeH0m111
                + (64. / 3 * 1. / x * (12 * x3 - 15 * x2 - 15 * x - 8)
                   - 384 * (x + 1) * H0)
                      * tildeH01m1m1
                + (-64. / 3 * 1. / x * (x - 1) * (4 * x2 + 7 * x + 4)
                   + 128 * (x + 1) * H0)
                      * tildeH01m11
                + (64. / 3 * 1. / x * (4 * x3 - 45 * x2 - 15 * x + 4)
                   - 128 * (x + 1) * H0)
                      * tildeH011m1
                + (-64 * (4 * x2 - 21 * x - 9) + 384 * (x + 1) * H0)
                      * tildeH0111
                - 384 * (x + 1) * tildeH0m1m1m11
                - 256 * (x + 1) * tildeH0m1m11m1 + 384 * (x + 1) * tildeH0m1m111
                - 128 * (x + 1) * tildeH0m11m1m1 - 384 * (x + 1) * tildeH0m111m1
                + 768 * (x + 1) * tildeH0m1111 - 384 * (x + 1) * tildeH01m1m11
                - 256 * (x + 1) * tildeH01m11m1 + 384 * (x + 1) * tildeH01m111
                - 128 * (x + 1) * tildeH011m1m1 - 384 * (x + 1) * tildeH0111m1
                + 768 * (x + 1) * tildeH01111
                + 64 * (x + 1)
                      * (tildeH0m1m1 - tildeH0m11 + tildeH01m1 - tildeH011)
                      * zeta2
                - 128 * (x + 1) * (tildeH0m1 + tildeH01) * zeta2 * ln2
                + (-128. / 3 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H1
                       * tildeH0m1
                   - 32. / 3 * 1. / x * (32 * x3 + 200 * x2 - 104 * x + 1)
                         * tildeH0m1
                   - 128. / 3 * (4 * x2 - 21 * x - 9) * H0 * tildeH0m1
                   + 128 * (x + 1) * H0_2 * tildeH0m1
                   + (-128. / 3 * (x - 1) / x * (4 * x2 + 7 * x + 4) * H1
                      - 32. / 3 * 1. / x * (32 * x3 + 200 * x2 - 104 * x + 1)
                      - 128. / 3 * (4 * x2 - 21 * x - 9) * H0
                      + 128 * (x + 1) * H0_2)
                         * tildeH01
                   + (256 * (x + 1) * tildeH0m1 + 256 * (x + 1) * tildeH01)
                         * H01
                   + (-128. / 3 * 1. / x * (8 * x3 + 18 * x2 - 3 * x - 10)
                      + 512 * (x + 1) * H0)
                         * tildeH0m1m1
                   + (-128. / 3 * 1. / x * (8 * x3 - 30 * x2 - 15 * x - 2)
                      + 512 * (x + 1) * H0)
                         * tildeH0m11
                   + (-128. / 3 * 1. / x * (8 * x3 - 6 * x2 - 9 * x - 6)
                      + 512 * (x + 1) * H0)
                         * tildeH01m1
                   + (-128. / 3 * 1. / x * (8 * x3 - 54 * x2 - 21 * x + 2)
                      + 512 * (x + 1) * H0)
                         * tildeH011
                   - 384 * (x + 1) * tildeH0m1m1m1
                   + 640 * (x + 1) * tildeH0m1m11 + 128 * (x + 1) * tildeH0m11m1
                   + 1152 * (x + 1) * tildeH0m111 - 384 * (x + 1) * tildeH01m1m1
                   + 640 * (x + 1) * tildeH01m11 + 128 * (x + 1) * tildeH011m1
                   + 1152 * (x + 1) * tildeH0111)
                      * ln2
                + (256 * (12 * x2 + 3 * x - 2) / (3 * x)
                       * (tildeH0m1 + tildeH01)
                   + 512 * (x + 1)
                         * (tildeH0m1m1 + tildeH0m11 + tildeH01m1 + tildeH011))
                      * ln2_2));

        return aQqPS30 + tildeaQqPS30;
    } else if (v == 1) {
        double L1 = log(1. - x);
        double L = log(x);
        return (
            (1. - x)
                * (232.9555 * L1 * L1 * L1 + 1309.528 * L1 * L1 - 31729.716 * x2
                   + 66638.193 * x + 2825.641 / x)
            + 41850.518 * x * L + 688.396 / x * L
        );
    } else if (v == -1) {
        double L1 = log(1. - x);
        double L = log(x);
        return (
            (1. - x)
                * (126.3546 * L1 * L1 + 353.8539 * L1 + 6787.608 * x
                   + 3780.192 / x)
            + 8571.165 * x * L - 2346.893 * L * L + 688.396 / x * L
        );
    } else {
        std::cout << "a_Qq_PS_30: Choose either v=0, v=1 or v=-1!!\nExiting!!\n"
                  << std::endl;
        exit(-1);
    }
}
