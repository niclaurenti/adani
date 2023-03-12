#include "adani/MasslessCoefficientFunctions.h"
#include "adani/Constants.h"
#include<cmath>


//==========================================================================================//
//  Massless gluon coefficient functions for F2 at O(alpha_s) for mu=Q.
//
//  Eq. (4.4) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double C2_g1_massless(double x, int nf) {

    if (x>=1 || x<0) return 0;

    return 4. * nf * TR * (
        - 8. * x * x + 8. * x - 1.
        + log((1. - x) / x) * (2. * x * x - 2. * x + 1.)
    ) ;

}

//==========================================================================================//
//  Massless gluon coefficient functions for FL at O(alpha_s) for mu=Q.
//
// Eq. (3) of Ref. [arXiv:hep-ph/0411112v2]
//------------------------------------------------------------------------------------------//

double CL_g1_massless(double x, int nf) {

    if (x>=1 || x<0) return 0;

    return 16. * nf * TR * x * (1. - x) ;

}

//==========================================================================================//
//  Massless gluon coefficient functions for F2 at O(alpha_s^2) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long) result
//
//  Eq. (4.10) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double C2_g2_massless(double x, int nf) {

    if (x>=1 || x<0) return 0;

    double x1 = 1. - x;
    double L0 = log(x);
    double L1 = log(x1);

    return nf * (
        58./9. * L1 * L1 * L1 - 24 * L1 * L1 - 34.88 * L1 + 30.586
        - (25.08 + 760.3 * x + 29.65 * L1 * L1 * L1) * x1 + 1204 * x * L0 * L0
        + L0 * L1 * (293.8 + 711.2 * x + 1043 * L0)
        + 115.6 * L0 - 7.109 * L0 * L0
        + 70./9. * L0 * L0 * L0 + 11.9033 * x1 / x
    ) ;

}

//==========================================================================================//
//  Massless quark coefficient functions for F2 at O(alpha_s^2) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long) result
//
//  Eq. (4.9) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double C2_ps2_massless(double x, int nf) {

    if (x>=1 || x<0) return 0;

    double x1 = 1 - x ;
    double L0 = log(x);
    double L02 = L0 * L0 ;
    double L03 = L02 * L0 ;

    double L1 = log(x1);
    double L12 = L1 * L1 ;
    double L13 = L12 * L1 ;

    return nf * (
        (8./3 * L12 - 32./3 * L1 + 9.8937) * x1
        + (9.57 - 13.41 * x + 0.08 * L13) * x1 * x1 + 5.667 * x * L03
        - L02 * L1 * (20.26 - 33.93 * x) + 43.36 * x1 * L0 - 1.053 * L02
        + 40./9 * L03 + 5.2903 / x * x1 * x1
    ) ;

}

//==========================================================================================//
//  Massless gluon coefficient functions for FL at O(alpha_s^2) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long) result
//
//  Eq. (6) of Ref. [arXiv:hep-ph/0411112v2]
//------------------------------------------------------------------------------------------//

double CL_g2_massless(double x, int nf) {

    if (x>=1 || x<0) return 0;

    double L0 = log(x);
    double L02 = L0 * L0;

    double x1 = 1. - x ;
    double L1 = log(x1) ;
    double L12 = L1 * L1 ;

    return nf * (
        (94.74 - 49.2 * x) * x1 * L12 + 864.8 * x1 * L1
        + 1161 * x * L1 * L0 + 60.06 * x * L02 + 39.66 * x1 * L0
        - 5.333 * (1./x - 1)
    );

}

//==========================================================================================//
//  Massless gluon coefficient functions for FL at O(alpha_s^2) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long) result
//
//  Eq. (5) of Ref. [arXiv:hep-ph/0411112v2]
//------------------------------------------------------------------------------------------//

double CL_ps2_massless(double x, int nf) {

    if (x>=1 || x<0) return 0;

    double L0 = log(x);
    double L02 = L0 * L0;

    double x1 = 1. - x;
    double L1 = log(x1);

    return nf * (
        (15.94 - 5.212 * x) * x1 * x1 * L1
        + (0.421 + 1.520 * x) * L02 + 28.09 * x1 * L0
        - (2.370/x - 19.27) * x1 * x1 * x1
    ) ;

}

//==========================================================================================//
//  Massless gluon coefficient functions for F2 at O(alpha_s^3) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long) result.
//  The term fl_g_11 is put to zero for the reason explained in page 15 of arXiv:1205.5727
//
//  Eq. (4.13) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

double C2_g3_massless(double x, int nf) {//remember that there is a delta(x1) that has been omitted

    if(x<0 || x>=1) return 0;

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

    double c_nf = (
        966./81. * L15 - 1871./18. * L14 + 89.31 * L13 + 979.2 * L12
        - 2405 * L1 + 1372 * x1 * L14 - 15729 - 310510 * x + 331570 * x2
        - 244150 * x * L02 - 253.3 * x * L05 + L0 * L1 * (138230 - 237010 * L0)
        - 11860 * L0 - 700.8 * L02 - 1440 * L03 + 4961./162. * L04
        - 134./9. * L05 - (6362.54 + 932.089 * L0)/x // there's a typo in the paper: -932.089 -> +932.089
    ) ; //there is + 0.625*delta(x1)

    double c_nf2 = (
        131./81. * L14 - 14.72 * L13 + 3.607 * L12 - 226.1 * L1
        + 4.762 - 190 * x - 818.4 * x2 - 4019 * x * L02
        - L0 * L1 * (791.5 + 4646 * L0) + 739.0 * L0 + 418.0 * L02
        + 104.3 * L03 + 809./81. * L04 + 12./9. * L05 + 84.423 / x
    ) ;

    /*double c_nf_fl = (
        3.211 * L12 + 19.04 * x * L1 + 0.623 * x1 * L13 - 64.47
        * x + 121.6 * x2 - 45.82 * x3 - x * L0 * L1 * (31.68 +
        37.24 * L0) + 11.27 * x2 * L03 - 82.40 * x * L0 - 16.08
        * x * L02 + 520./81. * x * L03 + 20./27. * x * L04
    );*/

    return (
        c_nf * nf
        + c_nf2 * nf * nf
        //+ c_nf_fl * nf * nf * fl_g_11
    );

}

//==========================================================================================//
//  Massless quark coefficient functions for F2 at O(alpha_s^3) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long) result.
//  The term fl_ps_11 is put to zero for the reason explained in page 15 of arXiv:1205.5727
//
//  Eq. (4.12) of Ref. [hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//


double C2_ps3_massless(double x, int nf) {//remember that there is a delta(x1) that has been omitted

    if(x<0 || x>=1) return 0;

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

    double c_nf = (
        (
            856./81 * L14 - 6032./81 * L13 + 130.57 * L12
            - 542 * L1 + 8501 - 4714 * x + 61.5 * x2
        ) * x1
        + L0 * L1 * (8831 * L0 + 4162 * x1)
        - 15.44 * x * L05 + 3333 * x * L02 + 1615 * L0 + 1208 * L02
        - 333.73 * L03 + 4244./81 * L04 - 40./9 * L05
        - 1./x * (2731.82 * x1 + 414.262 * L0)
    );

    double c_nf2 = (
        (
            -64./81 * L13 + 208./81 * L12 + 23.09 * L1 - 220.27
            + 59.80 * x - 177.6 * x2
        ) * x1
        - L0 * L1 * (160.3 * L0 + 135.4 * x1) - 24.14 * x * L03
        - 215.4 * x * L02 - 209.8 * L0 - 90.38 * L02 - 3568./243 * L03
        - 184./81 * L04 + 40.2426 * x1/x
    ) ;

    // double c_fl_nf = (
    //     (126.42 - 50.29 * x - 50.15 * x2) * x1 - 26.717
    //     - 9.075 * x * x1 * L1 - x * L02 * (101.8 + 34.79 * L0 + 3.070 * L02)
    //     + 59.59 * L0 - 320./81 * L02 * (5 + L0)
    // ) * x ;

    return (
        c_nf * nf
        + c_nf2 * nf * nf
        //+ c_fl_nf * fl_ps_11 * nf
    );

}

//==========================================================================================//
//  Massless gluon coefficient functions for FL at O(alpha_s^3) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long) result.
//  The term fl_ps_11 is put to zero for the reason explained in page 15 of arXiv:1205.5727
//
//  Eq. (10) of Ref. [arXiv:hep-ph/0411112v2]
//------------------------------------------------------------------------------------------//

double CL_g3_massless(double x, int nf) {//remember that there is a delta(x1) that has been omitted

    if(x<0 || x>=1) return 0;

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

    double c_nf = (
        (144. * L14 - 47024./27. * L13 + 6319. * L12 + 53160. * L1) * x1
        + 72549. * L0 * L1 + 88238. * L02 * L1
        + (3709. - 33514. * x - 9533. * x2) * x1 + 66773. * x * L02 - 1117. * L0
        + 45.37 * L02 - 5360./27. * L03 - (2044.70 * x1 + 409.506 * L0) / x
    ) ;

    double c_nf2 = (
        (32./3. * L13 - 1216./9. * L12 - 592.3 * L1 + 1511. * x * L1) * x1
        + 311.3 * L0 * L1 + 14.24 * L02 * L1 + (577.3 - 729. * x) * x1
        + 30.78 * x * L03 + 366. * L0 + 1000./9. * L02 + 160./9. * L03
        + 88.5037/x * x1
    ) ;

    /*double c_nf_fl = (
        (
            -0.0105 * L13 + 1.55 * L12 + 19.72 * x * L1 - 66.745 * x
            + 0.615 * x2
        ) * x1
        + 20./27. * x * L04 + (280./81. + 2.26 * x) * x * L03
        - (15.4 - 2.201 * x) * x * L02 - (71.66 - 0.121 * x) * x * L0
    ) ;*/

    return (
        c_nf * nf
        + c_nf2 * nf * nf
        //+ c_nf_fl * nf * nf * fl_g_11
    ) ;

}

//==========================================================================================//
//  Massless gluon coefficient functions for FL at O(alpha_s^3) for mu=Q.
//  Observe that this result is a parameterization of the exact (known but long) result.
//  The term fl_ps_11 is put to zero for the reason explained in page 15 of arXiv:1205.5727
//
//  Eq. (9) of Ref. [arXiv:hep-ph/0411112v2]
//------------------------------------------------------------------------------------------//

double CL_ps3_massless(double x, int nf) {//remember that there is a delta(x1) that has been omitted

    if(x<0 || x>=1) return 0;

    // double fl_ps_11 = fl11ps(nf) ;

    // double x2=x*x;

    double x1 = 1. - x;
    double L0 = log(x);
    double L1 = log(x1);

    double L02 = L0 * L0;
    double L03 = L02 * L0;

    double L12 = L1 * L1;
    double L13 = L12 * L1;

    double c_nf =  (
        (1568./27 * L13 - 3968./9 * L12 + 5124 * L1) * x1 * x1
        + (2184 * L0 + 6059 * x1) * L0 * L1
        - (795.6 + 1036 * x) * x1 * x1 - 143.6 * x1 * L0 + 2848./9 * L02
        - 1600./27 * L03 - (885.53 * x1 + 182.00 * L0)/x * x1
    );

    double c_nf2 = (
        (-32./9 * L12 + 29.52 * L1) * x1 * x1
        + (35.18 * L0 + 73.06 * x1) * L0 * L1 - 35.24 * x * L02
        - (14.16 - 69.84 * x) * x1 * x1 - 69.41 * x1 * L0
        - 128./9 * L02 + 40.239 /x * x1 * x1
    );

    /*double c_fl_nf = (
        ((107.0 + 321.05 * x - 54.62 * x2) * x1 - 26.717 + 9.773
        * L0 + (363.8 + 68.32 * L0) * x * L0 - 320./81 * L02 * (2
        + L0)) * x
    );*/

    return (
        c_nf * nf
        + c_nf2 * nf * nf
        //+ c_fl_nf * fl_ps_11 * nf
    ) ;

}

//==========================================================================================//
//  Charge factors for the flavour topologies entering up to three loops.
//
//  Table 2 of Ref. [arXiv:hep-ph/0504242v1]
//------------------------------------------------------------------------------------------//

//                    u       d       s      c       b       t
double charges[] = {2./3., -1./3., -1./3., 2./3., -1./3., 2./3.} ;

double fl11g(int nf) {

    double eavg = 0., e2avg = 0. ;

    for(int i=0; i<nf; i++) {
        eavg += charges[i] ;
        e2avg += charges[i] * charges[i] ;
    }
    eavg /= nf ;
    e2avg /= nf ;

    return eavg * eavg / e2avg ;
}

double fl11ps(int nf) {

    double eavg = 0., e2avg = 0. ;

    for(int i=0; i<nf; i++) {
        eavg += charges[i] ;
        e2avg += charges[i] * charges[i] ;
    }
    eavg /= nf ;
    e2avg /= nf ;

    return eavg * eavg / e2avg - 3 * eavg ;
}
