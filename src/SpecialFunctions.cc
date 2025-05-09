#include "adani/SpecialFunctions.h"
#include "adani/Constants.h"
#include <cmath>

#define Li2_1_2 0.5822405265
#define Li3_1_2 0.5372131936

//==========================================================================================//
//  Beta functions.
//------------------------------------------------------------------------------------------//

double beta0(int nf) { return (11. / 3 * CA - 2. / 3 * nf); }

double beta1(int nf) {
    return (34. / 3. * CA * CA)
           + ((-20. / 3. * CA * TR) + (-4. * CF * TR)) * nf;
}

//==========================================================================================//
//  Theta function.
//------------------------------------------------------------------------------------------//

double theta(double x) {

    if (x > 0)
        return 1.;
    else
        return 0.;
}

//==========================================================================================//
//  Dilogarithm. From Marco Bonvini.
//------------------------------------------------------------------------------------------//

double Li2(double x) {

    double x_0 = -0.30;
    double x_1 = 0.25;
    double x_2 = 0.51;

    if (x == 1.)
        return zeta2;

    if (x <= x_0) {
        double temp = log(fabs(1.0 - x));
        return -Li2(-x / (1.0 - x)) - temp * temp / 2;
    }

    else if (x < x_1) {
        double z = -log(1.0 - x);
        // clang-format off
        double temp = (
            z*(1.0-z/4.0*(1.0-z/9.0*(1.0-z*z/100.0
            *(1.0-5.0*z*z/294.0*(1.0-7.0*z*z/360.0
            *(1.0-5.0*z*z/242.0*(1.0-7601.0*z*z/354900.0
            *(1.0-91.0*z*z/4146.0*(1.0-3617.0*z*z/161840.0)
           ))))))))
        ) ;
        // clang-format on
        return temp;
    }

    else if (x < x_2)
        return -Li2(-x) + Li2(x * x) / 2.0;

    else {
        return zeta2 - Li2(1.0 - x) - log(fabs(x)) * log(fabs(1.0 - x));
    };
}

//==========================================================================================//
//  Polylogarithm of weight 3. From Marco Bonvini.
//------------------------------------------------------------------------------------------//

double Li3(double x) {

    double x_0 = -1.0;
    double x_1 = -0.85;
    double x_2 = 0.25;
    double x_3 = 0.63;
    double x_4 = 1.0;

    if (x == 1.)
        return zeta3;

    if (x == -1.)
        return -0.75 * zeta3;

    if (x <= x_0) {
        double lnx = log(-x);
        return Li3(1.0 / x) - zeta2 * lnx - lnx * lnx * lnx / 6.0;
    }

    else if (x < x_1) {
        return Li3(x * x) / 4.0 - Li3(-x);
    }

    else if (x < x_2) {
        double z = -log(1.0 - x);
        // clang-format off
        double temp = (
            z*(1.0-3.0*z/8.0*(1.0-17.0*z/81.0*(1.0-15*z/136.0
            *(1.0-28.0*z/1875.0*(1.0+5.0*z/8.0*(1.0-304.0*z/7203.0
            *(1.0+945.0*z/2432.0*(1.0-44.0*z/675.0*(1.0+7.0*z/24.0
            *(1.0-26104.0*z/307461.0*(1.0+1925.0*z/8023.0
            *(1.0-53598548.0*z/524808375.0
            *(1.0+22232925.0*z/107197096.0
           )))))))))))))
        ) ;
        // clang-format on
        return temp;
    }

    else if (x < x_3) {
        return Li3(x * x) / 4.0 - Li3(-x);
    }

    else if (x < x_4) {
        double ln1x = log(1.0 - x);
        return (
            -Li3(1.0 - x) - Li3(-x / (1.0 - x)) + zeta3 + zeta2 * ln1x
            - log(x) * ln1x * ln1x / 2.0 + ln1x * ln1x * ln1x / 6.0
        );
    }

    else {
        double lnx = log(x);
        return Li3(1. / x) + 2.0 * zeta2 * lnx - lnx * lnx * lnx / 6.0;
    }
}

//==========================================================================================//
//  Nielsen function. From Marco Bonvini.
//------------------------------------------------------------------------------------------//

double S12(double x) {

    if (x > 1) {
        return (
            -Li3(1. - x) + zeta3 + log(x - 1.) * Li2(1. - x)
            + 0.5 * log(x) * (log(x - 1.) * log(x - 1.) - M_PI * M_PI)
        );
    }

    else if (x == 1) {
        return zeta3;
    }

    else if ((0 < x) && (x < 1)) {
        double logxm = log(1. - x);
        return (
            -Li3(1. - x) + zeta3 + logxm * Li2(1. - x)
            + 0.5 * log(x) * logxm * logxm
        );
    }

    else if (x == 0) {
        return 0.;
    }

    else if (x < 0) {
        double c = 1. / (1. - x);
        double logc = log(c);
        return (
            -Li3(c) + zeta3 + logc * Li2(c) + 0.5 * logc * logc * log(1. - c)
            - 1. / 6. * logc * logc * logc
        );
    }

    else {
        return 0.;
    }
}

//==========================================================================================//
//  Harmonic polylogarithms of weight 1.
//
//  Eq. (1) of Ref. [arXiv:hep-ph/9905237v1]
//------------------------------------------------------------------------------------------//

double H_0(double x) { return log(x); }

//------------------------------------------------------------------------------------------//

double H_1(double x) { return -log(1. - x); }

//------------------------------------------------------------------------------------------//

double H_m1(double x) { return log(1. + x); }

//==========================================================================================//
//  Harmonic polylogarithms of weight 2.
//
//  Eq. (9,11) of Ref. [arXiv:hep-ph/9905237v1]
//------------------------------------------------------------------------------------------//

double H_m1m1(double x) {

    double Lp = log(1. + x);
    return 0.5 * Lp * Lp;
}

//------------------------------------------------------------------------------------------//

double H_m10(double x) { return log(x) * log(1. + x) + Li2(-x); }

//------------------------------------------------------------------------------------------//

double H_m11(double x) {

    double xp = 1. + x;
    return Li2(xp / 2) - ln2 * log(xp) - Li2_1_2;
}

//------------------------------------------------------------------------------------------//

double H_0m1(double x) { return -Li2(-x); }

//------------------------------------------------------------------------------------------//

double H_00(double x) {

    double L = log(x);
    return 0.5 * L * L;
}

//------------------------------------------------------------------------------------------//

double H_01(double x) { return Li2(x); }

//------------------------------------------------------------------------------------------//

double H_1m1(double x) {

    double xm = 1. - x;
    return Li2(xm / 2.) - ln2 * log(xm) - Li2_1_2;
}

double H_10(double x) { return -log(x) * log(1. - x) - Li2(x); }

//------------------------------------------------------------------------------------------//

double H_11(double x) {

    double Lm = log(1. - x);
    return 0.5 * Lm * Lm;
}

//==========================================================================================//
//  Harmonic polylogarithms of weight 3.
//
//  Eq. (78-89) of Ref. [arXiv:hep-ph/0507152v2]
//------------------------------------------------------------------------------------------//

double H_m1m1m1(double x) {

    double Lp = log(1. + x);
    return 1. / 6. * Lp * Lp * Lp;
}

//------------------------------------------------------------------------------------------//

double H_m1m10(double x) {

    return H_0(x) * H_m1m1(x) - H_0m1m1(x) - H_m10m1(x);
}

//------------------------------------------------------------------------------------------//

double H_m1m11(double x) { return 0.5 * H_m1(x) * H_m11(x) - 0.5 * H_m11m1(x); }

//------------------------------------------------------------------------------------------//

double H_m10m1(double x) { return H_0m1(x) * H_m1(x) - 2. * H_0m1m1(x); }

//------------------------------------------------------------------------------------------//

double H_m100(double x) {

    double L = log(x);
    return 0.5 * L * L * log(1. + x) + L * Li2(-x) - Li3(-x);
}

//------------------------------------------------------------------------------------------//

double H_m101(double x) { return H_m1(x) * H_01(x) - H_0m11(x) - H_01m1(x); }

//------------------------------------------------------------------------------------------//

double H_m11m1(double x) { return H_m1(x) * H_1m1(x) - 2. * H_1m1m1(x); }

//------------------------------------------------------------------------------------------//

double H_m110(double x) { return H_0(x) * H_m11(x) - H_0m11(x) - H_m101(x); }

//------------------------------------------------------------------------------------------//

double H_m111(double x) {

    double xm = 1. - x;
    double Lm = log(xm);
    return 0.5 * log((1. + x) / 2.) * Lm * Lm + Li3(0.5) + Lm * Li2(xm / 2.)
           - Li3(xm / 2);
}

//------------------------------------------------------------------------------------------//

double H_0m1m1(double x) { return S12(-x); }

//------------------------------------------------------------------------------------------//

double H_0m10(double x) { return -log(x) * Li2(-x) + 2. * Li3(-x); }

//------------------------------------------------------------------------------------------//

double H_0m11(double x) {

    double xm = 1. - x;
    double Lm = log(xm);
    double Lm2 = Lm * Lm;
    double Lm3 = Lm2 * Lm;

    return -S12(x) + Li3(-2. * x / xm) - Li3(xm / 2) - Li3(-x) + Li3_1_2
           + Li3(x) + Lm * Li2(-x) + Lm * Li2_1_2 - Lm * Li2(x)
           + 0.5 * ln2 * Lm2 - 1. / 6. * Lm3;
}

//------------------------------------------------------------------------------------------//

double H_00m1(double x) { return -Li3(-x); }

//------------------------------------------------------------------------------------------//

double H_000(double x) {

    double L = log(x);
    return 1. / 6. * L * L * L;
}

//------------------------------------------------------------------------------------------//

double H_001(double x) { return Li3(x); }

//------------------------------------------------------------------------------------------//

double H_01m1(double x) {

    double xp = 1. + x;
    double Lp = log(xp);
    return Li3(2. * x / xp) - Li3(x / xp) - Li3(xp / 2.) - Li3(x)
           + Lp * Li2(0.5) + Lp * Li2(x) + 0.5 * ln2 * Lp * Lp + Li3_1_2;
}

//------------------------------------------------------------------------------------------//

double H_010(double x) { return log(x) * Li2(x) - 2. * Li3(x); }

//------------------------------------------------------------------------------------------//

double H_011(double x) { return S12(x); }

//------------------------------------------------------------------------------------------//

double H_1m1m1(double x) {

    double xp = 1. + x;
    double Lp = log(xp);
    return -0.5 * log((1. - x) / 2.) * Lp * Lp - Li3_1_2 - Lp * Li2(xp / 2.)
           + Li3(xp / 2.);
}

//------------------------------------------------------------------------------------------//

double H_1m10(double x) { return H_0(x) * H_1m1(x) - H_01m1(x) - H_10m1(x); }

//------------------------------------------------------------------------------------------//

double H_1m11(double x) { return H_1(x) * H_m11(x) - 2. * H_m111(x); }

//------------------------------------------------------------------------------------------//

double H_10m1(double x) { return H_1(x) * H_0m1(x) - H_01m1(x) - H_0m11(x); }

//------------------------------------------------------------------------------------------//

double H_100(double x) {

    double L = log(x);
    return -1. / 2 * L * L * log(1. - x) - L * Li2(x) + Li3(x);
}

//------------------------------------------------------------------------------------------//

double H_101(double x) { return -2. * S12(x) - log(1. - x) * Li2(x); }

//------------------------------------------------------------------------------------------//

double H_11m1(double x) { return 0.5 * H_1(x) * H_1m1(x) - 0.5 * H_1m11(x); }

//------------------------------------------------------------------------------------------//

double H_110(double x) {

    double xm = 1. - x;
    return zeta2 * log(xm) - Li3(xm) + zeta3;
}

//------------------------------------------------------------------------------------------//

double H_111(double x) {

    double Lm = log(1. - x);
    return -1. / 6. * Lm * Lm * Lm;
}
