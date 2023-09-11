#include "adani/ExactCoefficientFunctions.h"
#include "adani/Constants.h"
#include "adani/Convolutions.h"
#include "adani/MatchingConditions.h"
#include "adani/SpecialFunctions.h"
#include "adani/SplittingFunctions.h"
#include <cmath>
#include <iostream>

//==========================================================================================//
//  Exact massive gluon coefficient functions for F2 at O(alpha_s)
//
//  Eq. (50) from Ref. [arXiv:1001.2312]
//------------------------------------------------------------------------------------------//

double C2_g1(double x, double m2Q2) { // m2Q2=m^2/Q^2

    double beta = sqrt(1. - 4. * m2Q2 * x / (1 - x));
    double x2 = x * x;
    double m4Q4 = m2Q2 * m2Q2;
    double L = log((1. + beta) / (1. - beta));

    return 4. * TR
           * (L
                  * (-8. * x2 * m4Q4 - 4. * x * m2Q2 * (3. * x - 1) + 2. * x2
                     - 2. * x + 1.)
              + beta * (4. * m2Q2 * x * (x - 1.) - (8. * x2 - 8. * x + 1.)));
}

//==========================================================================================//
//  Exact massive gluon coefficient functions for FL at O(alpha_s)
//
//  Eq. (2.9) from Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double CL_g1(double x, double m2Q2) {

    double beta = sqrt(1. - 4. * m2Q2 * x / (1. - x));
    double x2 = x * x;
    double L = log((1. + beta) / (1. - beta));

    return 16. * TR * (x * (1. - x) * beta - 2. * x2 * m2Q2 * L);
}

//==========================================================================================//
//  OBSERVATION: in the O(alpha_s^2) exact coefficeint functions the
//  mu-independent part is an interpolation in a certain grid. When this
//  function is called for a (x,Q) value outside this grid, the value 0 is
//  returned. The mu-dependent part is defined everywhere. However, also this
//  one is put to zero for values outside the grid otherwise we would have that
//  one term contributes while the other doesn't. Unfortunately in this way the
//  check if(eta > 1e6 || eta < 1e-6 || xi<1e-3 || xi>1e5) return 0. must be
//  performed twice (but this is not a big issue since in the integrals we have
//  only the interpolated result).
//------------------------------------------------------------------------------------------//

//==========================================================================================//
//  Exact massive gluon coefficient functions for F2 at O(alpha_s^2)
//
//  Exact (but numerical) result from [arXiv:hep-ph/9411431].
//  Taken from the Fortran code 'src/hqcoef.f'
//------------------------------------------------------------------------------------------//

double C2_g2(double x, double m2Q2, double m2mu2) {

    double xi = 1. / m2Q2;
    double eta = 0.25 * xi * (1 - x) / x - 1.;

    if (eta > 1e6 || eta < 1e-6 || xi < 1e-3 || xi > 1e5)
        return 0.;

    return C2_g20(x, m2Q2) + C2_g21(x, m2Q2) * log(1. / m2mu2);
}

//==========================================================================================//
//  Exact massive quark coefficient functions for F2 at O(alpha_s)
//
//  Exact (but numerical) result from [arXiv:hep-ph/9411431].
//  Taken from the Fortran code 'src/hqcoef.f'
//------------------------------------------------------------------------------------------//

double C2_ps2(double x, double m2Q2, double m2mu2) {

    double xi = 1. / m2Q2;
    double eta = 0.25 * xi * (1 - x) / x - 1.;

    if (eta > 1e6 || eta < 1e-6 || xi < 1e-3 || xi > 1e5)
        return 0.;

    return C2_ps20(x, m2Q2) + C2_ps21(x, m2Q2) * log(1. / m2mu2);
}

//==========================================================================================//
//  Exact massive gluon coefficient functions for FL at O(alpha_s)
//
//  Exact (but numerical) result from [arXiv:hep-ph/9411431].
//  Taken from the Fortran code 'src/hqcoef.f'
//------------------------------------------------------------------------------------------//

double CL_g2(double x, double m2Q2, double m2mu2) {

    double xi = 1. / m2Q2;
    double eta = 0.25 * xi * (1 - x) / x - 1.;

    if (eta > 1e6 || eta < 1e-6 || xi < 1e-3 || xi > 1e5)
        return 0.;

    return CL_g20(x, m2Q2) + CL_g21(x, m2Q2) * log(1. / m2mu2);
}

//==========================================================================================//
//  Exact massive quarkk coefficient functions for FL at O(alpha_s)
//
//  Exact (but numerical) result from [arXiv:hep-ph/9411431].
//  Taken from the Fortran code 'src/hqcoef.f'
//------------------------------------------------------------------------------------------//

double CL_ps2(double x, double m2Q2, double m2mu2) {

    double xi = 1. / m2Q2;
    double eta = 0.25 * xi * (1 - x) / x - 1.;

    if (eta > 1e6 || eta < 1e-6 || xi < 1e-3 || xi > 1e5)
        return 0.;

    return CL_ps20(x, m2Q2) + CL_ps21(x, m2Q2) * log(1. / m2mu2);
}

//==========================================================================================//
//  Exact massive quark coefficient functions for F2 at O(alpha_s^2):
//  mu independent term.
//
//  Eq. (4.4) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2_ps20(double x, double m2Q2) {

    double xi = 1. / m2Q2;
    double eta = 0.25 * xi * (1 - x) / x - 1.;

    if (eta > 1e6 || eta < 1e-6 || xi < 1e-3 || xi > 1e5)
        return 0.;

    return 16 * M_PI * xi * c2nloq_(&eta, &xi) / x;
}

//==========================================================================================//
//  Exact massive quark coefficient functions for F2 at O(alpha_s^2):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.1) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2_ps21(double x, double m2Q2) {

    return -C2_g1_x_Pgq0(x, m2Q2);
    // The minus sign comes from the fact that in [arXiv:1205.5727]
    // the expansion is performed in terms of log(m^2/mu^2)
    // (even if it says the opposite) but we are
    // expanding in terms of log(mu^2/m^2)
}

//==========================================================================================//
//  Exact massive quark coefficient functions for FL at O(alpha_s^2):
//  mu independent term.
//
//  Eq. (4.4) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double CL_ps20(double x, double m2Q2) {

    double xi = 1. / m2Q2;
    double eta = 0.25 * xi * (1 - x) / x - 1.;

    if (eta > 1e6 || eta < 1e-6 || xi < 1e-3 || xi > 1e5)
        return 0.;

    return 16 * M_PI * xi * clnloq_(&eta, &xi) / x;
}

//==========================================================================================//
//  Exact massive quark coefficient functions for FL at O(alpha_s^2):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.1) of Ref. [arXiv:1205.5727] for FL
//------------------------------------------------------------------------------------------//

double CL_ps21(double x, double m2Q2) { return -CL_g1_x_Pgq0(x, m2Q2); }

//==========================================================================================//
//  Exact massive gluon coefficient functions for F2 at O(alpha_s^2):
//  mu independent term.
//
//  Eq. (4.4) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2_g20(double x, double m2Q2) {

    double xi = 1. / m2Q2;
    double eta = 0.25 * xi * (1 - x) / x - 1.;

    if (eta > 1e6 || eta < 1e-6 || xi < 1e-3 || xi > 1e5)
        return 0.;

    return 16 * M_PI * xi * c2nlog_(&eta, &xi) / x;
}

//==========================================================================================//
//  Exact massive gluon coefficient functions for F2 at O(alpha_s^2):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.4) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2_g21(double x, double m2Q2) {

    int nf = 1;
    // Put nf to 1 since the nf contribution cancels for any value of nf

    return -(C2_g1_x_Pgg0(x, m2Q2, nf) - C2_g1(x, m2Q2) * beta(0, nf));
}

//==========================================================================================//
//  Exact massive gluon coefficient functions for FL at O(alpha_s^2):
//  mu independent term.
//
//  Eq. (4.4) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double CL_g20(double x, double m2Q2) {

    double xi = 1. / m2Q2;
    double eta = 0.25 * xi * (1 - x) / x - 1.;

    if (eta > 1e6 || eta < 1e-6 || xi < 1e-3 || xi > 1e5)
        return 0.;

    return 16 * M_PI * xi * clnlog_(&eta, &xi) / x;
}

//==========================================================================================//
//  Exact massive gluon coefficient functions for FL at O(alpha_s^2):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.4) of Ref. [arXiv:1205.5727] for FL
//------------------------------------------------------------------------------------------//

double CL_g21(double x, double m2Q2) {

    int nf = 1;
    // Put nf to 1 since the nf contribution cancels for any value of nf

    return -(CL_g1_x_Pgg0(x, m2Q2, nf) - CL_g1(x, m2Q2) * beta(0, nf));
}

//==========================================================================================//
//  Exact massive quark coefficient functions for F2 at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.2) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2_ps31(double x, double m2Q2, int nf) {

    return -(
        C2_g1_x_Pgq1(x, m2Q2, nf) + C2_g20_x_Pgq0(x, m2Q2)
        + C2_ps20_x_Pqq0(x, m2Q2, nf) - 2. * beta(0, nf) * C2_ps20(x, m2Q2)
    );
}

//==========================================================================================//
//  Exact massive quark coefficient functions for FL at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.2) of Ref. [arXiv:1205.5727] for FL
//------------------------------------------------------------------------------------------//

double CL_ps31(double x, double m2Q2, int nf) {

    return -(
        CL_g1_x_Pgq1(x, m2Q2, nf) + CL_g20_x_Pgq0(x, m2Q2)
        + CL_ps20_x_Pqq0(x, m2Q2, nf) - 2. * beta(0, nf) * CL_ps20(x, m2Q2)
    );
}

//==========================================================================================//
//  Exact massive quark coefficient functions for F2 at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)^2
//
//  Eq. (4.3) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2_ps32(double x, double m2Q2, int nf) {

    return 0.5
               * (C2_g1_x_Pgg0_x_Pgq0(x, m2Q2, nf)
                  + C2_g1_x_Pqq0_x_Pgq0(x, m2Q2, nf))
           - 3. / 2 * beta(0, nf) * C2_g1_x_Pgq0(x, m2Q2);
}

//==========================================================================================//
//  Exact massive quark coefficient functions for FL at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)^2
//
//  Eq. (4.3) of Ref. [arXiv:1205.5727] for FL
//------------------------------------------------------------------------------------------//

double CL_ps32(double x, double m2Q2, int nf) {

    return 0.5
               * (CL_g1_x_Pgg0_x_Pgq0(x, m2Q2, nf)
                  + CL_g1_x_Pqq0_x_Pgq0(x, m2Q2, nf))
           - 3. / 2 * beta(0, nf) * CL_g1_x_Pgq0(x, m2Q2);
}

//==========================================================================================//
//  Exact massive gluon coefficient functions for F2 at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.5) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2_g31(double x, double m2Q2, int nf) {

    return -(
        C2_g1_x_Pgg1(x, m2Q2, nf) - beta(1, nf) * C2_g1(x, m2Q2)
        + C2_ps20_x_Pqg0(x, m2Q2, nf) + C2_g20_x_Pgg0(x, m2Q2, nf)
        - 2. * beta(0, nf) * C2_g20(x, m2Q2)
    );
}

//==========================================================================================//
//  Exact massive gluon coefficient functions for FL at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.5) of Ref. [arXiv:1205.5727] for FL
//------------------------------------------------------------------------------------------//

double CL_g31(double x, double m2Q2, int nf) {

    return -(
        CL_g1_x_Pgg1(x, m2Q2, nf) - beta(1, nf) * CL_g1(x, m2Q2)
        + CL_ps20_x_Pqg0(x, m2Q2, nf) + CL_g20_x_Pgg0(x, m2Q2, nf)
        - 2 * beta(0, nf) * CL_g20(x, m2Q2)
    );
}

//==========================================================================================//
//  Exact massive gluon coefficient functions for F2 at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)^2
//
//  Eq. (4.6) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2_g32(double x, double m2Q2, int nf, int method_flag) {

    double C2_g1xPgg0xPgg0;

    double beta0 = beta(0, nf);

    if (method_flag == 0)
        C2_g1xPgg0xPgg0 = C2_g1_x_Pgg0_x_Pgg0(x, m2Q2, nf);
    else if (method_flag == 1)
        C2_g1xPgg0xPgg0 = C2_g1_x_Pgg0_x_Pgg0_MC(x, m2Q2, nf);
    else {
        std::cout << "C2_g32: Choose either method_flag = 0 or method_flag = 1"
                  << std::endl;
        exit(-1);
    }

    return 0.5 * C2_g1xPgg0xPgg0 + 0.5 * C2_g1_x_Pqg0_x_Pgq0(x, m2Q2, nf)
           - 3. / 2 * beta0 * C2_g1_x_Pgg0(x, m2Q2, nf)
           + beta0 * beta0 * C2_g1(x, m2Q2);
}

//==========================================================================================//
//  Exact massive gluon coefficient functions for FL at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)^2
//
//  Eq. (4.6) of Ref. [arXiv:1205.5727] for FL
//------------------------------------------------------------------------------------------//

double CL_g32(double x, double m2Q2, int nf, int method_flag) {

    double CL_g1xPgg0xPgg0;

    double beta0 = beta(0, nf);

    if (method_flag == 0)
        CL_g1xPgg0xPgg0 = CL_g1_x_Pgg0_x_Pgg0(x, m2Q2, nf);
    else if (method_flag == 1)
        CL_g1xPgg0xPgg0 = CL_g1_x_Pgg0_x_Pgg0_MC(x, m2Q2, nf);
    else {
        std::cout << "C2_g32: Choose either method_flag = 0 or method_flag = 1"
                  << std::endl;
        exit(-1);
    }

    return 0.5 * CL_g1xPgg0xPgg0 + 0.5 * CL_g1_x_Pqg0_x_Pgq0(x, m2Q2, nf)
           - 3. / 2 * beta0 * CL_g1_x_Pgg0(x, m2Q2, nf)
           + beta0 * beta0 * CL_g1(x, m2Q2);
}
