#include "adani/ExactCoefficientFunctions.h"
#include "adani/MatchingConditions.h"
#include "adani/Convolutions.h"
#include "adani/ColorFactors.h"
#include "adani/Convolutions.h"
#include "adani/SplittingFunctions.h"
#include "adani/SpecialFunctions.h"
#include <cmath>
#include <iostream>

//==========================================================================================//
//  Exact massive gluon coefficient functions for F2 at O(alpha_s)
//
//  Eq. (50) from Ref. [arXiv:1001.2312]
//------------------------------------------------------------------------------------------//

double C2m_g1(double x, double mQ) { //mQ=m^2/Q^2

    double x_max = 1. / (1. + 4. * mQ) ;

    if (x>=x_max || x<0) return 0;

    double beta = sqrt(1. - 4. * mQ * x / (1 - x)) ;
    double x2 = x * x ;
    double mQ2 = mQ * mQ ;
    double L = log((1. + beta) / (1. - beta)) ;

    return 4 * TR * (
        L * (-8 * x2 * mQ2 - 4 * x * mQ * (3 * x - 1) + 2 * x2 - 2 * x + 1)
        + beta * (4 * mQ * x * (x - 1) - (8 * x2 - 8 * x + 1))
    ) / 4. / M_PI ;

}

//==========================================================================================//
//  Exact massive gluon coefficient functions for FL at O(alpha_s)
//
//  Eq. (2.9) from Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double CLm_g1(double x, double mQ) {

    double x_max = 1. / (1. + 4. * mQ) ;

    if (x>=x_max || x<0) return 0;

    double beta = sqrt(1. - 4. * mQ * x / (1. - x)) ;
    double x2 = x * x ;
    double L = log((1. + beta) / (1. - beta)) ;

    return 16. * TR * (
        x * (1. - x) * beta - 2. * x2 * mQ * L
    ) / 4. / M_PI ;

}


// double D2m_g1_4(double x, double mQ) {

//     return C2m_g1(x,mQ);

// }


// double DLm_g1_4(double x, double mQ) {

//     return CLm_g1(x,mQ);

// }


// double D2m_g1_5(double x, double mQ) {

//     return D2m_g1_4(x,mQ) - 2 * K_Qg1(x,mQ);

// }


// double DLm_g1_5(double x, double mQ) {

//     return DLm_g1_4(x,mQ);

// }


//==========================================================================================//
//  Exact massive gluon coefficient functions for F2 at O(alpha_s^2)
//
//  Exact (but numerical) result from [arXiv:hep-ph/9411431].
//  Taken from the Fortran code 'src/hqcoef.f'
//------------------------------------------------------------------------------------------//

double C2m_g2(double x, double mQ, double mMu) {

    double xi = 1. / mQ;
    double eta = 0.25 * xi * (1 - x) / x - 1. ;

    double x_max = 1. / (1 + 4 * mQ) ;

    if (x>=x_max || x<0) return 0;

    //if(eta > 1e6 || eta < 1e-6 || xi<1e-3 || xi>1e5) return __builtin_nan("");
    if(eta > 1e6 || eta < 1e-6 || xi<1e-3 || xi>1e5) return 0.;

    return (
        xi * c2nlog_(&eta, &xi) / x / M_PI
        + C2m_g21(x, mQ) * log(1. / mMu)
    ) ;

}

//==========================================================================================//
//  Exact massive quark coefficient functions for F2 at O(alpha_s)
//
//  Exact (but numerical) result from [arXiv:hep-ph/9411431].
//  Taken from the Fortran code 'src/hqcoef.f'
//------------------------------------------------------------------------------------------//

double C2m_ps2(double x, double mQ, double mMu) {

    double xi = 1. / mQ;
    double eta = 0.25 * xi * (1 - x) / x - 1. ;

    double x_max = 1. / (1 + 4 * mQ) ;

    if (x>=x_max || x<0) return 0;

    //if(eta > 1e6 || eta < 1e-6 || xi<1e-3 || xi>1e5) return __builtin_nan("");
    if(eta > 1e6 || eta < 1e-6 || xi<1e-3 || xi>1e5) return 0.;

    return (
        xi * c2nloq_(&eta, &xi) / x / M_PI
        + C2m_ps21(x, mQ) * log(1. / mMu)
    ) ;

}

//==========================================================================================//
//  Exact massive gluon coefficient functions for FL at O(alpha_s)
//
//  Exact (but numerical) result from [arXiv:hep-ph/9411431].
//  Taken from the Fortran code 'src/hqcoef.f'
//------------------------------------------------------------------------------------------//

double CLm_g2(double x, double mQ, double mMu) {

    double xi = 1. / mQ;
    double eta = 0.25 * xi * (1 - x) / x - 1. ;

    double x_max = 1. / (1 + 4 * mQ) ;

    if (x>=x_max || x<0) return 0;

    //if(eta > 1e6 || eta < 1e-6 || xi<1e-3 || xi>1e5) return __builtin_nan("");
    if(eta > 1e6 || eta < 1e-6 || xi<1e-3 || xi>1e5) return 0.;

    return (
        xi * clnlog_(&eta, &xi) / x / M_PI
        + CLm_g21(x, mQ) * log(1. / mMu)
    ) ;

}

//==========================================================================================//
//  Exact massive quarkk coefficient functions for FL at O(alpha_s)
//
//  Exact (but numerical) result from [arXiv:hep-ph/9411431].
//  Taken from the Fortran code 'src/hqcoef.f'
//------------------------------------------------------------------------------------------//

double CLm_ps2(double x, double mQ, double mMu) {

    double xi = 1. / mQ;
    double eta = 0.25 * xi * (1 - x) / x - 1. ;

    double x_max = 1. / (1 + 4 * mQ) ;

    if (x>=x_max || x<0) return 0;

    //if(eta > 1e6 || eta < 1e-6 || xi<1e-3 || xi>1e5) return __builtin_nan("");
    if(eta > 1e6 || eta < 1e-6 || xi<1e-3 || xi>1e5) return 0.;

    return (
        xi * clnloq_(&eta, &xi) / x / M_PI
        + CLm_ps21(x, mQ) * log(1. / mMu)
    ) ;

}


// double D2m_g2_4(double x, double mQ, double mMu) {

//     double Lmu=log(1./mMu);

//     return C2m_g2(x,mQ,mMu) - Lmu / 6. / M_PI * C2m_g1(x,mQ);

// }

// double DLm_g2_4(double x, double mQ, double mMu) {

//     double Lmu=log(1./mMu);

//     return CLm_g2(x,mQ,mMu) - Lmu / 6. / M_PI * CLm_g1(x,mQ);

// }

// double D2m_g2_5(double x, double mQ, double mMu) {

//     return (
//         D2m_g2_4(x, mQ, mMu)
//         -D2m_g1_4(x,mQ)*K_gg1_local(mMu)
//         -2*(K_Qg2(x,mMu)-K_Qg1(x,mMu)*K_gg1_local(mMu))
//         -2*C2_b1_x_K_Qg1(x,mQ)
//    ) ;

// }


// double DLm_g2_5(double x, double mQ, double mMu) {

//     return (
//         DLm_g2_4(x, mQ, mMu)
//         -DLm_g1_4(x,mQ)*K_gg1_local(mMu)
//         -2*CL_b1_x_K_Qg1(x,mQ)
//    ) ;

// }

//==========================================================================================//
//  Exact massive quark coefficient functions for F2 at O(alpha_s^2):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.1) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2m_ps21(double x, double mQ) {

    double x_max = 1. / (1. + 4. * mQ) ;

    if (x>=x_max || x<0) return 0;

    return - C2m_g1_x_Pgq0(x, mQ);
    // The minus sign comes from the fact that in [arXiv:1205.5727]
    // the expansion is performed in terms of log(m^2/mu^2)
    // (even if it says the opposite) but we are
    // expanding in terms of log(mu^2/m^2)

}

//==========================================================================================//
//  Exact massive quark coefficient functions for FL at O(alpha_s^2):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.1) of Ref. [arXiv:1205.5727] for FL
//------------------------------------------------------------------------------------------//

double CLm_ps21(double x, double mQ) {

    double x_max = 1. / (1. + 4. * mQ) ;

    if (x>=x_max || x<0) return 0;

    return - CLm_g1_x_Pgq0(x, mQ);

}

//==========================================================================================//
//  Exact massive gluon coefficient functions for F2 at O(alpha_s^2):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.4) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2m_g21(double x, double mQ) {

    double x_max = 1. / (1. + 4. * mQ) ;

    if (x>=x_max || x<0) return 0;

    int nf = 1 ;
    //Put nf to 1 since the nf contribution cancels for any value of nf

    return - (C2m_g1_x_Pgg0(x, mQ, nf) - C2m_g1(x, mQ) * beta(0,nf));

}

//==========================================================================================//
//  Exact massive gluon coefficient functions for FL at O(alpha_s^2):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.4) of Ref. [arXiv:1205.5727] for FL
//------------------------------------------------------------------------------------------//

double CLm_g21(double x, double mQ) {

    double x_max = 1. / (1. + 4. * mQ) ;

    if (x>=x_max || x<0) return 0;

    int nf = 1 ;
    //Put nf to 1 since the nf contribution cancels for any value of nf

    return - (CLm_g1_x_Pgg0(x, mQ, nf) - CLm_g1(x, mQ) * beta(0,nf));

}

//==========================================================================================//
//  Exact massive quark coefficient functions for F2 at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.2) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2m_ps31(double x, double mQ, int nf) {

    double x_max = 1. / (1. + 4. * mQ) ;

    if (x>=x_max || x<0) return 0;

    return - (
        C2m_g1_x_Pgq1(x, mQ, nf)
        + C2m_g20_x_Pgq0(x, mQ)
        + C2m_ps20_x_Pqq0(x, mQ, nf)
        - 2. * beta(0, nf) * C2m_ps2(x, mQ, 1)
   ) ;

}

//==========================================================================================//
//  Exact massive quark coefficient functions for FL at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.2) of Ref. [arXiv:1205.5727] for FL
//------------------------------------------------------------------------------------------//

double CLm_ps31(double x, double mQ, int nf) {

    double x_max = 1. / (1. + 4. * mQ) ;

    if (x>=x_max || x<0) return 0;

    return - (
        CLm_g1_x_Pgq1(x, mQ, nf)
        + CLm_g20_x_Pgq0(x, mQ)
        + CLm_ps20_x_Pqq0(x, mQ, nf)
        - 2. * beta(0, nf) * CLm_ps2(x, mQ, 1)
   ) ;

}

//==========================================================================================//
//  Exact massive quark coefficient functions for F2 at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)^2
//
//  Eq. (4.3) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2m_ps32(double x, double mQ, int nf) {

    double x_max = 1. / (1. + 4. * mQ) ;

    if (x>=x_max || x<0) return 0;

    return (
        0.5 * (C2m_g1_x_Pgg0_x_Pgq0(x, mQ, nf) + C2m_g1_x_Pqq0_x_Pgq0(x, mQ, nf))
        - 3. / 2 * beta(0, nf) * C2m_g1_x_Pgq0(x, mQ)
   ) ;

}

//==========================================================================================//
//  Exact massive quark coefficient functions for FL at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)^2
//
//  Eq. (4.3) of Ref. [arXiv:1205.5727] for FL
//------------------------------------------------------------------------------------------//

double CLm_ps32(double x, double mQ, int nf) {

    double x_max = 1. / (1. + 4. * mQ) ;

    if (x>=x_max || x<0) return 0;

    return (
        0.5 * (CLm_g1_x_Pgg0_x_Pgq0(x, mQ, nf) + CLm_g1_x_Pqq0_x_Pgq0(x, mQ, nf))
        - 3. / 2 * beta(0, nf) * CLm_g1_x_Pgq0(x, mQ)
   ) ;

}

//==========================================================================================//
//  Exact massive gluon coefficient functions for F2 at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.5) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2m_g31(double x, double mQ, int nf) {

    double x_max = 1. / (1. + 4. * mQ) ;

    if (x>=x_max || x<0) return 0. ;

    return -(
        C2m_g1_x_Pgg1(x, mQ, nf)
        - beta(1, nf) * C2m_g1(x, mQ)
        + C2m_ps20_x_Pqg0(x, mQ, nf)
        + C2m_g20_x_Pgg0(x, mQ, nf)
        - 2. * beta(0,nf) * C2m_g2(x, mQ, 1)
   );

}

//==========================================================================================//
//  Exact massive gluon coefficient functions for FL at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.5) of Ref. [arXiv:1205.5727] for FL
//------------------------------------------------------------------------------------------//

double CLm_g31(double x, double mQ, int nf) {

    double x_max = 1. / (1. + 4. * mQ) ;

    if (x>=x_max || x<0) return 0. ;

    return -(
        CLm_g1_x_Pgg1(x, mQ, nf)
        - beta(1, nf) * CLm_g1(x, mQ)
        + CLm_ps20_x_Pqg0(x, mQ, nf)
        + CLm_g20_x_Pgg0(x, mQ, nf)
        - 2 * beta(0,nf) * CLm_g2(x, mQ, 1)
   );

}

//==========================================================================================//
//  Exact massive gluon coefficient functions for F2 at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)^2
//
//  Eq. (4.6) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double C2m_g32(double x, double mQ, int nf, int method_flag, int calls) {

    double x_max = 1. / (1. + 4. * mQ) ;

    if (x>=x_max || x<0) return 0;

    double C2m_g1xPgg0xPgg0 ;

    if(method_flag == 0) C2m_g1xPgg0xPgg0 = C2m_g1_x_Pgg0_x_Pgg0(x, mQ, nf) ;
    else if(method_flag == 1) C2m_g1xPgg0xPgg0 = C2m_g1_x_Pgg0_x_Pgg0_MC(x, mQ, nf, calls) ;
    else {
        std::cout << "Choose either method_flag = 0 or method_flag = 1" << std::endl ;
        exit(-1);
    }

    return (
        0.5 * C2m_g1xPgg0xPgg0
        + 0.5 * C2m_g1_x_Pqg0_x_Pgq0(x, mQ, nf)
        - 3. / 2 * beta(0, nf) * C2m_g1_x_Pgg0(x, mQ, nf)
        + beta(0,nf) * beta(0,nf) * C2m_g1(x,mQ)
   );

}

//==========================================================================================//
//  Exact massive gluon coefficient functions for FL at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)^2
//
//  Eq. (4.6) of Ref. [arXiv:1205.5727] for FL
//------------------------------------------------------------------------------------------//

double CLm_g32(double x, double mQ, int nf, int method_flag, int calls) {

    double x_max = 1. / (1. + 4. * mQ) ;

    if (x>=x_max || x<0) return 0;

    double CLm_g1xPgg0xPgg0 ;

    if(method_flag == 0) CLm_g1xPgg0xPgg0 = CLm_g1_x_Pgg0_x_Pgg0(x, mQ, nf) ;
    else if(method_flag == 1) CLm_g1xPgg0xPgg0 = CLm_g1_x_Pgg0_x_Pgg0_MC(x, mQ, nf, calls) ;
    else {
        std::cout << "Choose either method_flag = 0 or method_flag = 1" << std::endl ;
        exit(-1);
    }

    return (
        0.5 * CLm_g1xPgg0xPgg0
        + 0.5 * CLm_g1_x_Pqg0_x_Pgq0(x, mQ, nf)
        - 3. / 2 * beta(0, nf) * CLm_g1_x_Pgg0(x, mQ, nf)
        + beta(0,nf) * beta(0,nf) * CLm_g1(x,mQ)
   );

}