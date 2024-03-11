#include "adani/ExactCoefficientFunctions.h"
#include "adani/Constants.h"
#include "adani/MatchingConditions.h"
#include "adani/SpecialFunctions.h"
#include "adani/SplittingFunctions.h"

#include <cmath>
#include <iostream>

using std::cout ;
using std::endl ;

ExactCoefficientFunction::ExactCoefficientFunction(const int& order, const char& kind, const char& channel, const int& method_flag, const double& abserr, const double& relerr, const int& MCcalls, const int& dim) : CoefficientFunction(order, kind, channel) {

    // if (GetOrder() == 3) {
    //     cout <<  "Error: exact coefficient functions are not known at O(as^3)!" << endl;
    //     exit(-1);
    // }

    SetMethodFlag(method_flag);
    SetAbserr(abserr);
    SetRelerr(relerr);
    SetMCcalls(MCcalls);
    SetDim(dim);

    if (GetOrder() == 1) {
        convolution_ = NULL;
    } else if (GetOrder() == 2) {
        if (GetChannel() == 'q') {
            convolutions_.push_back( Convolution(gluon_leadingorder_, &SplittingFunction(0, 'g', 'q')) );
        } else {
            convolutions_.push_back( Convolution(gluon_leadingorder_, &SplittingFunction(0, 'g', 'g')) );
        }
    } else if (GetOrder() == 3) {
        if (GetChannel() == 'q') {
            convolutions_.push_back( Convolution(gluon_leadingorder_, &SplittingFunction(1, 'g', 'q')) );
            convolutions_.push_back( Convolution(ExactCoefficientFunction(2, GetKind(), 'g'), &SplittingFunction(0, 'g', 'q')) );
            convolutions_.push_back( Convolution(ExactCoefficientFunction(2, GetKind(), 'q'), &SplittingFunction(0, 'q', 'q')) );
            convolutions_.push_back( Convolution(gluon_leadingorder_, &ConvolutedSplittingFunctions(1, 'g', 'q', 'q')) );
            convolutions_.push_back( Convolution(gluon_leadingorder_, SplittingFunction(0, 'g', 'q')) );
        } else {
            convolutions_.push_back( Convolution(gluon_leadingorder_, SplittingFunction(1, 'g', 'g')) );
            convolutions_.push_back( Convolution(ExactCoefficientFunction(2, GetKind(), 'q'), &SplittingFunction(0, 'q', 'g')) );
            convolutions_.push_back( Convolution(ExactCoefficientFunction(2, GetKind(), 'g'), &SplittingFunction(0, 'g', 'g')) );
            // MC integral
            convolutions_.push_back( Convolution(gluon_leadingorder_, &ConvolutedSplittingFunctions(0, 'g', 'q', 'g')) );
            convolutions_.push_back( Convolution(gluon_leadingorder_, &SplittingFunction(0, 'g', 'g')) );
        }
    }

}

double ExactCoefficientFunction::fx(const double x, const double m2Q2, const double m2mu2, const int nf) const {

    // O(as)
    if (GetOrder() == 1 && GetKind() == '2') return C2_g1(x, m2Q2) ;
    if (GetOrder() == 1 && GetKind() == 'L') return CL_g1(x, m2Q2) ;

    // O(as^2)
    if (GetOrder() == 2 && GetKind() == '2' && GetChannel() == 'g') return C2_g2(x, m2Q2, m2mu2) ;
    if (GetOrder() == 2 && GetKind() == '2' && GetChannel() == 'q') return C2_ps2(x, m2Q2, m2mu2) ;
    if (GetOrder() == 2 && GetKind() == 'L' && GetChannel() == 'g') return CL_g2(x, m2Q2, m2mu2) ;
    if (GetOrder() == 2 && GetKind() == 'L' && GetChannel() == 'q') return CL_ps2(x, m2Q2, m2mu2) ;

    else {
        cout << "Error: something has gone wrong!" << endl;
    }
}

double ExactCoefficientFunction::MuIndependentTerms(const double x, const double m2Q2, const int nf) const {
    if (GetOrder() == 1) {
        if (GetKind() == '2') return C2_g1(x, m2Q2) ;
        else if (GetKind() == 'L') return CL_g1(x, m2Q2) ;
    } else if (GetOrder() == 2) {
        if (GetKind() == '2') {
            if (GetChannel() == 'g') return return C2_g20(x, m2Q2) ;
            else if (GetChannel() == 'q') return return C2_ps20(x, m2Q2) ;
        } else if (GetKind() == 'L') {
            if (GetChannel() == 'g') return return CL_g20(x, m2Q2) ;
            else if (GetChannel() == 'q') return return CL_ps20(x, m2Q2) ;
        }
    } else if (GetOrder() == 3){
        cout << "Error: mu independent terms are not known at O(as^3)!" << endl;
        exit(-1);
    }
}

double ExactCoefficientFunction::MuDependentTerms(const double x, const double m2Q2, const double m2mu2, const int nf) const {
    if (GetOrder() == 1) return 0.;
    else if (GetOrder() == 2) {
        if (GetChannel() == 'q') return convolutions_[0].convolute(x, m2Q2, nf) * log(1. / m2mu2) ;
        else return (convolutions_[0].convolute(x, m2Q2, nf) - beta(0, nf) * gluon_leadingorder_ -> fx(x, m2Q2, static_cast<double>(nan("")), nf)) * log(1. / m2mu2) ;
    } else if (GetOrder() == 3) {
        if (GetChannel() == 'q') {

        }
    }
}

// void ExactCoefficientFunction::Set_fx() {

//     // O(as)
//     if (GetOrder() == 1 && GetKind() == '2') fx_ = C2_g1 ;
//     if (GetOrder() == 1 && GetKind() == 'L') fx_ = CL_g1 ;

//     // O(as^2)
//     if (GetOrder() == 2 && GetKind() == '2' && GetChannel() == 'g') fx_ = C2_g2 ;
//     if (GetOrder() == 2 && GetKind() == '2' && GetChannel() == 'q') fx_ = C2_q2 ;
//     if (GetOrder() == 2 && GetKind() == 'L' && GetChannel() == 'g') fx_ = CL_g2 ;
//     if (GetOrder() == 2 && GetKind() == 'L' && GetChannel() == 'q') fx_ = CL_q2 ;

//     else {
//         cout << "Error: something has gone wrong!" << endl;
//     }
// }

void ExactCoefficientFunction::SetMethodFlag(const int& method_flag) {
    // check method_flag
    if (method_flag != 0 && method_flag != 1) {
        cout << "Error: method_flag must be 0 or 1. Got" << method_flag << endl ;
        exit(-1) ;
    }
    method_flag_ = method_flag ;

}

void ExactCoefficientFunction::SetAbserr(const double& abserr) {
    // check abserr
    if (abserr <= 0) {
        cout << "Error: abserr must be positive. Got " << abserr << endl;
        exit(-1);
    }
    abserr_ = abserr;
}

void ExactCoefficientFunction::SetRelerr(const double& relerr) {
    // check relerr
    if (relerr <= 0) {
        cout << "Error: relerr must be positive. Got " << relerr << endl;
        exit(-1);
    }
    relerr_ = relerr;
}

void ExactCoefficientFunction::SetMCcalls(const int& MCcalls) {
    // check MCcalls
    if (MCcalls <= 0) {
        cout << "Error: MCcalls must be positive. Got " << MCcalls << endl;
        exit(-1);
    }
    MCcalls_ = MCcalls;
}

void ExactCoefficientFunction::SetDim(const int& dim) {
    // check dim
    if (dim <= 0) {
        cout << "Error: MCcalls must be positive. Got " << dim << endl;
        exit(-1);
    }
    dim_ = dim;
}

//==========================================================================================//
//  Exact massive gluon coefficient functions for F2 at O(alpha_s)
//
//  Eq. (50) from Ref. [arXiv:1001.2312]
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C2_g1(const double x, const double m2Q2) const { // m2Q2=m^2/Q^2

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

double ExactCoefficientFunction::CL_g1(const double x, const double m2Q2) const {

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

/// @cond UNNECESSARY
/**
 * @name Fortran massive coefficient functions
 * Fortran functions for the O(alpha_s^2)
 * coefficient functions from 'src/hqcoef.f'.
 */
///@{
extern "C" {
    double c2log_(double *wr,double *xi);
    double c2nlog_(double *wr, double *xi);
    double clnlog_(double *wr, double *xi);
    double c2nloq_(double *wr, double *xi);
    double clnloq_(double *wr, double *xi);
    // double c2nlobarg_(double *wr,double *xi);
    // double clnlobarg_(double *wr,double *xi);
    // double c2nlobarq_(double *wr,double *xi);
    // double clnlobarq_(double *wr,double *xi);
}
///@}
/// \endcond

//==========================================================================================//
//  Exact massive gluon coefficient functions for F2 at O(alpha_s^2)
//
//  Exact (but numerical) result from [arXiv:hep-ph/9411431].
//  Taken from the Fortran code 'src/hqcoef.f'
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C2_g2(const double x, const double m2Q2, const double m2mu2) const {

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

double ExactCoefficientFunction::C2_ps2(const double x, const double m2Q2, const double m2mu2) const {

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

double ExactCoefficientFunction::CL_g2(const double x, const double m2Q2, const double m2mu2) const {

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

double ExactCoefficientFunction::CL_ps2(const double x, const double m2Q2, const double m2mu2) const {

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

double ExactCoefficientFunction::C2_ps20(const double x, const double m2Q2) const {

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

double ExactCoefficientFunction::C2_ps21(const double x, const double m2Q2) const {

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

double ExactCoefficientFunction::CL_ps20(const double x, const double m2Q2) const {

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

double ExactCoefficientFunction::CL_ps21(const double x, const double m2Q2) const { return -CL_g1_x_Pgq0(x, m2Q2); }

//==========================================================================================//
//  Exact massive gluon coefficient functions for F2 at O(alpha_s^2):
//  mu independent term.
//
//  Eq. (4.4) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C2_g20(const double x, const double m2Q2) const {

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

double ExactCoefficientFunction::C2_g21(const double x, const double m2Q2) const {

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

double ExactCoefficientFunction::CL_g20(const double x, const double m2Q2) const {

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

double ExactCoefficientFunction::CL_g21(const double x, const double m2Q2) const {

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

double ExactCoefficientFunction::C2_ps31(const double x, const double m2Q2, const int nf) const {

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

double ExactCoefficientFunction::CL_ps31(const double x, const double m2Q2, const int nf) const {

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

double ExactCoefficientFunction::C2_ps32(const double x, const double m2Q2, const int nf) const {

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

double ExactCoefficientFunction::CL_ps32(const double x, const double m2Q2, const int nf) const {

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

double ExactCoefficientFunction::C2_g31(const double x, const double m2Q2, const int nf) const {

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

double ExactCoefficientFunction::CL_g31(const double x, const double m2Q2, const int nf) const {

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

double ExactCoefficientFunction::C2_g32(const double x, const double m2Q2, const int nf) const {

    double C2_g1xPgg0xPgg0;

    double beta0 = beta(0, nf);

    if (method_flag_ == 0)
        C2_g1xPgg0xPgg0 = C2_g1_x_Pgg0_x_Pgg0(x, m2Q2, nf);
    else if (method_flag_ == 1)
        C2_g1xPgg0xPgg0 = C2_g1_x_Pgg0_x_Pgg0_MC(x, m2Q2, nf);
    else {
        cout << "C2_g32: Choose either method_flag = 0 or method_flag = 1"
                  << endl;
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

double ExactCoefficientFunction::CL_g32(const double x, const double m2Q2, const int nf) const {

    double CL_g1xPgg0xPgg0;

    double beta0 = beta(0, nf);

    if (method_flag_ == 0)
        CL_g1xPgg0xPgg0 = CL_g1_x_Pgg0_x_Pgg0(x, m2Q2, nf);
    else if (method_flag_ == 1)
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
