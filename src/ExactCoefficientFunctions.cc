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

    SetMethodFlag(method_flag);
    SetAbserr(abserr);
    SetRelerr(relerr);
    SetMCcalls(MCcalls);
    SetDim(dim);

    gluon_lo_ = nullptr;
    gluon_nlo_ = nullptr;
    quark_nlo_ = nullptr;

    if (GetOrder() > 1) gluon_lo_ = new ExactCoefficientFunction(1, GetKind(), 'g');
    if (GetOrder() > 2) {
        gluon_nlo_ = new ExactCoefficientFunction(2, GetKind(), 'g');
        quark_nlo_ = new ExactCoefficientFunction(2, GetKind(), 'q');
    }

    if (GetOrder() == 2) {
        if (GetChannel() == 'q') {
            convolutions_lmu1.push_back( &Convolution(gluon_lo_, &SplittingFunction(0, 'g', 'q'), abserr, relerr, dim) );
        } else if (GetChannel() == 'g'){
            convolutions_lmu1.push_back( new Convolution(gluon_lo_, &SplittingFunction(0, 'g', 'g'), abserr, relerr, dim) );
            convolutions_lmu1.push_back( new Convolution(gluon_lo_, &((-1)*Delta())) );
        }
    } else if (GetOrder() == 3) {
        if (GetChannel() == 'q') {
            convolutions_lmu1.push_back( new Convolution(gluon_lo_, &SplittingFunction(1, 'g', 'q'), abserr, relerr, dim) );
            convolutions_lmu1.push_back( new Convolution(gluon_nlo_, &SplittingFunction(0, 'g', 'q'), abserr, relerr, dim) );
            convolutions_lmu1.push_back( new Convolution(quark_nlo_, &SplittingFunction(0, 'q', 'q'), abserr, relerr, dim) );
            convolutions_lmu1.push_back( new Convolution(quark_nlo_, &(-2*Delta())) );

            convolutions_lmu2.push_back( new Convolution(gluon_lo_, &(0.5*ConvolutedSplittingFunctions(1, 'g', 'g', 'g', 'q')), abserr, relerr, dim) );
            convolutions_lmu2.push_back( new Convolution(gluon_lo_, &(0.5*ConvolutedSplittingFunctions(1, 'q', 'q', 'g', 'q')), abserr, relerr, dim) );
            convolutions_lmu2.push_back( new Convolution(gluon_lo_, &(-3./2*SplittingFunction(0, 'g', 'q')), abserr, relerr, dim) );
        } else {
            convolutions_lmu1.push_back( new Convolution(gluon_lo_, &SplittingFunction(1, 'g', 'g'), abserr, relerr, dim) );
            convolutions_lmu1.push_back( new Convolution(gluon_lo_, &Delta()) );
            convolutions_lmu1.push_back( new Convolution(quark_nlo_, &SplittingFunction(0, 'q', 'g'), abserr, relerr, dim) );
            convolutions_lmu1.push_back( new Convolution(gluon_nlo_, &SplittingFunction(0, 'g', 'g'), abserr, relerr, dim) );
            convolutions_lmu1.push_back( new Convolution(gluon_nlo_, &Delta()) );

            convolutions_lmu2.push_back( new MonteCarloDoubleConvolution(gluon_lo_, &SplittingFunction(0, 'g', 'g'), abserr, relerr, dim, MCcalls) );
            convolutions_lmu2.push_back( new Convolution(gluon_lo_, &ConvolutedSplittingFunctions(0, 'g', 'q', 'q', 'g'), abserr, relerr, dim) );
            convolutions_lmu2.push_back( new Convolution(gluon_lo_, &SplittingFunction(0, 'g', 'g'), abserr, relerr, dim) );
            convolutions_lmu1.push_back( new Convolution(gluon_lo_, &Delta()) );
        }
    }

    SetFunctions();

}

double ExactCoefficientFunction::fx(const double x, const double m2Q2, const double m2mu2, const int nf) const {
    return MuIndependentTerms(x, m2Q2, nf) + MuDependentTerms(x, m2Q2, m2mu2, nf);
}

double ExactCoefficientFunction::MuIndependentTerms(const double x, const double m2Q2, const int nf) const {
    return (this->*mu_indep_)(x, m2Q2, nf);
}

double ExactCoefficientFunction::MuDependentTerms(const double x, const double m2Q2, const double m2mu2, const int nf) const {
    return (this->*mu_dep_)(x, m2Q2, m2mu2, nf);
}

void ExactCoefficientFunction::SetFunctions() {
    if (GetOrder() == 1) {
        if (GetKind() == '2') {
            mu_indep_ = &ExactCoefficientFunction::C2_g1 ;
        } else if (GetKind() =='L') {
            mu_indep_ = &ExactCoefficientFunction::CL_g1 ;
        }
        mu_dep_ = &ExactCoefficientFunction::ZeroFunction ;
    } else if (GetOrder() == 2) {
        if (GetChannel() == 'q') {
            if (GetKind() == '2') {
                mu_indep_ = &ExactCoefficientFunction::C2_ps20 ;
            } else if (GetKind() == 'L') {
                mu_indep_ = &ExactCoefficientFunction::CL_ps20 ;
            }
            mu_dep_ = &ExactCoefficientFunction::C_ps2_MuDep ;
        } else if (GetChannel() == 'g') {
            if (GetKind() == '2') {
                mu_indep_ = &ExactCoefficientFunction::C2_g20 ;
            } else if (GetKind() == 'L') {
                mu_indep_ = &ExactCoefficientFunction::CL_g20 ;
            }
            mu_dep_ = &ExactCoefficientFunction::C_g2_MuDep ;
        }
    } else if (GetOrder() == 3) {
        if (GetChannel() == 'q') {
            mu_dep_ = &ExactCoefficientFunction::C_ps3_MuDep ;
        } else if (GetChannel() == 'g') {
            mu_dep_ = &ExactCoefficientFunction::C_g3_MuDep ;
        }
        mu_indep_ = nullptr;
    }
}

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

// double ExactCoefficientFunction::C2_g2(const double x, const double m2Q2, const double m2mu2) const {

//     double xi = 1. / m2Q2;
//     double eta = 0.25 * xi * (1 - x) / x - 1.;

//     if (eta > 1e6 || eta < 1e-6 || xi < 1e-3 || xi > 1e5)
//         return 0.;

//     return C2_g20(x, m2Q2) + C2_g21(x, m2Q2) * log(1. / m2mu2);
// }

//==========================================================================================//
//  Exact massive quark coefficient functions for F2 at O(alpha_s)
//
//  Exact (but numerical) result from [arXiv:hep-ph/9411431].
//  Taken from the Fortran code 'src/hqcoef.f'
//------------------------------------------------------------------------------------------//

// double ExactCoefficientFunction::C2_ps2(const double x, const double m2Q2, const double m2mu2) const {

//     double xi = 1. / m2Q2;
//     double eta = 0.25 * xi * (1 - x) / x - 1.;

//     if (eta > 1e6 || eta < 1e-6 || xi < 1e-3 || xi > 1e5)
//         return 0.;

//     return C2_ps20(x, m2Q2) + C2_ps21(x, m2Q2) * log(1. / m2mu2);
// }

//==========================================================================================//
//  Exact massive gluon coefficient functions for FL at O(alpha_s)
//
//  Exact (but numerical) result from [arXiv:hep-ph/9411431].
//  Taken from the Fortran code 'src/hqcoef.f'
//------------------------------------------------------------------------------------------//

// double ExactCoefficientFunction::CL_g2(const double x, const double m2Q2, const double m2mu2) const {

//     double xi = 1. / m2Q2;
//     double eta = 0.25 * xi * (1 - x) / x - 1.;

//     if (eta > 1e6 || eta < 1e-6 || xi < 1e-3 || xi > 1e5)
//         return 0.;

//     return CL_g20(x, m2Q2) + CL_g21(x, m2Q2) * log(1. / m2mu2);
// }

//==========================================================================================//
//  Exact massive quarkk coefficient functions for FL at O(alpha_s)
//
//  Exact (but numerical) result from [arXiv:hep-ph/9411431].
//  Taken from the Fortran code 'src/hqcoef.f'
//------------------------------------------------------------------------------------------//

// double ExactCoefficientFunction::CL_ps2(const double x, const double m2Q2, const double m2mu2) const {

//     double xi = 1. / m2Q2;
//     double eta = 0.25 * xi * (1 - x) / x - 1.;

//     if (eta > 1e6 || eta < 1e-6 || xi < 1e-3 || xi > 1e5)
//         return 0.;

//     return CL_ps20(x, m2Q2) + CL_ps21(x, m2Q2) * log(1. / m2mu2);
// }

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

double ExactCoefficientFunction::C_ps21(const double x, const double m2Q2) const {

    int nf = static_cast<int>(nan(""));

    return - convolutions_lmu1[0] -> Convolute(x, m2Q2, nf);
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

// double ExactCoefficientFunction::CL_ps21(const double x, const double m2Q2) const { return -CL_g1_x_Pgq0(x, m2Q2); }

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

double ExactCoefficientFunction::C_g21(const double x, const double m2Q2) const {

    int nf = 1;
    // Put nf to 1 since the nf contribution cancels for any value of nf

    return -(convolutions_lmu1[0] -> Convolute(x, m2Q2, nf) + convolutions_lmu1[1] -> Convolute(x, m2Q2, nf) * beta(0, nf));
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

// double ExactCoefficientFunction::CL_g21(const double x, const double m2Q2) const {

//     int nf = 1;
//     // Put nf to 1 since the nf contribution cancels for any value of nf

//     return -(CL_g1_x_Pgg0(x, m2Q2, nf) - CL_g1(x, m2Q2) * beta(0, nf));
// }

double ExactCoefficientFunction::C_g2_MuDep(const double x, const double m2Q2, const double m2mu2) const {

    return C_g21(x, m2Q2) * log(1./m2mu2) ;
}

double ExactCoefficientFunction::C_ps2_MuDep(const double x, const double m2Q2, const double m2mu2) const {

    return C_ps21(x, m2Q2) * log(1./m2mu2) ;
}

//==========================================================================================//
//  Exact massive quark coefficient functions for F2 at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.2) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C_ps31(const double x, const double m2Q2, const int nf) const {

    return -(
        convolutions_lmu1[0] -> Convolute(x, m2Q2, nf) + convolutions_lmu1[1] -> Convolute(x, m2Q2, nf)
        + convolutions_lmu1[2] -> Convolute(x, m2Q2, nf) + beta(0, nf) * convolutions_lmu1[3] -> Convolute(x, m2Q2, nf)
    );
}

//==========================================================================================//
//  Exact massive quark coefficient functions for FL at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.2) of Ref. [arXiv:1205.5727] for FL
//------------------------------------------------------------------------------------------//

// double ExactCoefficientFunction::CL_ps31(const double x, const double m2Q2, const int nf) const {

//     return -(
//         CL_g1_x_Pgq1(x, m2Q2, nf) + CL_g20_x_Pgq0(x, m2Q2)
//         + CL_ps20_x_Pqq0(x, m2Q2, nf) - 2. * beta(0, nf) * CL_ps20(x, m2Q2)
//     );
// }

//==========================================================================================//
//  Exact massive quark coefficient functions for F2 at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)^2
//
//  Eq. (4.3) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C_ps32(const double x, const double m2Q2, const int nf) const {

    return convolutions_lmu2[0] -> Convolute(x, m2Q2, nf) + convolutions_lmu2[1] -> Convolute(x, m2Q2, nf)
           + beta(0, nf) * convolutions_lmu2[2] -> Convolute(x, m2Q2, nf);
}

//==========================================================================================//
//  Exact massive quark coefficient functions for FL at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)^2
//
//  Eq. (4.3) of Ref. [arXiv:1205.5727] for FL
//------------------------------------------------------------------------------------------//

// double ExactCoefficientFunction::CL_ps32(const double x, const double m2Q2, const int nf) const {

//     return 0.5
//                * (CL_g1_x_Pgg0_x_Pgq0(x, m2Q2, nf)
//                   + CL_g1_x_Pqq0_x_Pgq0(x, m2Q2, nf))
//            - 3. / 2 * beta(0, nf) * CL_g1_x_Pgq0(x, m2Q2);
// }

//==========================================================================================//
//  Exact massive gluon coefficient functions for F2 at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.5) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C_g31(const double x, const double m2Q2, const int nf) const {

    return -(
        convolutions_lmu1[0] -> Convolute(x, m2Q2, nf) + beta(1, nf) * convolutions_lmu1[1] -> Convolute(x, m2Q2, nf)
        + convolutions_lmu1[2] -> Convolute(x, m2Q2, nf) + convolutions_lmu1[3] -> Convolute(x, m2Q2, nf)
        + beta(0, nf) * convolutions_lmu1[4] -> Convolute(x, m2Q2, nf)
    );
}

//==========================================================================================//
//  Exact massive gluon coefficient functions for FL at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.5) of Ref. [arXiv:1205.5727] for FL
//------------------------------------------------------------------------------------------//

// double ExactCoefficientFunction::CL_g31(const double x, const double m2Q2, const int nf) const {

//     return -(
//         CL_g1_x_Pgg1(x, m2Q2, nf) - beta(1, nf) * CL_g1(x, m2Q2)
//         + CL_ps20_x_Pqg0(x, m2Q2, nf) + CL_g20_x_Pgg0(x, m2Q2, nf)
//         - 2 * beta(0, nf) * CL_g20(x, m2Q2)
//     );
// }

//==========================================================================================//
//  Exact massive gluon coefficient functions for F2 at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)^2
//
//  Eq. (4.6) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C_g32(const double x, const double m2Q2, const int nf) const {

    // double C2_g1xPgg0xPgg0;

    double beta0 = beta(0, nf);

    // if (method_flag_ == 0)
    //     C2_g1xPgg0xPgg0 = C2_g1_x_Pgg0_x_Pgg0(x, m2Q2, nf);
    // else if (method_flag_ == 1)
    //     C2_g1xPgg0xPgg0 = C2_g1_x_Pgg0_x_Pgg0_MC(x, m2Q2, nf);
    // else {
    //     cout << "C2_g32: Choose either method_flag = 0 or method_flag = 1"
    //               << endl;
    //     exit(-1);
    // }

    return convolutions_lmu2[0] -> Convolute(x, m2Q2, nf) + convolutions_lmu2[1] -> Convolute(x, m2Q2, nf)
           + beta0 * convolutions_lmu2[2] -> Convolute(x, m2Q2, nf)
           + beta0 * beta0 * convolutions_lmu2[3] -> Convolute(x, m2Q2, nf);
}

//==========================================================================================//
//  Exact massive gluon coefficient functions for FL at O(alpha_s^3):
//  Term proportional to log(mu^2/m^2)^2
//
//  Eq. (4.6) of Ref. [arXiv:1205.5727] for FL
//------------------------------------------------------------------------------------------//

// double ExactCoefficientFunction::CL_g32(const double x, const double m2Q2, const int nf) const {

//     double CL_g1xPgg0xPgg0;

//     double beta0 = beta(0, nf);

//     if (method_flag_ == 0)
//         CL_g1xPgg0xPgg0 = CL_g1_x_Pgg0_x_Pgg0(x, m2Q2, nf);
//     else if (method_flag_ == 1)
//         CL_g1xPgg0xPgg0 = CL_g1_x_Pgg0_x_Pgg0_MC(x, m2Q2, nf);
//     else {
//         std::cout << "C2_g32: Choose either method_flag = 0 or method_flag = 1"
//                   << std::endl;
//         exit(-1);
//     }

//     return 0.5 * CL_g1xPgg0xPgg0 + 0.5 * CL_g1_x_Pqg0_x_Pgq0(x, m2Q2, nf)
//            - 3. / 2 * beta0 * CL_g1_x_Pgg0(x, m2Q2, nf)
//            + beta0 * beta0 * CL_g1(x, m2Q2);
// }
