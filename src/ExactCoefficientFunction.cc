#include "adani/ExactCoefficientFunction.h"
#include "adani/Constants.h"
#include "adani/SpecialFunctions.h"

#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

//==========================================================================================//
//  ExactCoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

ExactCoefficientFunction::ExactCoefficientFunction(
    const int &order, const char &kind, const char &channel,
    const double &abserr, const double &relerr, const int &dim,
    const bool &MCintegral, const int &MCcalls
)
    : CoefficientFunction(order, kind, channel) {

    gluon_as1_ = nullptr;
    gluon_as2_ = nullptr;
    quark_as2_ = nullptr;

    Pgq0_ = nullptr;
    Pgg0_ = nullptr;
    Pgq1_ = nullptr;
    Pqq0_ = nullptr;
    Pgg0Pgq0_ = nullptr;
    Pqq0Pgq0_ = nullptr;
    Pgg1_ = nullptr;
    Pqg0_ = nullptr;
    Pgq0Pqg0_ = nullptr;

    delta_ = nullptr;

    if (GetOrder() > 1) {
        gluon_as1_ = new ExactCoefficientFunction(1, GetKind(), 'g');

        Pgq0_ = new SplittingFunction(0, 'g', 'q');
        Pgg0_ = new SplittingFunction(0, 'g', 'g');
        delta_ = new Delta();
    }
    if (GetOrder() > 2) {
        gluon_as2_ = new ExactCoefficientFunction(2, GetKind(), 'g');
        quark_as2_ = new ExactCoefficientFunction(2, GetKind(), 'q');

        Pgq1_ = new SplittingFunction(1, 'g', 'q');
        Pqq0_ = new SplittingFunction(0, 'q', 'q');
        Pgg0Pgq0_ = new ConvolutedSplittingFunctions(0, 'g', 'g', 0, 'g', 'q');
        Pqq0Pgq0_ = new ConvolutedSplittingFunctions(0, 'q', 'q', 0, 'g', 'q');
        Pgg1_ = new SplittingFunction(1, 'g', 'g');
        Pqg0_ = new SplittingFunction(0, 'q', 'g');
        Pgq0Pqg0_ = new ConvolutedSplittingFunctions(0, 'g', 'q', 0, 'q', 'g');
    }

    if (GetOrder() == 2) {
        if (GetChannel() == 'q') {
            convolutions_lmu1_.push_back(
                new Convolution(gluon_as1_, Pgq0_, abserr, relerr, dim)
            );
        } else if (GetChannel() == 'g') {
            convolutions_lmu1_.push_back(
                new Convolution(gluon_as1_, Pgg0_, abserr, relerr, dim)
            );
            convolutions_lmu1_.push_back(new Convolution(gluon_as1_, delta_));
        }
    } else if (GetOrder() == 3) {
        if (GetChannel() == 'q') {
            convolutions_lmu1_.push_back(
                new Convolution(gluon_as1_, Pgq1_, abserr, relerr, dim)
            );
            convolutions_lmu1_.push_back(
                new Convolution(gluon_as2_, Pgq0_, abserr, relerr, dim)
            );
            convolutions_lmu1_.push_back(
                new Convolution(quark_as2_, Pqq0_, abserr, relerr, dim)
            );
            convolutions_lmu1_.push_back(new Convolution(quark_as2_, delta_));

            convolutions_lmu2_.push_back(
                new Convolution(gluon_as1_, Pgg0Pgq0_, abserr, relerr, dim)
            );
            convolutions_lmu2_.push_back(
                new Convolution(gluon_as1_, Pqq0Pgq0_, abserr, relerr, dim)
            );
            convolutions_lmu2_.push_back(
                new Convolution(gluon_as1_, Pgq0_, abserr, relerr, dim)
            );
        } else {
            convolutions_lmu1_.push_back(
                new Convolution(gluon_as1_, Pgg1_, abserr, relerr, dim)
            );
            convolutions_lmu1_.push_back(new Convolution(gluon_as1_, delta_));
            convolutions_lmu1_.push_back(
                new Convolution(quark_as2_, Pqg0_, abserr, relerr, dim)
            );
            convolutions_lmu1_.push_back(
                new Convolution(gluon_as2_, Pgg0_, abserr, relerr, dim)
            );
            convolutions_lmu1_.push_back(new Convolution(gluon_as2_, delta_));

            convolutions_lmu2_.push_back(new DoubleConvolution(
                gluon_as1_, Pgg0_, abserr, relerr, dim, MCintegral, MCcalls
            ));
            convolutions_lmu2_.push_back(
                new Convolution(gluon_as1_, Pgq0Pqg0_, abserr, relerr, dim)
            );
            convolutions_lmu2_.push_back(
                new Convolution(gluon_as1_, Pgg0_, abserr, relerr, dim)
            );
            convolutions_lmu2_.push_back(new Convolution(gluon_as1_, delta_));
        }
    }

    SetFunctions();
}

//==========================================================================================//
//  ExactCoefficientFunction: destructor
//------------------------------------------------------------------------------------------//

ExactCoefficientFunction::~ExactCoefficientFunction() {

    delete gluon_as1_;
    delete gluon_as2_;
    delete quark_as2_;

    for (long unsigned int i = 0; i < convolutions_lmu1_.size(); i++) {
        delete convolutions_lmu1_[i];
    }

    for (long unsigned int i = 0; i < convolutions_lmu2_.size(); i++) {
        delete convolutions_lmu2_[i];
    }

    delete Pgq0_;
    delete Pgg0_;
    delete Pgq1_;
    delete Pqq0_;
    delete Pgg0Pgq0_;
    delete Pqq0Pgq0_;
    delete Pgg1_;
    delete Pqg0_;
    delete Pgq0Pqg0_;

    delete delta_;
}

//==========================================================================================//
//  ExactCoefficientFunction: central value of the exact coefficient function
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::fx(
    double x, double m2Q2, double m2mu2, int nf
) const {
    return MuIndependentTerms(x, m2Q2, nf)
           + MuDependentTerms(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  ExactCoefficientFunction: mu-independent terms
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::MuIndependentTerms(
    double x, double m2Q2, int nf
) const {
    return (this->*mu_indep_)(x, m2Q2, nf);
}

//==========================================================================================//
//  ExactCoefficientFunction: mu dependent terms
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::MuDependentTerms(
    double x, double m2Q2, double m2mu2, int nf
) const {

    double x_max = 1. / (1 + 4 * m2Q2);
    if (x <= 0 || x > x_max)
        return 0.;

    return (this->*mu_dep_)(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  ExactCoefficientFunction: implemented only because it is pure virtual in
//  base class. Returning three identical values
//------------------------------------------------------------------------------------------//

Value ExactCoefficientFunction::fxBand(
    double x, double m2Q2, double m2mu2, int nf
) const {
    return Value(fx(x, m2Q2, m2mu2, nf));
}

//==========================================================================================//
//  ExactCoefficientFunction: function that sets the pointer for mu_indep_ and
//  mu_dep_
//------------------------------------------------------------------------------------------//

void ExactCoefficientFunction::SetFunctions() {
    if (GetOrder() == 1) {
        if (GetKind() == '2') {
            mu_indep_ = &ExactCoefficientFunction::C2_g1;
        } else if (GetKind() == 'L') {
            mu_indep_ = &ExactCoefficientFunction::CL_g1;
        }
        mu_dep_ = &ExactCoefficientFunction::ZeroFunction;
    } else if (GetOrder() == 2) {
        if (GetChannel() == 'q') {
            if (GetKind() == '2') {
                mu_indep_ = &ExactCoefficientFunction::C2_ps20;
            } else if (GetKind() == 'L') {
                mu_indep_ = &ExactCoefficientFunction::CL_ps20;
            }
            mu_dep_ = &ExactCoefficientFunction::C_ps2_MuDep;
        } else if (GetChannel() == 'g') {
            if (GetKind() == '2') {
                mu_indep_ = &ExactCoefficientFunction::C2_g20;
            } else if (GetKind() == 'L') {
                mu_indep_ = &ExactCoefficientFunction::CL_g20;
            }
            mu_dep_ = &ExactCoefficientFunction::C_g2_MuDep;
        }
    } else if (GetOrder() == 3) {
        if (GetChannel() == 'q') {
            mu_dep_ = &ExactCoefficientFunction::C_ps3_MuDep;
        } else if (GetChannel() == 'g') {
            mu_dep_ = &ExactCoefficientFunction::C_g3_MuDep;
        }
        mu_indep_ = &ExactCoefficientFunction::WarningFunc;
    }
}

//==========================================================================================//
//  Exact massive gluon coefficient functions for F2 at O(as)
//
//  Eq. (50) from Ref. [arXiv:1001.2312]
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C2_g1(double x, double m2Q2, int /*nf*/)
    const { // m2Q2=m^2/Q^2

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
//  Exact massive gluon coefficient functions for FL at O(as)
//
//  Eq. (2.9) from Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double
ExactCoefficientFunction::CL_g1(double x, double m2Q2, int /*nf*/) const {

    double beta = sqrt(1. - 4. * m2Q2 * x / (1. - x));
    double x2 = x * x;
    double L = log((1. + beta) / (1. - beta));

    return 16. * TR * (x * (1. - x) * beta - 2. * x2 * m2Q2 * L);
}

//==========================================================================================//
//  OBSERVATION: in the O(as^2) exact coefficeint functions the
//  mu-independent part is an interpolation in a certain grid. When this
//  function is called for a (x,Q) value outside this grid, the value 0 is
//  returned. The mu-dependent part is defined everywhere. This means that
//  outside the grid one contribution is set to zero while the other is not.
//------------------------------------------------------------------------------------------//

/// @cond UNNECESSARY
/**
 * @name Fortran massive coefficient functions
 * Fortran functions for the O(as^2)
 * coefficient functions from 'src/hqcoef.f'.
 */
///@{
extern "C" {
    // double c2log_(double *wr,double *xi);
    double c2nlog_(double *wr, double *xi);
    double clnlog_(double *wr, double *xi);
    double c2nloq_(double *wr, double *xi);
    double clnloq_(double *wr, double *xi);
    // double c2nlobarg_(double *wr, double *xi);
    // double clnlobarg_(double *wr, double *xi);
    // double c2nlobarq_(double *wr, double *xi);
    // double clnlobarq_(double *wr, double *xi);
}
///@}
/// \endcond

//==========================================================================================//
//  mu independent part of the exact massive quark coefficient functions for F2
//  at O(as^2)
//
//  Exact (but numerical) result from [arXiv:hep-ph/9411431].
//  Taken from the Fortran code 'src/hqcoef.f'
//------------------------------------------------------------------------------------------//

double
ExactCoefficientFunction::C2_ps20(double x, double m2Q2, int /*nf*/) const {

    double xi = 1. / m2Q2;
    double eta = 0.25 * xi * (1 - x) / x - 1.;

    if (eta > 1e6 || eta < 1e-6 || xi < 1e-3 || xi > 1e5)
        return 0.;

    return 16 * M_PI * xi * c2nloq_(&eta, &xi) / x;
}

//==========================================================================================//
//  Exact massive quark coefficient functions at O(as^2):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.1) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C_ps21(double x, double m2Q2) const {

    int nf = static_cast<int>(nan(""));

    return -convolutions_lmu1_[0]->Convolute(x, m2Q2, nf);
    // The minus sign comes from the fact that in [arXiv:1205.5727]
    // the expansion is performed in terms of log(m^2/mu^2)
    // (even if it says the opposite) but we are
    // expanding in terms of log(mu^2/m^2)
}

//==========================================================================================//
//  mu independent part of the exact massive quark coefficient functions for FL
//  at O(as)
//
//  Exact (but numerical) result from [arXiv:hep-ph/9411431].
//  Taken from the Fortran code 'src/hqcoef.f'
//------------------------------------------------------------------------------------------//

double
ExactCoefficientFunction::CL_ps20(double x, double m2Q2, int /*nf*/) const {

    double xi = 1. / m2Q2;
    double eta = 0.25 * xi * (1 - x) / x - 1.;

    if (eta > 1e6 || eta < 1e-6 || xi < 1e-3 || xi > 1e5)
        return 0.;

    return 16 * M_PI * xi * clnloq_(&eta, &xi) / x;
}

//==========================================================================================//
//  mu independent part of the exact massive gluon coefficient functions for F2
//  at O(as^2)
//
//  Exact (but numerical) result from [arXiv:hep-ph/9411431].
//  Taken from the Fortran code 'src/hqcoef.f'
//------------------------------------------------------------------------------------------//

double
ExactCoefficientFunction::C2_g20(double x, double m2Q2, int /*nf*/) const {

    double xi = 1. / m2Q2;
    double eta = 0.25 * xi * (1 - x) / x - 1.;

    if (eta > 1e6 || eta < 1e-6 || xi < 1e-3 || xi > 1e5)
        return 0.;

    return 16 * M_PI * xi * c2nlog_(&eta, &xi) / x;
}

//==========================================================================================//
//  Exact massive gluon coefficient functions at O(as^2):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.4) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C_g21(double x, double m2Q2) const {

    int nf_one = 1;
    // Put nf to 1 since the nf contribution cancels for any value of nf

    return -(
        convolutions_lmu1_[0]->Convolute(x, m2Q2, nf_one)
        - convolutions_lmu1_[1]->Convolute(x, m2Q2, nf_one) * beta(0, nf_one)
    );
}

//==========================================================================================//
//  mu independent part of the exact massive gluon coefficient functions for FL
//  at O(as^2)
//
//  Exact (but numerical) result from [arXiv:hep-ph/9411431].
//  Taken from the Fortran code 'src/hqcoef.f'
//------------------------------------------------------------------------------------------//

double
ExactCoefficientFunction::CL_g20(double x, double m2Q2, int /*nf*/) const {

    double xi = 1. / m2Q2;
    double eta = 0.25 * xi * (1 - x) / x - 1.;

    if (eta > 1e6 || eta < 1e-6 || xi < 1e-3 || xi > 1e5)
        return 0.;

    return 16 * M_PI * xi * clnlog_(&eta, &xi) / x;
}

//==========================================================================================//
//  ExactCoefficientFunction: mu dependent part of the exact massive gluon
//  coefficient function at O(as^2):
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::
    C_g2_MuDep(double x, double m2Q2, double m2mu2, int /*nf*/) const {

    return C_g21(x, m2Q2) * log(1. / m2mu2);
}

//==========================================================================================//
//  ExactCoefficientFunction: mu dependent part of the exact massive quark
//  coefficient function at O(as^2):
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::
    C_ps2_MuDep(double x, double m2Q2, double m2mu2, int /*nf*/) const {

    return C_ps21(x, m2Q2) * log(1. / m2mu2);
}

//==========================================================================================//
//  Exact massive quark coefficient functions at O(as^3):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.2) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C_ps31(double x, double m2Q2, int nf) const {

    return -(
        convolutions_lmu1_[0]->Convolute(x, m2Q2, nf)
        + convolutions_lmu1_[1]->Convolute(x, m2Q2, nf)
        + convolutions_lmu1_[2]->Convolute(x, m2Q2, nf)
        - 2. * beta(0, nf) * convolutions_lmu1_[3]->Convolute(x, m2Q2, nf)
    );
}

//==========================================================================================//
//  Exact massive quark coefficient functions at O(as^3):
//  Term proportional to log(mu^2/m^2)^2
//
//  Eq. (4.3) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C_ps32(double x, double m2Q2, int nf) const {

    return 0.5
               * (convolutions_lmu2_[0]->Convolute(x, m2Q2, nf)
                  + convolutions_lmu2_[1]->Convolute(x, m2Q2, nf))
           - 3. / 2 * beta(0, nf)
                 * convolutions_lmu2_[2]->Convolute(x, m2Q2, nf);
}

//==========================================================================================//
//  Exact massive gluon coefficient functions for F2 at O(as^3):
//  Term proportional to log(mu^2/m^2)
//
//  Eq. (4.5) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C_g31(double x, double m2Q2, int nf) const {

    return -(
        convolutions_lmu1_[0]->Convolute(x, m2Q2, nf)
        - beta(1, nf) * convolutions_lmu1_[1]->Convolute(x, m2Q2, nf)
        + convolutions_lmu1_[2]->Convolute(x, m2Q2, nf)
        + convolutions_lmu1_[3]->Convolute(x, m2Q2, nf)
        - 2. * beta(0, nf) * convolutions_lmu1_[4]->Convolute(x, m2Q2, nf)
    );
}

//==========================================================================================//
//  Exact massive gluon coefficient functions at O(as^3):
//  Term proportional to log(mu^2/m^2)^2
//
//  Eq. (4.6) of Ref. [arXiv:1205.5727]
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C_g32(double x, double m2Q2, int nf) const {

    double beta0 = beta(0, nf);

    return 0.5
               * (convolutions_lmu2_[0]->Convolute(x, m2Q2, nf)
                  + convolutions_lmu2_[1]->Convolute(x, m2Q2, nf))
           - 3. / 2 * beta0 * convolutions_lmu2_[2]->Convolute(x, m2Q2, nf)
           + beta0 * beta0 * convolutions_lmu2_[3]->Convolute(x, m2Q2, nf);
}

//==========================================================================================//
//  ExactCoefficientFunction: mu dependent part of the exact massive gluon
//  coefficient function at O(as^3):
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C_g3_MuDep(
    double x, double m2Q2, double m2mu2, int nf
) const {

    double lmu = log(1. / m2mu2);

    return C_g31(x, m2Q2, nf) * lmu + C_g32(x, m2Q2, nf) * lmu * lmu;
}

//==========================================================================================//
//  ExactCoefficientFunction: mu dependent part of the exact massive quark
//  coefficient function at O(as^3):
//------------------------------------------------------------------------------------------//

double ExactCoefficientFunction::C_ps3_MuDep(
    double x, double m2Q2, double m2mu2, int nf
) const {

    double lmu = log(1. / m2mu2);

    return C_ps31(x, m2Q2, nf) * lmu + C_ps32(x, m2Q2, nf) * lmu * lmu;
}

//==========================================================================================//
//  ExactCoefficientFunction: print warning for mu independent part of O(as^3)
//------------------------------------------------------------------------------------------//

double
ExactCoefficientFunction::WarningFunc(double /*x*/, double /*m2Q2*/, int /*nf*/)
    const {
    cout << "Error: mu-independent terms of the exact coefficient function at "
            "O(a_s^3) are not known!"
         << endl;
    exit(-1);
}
