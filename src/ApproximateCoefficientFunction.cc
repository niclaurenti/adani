#include "adani/ApproximateCoefficientFunction.h"

#include "adani/Constants.h"
#include "adani/HighEnergyCoefficientFunction.h"
#include "adani/HighScaleCoefficientFunction.h"
#include "adani/SpecialFunctions.h"

#include <cmath>
#include <iostream>
#include <vector>

using std::cout;
using std::endl;
using std::vector;

//==========================================================================================//
//  AbstractApproximate: constructor
//------------------------------------------------------------------------------------------//

AbstractApproximate::AbstractApproximate(
    const int &order, const char &kind, const char &channel,
    const double &abserr, const double &relerr, const int &dim,
    const bool &MCintegral, const int &MCcalls
)
    : CoefficientFunction(order, kind, channel) {

    muterms_ = new ExactCoefficientFunction(
        order, kind, channel, abserr, relerr, dim, MCintegral, MCcalls
    );
}

//==========================================================================================//
//  AbstractApproximate: destructor
//------------------------------------------------------------------------------------------//

AbstractApproximate::~AbstractApproximate() { delete muterms_; }

double
AbstractApproximate::MuIndependentTerms(double x, double m2Q2, int nf) const {
    return MuIndependentTermsBand(x, m2Q2, nf).GetCentral();
}

//==========================================================================================//
//  AbstractApproximate: band of the total result
//------------------------------------------------------------------------------------------//

Value AbstractApproximate::fxBand(
    double x, double m2Q2, double m2mu2, int nf
) const {
    double x_max = 1. / (1. + 4 * m2Q2);
    if (x >= x_max || x <= 0)
        return 0;

    return MuIndependentTermsBand(x, m2Q2, nf)
           + MuDependentTerms(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  AbstractApproximate: mu dependent terms
//------------------------------------------------------------------------------------------//

double AbstractApproximate::MuDependentTerms(
    double x, double m2Q2, double m2mu2, int nf
) const {
    return muterms_->MuDependentTerms(x, m2Q2, m2mu2, nf);
}

//==========================================================================================//
//  ApproximateCoefficientFunction: parameters of the approximation
//------------------------------------------------------------------------------------------//

#define a 2.5
#define b 5.

struct approximation_parameters C2_g1_params = { 0.2, 2.5, 2.5, 1.2 };
struct approximation_parameters CL_g1_params = { 20., 11., 3., 2. };

struct approximation_parameters C2_g2_params = { 1.7, 2.5, 2.5, 1.2 };
struct approximation_parameters CL_g2_params = { 20., 11., 3., 2. };
struct approximation_parameters C2_ps2_params = { 1.7, 2.5, 2.5, 1.2 };
struct approximation_parameters CL_ps2_params = { 20., 11., 3., 2. };

struct approximation_parameters C2_g3_params = { 0.3, 2.5, 2.5, 1.2 };
struct approximation_parameters CL_g3_params = { 10., 11., 3., 2. };
struct approximation_parameters C2_ps3_params = { 0.3, 2.5, 2.5, 1.2 };
struct approximation_parameters CL_ps3_params = { 20., 11., 3., 2. };

struct variation_parameters C2_var = { 0.3, 3. };
struct variation_parameters CL_var = { 0.2, 2. };

//==========================================================================================//
//  ApproximateCoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

ApproximateCoefficientFunction::ApproximateCoefficientFunction(
    const int &order, const char &kind, const char &channel, const bool &NLL,
    const string &highscale_version, const double &abserr, const double &relerr,
    const int &dim, const bool &MCintegral, const int &MCcalls
)
    : AbstractApproximate(
        order, kind, channel, abserr, relerr, dim, MCintegral, MCcalls
    ) {

    threshold_ = new ThresholdCoefficientFunction(order, kind, channel);
    asymptotic_ = new AsymptoticCoefficientFunction(
        order, kind, channel, NLL, highscale_version
    );

    if (order == 1) {
        if (kind == '2') {
            if (channel == 'g')
                approximation_ = C2_g1_params;
            variation_ = C2_var;
        } else if (kind == 'L') {
            if (channel == 'g')
                approximation_ = CL_g2_params;
            variation_ = CL_var;
        }
    }
    if (order == 2) {
        if (kind == '2') {
            if (channel == 'g')
                approximation_ = C2_g2_params;
            else if (channel == 'q')
                approximation_ = C2_ps2_params;
            variation_ = C2_var;
        } else if (kind == 'L') {
            if (channel == 'g')
                approximation_ = CL_g2_params;
            else if (channel == 'q')
                approximation_ = CL_ps2_params;
            variation_ = CL_var;
        }
    } else if (order == 3) {
        if (kind == '2') {
            if (channel == 'g')
                approximation_ = C2_g3_params;
            else if (channel == 'q')
                approximation_ = C2_ps3_params;
            else {
                cout << "Error" << endl;
                exit(-1);
            }
            variation_ = C2_var;
        } else if (kind == 'L') {
            if (channel == 'g')
                approximation_ = CL_g3_params;
            else if (channel == 'q')
                approximation_ = CL_ps3_params;
            variation_ = CL_var;
        }
    }
}

//==========================================================================================//
//  ApproximateCoefficientFunction: destructor
//------------------------------------------------------------------------------------------//

ApproximateCoefficientFunction::~ApproximateCoefficientFunction() {
    delete threshold_;
    delete asymptotic_;
}

//==========================================================================================//
//  ApproximateCoefficientFunction: band of the approximate mu independent terms
//------------------------------------------------------------------------------------------//

Value ApproximateCoefficientFunction::MuIndependentTermsBand(
    double x, double m2Q2, int nf
) const {

    double x_max = 1. / (1 + 4 * m2Q2);
    if (x <= 0 || x > x_max)
        return 0.;

    double A = approximation_.A, B = approximation_.B, C = approximation_.C,
           D = approximation_.D;
    double var = variation_.var, fact = variation_.fact;

    double Amax = fact * A, Amin = A / fact, Bmax = B * fact, Bmin = B / fact;
    double Cmax = (1. + var) * C, Cmin = (1. - var) * C, Dmax = (1. + var) * D,
           Dmin = (1. - var) * D;

    double Avec[3] = { A, Amax, Amin };
    double Bvec[3] = { B, Bmax, Bmin };
    double Cvec[3] = { C, Cmax, Cmin };
    double Dvec[3] = { D, Dmax, Dmin };

    vector<double> asy =
        (asymptotic_->MuIndependentTermsBand(x, m2Q2, nf)).ToVect();
    vector<double> thresh =
        (threshold_->MuIndependentTermsBand(x, m2Q2, nf)).ToVect();

    double central = Approximation(x, m2Q2, asy[0], thresh[0], A, B, C, D);
    double higher = central, lower = central, tmp;

    for (int i = 0; i < int(asy.size()); i++) {
        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    for (int m = 0; m < 3; m++) {
                        for (int n = 0; n < 3; n++) {
                            tmp = Approximation(
                                x, m2Q2, asy[i], thresh[j], Avec[k], Bvec[l],
                                Cvec[m], Dvec[n]
                            );
                            if (tmp > higher)
                                higher = tmp;
                            if (tmp < lower)
                                lower = tmp;
                        }
                    }
                }
            }
        }
    }

    return Value(central, higher, lower);
}

//==========================================================================================//
//  ApproximateCoefficientFunction: functional form of the approximation
//------------------------------------------------------------------------------------------//

double ApproximateCoefficientFunction::Approximation(
    double x, double m2Q2, double asy, double thresh, double A, double B,
    double C, double D
) const {

    double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;
    double xi = 1. / m2Q2;

    double h = A + (B - A) / (1. + exp(a * (log(xi) - b)));
    double k = C + (D - C) / (1. + exp(a * (log(xi) - b)));

    double damp_thr = 1. / (1. + pow(eta / h, k));
    double damp_asy = 1. - damp_thr;

    return asy * damp_asy + thresh * damp_thr;
}

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: parameters of the approximation
//------------------------------------------------------------------------------------------//

struct klmv_params klmv_C2g2A = { 1., 42.5, 0, 0, 0 };
struct klmv_params klmv_C2g2B = { 0.8, 19.4, 0, 0, 0 };
struct klmv_params klmv_C2q2A = { 1., 42.5, 0, 0, 0 };
struct klmv_params klmv_C2q2B = { 0.8, 19.4, 0, 0, 0 };

struct klmv_params klmv_C2g3A = { 1., 20., 0.007, 4, 0.28 };
struct klmv_params klmv_C2g3B = { 0.8, 10.7, 0.055, 2, 0.423 };
struct klmv_params klmv_C2q3A = { 1., 20., 0.004, 4, 0.125 };
struct klmv_params klmv_C2q3B = { 0.8, 10.7, 0.0245, 2, 0.17 };
struct klmv_params klmv_C2g3B_lowxi = { 0.8, 10.7, 0.055125, 2, 0.3825 };

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: constructor
//------------------------------------------------------------------------------------------//

ApproximateCoefficientFunctionKLMV::ApproximateCoefficientFunctionKLMV(
    const int &order, const char &kind, const char &channel,
    const string &highscale_version, const bool &lowxi, const double &abserr,
    const double &relerr, const int &dim, const bool &MCintegral,
    const int &MCcalls
)
    : AbstractApproximate(
        order, kind, channel, abserr, relerr, dim, MCintegral, MCcalls
    ) {
    if (GetOrder() == 1) {
        cout << "Error: KLMV approximation is not implemented at O(as)!"
             << endl;
        exit(-1);
    }
    if (GetKind() == 'L') {
        cout << "Error: KLMV approximation is not implemented for kind = 'L'!"
             << endl;
        exit(-1);
    }

    if (GetOrder() == 2) {
        if (GetChannel() == 'g') {
            params_A_ = klmv_C2g2A;
            params_B_ = klmv_C2g2B;
        } else if (GetChannel() == 'q') {
            params_A_ = klmv_C2q2A;
            params_B_ = klmv_C2q2B;
        }
    } else if (GetOrder() == 3) {
        if (GetChannel() == 'g') {
            params_A_ = klmv_C2g3A;
            if (lowxi)
                params_B_ = klmv_C2g3B_lowxi;
            else
                params_B_ = klmv_C2g3B;
        } else if (GetChannel() == 'q') {
            params_A_ = klmv_C2q3A;
            params_B_ = klmv_C2q3B;
        }
    }

    threshold_ = new ThresholdCoefficientFunction(order, kind, channel);

    highscale_ = new HighScaleCoefficientFunction(
        order, kind, channel, highscale_version
    );

    highenergy_ =
        new HighEnergyCoefficientFunction(order, kind, channel, false);
}

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: destructor
//------------------------------------------------------------------------------------------//

ApproximateCoefficientFunctionKLMV::~ApproximateCoefficientFunctionKLMV() {
    delete threshold_;
    delete highscale_;
    delete highenergy_;
}

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: band of the mu independent terms
//------------------------------------------------------------------------------------------//

Value ApproximateCoefficientFunctionKLMV::MuIndependentTermsBand(
    double x, double m2Q2, int nf
) const {

    double xmax = 1. / (1. + 4 * m2Q2);
    if (x <= 0 || x >= xmax)
        return 0.;

    double thr = threshold_->MuIndependentTerms(x, m2Q2, nf);
    double thr_const = threshold_->BetaIndependentTerms(x, m2Q2, 1.);

    double he_ll = highenergy_->MuIndependentTerms(x, m2Q2, nf);

    vector<double> hs = highscale_->fxBand_NotOrdered(x, m2Q2, 1., nf);

    double gamma = params_A_.gamma, C = params_A_.C, delta = params_B_.gamma,
           D = params_B_.C;

    double res_A, res_B;

    if (GetOrder() == 2) {
        res_A =
            ApproximationA(x, m2Q2, 0., he_ll, hs[0], thr, thr_const, gamma, C);
        res_B = ApproximationB(x, m2Q2, 0., he_ll, hs[0], thr, 0., delta, D);
    } else if (GetOrder() == 3) {
        Value he_nll = ApproximateNLL(x, m2Q2);
        res_A = ApproximationA(
            x, m2Q2, he_ll, he_nll.GetHigher(), hs[1], thr, thr_const, gamma, C
        );
        res_B = ApproximationB(
            x, m2Q2, he_ll, he_nll.GetLower(), hs[2], thr, thr_const, delta, D
        );
    }

    if (res_A > res_B)
        return Value(res_A, res_B);
    else
        return Value(res_B, res_A);
}

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: functional form of the approximation A
//------------------------------------------------------------------------------------------//

double ApproximateCoefficientFunctionKLMV::ApproximationA(
    double x, double m2Q2, double he_logx, double he_const, double hs,
    double thr, double thr_const, double gamma, double C
) const {
    double beta = sqrt(1. - 4 * m2Q2 * x / (1. - x));
    double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

    double xi = 1. / m2Q2;
    double f = 1. / (1. + exp(2. * (xi - 4.)));

    double eta_gamma = pow(eta, gamma);

    double beta3 = beta * beta * beta;

    return thr - thr_const + (1. - f) * beta * hs
           + f * beta3
                 * (-log(eta) / log(x) * he_logx
                    + he_const * eta_gamma / (C + eta_gamma));
}

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: functional form of the approximation B
//------------------------------------------------------------------------------------------//

double ApproximateCoefficientFunctionKLMV::ApproximationB(
    double x, double m2Q2, double he_logx, double he_const, double hs,
    double thr, double thr_const, double delta, double D
) const {

    double beta = sqrt(1. - 4 * m2Q2 * x / (1. - x));
    double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

    double xi = 1. / m2Q2;
    double f = 1. / (1. + exp(2. * (xi - 4.)));

    double eta_gamma = pow(eta, delta);

    double beta3 = beta * beta * beta;

    return thr - thr_const + 2 * f * thr_const + (1. - f) * beta3 * hs
           + f * beta3
                 * (-log(eta) / log(x) * he_logx
                    + he_const * eta_gamma / (D + eta_gamma));
}

//==========================================================================================//
//  ApproximateCoefficientFunctionKLMV: NLL coefficient of the small-x
//  approximation
//------------------------------------------------------------------------------------------//

Value ApproximateCoefficientFunctionKLMV::ApproximateNLL(
    double x, double m2Q2
) const {
    double pi3 = M_PI * M_PI * M_PI;
    double tmp_A =
        (64. * pi3)
        * (params_A_.log_coeff * pow(log(1. / m2Q2) / log(5), params_A_.log_pow)
           - params_A_.const_coeff)
        * 4. / m2Q2 / x;
    double tmp_B =
        (64. * pi3)
        * (params_B_.log_coeff * pow(log(1. / m2Q2) / log(5), params_B_.log_pow)
           - params_B_.const_coeff)
        * 4. / m2Q2 / x;
    return Value(tmp_A, tmp_B);
}
