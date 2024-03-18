#include "adani/ApproximateCoefficientFunctions.h"
#include "adani/Constants.h"

#include "adani/HighEnergyCoefficientFunctions.h"
#include "adani/HighScaleCoefficientFunctions.h"
#include "adani/SpecialFunctions.h"

#include <cmath>
#include <iostream>
#include <vector>

using std::cout;
using std::endl;
using std::vector;

AbstractApproximate::AbstractApproximate(const int& order, const char& kind, const char& channel, const double& abserr, const double& relerr, const int& dim, const int& method_flag, const int& MCcalls) : CoefficientFunction(order, kind, channel) {
    
    muterms_ = new ExactCoefficientFunction(order, kind, channel, abserr, relerr, dim, method_flag, MCcalls);
}

AbstractApproximate::~AbstractApproximate() {
    delete muterms_;
}

double AbstractApproximate::MuIndependentTerms(double x, double m2Q2, int nf) const {
    return MuIndependentTermsBand(x, m2Q2, nf).GetCentral();
}

Value AbstractApproximate::fxBand(double x, double m2Q2, double m2mu2, int nf) const  {
    double x_max = 1. / (1. + 4 * m2Q2);
    if (x >= x_max || x < 0)
        return 0;

    return MuIndependentTermsBand(x, m2Q2, nf) + MuDependentTerms(x, m2Q2, m2mu2, nf);
}

double AbstractApproximate::MuDependentTerms(double x, double m2Q2, double m2mu2, int nf) const {
    return muterms_ -> MuDependentTerms(x, m2Q2, m2mu2, nf);
}

#define a 2.5
#define b 5.

struct approximation_parameters C2_g1_params = {0.2, 2.5, 2.5, 1.2};
struct approximation_parameters CL_g1_params = {20., 11., 3., 2.};

struct approximation_parameters C2_g2_params = {1.7, 2.5, 2.5, 1.2};
struct approximation_parameters CL_g2_params = {20., 11., 3., 2.};
struct approximation_parameters C2_ps2_params = {1.7, 2.5, 2.5, 1.2};
struct approximation_parameters CL_ps2_params = {20., 11., 3., 2.};

struct approximation_parameters C2_g3_params = {0.3, 2.5, 2.5, 1.2};
struct approximation_parameters CL_g3_params = {10., 11., 3., 2.};
struct approximation_parameters C2_ps3_params = {0.3, 2.5, 2.5, 1.2};
struct approximation_parameters CL_ps3_params = {20., 11., 3., 2.};

struct variation_parameters C2_var = {0.3, 3.};
struct variation_parameters CL_var = {0.2, 2.};

ApproximateCoefficientFunction::ApproximateCoefficientFunction(const int& order, const char& kind, const char& channel, const bool& NLL, const bool& exact_highscale, const bool& revised_approx_highscale, const double& abserr, const double& relerr, const int& dim, const int& method_flag, const int& MCcalls) : AbstractApproximate(order, kind, channel, abserr, relerr, dim, method_flag, MCcalls) {
    
    threshold_ = new ThresholdCoefficientFunction(order, kind, channel);
    asymptotic_ = new AsymptoticCoefficientFunction(order, kind, channel, NLL, exact_highscale, revised_approx_highscale);

    if (order == 1) {
        if (kind == '2'){
            if (channel == 'g') approximation_ = C2_g1_params;
            variation_ = C2_var ;
        } else if (kind == 'L'){
            if (channel == 'g') approximation_ = CL_g2_params;
            variation_ = CL_var ;
        } 
    }
    if (order == 2) {
        if (kind == '2'){
            if (channel == 'g') approximation_ = C2_g2_params;
            else if (channel == 'q') approximation_ = C2_ps2_params;
            variation_ = C2_var ;
        } else if (kind == 'L'){
            if (channel == 'g') approximation_ = CL_g2_params;
            else if (channel == 'q') approximation_ = CL_ps2_params;
            variation_ = CL_var ;
        } 
    } else if (order == 3) {
        if (kind == '2') {
            if (channel == 'g') approximation_ = C2_g3_params;
            else if (channel == 'q') approximation_ = C2_ps3_params;
            else {
                cout << "Error" << endl;
                exit(-1);
            }
            variation_ = C2_var ;
        } else if (kind == 'L') {
            if (channel == 'g') approximation_ = CL_g3_params;
            else if (channel == 'q') approximation_ = CL_ps3_params;
            variation_ = CL_var ;
        } 
    }
}

ApproximateCoefficientFunction::~ApproximateCoefficientFunction() {
    delete threshold_;
    delete asymptotic_;
}

Value ApproximateCoefficientFunction::MuIndependentTermsBand(double x, double m2Q2, int nf) const {

    double A = approximation_.A, B = approximation_.B, C = approximation_.C, D = approximation_.D; 
    double var = variation_.var, fact = variation_.fact;

    double Amax = fact * A, Amin = A / fact, Bmax = B * fact, Bmin = B / fact;
    double Cmax = (1. + var) * C, Cmin = (1. - var) * C, Dmax = (1. + var) * D, Dmin = (1. - var) * D;

    double Avec[3] = { A, Amin, Amax };
    double Bvec[3] = { B, Bmax, Bmin };
    double Cvec[3] = { C, Cmax, Cmin };
    double Dvec[3] = { D, Dmax, Dmin };

    vector<double> asy = (asymptotic_ -> MuIndependentTermsBand(x, m2Q2, nf)).ToVect();
    vector<double> thresh  = (threshold_ -> MuIndependentTermsBand(x, m2Q2, nf)).ToVect();

    double central = Approximation(x, m2Q2, asy[0], thresh[0], A, B, C, D);
    double higher = central, lower = central, tmp;

    for(int i = 0; i<3; i++) {
        for(int j = 0; j<3; j++) {
            for(int k = 0; k<3; k++) {
                for(int l = 0; l<3; l++) {
                    for(int m = 0; m<3; m++) {
                        for(int n = 0; n<3; n++) {
                            tmp = Approximation(x, m2Q2, asy[i], thresh[j], Avec[k], Bvec[l], Cvec[m], Dvec[n]) ;
                            if (tmp > higher) higher = tmp;
                            if (tmp < lower) lower = tmp;
                        }
                    }
                }
            }
        }
    }

    return Value(central, higher, lower);

}

double ApproximateCoefficientFunction::Approximation(double x, double m2Q2, double asy, double thresh, double A, double B, double C, double D) const {
    
    double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;
    double xi = 1. / m2Q2;

    double h = A + (B - A) / (1. + exp(a * (log(xi) - b)));
    double k = C + (D - C) / (1. + exp(a * (log(xi) - b)));

    double damp_thr = 1. / (1. + pow(eta / h, k));
    double damp_asy = 1. - damp_thr;

    return asy * damp_asy + thresh * damp_thr;
}

struct klmv_params klmv_C2g2A = {1., 42.5, 0, 0, 0};
struct klmv_params klmv_C2g2B = {0.8, 19.4, 0, 0, 0};
struct klmv_params klmv_C2q2A = {1., 42.5, 0, 0, 0};
struct klmv_params klmv_C2q2B = {0.8, 19.4, 0, 0, 0};

struct klmv_params klmv_C2g3A = {1., 20., 0.007, 4, 0.28};
struct klmv_params klmv_C2g3B = {0.8, 10.7, 0.055, 2, 0.423};
struct klmv_params klmv_C2q3A = {1., 20., 0.004, 4, 0.125};
struct klmv_params klmv_C2q3B = {0.8, 10.7, 0.0245, 2, 0.17};
// struct klmv_params klmv_C2gB_lowxi = {0.8, 10.7, 0.055, 2, 0.423};

ApproximateCoefficientFunctionKLMV::ApproximateCoefficientFunctionKLMV(const int& order, const char& kind, const char& channel, const bool& revised_approx_highscale, const double& abserr, const double& relerr, const int& dim, const int& method_flag, const int& MCcalls) : AbstractApproximate(order, kind, channel, abserr, relerr, dim, method_flag, MCcalls) {
    if (GetOrder() == 1) {
        cout << "Error: KLMV approximation is not implemented at O(as)!" << endl;
        exit(-1);
    }
    if (GetKind() == 'L') {
        cout << "Error: KLMV approximation is not implemented for kind = 'L'!" << endl;
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
            params_B_ = klmv_C2g3B;
        } else if (GetChannel() == 'q') {
            params_A_ = klmv_C2q3A;
            params_B_ = klmv_C2q3B;
        }
    }

    threshold_ = new ThresholdCoefficientFunction(order, kind, channel);
    highscale_ = new HighScaleCoefficientFunction(order, kind, channel, false, revised_approx_highscale);
    highenergy_ = new HighEnergyCoefficientFunction(order, kind, channel, false);

}

ApproximateCoefficientFunctionKLMV::~ApproximateCoefficientFunctionKLMV() {
    delete threshold_;
    delete highscale_;
    delete threshold_;
}

Value ApproximateCoefficientFunctionKLMV::MuIndependentTermsBand(double x, double m2Q2, int nf) const {
    
    double thr = threshold_ -> MuIndependentTerms(x, m2Q2, nf);
    double thr_const = threshold_ -> BetaIndependentTerms(x, m2Q2, 1.);
    
    double he_ll = highenergy_ -> MuIndependentTerms(x, m2Q2, nf);

    Value hs = highenergy_ -> MuIndependentTermsBand(x, m2Q2, nf);

    double gamma = params_A_.gamma, C = params_A_.C, delta = params_B_.gamma, D = params_B_.C;

    double res_A, res_B;

    if (GetOrder() == 2) {
        res_A = ApproximationA(x, m2Q2, 0., he_ll, hs.GetCentral(), thr, thr_const, gamma, C);
        res_B = ApproximationB(x, m2Q2, 0., he_ll, hs.GetCentral(), thr, 0., delta, D);
    } else if (GetOrder() == 3) {
        Value he_nll = ApproximateNLL(x, m2Q2);
        res_A = ApproximationA(x, m2Q2, he_ll, he_nll.GetHigher(), hs.GetHigher(), thr, thr_const, gamma, C);
        res_B = ApproximationB(x, m2Q2, he_ll, he_nll.GetLower(), hs.GetLower(), thr, thr_const, delta, D);
    }

    if (res_A > res_B) return Value(res_A, res_B);
    else return Value(res_B, res_A);
}

double ApproximateCoefficientFunctionKLMV::ApproximationA(double x, double m2Q2, double he_logx, double he_const, double hs, double thr, double thr_const, double gamma, double C) const {
    double beta = sqrt(1. - 4 * m2Q2 * x / (1. - x));
    double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

    double xi = 1. / m2Q2;
    double f = 1. / (1. + exp(2. * (xi - 4.)));

    double eta_gamma = pow(eta, gamma);
    
    double beta3 = beta * beta * beta;

    return thr - thr_const + (1. - f) * beta * hs
            + f * beta3 * (-log(eta) / log(x) * he_logx
            + he_const * eta_gamma / (C + eta_gamma));

}

double ApproximateCoefficientFunctionKLMV::ApproximationB(double x, double m2Q2, double he_logx, double he_const, double hs, double thr, double thr_const, double delta, double D) const {
    
    double beta = sqrt(1. - 4 * m2Q2 * x / (1. - x));
    double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

    double xi = 1. / m2Q2;
    double f = 1. / (1. + exp(2. * (xi - 4.)));

    double eta_gamma = pow(eta, delta);
    
    double beta3 = beta * beta * beta;

    return thr - thr_const + 2 * f * thr_const +
            (1. - f) * beta3 * hs
            + f * beta3 * (-log(eta) / log(x) * he_logx
            + he_const * eta_gamma / (D + eta_gamma));

}

Value ApproximateCoefficientFunctionKLMV::ApproximateNLL(double x, double m2Q2) const {
    double pi3 = M_PI * M_PI * M_PI;
    double tmp_A = (64. * pi3) * (params_A_.log_coeff * pow(log(1. / m2Q2) / log(5), params_A_.log_pow) - params_A_.const_coeff) * 4./ m2Q2 / x;
    double tmp_B = (64. * pi3) * (params_B_.log_coeff * pow(log(1. / m2Q2) / log(5), params_B_.log_pow) - params_B_.const_coeff) * 4./ m2Q2 / x;
    return Value(tmp_A, tmp_B);
}

// //==========================================================================================//
// //  Approximate gluon coefficient function for F2 at O(alpha_s) with implicit
// //  parameters.
// //------------------------------------------------------------------------------------------//

// double C2_g1_approximation_implicit(double x, double m2Q2, double k, double h) {

//     double eta;

//     eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

//     double damp_thr = 1. / (1. + pow(eta / h, k));
//     double damp_asy = 1. - damp_thr;

//     double C_const = C2_g1_asymptotic(x, m2Q2) * damp_asy
//                      + C2_g1_threshold(x, m2Q2) * damp_thr;

//     return C_const;
// }

// //==========================================================================================//
// //  Approximate gluon coefficient function for F2 at O(alpha_s).
// //  Used for validation of the method.
// //------------------------------------------------------------------------------------------//

// double C2_g1_approximation(double x, double m2Q2) {

//     double xmax = 1. / (1. + 4. * m2Q2);

//     if (x >= xmax || x <= 0)
//         return 0.;

//     double xi = 1. / m2Q2;

//     double eta;

//     eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

//     double h = 0.2 + 2.3 / (1. + exp(2.5 * (log(xi) - 5.)));
//     double k = 2.5 - 1.3 / (1. + exp(2.5 * (log(xi) - 5.)));

//     double damp_thr = 1. / (1 + pow(eta / h, k));
//     double damp_asy = 1. - damp_thr;

//     double C_const = C2_g1_asymptotic(x, m2Q2) * damp_asy
//                      + C2_g1_threshold(x, m2Q2) * damp_thr;

//     return C_const;
// }

// //==========================================================================================//
// //  Approximate gluon coefficient function for FL at O(alpha_s).
// //  Used for validation of the method.
// //------------------------------------------------------------------------------------------//

// double CL_g1_approximation(double x, double m2Q2) {

//     double xmax = 1. / (1. + 4. * m2Q2);

//     if (x >= xmax || x <= 0)
//         return 0.;

//     double xi = 1. / m2Q2;

//     double eta;

//     eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

//     double h = 20. + (11. - 20.) / (1. + exp(2.5 * (log(xi) - 5.))); //the 20 must be retuned
//     double k = 3. - (3. - 2.) / (1. + exp(2.5 * (log(xi) - 5.)));

//     double damp_thr = 1. / (1. + pow(eta / h, k));
//     double damp_asy = 1. - damp_thr;

//     double C_const = CL_g1_highscale(x) * damp_asy;
//     // CL_g1_threshold=0

//     return C_const;
// }

// //==========================================================================================//
// //  Approximate gluon coefficient function for F2 at O(alpha_s^2).
// //  Used for validation of the method.
// //------------------------------------------------------------------------------------------//

// double C2_g2_approximation(double x, double m2Q2, double m2mu2, int v) {

//     double xmax = 1. / (1. + 4. * m2Q2);

//     if (x >= xmax || x <= 0)
//         return 0.;

//     double C_const = C2_g20_approximation_BAND(x, m2Q2, v, var_C2, fact_C2);

//     double C_log = C2_g21(x, m2Q2);

//     return C_const + C_log * log(1. / m2mu2);
// }

// //==========================================================================================//
// //  Approximate mu-independent term of the gluon coefficient function for F2 at
// //  O(alpha_s^2).
// //------------------------------------------------------------------------------------------//

// double C2_g20_approximation(double x, double m2Q2) {

//     double A = C2_g2_params.A;
//     double B = C2_g2_params.B;
//     double C = C2_g2_params.C;
//     double D = C2_g2_params.D;
//     double a = C2_g2_params.a;
//     double b = C2_g2_params.b;

//     return C2_g20_approximation_implicit(x, m2Q2, A, B, C, D, a, b);
// }

// //==========================================================================================//
// //  Approximate mu-independent term of the gluon coefficient function for F2 at
// //  O(alpha_s) with implicit parameters.
// //------------------------------------------------------------------------------------------//

// double C2_g20_approximation_implicit(
//     double x, double m2Q2, double A, double B, double C, double D, double a,
//     double b
// ) {

//     double xi = 1. / m2Q2;

//     double eta;

//     eta = 0.25 / m2Q2 * (1 - x) / x - 1;

//     double h = A + (B - A) / (1. + exp(a * (log(xi) - b)));
//     double k = C + (D - C) / (1. + exp(a * (log(xi) - b)));

//     double damp_thr = 1. / (1. + pow(eta / h, k));
//     double damp_asy = 1. - damp_thr;

//     return C2_g2_asymptotic(x, m2Q2, 1.) * damp_asy
//            + C2_g2_threshold(x, m2Q2, 1.) * damp_thr;
// }

// //==========================================================================================//
// //  Band of the approximate mu-independent term of the gluon coefficient
// //  function for F2 at O(alpha_s^2).
// //------------------------------------------------------------------------------------------//

// double C2_g20_approximation_BAND(
//     double x, double m2Q2, int v, double var, double fact
// ) {

//     double A = C2_g2_params.A;
//     double B = C2_g2_params.B;
//     double C = C2_g2_params.C;
//     double D = C2_g2_params.D;
//     double a = C2_g2_params.a;
//     double b = C2_g2_params.b;

//     double Cavg = C2_g20_approximation(x, m2Q2);

//     if (v == 0)
//         return Cavg;

//     double Amax = fact * A, Amin = A / fact, Bmax = B * fact, Bmin = B / fact;
//     double Cmax = (1. + var) * C, Cmin = (1. - var) * C, Dmax = (1. + var) * D,
//            Dmin = (1. - var) * D;

//     double Avec[3] = { A, Amin, Amax };
//     double Bvec[3] = { B, Bmax, Bmin };
//     double Cvec[3] = { C, Cmax, Cmin };
//     double Dvec[3] = { D, Dmax, Dmin };

//     vector<double> replicas;

//     for (int i = 0; i < 3; i++) {
//         for (int j = 0; j < 3; j++) {
//             for (int h = 0; h < 3; h++) {
//                 for (int k = 0; k < 3; k++) {
//                     replicas.push_back(C2_g20_approximation_implicit(
//                         x, m2Q2, Avec[i], Bvec[j], Cvec[h], Dvec[k], a, b
//                     ));
//                 }
//             }
//         }
//     }

//     double min = Cavg, max = Cavg;

//     for (double rep : replicas) {
//         if (rep > max)
//             max = rep;
//         if (rep < min)
//             min = rep;
//     }

//     if (v == 1)
//         return max;
//     if (v == -1)
//         return min;
//     else {
//         cout << "C2_g20_approximation_BAND: Choose either v=0 or v=1 or "
//                 "v=2!!\nExiting!!\n"
//              << endl;
//         exit(-1);
//     }
// }

// //==========================================================================================//
// //  Approximate quark coefficient function for F2 at O(alpha_s^2).
// //  Used for validation of the method.
// //------------------------------------------------------------------------------------------//

// double C2_ps2_approximation(double x, double m2Q2, double m2mu2, int v) {

//     double xmax = 1. / (1. + 4. * m2Q2);

//     if (x >= xmax || x <= 0)
//         return 0.;

//     double C_const = C2_ps20_approximation_BAND(x, m2Q2, v, var_C2, fact_C2);

//     double C_log = C2_ps21(x, m2Q2);

//     return C_const + C_log * log(1. / m2mu2);
// }

// //==========================================================================================//
// //  Approximate mu-independent term of the quark coefficient function for F2 at
// //  O(alpha_s^2).
// //------------------------------------------------------------------------------------------//

// double C2_ps20_approximation(double x, double m2Q2) {

//     double A = C2_ps2_params.A;
//     double B = C2_ps2_params.B;
//     double C = C2_ps2_params.C;
//     double D = C2_ps2_params.D;
//     double a = C2_ps2_params.a;
//     double b = C2_ps2_params.b;

//     return C2_ps20_approximation_implicit(x, m2Q2, A, B, C, D, a, b);
// }

// //==========================================================================================//
// //  Approximate mu-independent term of the quark coefficient function for F2 at
// //  O(alpha_s^2) with implicit parameters.
// //------------------------------------------------------------------------------------------//

// double C2_ps20_approximation_implicit(
//     double x, double m2Q2, double A, double B, double C, double D, double a,
//     double b
// ) {

//     double xi = 1. / m2Q2;

//     double eta;

//     eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

//     double h = A + (B - A) / (1. + exp(a * (log(xi) - b)));
//     double k = C + (D - C) / (1. + exp(a * (log(xi) - b)));

//     double damp_thr = 1. / (1. + pow(eta / h, k));
//     double damp_asy = 1. - damp_thr;

//     return C2_ps2_asymptotic(x, m2Q2, 1.) * damp_asy;
//     // threshPS = 0
// }

// //==========================================================================================//
// //  Band of the approximate mu-independent term of the quark coefficient
// //  function for F2 at O(alpha_s^2).
// //------------------------------------------------------------------------------------------//

// double C2_ps20_approximation_BAND(
//     double x, double m2Q2, int v, double var, double fact
// ) {

//     double A = C2_ps2_params.A;
//     double B = C2_ps2_params.B;
//     double C = C2_ps2_params.C;
//     double D = C2_ps2_params.D;
//     double a = C2_ps2_params.a;
//     double b = C2_ps2_params.b;

//     double Cavg = C2_ps20_approximation(x, m2Q2);

//     if (v == 0)
//         return Cavg;

//     double Amax = fact * A, Amin = A / fact, Bmax = B * fact, Bmin = B / fact;
//     double Cmax = (1. + var) * C, Cmin = (1. - var) * C, Dmax = (1. + var) * D,
//            Dmin = (1. - var) * D;

//     double Avec[3] = { A, Amin, Amax };
//     double Bvec[3] = { B, Bmax, Bmin };
//     double Cvec[3] = { C, Cmax, Cmin };
//     double Dvec[3] = { D, Dmax, Dmin };

//     vector<double> replicas;

//     for (int i = 0; i < 3; i++) {
//         for (int j = 0; j < 3; j++) {
//             for (int h = 0; h < 3; h++) {
//                 for (int k = 0; k < 3; k++) {
//                     replicas.push_back(C2_ps20_approximation_implicit(
//                         x, m2Q2, Avec[i], Bvec[j], Cvec[h], Dvec[k], a, b
//                     ));
//                 }
//             }
//         }
//     }

//     double min = Cavg, max = Cavg;

//     for (double rep : replicas) {
//         if (rep > max)
//             max = rep;
//         if (rep < min)
//             min = rep;
//     }

//     if (v == 1)
//         return max;
//     if (v == -1)
//         return min;
//     else {
//         cout << "C2_g20_approximation_BAND: Choose either v=0 or v=1 or "
//                 "v=2!!\nExiting!!\n"
//              << endl;
//         exit(-1);
//     }
// }

// //==========================================================================================//
// //  Approximate gluon coefficient function for FL at O(alpha_s^2).
// //  Used for validation of the method.
// //------------------------------------------------------------------------------------------//

// double CL_g2_approximation(double x, double m2Q2, double m2mu2, int v) {

//     double xmax = 1. / (1. + 4 * m2Q2);

//     if (x >= xmax || x <= 0)
//         return 0.;

//     double C_const = CL_g20_approximation_BAND(x, m2Q2, v, var_CL, fact_CL);

//     double C_log = CL_g21(x, m2Q2);

//     return C_const + C_log * log(1. / m2mu2);
// }

// //==========================================================================================//
// //  Approximate mu-independent term of the gluon coefficient function for FL at
// //  O(alpha_s^2).
// //------------------------------------------------------------------------------------------//

// double CL_g20_approximation(double x, double m2Q2) {

//     double A = CL_g2_params.A;
//     double B = CL_g2_params.B;
//     double C = CL_g2_params.C;
//     double D = CL_g2_params.D;
//     double a = CL_g2_params.a;
//     double b = CL_g2_params.b;

//     return CL_g20_approximation_implicit(x, m2Q2, A, B, C, D, a, b);
// }

// //==========================================================================================//
// //  Approximate mu-independent term of the gluon coefficient function for FL at
// //  O(alpha_s) with implicit parameters.
// //------------------------------------------------------------------------------------------//

// double CL_g20_approximation_implicit(
//     double x, double m2Q2, double A, double B, double C, double D, double a,
//     double b
// ) {

//     double xi = 1. / m2Q2;

//     double eta;

//     eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

//     double h = A + (B - A) / (1. + exp(a * (log(xi) - b)));
//     double k = C + (D - C) / (1. + exp(a * (log(xi) - b)));

//     double damp_thr = 1. / (1. + pow(eta / h, k));
//     double damp_asy = 1. - damp_thr;

//     return CL_g2_asymptotic(x, m2Q2, 1.) * damp_asy
//            + CL_g2_threshold(x, m2Q2, 1.) * damp_thr;
// }

// //==========================================================================================//
// //  Band of the approximate mu-independent term of the gluon coefficient
// //  function for FL at O(alpha_s^2).
// //------------------------------------------------------------------------------------------//

// double CL_g20_approximation_BAND(
//     double x, double m2Q2, int v, double var, double fact
// ) {

//     double A = CL_g2_params.A;
//     double B = CL_g2_params.B;
//     double C = CL_g2_params.C;
//     double D = CL_g2_params.D;
//     double a = CL_g2_params.a;
//     double b = CL_g2_params.b;

//     double Cavg = CL_g20_approximation(x, m2Q2);

//     if (v == 0)
//         return Cavg;

//     double Amax = fact * A, Amin = A / fact, Bmax = B * fact, Bmin = B / fact;
//     double Cmax = (1. + var) * C, Cmin = (1. - var) * C, Dmax = (1. + var) * D,
//            Dmin = (1. - var) * D;

//     double Avec[3] = { A, Amin, Amax };
//     double Bvec[3] = { B, Bmax, Bmin };
//     double Cvec[3] = { C, Cmax, Cmin };
//     double Dvec[3] = { D, Dmax, Dmin };

//     vector<double> replicas;

//     for (int i = 0; i < 3; i++) {
//         for (int j = 0; j < 3; j++) {
//             for (int h = 0; h < 3; h++) {
//                 for (int k = 0; k < 3; k++) {
//                     replicas.push_back(CL_g20_approximation_implicit(
//                         x, m2Q2, Avec[i], Bvec[j], Cvec[h], Dvec[k], a, b
//                     ));
//                 }
//             }
//         }
//     }

//     double min = Cavg, max = Cavg;

//     for (double rep : replicas) {
//         if (rep > max)
//             max = rep;
//         if (rep < min)
//             min = rep;
//     }

//     if (v == 1)
//         return max;
//     if (v == -1)
//         return min;
//     else {
//         cout << "CL_g20_approximation_BAND: Choose either v=0 or v=1 or "
//                 "v=2!!\nExiting!!\n"
//              << endl;
//         exit(-1);
//     }
// }

// //==========================================================================================//
// //  Approximate quark coefficient function for FL at O(alpha_s^2).
// //  Used for validation of the method.
// //------------------------------------------------------------------------------------------//

// double CL_ps2_approximation(double x, double m2Q2, double m2mu2, int v) {

//     double xmax = 1. / (1. + 4. * m2Q2);

//     if (x >= xmax || x <= 0)
//         return 0.;

//     double C_const = CL_ps20_approximation_BAND(x, m2Q2, v, var_CL, fact_CL);

//     double C_log = CL_ps21(x, m2Q2);

//     return C_const + C_log * log(1. / m2mu2);
// }

// //==========================================================================================//
// //  Approximate mu-independent term of the quark coefficient function for FL at
// //  O(alpha_s^2).
// //------------------------------------------------------------------------------------------//

// double CL_ps20_approximation(double x, double m2Q2) {

//     double A = CL_ps2_params.A;
//     double B = CL_ps2_params.B;
//     double C = CL_ps2_params.C;
//     double D = CL_ps2_params.D;
//     double a = CL_ps2_params.a;
//     double b = CL_ps2_params.b;

//     return CL_ps20_approximation_implicit(x, m2Q2, A, B, C, D, a, b);
// }

// //==========================================================================================//
// //  Approximate mu-independent term of the quark coefficient function for FL at
// //  O(alpha_s^2) with implicit parameters.
// //------------------------------------------------------------------------------------------//

// double CL_ps20_approximation_implicit(
//     double x, double m2Q2, double A, double B, double C, double D, double a,
//     double b
// ) {

//     double xi = 1. / m2Q2;

//     double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

//     double h = A + (B - A) / (1. + exp(a * (log(xi) - b)));
//     double k = C + (D - C) / (1. + exp(a * (log(xi) - b)));

//     double damp_thr = 1. / (1. + pow(eta / h, k));
//     double damp_asy = 1. - damp_thr;

//     return CL_ps2_asymptotic(x, m2Q2, 1.) * damp_asy;
//     // threshPS = 0
// }

// //==========================================================================================//
// //  Band of the approximate mu-independent term of the quark coefficient
// //  function for FL at O(alpha_s^2).
// //------------------------------------------------------------------------------------------//

// double CL_ps20_approximation_BAND(
//     double x, double m2Q2, int v, double var, double fact
// ) {

//     double A = CL_ps2_params.A;
//     double B = CL_ps2_params.B;
//     double C = CL_ps2_params.C;
//     double D = CL_ps2_params.D;
//     double a = CL_ps2_params.a;
//     double b = CL_ps2_params.b;

//     double Cavg = CL_ps20_approximation(x, m2Q2);

//     if (v == 0)
//         return Cavg;

//     double Amax = fact * A, Amin = A / fact, Bmax = B * fact, Bmin = B / fact;
//     double Cmax = (1. + var) * C, Cmin = (1. - var) * C, Dmax = (1. + var) * D,
//            Dmin = (1. - var) * D;

//     double Avec[3] = { A, Amin, Amax };
//     double Bvec[3] = { B, Bmax, Bmin };
//     double Cvec[3] = { C, Cmax, Cmin };
//     double Dvec[3] = { D, Dmax, Dmin };

//     vector<double> replicas;

//     for (int i = 0; i < 3; i++) {
//         for (int j = 0; j < 3; j++) {
//             for (int h = 0; h < 3; h++) {
//                 for (int k = 0; k < 3; k++) {
//                     replicas.push_back(CL_ps20_approximation_implicit(
//                         x, m2Q2, Avec[i], Bvec[j], Cvec[h], Dvec[k], a, b
//                     ));
//                 }
//             }
//         }
//     }

//     double min = Cavg, max = Cavg;

//     for (double rep : replicas) {
//         if (rep > max)
//             max = rep;
//         if (rep < min)
//             min = rep;
//     }

//     if (v == 1)
//         return max;
//     if (v == -1)
//         return min;
//     else {
//         cout << "C2_g20_approximation_BAND: Choose either v=0 or v=1 or "
//                 "v=2!!\nExiting!!\n"
//              << endl;
//         exit(-1);
//     }
// }
// //==========================================================================================//
// //  Approximate gluon coefficient function for F2 at O(alpha_s^3).
// //------------------------------------------------------------------------------------------//

// double C2_g3_approximation(
//     double x, double m2Q2, double m2mu2, int nf, int v, int method_flag
// ) {

//     double xmax = 1. / (1. + 4. * m2Q2);

//     if (x >= xmax || x <= 0)
//         return 0.;

//     double C30 = C2_g30_approximation_BAND(x, m2Q2, nf, v, var_C2, fact_C2);

//     if (m2mu2 == 1.)
//         return C30;

//     double Lmu = -log(m2mu2);
//     double Lmu2 = Lmu * Lmu;

//     return C30 + C2_g31(x, m2Q2, nf) * Lmu
//            + C2_g32(x, m2Q2, nf, method_flag) * Lmu2;
// }

// //==========================================================================================//
// //  Approximate mu-independent term of the gluon coefficient function for F2 at
// //  O(alpha_s^3).
// //------------------------------------------------------------------------------------------//

// double C2_g30_approximation(double x, double m2Q2, int nf) {

//     double A = C2_g3_params.A;
//     double B = C2_g3_params.B;
//     double C = C2_g3_params.C;
//     double D = C2_g3_params.D;
//     double a = C2_g3_params.a;
//     double b = C2_g3_params.b;
//     int v1 = 0, v2 = 0;

//     return C2_g30_approximation_implicit(x, m2Q2, nf, A, B, C, D, a, b, v1, v2);
// }

// //==========================================================================================//
// //  Approximate mu-independent term of the gluon coefficient function for F2 at
// //  O(alpha_s^3) with implicit parameters.
// //------------------------------------------------------------------------------------------//

// double C2_g30_approximation_implicit(
//     double x, double m2Q2, int nf, double A, double B, double C, double D,
//     double a, double b, int v1, int v2
// ) {

//     double xi = 1. / m2Q2;

//     double eta = 0.25 / m2Q2 * (1 - x) / x - 1.;

//     double h = A + (B - A) / (1. + exp(a * (log(xi) - b)));
//     double k = C + (D - C) / (1. + exp(a * (log(xi) - b)));

//     double damp_thr = 1. / (1. + pow(eta / h, k));
//     double damp_asy = 1. - damp_thr;

//     return C2_g3_asymptotic(x, m2Q2, 1., nf, v1, v2) * damp_asy
//            + C2_g3_threshold(x, m2Q2, 1., nf) * damp_thr;
// }

// //==========================================================================================//
// //  Band of the approximate mu-independent term of the gluon coefficient
// //  function for F2 at O(alpha_s^3).
// //------------------------------------------------------------------------------------------//

// double C2_g30_approximation_BAND(
//     double x, double m2Q2, int nf, int v, double var, double fact
// ) {

//     double A = C2_g3_params.A;
//     double B = C2_g3_params.B;
//     double C = C2_g3_params.C;
//     double D = C2_g3_params.D;
//     double a = C2_g3_params.a;
//     double b = C2_g3_params.b;

//     double Cavg = C2_g30_approximation(x, m2Q2, nf);

//     if (v == 0)
//         return Cavg;

//     double Amax = fact * A, Amin = A / fact, Bmax = B * fact, Bmin = B / fact;
//     double Cmax = (1. + var) * C, Cmin = (1. - var) * C, Dmax = (1. + var) * D,
//            Dmin = (1. - var) * D;

//     double Avec[3] = { A, Amin, Amax };
//     double Bvec[3] = { B, Bmax, Bmin };
//     double Cvec[3] = { C, Cmax, Cmin };
//     double Dvec[3] = { D, Dmax, Dmin };

//     vector<double> replicas;

//     for (int i = 0; i < 3; i++) {
//         for (int j = 0; j < 3; j++) {
//             for (int h = 0; h < 3; h++) {
//                 for (int k = 0; k < 3; k++) {
//                     for (int v1 = -1; v1 <= 1; v1++) {
//                         for (int v2 = -1; v2 <= 1; v2++) {
//                             replicas.push_back(C2_g30_approximation_implicit(
//                                 x, m2Q2, nf, Avec[i], Bvec[j], Cvec[h], Dvec[k],
//                                 a, b, v1, v2
//                             ));
//                         }
//                     }
//                 }
//             }
//         }
//     }

//     double min = Cavg, max = Cavg;

//     for (double rep : replicas) {
//         if (rep > max)
//             max = rep;
//         if (rep < min)
//             min = rep;
//     }

//     if (v == 1)
//         return max;
//     if (v == -1)
//         return min;

//     else {
//         cout << "C2_g30_approximation_BAND: Choose either v=0 or v=1 or "
//                 "v=-1!!\nExiting!!\n"
//              << endl;
//         exit(-1);
//     }
// }

// //==========================================================================================//
// //  Approximate quark coefficient function for F2 at O(alpha_s^3).
// //------------------------------------------------------------------------------------------//

// double
// C2_ps3_approximation(double x, double m2Q2, double m2mu2, int nf, int v) {

//     double xmax = 1. / (1. + 4. * m2Q2);

//     if (x >= xmax || x <= 0)
//         return 0.;

//     double C30 = C2_ps30_approximation_BAND(x, m2Q2, nf, v, var_C2, fact_C2);

//     if (m2mu2 == 1.)
//         return C30;

//     double Lmu = -log(m2mu2);
//     double Lmu2 = Lmu * Lmu;

//     return C30 + C2_ps31(x, m2Q2, nf) * Lmu + C2_ps32(x, m2Q2, nf) * Lmu2;
// }

// //==========================================================================================//
// //  Approximate mu-independent term of the quark coefficient function for F2 at
// //  O(alpha_s^3).
// //------------------------------------------------------------------------------------------//

// double C2_ps30_approximation(double x, double m2Q2, int nf) {

//     double A = C2_ps3_params.A;
//     double B = C2_ps3_params.B;
//     double C = C2_ps3_params.C;
//     double D = C2_ps3_params.D;
//     double a = C2_ps3_params.a;
//     double b = C2_ps3_params.b;

//     return C2_ps30_approximation_implicit(x, m2Q2, nf, A, B, C, D, a, b, 0);
// }

// //==========================================================================================//
// //  Approximate mu-independent term of the quark coefficient function for F2 at
// //  O(alpha_s^3) with implicit parameters.
// //------------------------------------------------------------------------------------------//

// double C2_ps30_approximation_implicit(
//     double x, double m2Q2, int nf, double A, double B, double C, double D,
//     double a, double b, int v
// ) {

//     double xi = 1. / m2Q2;

//     double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

//     double h = A + (B - A) / (1. + exp(a * (log(xi) - b)));
//     double k = C + (D - C) / (1. + exp(a * (log(xi) - b)));

//     double damp_thr = 1. / (1. + pow(eta / h, k));
//     double damp_asy = 1. - damp_thr;

//     return C2_ps3_asymptotic(x, m2Q2, 1., nf, v) * damp_asy;
// }

// //==========================================================================================//
// //  Band of the approximate mu-independent term of the quark coefficient
// //  function for F2 at O(alpha_s^3).
// //------------------------------------------------------------------------------------------//

// double C2_ps30_approximation_BAND(
//     double x, double m2Q2, int nf, int v, double var, double fact
// ) {

//     double A = C2_ps3_params.A;
//     double B = C2_ps3_params.B;
//     double C = C2_ps3_params.C;
//     double D = C2_ps3_params.D;
//     double a = C2_ps3_params.a;
//     double b = C2_ps3_params.b;

//     double Cavg = C2_ps30_approximation(x, m2Q2, nf);

//     if (v == 0)
//         return Cavg;

//     double Amax = fact * A, Amin = A / fact, Bmax = B * fact, Bmin = B / fact;
//     double Cmax = (1. + var) * C, Cmin = (1. - var) * C, Dmax = (1. + var) * D,
//            Dmin = (1. - var) * D;

//     double Avec[3] = { A, Amin, Amax };
//     double Bvec[3] = { B, Bmax, Bmin };
//     double Cvec[3] = { C, Cmax, Cmin };
//     double Dvec[3] = { D, Dmax, Dmin };

//     vector<double> replicas;

//     for (int i = 0; i < 3; i++) {
//         for (int j = 0; j < 3; j++) {
//             for (int h = 0; h < 3; h++) {
//                 for (int k = 0; k < 3; k++) {
//                     for (int v1 = -1; v1 <= 1; v1++) {
//                         replicas.push_back(C2_ps30_approximation_implicit(
//                             x, m2Q2, nf, Avec[i], Bvec[j], Cvec[h], Dvec[k], a,
//                             b, v1
//                         ));
//                     }
//                 }
//             }
//         }
//     }

//     double min = Cavg, max = Cavg;

//     for (double rep : replicas) {
//         if (rep > max)
//             max = rep;
//         if (rep < min)
//             min = rep;
//     }

//     if (v == 1)
//         return max;
//     if (v == -1)
//         return min;

//     else {
//         cout << "C2_ps30_approximation_BAND: Choose either v=0 or v=1 or "
//                 "v=-1!!\nExiting!!\n"
//              << endl;
//         exit(-1);
//     }
// }

// //==========================================================================================//
// //  Approximate gluon coefficient function for FL at O(alpha_s^3).
// //------------------------------------------------------------------------------------------//

// double CL_g3_approximation(
//     double x, double m2Q2, double m2mu2, int nf, int v, int method_flag
// ) {

//     double xmax = 1. / (1. + 4. * m2Q2);

//     if (x >= xmax || x <= 0)
//         return 0.;

//     double C30 = CL_g30_approximation_BAND(x, m2Q2, nf, v, var_CL, fact_CL);

//     if (m2mu2 == 1.)
//         return C30;

//     double Lmu = -log(m2mu2);
//     double Lmu2 = Lmu * Lmu;

//     return C30 + CL_g31(x, m2Q2, nf) * Lmu
//            + CL_g32(x, m2Q2, nf, method_flag) * Lmu2;
// }

// //==========================================================================================//
// //  Approximate mu-independent term of the gluon coefficient function for F2 at
// //  O(alpha_s^3).
// //------------------------------------------------------------------------------------------//

// double CL_g30_approximation(double x, double m2Q2, int nf) {

//     double A = CL_g3_params.A;
//     double B = CL_g3_params.B;
//     double C = CL_g3_params.C;
//     double D = CL_g3_params.D;
//     double a = CL_g3_params.a;
//     double b = CL_g3_params.b;

//     return CL_g30_approximation_implicit(x, m2Q2, nf, A, B, C, D, a, b, 0);
// }

// //==========================================================================================//
// //  Approximate mu-independent term of the gluon coefficient function for FL at
// //  O(alpha_s^3) with implicit parameters.
// //------------------------------------------------------------------------------------------//

// double CL_g30_approximation_implicit(
//     double x, double m2Q2, int nf, double A, double B, double C, double D,
//     double a, double b, int v
// ) {

//     double xi = 1. / m2Q2;

//     double eta = 0.25 / m2Q2 * (1 - x) / x - 1.;

//     double h = A + (B - A) / (1. + exp(a * (log(xi) - b)));
//     double k = C + (D - C) / (1. + exp(a * (log(xi) - b)));

//     double damp_thr = 1. / (1. + pow(eta / h, k));
//     double damp_asy = 1. - damp_thr;

//     return CL_g3_asymptotic(x, m2Q2, 1., nf, v) * damp_asy
//            + CL_g3_threshold(x, m2Q2, 1., nf) * damp_thr;
// }

// //==========================================================================================//
// //  Band of the approximate mu-independent term of the gluon coefficient
// //  function for FL at O(alpha_s^3).
// //------------------------------------------------------------------------------------------//

// double CL_g30_approximation_BAND(
//     double x, double m2Q2, int nf, int v, double var, double fact
// ) {

//     double A = CL_g3_params.A;
//     double B = CL_g3_params.B;
//     double C = CL_g3_params.C;
//     double D = CL_g3_params.D;
//     double a = CL_g3_params.a;
//     double b = CL_g3_params.b;

//     double Cavg = CL_g30_approximation(x, m2Q2, nf);

//     if (v == 0)
//         return Cavg;

//     double Amax = fact * A, Amin = A / fact, Bmax = B * fact, Bmin = B / fact;
//     double Cmax = (1. + var) * C, Cmin = (1. - var) * C, Dmax = (1. + var) * D,
//            Dmin = (1. - var) * D;

//     double Avec[3] = { A, Amin, Amax };
//     double Bvec[3] = { B, Bmax, Bmin };
//     double Cvec[3] = { C, Cmax, Cmin };
//     double Dvec[3] = { D, Dmax, Dmin };

//     vector<double> replicas;

//     for (int i = 0; i < 3; i++) {
//         for (int j = 0; j < 3; j++) {
//             for (int h = 0; h < 3; h++) {
//                 for (int k = 0; k < 3; k++) {
//                     for (int v1 = -1; v1 <= 1; v1++) {
//                         replicas.push_back(CL_g30_approximation_implicit(
//                             x, m2Q2, nf, Avec[i], Bvec[j], Cvec[h], Dvec[k], a,
//                             b, v1
//                         ));
//                     }
//                 }
//             }
//         }
//     }

//     double min = Cavg, max = Cavg;

//     for (double rep : replicas) {
//         if (rep > max)
//             max = rep;
//         if (rep < min)
//             min = rep;
//     }

//     if (v == 1)
//         return max;
//     if (v == -1)
//         return min;

//     else {
//         cout << "CL_g30_approximation_BAND: Choose either v=0 or v=1 or "
//                 "v=-1!!\nExiting!!\n"
//              << endl;
//         exit(-1);
//     }
// }

// //==========================================================================================//
// //  Approximate quark coefficient funtcions for FL at O(alpha_s^3).
// //------------------------------------------------------------------------------------------//

// double
// CL_ps3_approximation(double x, double m2Q2, double m2mu2, int nf, int v) {

//     double xmax = 1. / (1. + 4. * m2Q2);

//     if (x >= xmax || x <= 0)
//         return 0.;

//     double C30 = CL_ps30_approximation_BAND(x, m2Q2, nf, v, var_CL, fact_CL);

//     if (m2mu2 == 1.)
//         return C30;

//     double Lmu = -log(m2mu2);
//     double Lmu2 = Lmu * Lmu;

//     return C30 + CL_ps31(x, m2Q2, nf) * Lmu + CL_ps32(x, m2Q2, nf) * Lmu2;
// }

// //==========================================================================================//
// //  Approximate mu-independent term of the quark coefficient function for F2 at
// //  O(alpha_s^3).
// //------------------------------------------------------------------------------------------//

// double CL_ps30_approximation(double x, double m2Q2, int nf) {

//     double A = CL_ps3_params.A;
//     double B = CL_ps3_params.B;
//     double C = CL_ps3_params.C;
//     double D = CL_ps3_params.D;
//     double a = CL_ps3_params.a;
//     double b = CL_ps3_params.b;

//     return CL_ps30_approximation_implicit(x, m2Q2, nf, A, B, C, D, a, b, 0);
// }

// //==========================================================================================//
// //  Approximate mu-independent term of the quark coefficient function for FL at
// //  O(alpha_s^3) with implicit parameters.
// //------------------------------------------------------------------------------------------//

// double CL_ps30_approximation_implicit(
//     double x, double m2Q2, int nf, double A, double B, double C, double D,
//     double a, double b, int v
// ) {

//     double xi = 1. / m2Q2;

//     double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

//     double h = A + (B - A) / (1. + exp(a * (log(xi) - b)));
//     double k = C + (D - C) / (1. + exp(a * (log(xi) - b)));

//     double damp_thr = 1. / (1. + pow(eta / h, k));
//     double damp_asy = 1. - damp_thr;

//     return CL_ps3_asymptotic(x, m2Q2, 1., nf, v) * damp_asy;
// }

// //==========================================================================================//
// //  Band of the approximate mu-independent term of the quark coefficient
// //  function for FL at O(alpha_s^3).
// //------------------------------------------------------------------------------------------//

// double CL_ps30_approximation_BAND(
//     double x, double m2Q2, int nf, int v, double var, double fact
// ) {

//     double A = CL_ps3_params.A;
//     double B = CL_ps3_params.B;
//     double C = CL_ps3_params.C;
//     double D = CL_ps3_params.D;
//     double a = CL_ps3_params.a;
//     double b = CL_ps3_params.b;

//     double Cavg = CL_ps30_approximation(x, m2Q2, nf);

//     if (v == 0)
//         return Cavg;

//     double Amax = fact * A, Amin = A / fact, Bmax = B * fact, Bmin = B / fact;
//     double Cmax = (1. + var) * C, Cmin = (1. - var) * C, Dmax = (1. + var) * D,
//            Dmin = (1. - var) * D;

//     double Avec[3] = { A, Amin, Amax };
//     double Bvec[3] = { B, Bmax, Bmin };
//     double Cvec[3] = { C, Cmax, Cmin };
//     double Dvec[3] = { D, Dmax, Dmin };

//     vector<double> replicas;

//     for (int i = 0; i < 3; i++) {
//         for (int j = 0; j < 3; j++) {
//             for (int h = 0; h < 3; h++) {
//                 for (int k = 0; k < 3; k++) {
//                     for (int v1 = -1; v1 <= 1; v1++) {
//                         replicas.push_back(CL_ps30_approximation_implicit(
//                             x, m2Q2, nf, Avec[i], Bvec[j], Cvec[h], Dvec[k], a,
//                             b, v1
//                         ));
//                     }
//                 }
//             }
//         }
//     }

//     double min = Cavg, max = Cavg;

//     for (double rep : replicas) {
//         if (rep > max)
//             max = rep;
//         if (rep < min)
//             min = rep;
//     }

//     if (v == 1)
//         return max;
//     if (v == -1)
//         return min;

//     else {
//         cout << "CL_ps30_approximation_BAND: Choose either v=0 or v=1 or "
//                 "v=-1!!\nExiting!!\n"
//              << endl;
//         exit(-1);
//     }
// }

// //==========================================================================================//
// //  Approximate gluon coefficient funtcions for F2 at O(alpha_s^2) from
// //  [arXiv:1205.5727]. klmv = Kawamura, Lo Presti, Moch, Vogt. Upper band (i.e.
// //  approximation A)
// //
// //  Eq. (4.10) of Ref. [arXiv:1205.5727].
// //------------------------------------------------------------------------------------------//

// double C2_g2_approximationA_klmv(double x, double m2Q2, double m2mu2) {

//     double x_max = 1. / (1. + 4 * m2Q2);

//     if (x >= x_max || x < 0)
//         return 0;

//     double beta = sqrt(1. - 4 * m2Q2 * x / (1. - x));

//     double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

//     double xi = 1. / m2Q2;

//     double f = 1. / (1. + exp(2. * (xi - 4.)));

//     double gamma = 1.0, C = 42.5;

//     double eta_gamma = pow(eta, gamma);

//     double beta3 = beta * beta * beta;

//     double C_const = C2_g2_threshold(x, m2Q2, 1)
//                      - C2_g2_threshold_const(x, m2Q2, 1.)
//                      + (1. - f) * beta * C2_g2_highscale(x, m2Q2, 1.)
//                      + f * beta3 * C2_g2_highenergy(x, m2Q2, 1.) * eta_gamma
//                            / (C + eta_gamma);

//     double C_log = C2_g21(x, m2Q2);

//     return C_const + C_log * log(1. / m2mu2);
// }

// //==========================================================================================//
// //  Approximate gluon coefficient funtcions for F2 at O(alpha_s^2) from
// //  [arXiv:1205.5727]. klmv = Kawamura, Lo Presti, Moch, Vogt Lower band (i.e.
// //  approximation B)
// //
// //  Eq. (4.11) of Ref. [arXiv:1205.5727].
// //------------------------------------------------------------------------------------------//

// double C2_g2_approximationB_klmv(double x, double m2Q2, double m2mu2) {

//     double x_max = 1. / (1. + 4. * m2Q2);

//     if (x >= x_max || x < 0)
//         return 0;

//     double beta = sqrt(1. - 4. * m2Q2 * x / (1. - x));

//     double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

//     double xi = 1. / m2Q2;

//     double f = 1. / (1. + exp(2. * (xi - 4.)));

//     double delta = 0.8, D = 19.4;

//     double eta_delta = pow(eta, delta);

//     double beta3 = beta * beta * beta;

//     double C_const = C2_g2_threshold(x, m2Q2, 1.)
//                      + (1. - f) * beta3 * C2_g2_highscale(x, m2Q2, 1.)
//                      + f * beta3 * C2_g2_highenergy(x, m2Q2, 1.) * eta_delta
//                            / (D + eta_delta);

//     double C_log = C2_g21(x, m2Q2);

//     return C_const + C_log * log(1. / m2mu2);
// }

// //==========================================================================================//
// //  Approximate quark coefficient funtcions for F2 at O(alpha_s^2) from
// //  [arXiv:1205.5727]. klmv = Kawamura, Lo Presti, Moch, Vogt. Upper band (i.e.
// //  approximation A)
// //
// //  Eq. (4.10) of Ref. [arXiv:1205.5727] but for the quark coefficient function.
// //------------------------------------------------------------------------------------------//

// double C2_ps2_approximationA_klmv(double x, double m2Q2, double m2mu2) {

//     double x_max = 1. / (1. + 4 * m2Q2);

//     if (x >= x_max || x < 0)
//         return 0;

//     double beta = sqrt(1. - 4. * m2Q2 * x / (1. - x));

//     double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

//     double xi = 1. / m2Q2;

//     double f = 1. / (1. + exp(2. * (xi - 4.)));

//     double gamma = 1.0, C = 42.5;

//     double eta_gamma = pow(eta, gamma);

//     double beta3 = beta * beta * beta;

//     double C_const = (1. - f) * beta * C2_ps2_highscale(x, m2Q2, 1.)
//                      + f * beta3 * C2_ps2_highenergy(x, m2Q2, 1.) * eta_gamma
//                            / (C + eta_gamma);

//     double C_log = C2_g21(x, m2Q2);

//     return C_const + C_log * log(1. / m2mu2);
// }

// //==========================================================================================//
// //  Approximate quark coefficient funtcions for F2 at O(alpha_s^2) from
// //  [arXiv:1205.5727]. klmv = Kawamura, Lo Presti, Moch, Vogt. Lower band (i.e.
// //  approximation B)
// //
// //  Eq. (4.11) of Ref. [arXiv:1205.5727] but for the quark coefficient function.
// //------------------------------------------------------------------------------------------//

// double C2_ps2_approximationB_klmv(double x, double m2Q2, double m2mu2) {

//     double x_max = 1. / (1. + 4 * m2Q2);

//     if (x >= x_max || x < 0)
//         return 0;

//     double beta = sqrt(1. - 4 * m2Q2 * x / (1. - x));

//     double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

//     double xi = 1. / m2Q2;

//     double f = 1. / (1. + exp(2. * (xi - 4.)));

//     double delta = 0.8, D = 19.4;

//     double eta_delta = pow(eta, delta);

//     double beta3 = beta * beta * beta;

//     double C_const = (1. - f) * beta3 * C2_ps2_highscale(x, m2Q2, 1.)
//                      + f * beta3 * C2_ps2_highenergy(x, m2Q2, 1.) * eta_delta
//                            / (D + eta_delta);

//     double C_log = C2_g21(x, m2Q2);

//     return C_const + C_log * log(1. / m2mu2);
// }

// //==========================================================================================//
// //  Approximate gluon coefficient funtcions for F2 at O(alpha_s^3) from
// //  [arXiv:1205.5727]. klmv = Kawamura, Lo Presti, Moch, Vogt. Upper band (i.e.
// //  approximation A)
// //
// //  Eq. (4.17) of Ref. [arXiv:1205.5727].
// //  This equation uses the approximate form of aQg30 given in Eq. (3.49) of
// //  [arXiv:1205.5727].
// //------------------------------------------------------------------------------------------//

// double C2_g3_approximationA_klmv(
//     double x, double m2Q2, double m2mu2, int nf, int method_flag
// ) {

//     double x_max = 1. / (1. + 4 * m2Q2);

//     if (x >= x_max || x < 0)
//         return 0;

//     double pi3 = M_PI * M_PI * M_PI;

//     double beta = sqrt(1. - 4 * m2Q2 * x / (1. - x));

//     double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

//     double xi = 1. / m2Q2;

//     double f = 1. / (1. + exp(2. * (xi - 4)));

//     double gamma = 1.0, C = 20.0;

//     double eta_gamma = pow(eta, gamma);

//     double beta3 = beta * beta * beta;

//     double C2_g3_highenergyA =
//         (64. * pi3) * (0.007 * pow(log(1. / m2Q2) / log(5), 4) - 0.28) * 4.
//         / m2Q2 / x;

//     double C30 = C2_g3_threshold(x, m2Q2, 1., nf)
//                  - C2_g3_threshold_const(x, m2Q2, 1.)
//                  + (1. - f) * beta * C2_g3_highscale(x, m2Q2, 1., nf, 1)
//                  + f * beta3
//                        * (-log(eta) / log(x) * C2_g3_highenergyLL(x, m2Q2, 1.)
//                           + C2_g3_highenergyA * eta_gamma / (C + eta_gamma));

//     if (m2mu2 == 1.)
//         return C30;

//     double Lmu = -log(m2mu2);
//     double Lmu2 = Lmu * Lmu;

//     return C30 + C2_g31(x, m2Q2, nf) * Lmu
//            + C2_g32(x, m2Q2, nf, method_flag) * Lmu2;
// }

// //==========================================================================================//
// //  Approximate gluon coefficient funtcions for F2 at O(alpha_s^3) from
// //  [arXiv:1205.5727]. klmv = Kawamura, Lo Presti, Moch, Vogt. Lower band (i.e.
// //  approximation B)
// //
// //  Eq. (4.18) of Ref. [arXiv:1205.5727].
// //  This equation uses the approximate form of aQg30 given in Eq. (16) of Ref.
// //  [arXiv:1701.05838].
// //------------------------------------------------------------------------------------------//

// double C2_g3_approximationB_klmv(
//     double x, double m2Q2, double m2mu2, int nf, int method_flag
// ) {

//     double x_max = 1. / (1. + 4 * m2Q2);

//     if (x >= x_max || x < 0)
//         return 0;

//     double pi3 = M_PI * M_PI * M_PI;

//     double beta = sqrt(1. - 4 * m2Q2 * x / (1. - x));

//     double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

//     double xi = 1. / m2Q2;

//     double f = 1. / (1. + exp(2. * (xi - 4)));

//     double delta = 0.8, D = 10.7;

//     double eta_delta = pow(eta, delta);

//     double beta3 = beta * beta * beta;

//     double C2_g3_highenergyB =
//         (64. * pi3) * (0.055 * pow(log(1. / m2Q2) / log(5), 2) - 0.423) * 4.
//         / m2Q2 / x;

//     double C30 =
//         (C2_g3_threshold(x, m2Q2, 1., nf) - C2_g3_threshold_const(x, m2Q2, 1.))
//         + f * 2. * C2_g3_threshold_const(x, m2Q2, 1.)
//         + (1. - f) * beta3 * C2_g3_highscale(x, m2Q2, 1., nf, -1)
//         + f * beta3
//               * (-log(eta) / log(x) * C2_g3_highenergyLL(x, m2Q2, 1.)
//                  + C2_g3_highenergyB * eta_delta / (D + eta_delta));

//     if (m2mu2 == 1.)
//         return C30;

//     double Lmu = -log(m2mu2);
//     double Lmu2 = Lmu * Lmu;

//     return C30 + C2_g31(x, m2Q2, nf) * Lmu
//            + C2_g32(x, m2Q2, nf, method_flag) * Lmu2;
// }

// //==========================================================================================//
// //  Approximate gluon coefficient funtcions for F2 at O(alpha_s^3) from
// //  [arXiv:1205.5727]. klmv = Kawamura, Lo Presti, Moch, Vogt. Lower band (i.e.
// //  approximation B).
// //
// //  Eq. (4.18) of Ref. [arXiv:1205.5727].
// //  This equation uses the approximate form of aQg30 given in Eq. (3.50) of
// //  [arXiv:1205.5727] instead of the one given in Eq. (16) of Ref.
// //  [arXiv:1701.05838].
// //------------------------------------------------------------------------------------------//

// double C2_g3_approximationB_klmv_paper(
//     double x, double m2Q2, double m2mu2, int nf, int method_flag
// ) {

//     double x_max = 1. / (1. + 4 * m2Q2);

//     if (x >= x_max || x < 0)
//         return 0;

//     double pi3 = M_PI * M_PI * M_PI;

//     double beta = sqrt(1. - 4 * m2Q2 * x / (1. - x));

//     double eta = 0.25 / m2Q2 * (1. - x) / x - 1;

//     double xi = 1. / m2Q2;

//     double f = 1. / (1. + exp(2. * (xi - 4.)));

//     double delta = 0.8, D = 10.7;

//     double eta_delta = pow(eta, delta);

//     double beta3 = beta * beta * beta;

//     double C2_g3_highenergyB =
//         (64. * pi3) * (0.055 * pow(log(1. / m2Q2) / log(5), 2) - 0.423) * 4.
//         / m2Q2 / x;

//     double C30 =
//         (C2_g3_threshold(x, m2Q2, 1., nf) - C2_g3_threshold_const(x, m2Q2, 1.))
//         + f * 2. * C2_g3_threshold_const(x, m2Q2, 1.)
//         + (1. - f) * beta3 * C2_g3_highscale(x, m2Q2, 1., nf, -12)
//         + f * beta3
//               * (-log(eta) / log(x) * C2_g3_highenergyLL(x, m2Q2, 1.)
//                  + C2_g3_highenergyB * eta_delta / (D + eta_delta));

//     if (m2mu2 == 1.)
//         return C30;

//     double Lmu = -log(m2mu2);
//     double Lmu2 = Lmu * Lmu;

//     return C30 + C2_g31(x, m2Q2, nf) * Lmu
//            + C2_g32(x, m2Q2, nf, method_flag) * Lmu2;
// }

// //==========================================================================================//
// //  Approximate gluon coefficient funtcions for F2 at O(alpha_s^3) from
// //  [arXiv:1205.5727]. klmv = Kawamura, Lo Presti, Moch, Vogt. Lower band (i.e.
// //  approximation B) with the low xi limit.
// //
// //  Eq. (4.21) of Ref. [arXiv:1205.5727].
// //  This equation uses the exact expression of aQqPS30 given in Eq.
// //  (5.41, 5.42, 5.45) of Ref. [arXiv:1409.1135], and as approximation B for the
// //  small x limit at NLL it uses Eq. (4.25) of [arXiv:1205.5727] instead of Eq.
// //  (4.20)
// //------------------------------------------------------------------------------------------//

// double C2_g3_approximationBlowxi_klmv(
//     double x, double m2Q2, double m2mu2, int nf, int method_flag
// ) {

//     double x_max = 1. / (1. + 4 * m2Q2);

//     if (x >= x_max || x < 0)
//         return 0;

//     double pi3 = M_PI * M_PI * M_PI;

//     double beta = sqrt(1. - 4 * m2Q2 * x / (1. - x));

//     double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

//     double xi = 1. / m2Q2;

//     double f = 1. / (1. + exp(2. * (xi - 4.)));

//     double delta = 0.8, D = 10.7;

//     double eta_delta = pow(eta, delta);

//     double beta3 = beta * beta * beta;

//     double C2_g3_highenergyB =
//         CA / CF * (64. * pi3)
//         * (0.0245 * pow(log(1. / m2Q2) / log(5), 2) - 0.17) * 4. / m2Q2 / x;

//     double C30 =
//         (C2_g3_threshold(x, m2Q2, 1., nf) - C2_g3_threshold_const(x, m2Q2, 1.))
//         + f * 2. * C2_g3_threshold_const(x, m2Q2, 1.)
//         + (1. - f) * beta3 * C2_g3_highscale(x, m2Q2, 1., nf, -12)
//         + f * beta3
//               * (-log(eta) / log(x) * C2_g3_highenergyLL(x, m2Q2, 1.)
//                  + C2_g3_highenergyB * eta_delta / (D + eta_delta));

//     if (m2mu2 == 1.)
//         return C30;

//     double Lmu = -log(m2mu2);
//     double Lmu2 = Lmu * Lmu;

//     return C30 + C2_g31(x, m2Q2, nf) * Lmu
//            + C2_g32(x, m2Q2, nf, method_flag) * Lmu2;
// }

// //==========================================================================================//
// //  Approximate quark coefficient funtcions for F2 at O(alpha_s^3) from
// //  [arXiv:1205.5727]. klmv = Kawamura, Lo Presti, Moch, Vogt. Upper band (i.e.
// //  approximation A).
// //
// //  Eq. (4.21) of Ref. [arXiv:1205.5727].
// //  This equation uses the exact expression of aQqPS30 given in Eq.
// //  (5.41, 5.42, 5.45) of Ref. [arXiv:1409.1135]
// //------------------------------------------------------------------------------------------//

// double C2_ps3_approximationA_klmv(double x, double m2Q2, double m2mu2, int nf) {

//     double x_max = 1. / (1. + 4 * m2Q2);

//     if (x >= x_max || x < 0)
//         return 0;

//     double pi3 = M_PI * M_PI * M_PI;

//     double beta = sqrt(1. - 4. * m2Q2 * x / (1. - x));

//     double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

//     double xi = 1. / m2Q2;

//     double f = 1. / (1. + exp(2. * (xi - 4)));

//     double gamma = 1.0, C = 20.0;

//     double eta_gamma = pow(eta, gamma);

//     double beta3 = beta * beta * beta;

//     double C2_ps30_highenergyNLLA =
//         (64. * pi3) * (0.004 * pow(log(1. / m2Q2) / log(5), 4) - 0.125) * 4.
//         / m2Q2 / x;

//     double C30 =
//         (1. - f) * beta * C2_ps3_highscale(x, m2Q2, 1., nf)
//         + f * beta3
//               * (-log(eta) / log(x) * C2_ps3_highenergyLL(x, m2Q2, 1.)
//                  + C2_ps30_highenergyNLLA * eta_gamma / (C + eta_gamma));

//     if (m2mu2 == 1.)
//         return C30;

//     double Lmu = -log(m2mu2);
//     double Lmu2 = Lmu * Lmu;

//     return C30 + C2_ps31(x, m2Q2, nf) * Lmu + C2_ps32(x, m2Q2, nf) * Lmu2;
// }

// //==========================================================================================//
// //  Approximate quark coefficient funtcions for F2 at O(alpha_s^3) from
// //  [arXiv:1205.5727]. klmv = Kawamura, Lo Presti, Moch, Vogt. Lower band (i.e.
// //  approximation B).
// //
// //  Eq. (4.22) of Ref. [arXiv:1205.5727].
// //  This equation uses the exact expression of aQqPS30 given in Eq.
// //  (5.41, 5.42, 5.45) of Ref. [arXiv:1409.1135]
// //------------------------------------------------------------------------------------------//

// double C2_ps3_approximationB_klmv(double x, double m2Q2, double m2mu2, int nf) {

//     double x_max = 1. / (1. + 4. * m2Q2);

//     if (x >= x_max || x < 0)
//         return 0;

//     double pi3 = M_PI * M_PI * M_PI;

//     double beta = sqrt(1. - 4. * m2Q2 * x / (1. - x));

//     double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

//     double xi = 1. / m2Q2;

//     double f = 1. / (1. + exp(2. * (xi - 4)));

//     double delta = 0.8, D = 10.7;

//     double eta_delta = pow(eta, delta);

//     double beta3 = beta * beta * beta;

//     double C2_ps30_highenergyNLLB =
//         (64. * pi3) * (0.0245 * pow(log(1. / m2Q2) / log(5), 2) - 0.17) * 4.
//         / m2Q2 / x;

//     double C30 =
//         (1. - f) * beta3 * C2_ps3_highscale(x, m2Q2, 1., nf)
//         + f * beta3
//               * (-log(eta) / log(x) * C2_ps3_highenergyLL(x, m2Q2, 1.)
//                  + C2_ps30_highenergyNLLB * eta_delta / (D + eta_delta));

//     if (m2mu2 == 1.)
//         return C30;

//     double Lmu = -log(m2mu2);
//     double Lmu2 = Lmu * Lmu;

//     return C30 + C2_ps31(x, m2Q2, nf) * Lmu + C2_ps32(x, m2Q2, nf) * Lmu2;
// }

// //==========================================================================================//
// //  Approximate quark coefficient funtcions for F2 at O(alpha_s^3) from
// //  [arXiv:1205.5727]. klmv = Kawamura, Lo Presti, Moch, Vogt. Upper band (i.e.
// //  approximation A).
// //
// //  Eq. (4.21) of Ref. [arXiv:1205.5727].
// //  This equation uses the approximate form of aQqPS30 given in Eq. (3.52) of
// //  [arXiv:1205.5727] instead of the exact expression given in Eq.
// //  (5.41, 5.42, 5.45) of Ref. [arXiv:1409.1135] Used only as a benchmark
// //  against the plots of the paper [arXiv:1205.5727].
// //------------------------------------------------------------------------------------------//

// double
// C2_ps3_approximationA_klmv_paper(double x, double m2Q2, double m2mu2, int nf) {

//     double x_max = 1. / (1. + 4. * m2Q2);

//     if (x >= x_max || x < 0)
//         return 0;

//     double pi3 = M_PI * M_PI * M_PI;

//     double beta = sqrt(1. - 4. * m2Q2 * x / (1. - x));

//     double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

//     double xi = 1. / m2Q2;

//     double f = 1. / (1. + exp(2. * (xi - 4.)));

//     double gamma = 1.0, C = 20.0;

//     double eta_gamma = pow(eta, gamma);

//     double beta3 = beta * beta * beta;

//     double C2_ps30_highenergyNLLA =
//         (64. * pi3) * (0.004 * pow(log(1. / m2Q2) / log(5.), 4) - 0.125) * 4.
//         / m2Q2 / x;

//     double C30 =
//         (1. - f) * beta * C2_ps3_highscale_klmv_paper(x, m2Q2, 1., nf, 1)
//         + f * beta3
//               * (-log(eta) / log(x) * C2_ps3_highenergyLL(x, m2Q2, 1.)
//                  + C2_ps30_highenergyNLLA * eta_gamma / (C + eta_gamma));

//     if (m2mu2 == 1.)
//         return C30;

//     double Lmu = -log(m2mu2);
//     double Lmu2 = Lmu * Lmu;

//     return C30 + C2_ps31(x, m2Q2, nf) * Lmu + C2_ps32(x, m2Q2, nf) * Lmu2;
// }

// //==========================================================================================//
// //  Approximate quark coefficient funtcions for F2 at O(alpha_s^3) from
// //  [arXiv:1205.5727]. klmv = Kawamura, Lo Presti, Moch, Vogt. Lower band (i.e.
// //  approximation B).
// //
// //  Eq. (4.22) of Ref. [arXiv:1205.5727].
// //  This equation uses the approximate form of aQqPS30 given in Eq. (3.53) of
// //  [arXiv:1205.5727] instead of the exact expression given in Eq.
// //  (5.41, 5.42, 5.45) of Ref. [arXiv:1409.1135] Used only as a benchmark
// //  against the plots of the paper [arXiv:1205.5727].
// //------------------------------------------------------------------------------------------//

// double
// C2_ps3_approximationB_klmv_paper(double x, double m2Q2, double m2mu2, int nf) {

//     double x_max = 1. / (1. + 4. * m2Q2);

//     if (x >= x_max || x < 0)
//         return 0;

//     double pi3 = M_PI * M_PI * M_PI;

//     double beta = sqrt(1. - 4. * m2Q2 * x / (1. - x));

//     double eta = 0.25 / m2Q2 * (1. - x) / x - 1.;

//     double xi = 1. / m2Q2;

//     double f = 1. / (1. + exp(2. * (xi - 4.)));

//     double delta = 0.8, D = 10.7;

//     double eta_delta = pow(eta, delta);

//     double beta3 = beta * beta * beta;

//     double C2_ps30_highenergyNLLB =
//         (64. * pi3) * (0.0245 * pow(log(1. / m2Q2) / log(5.), 2) - 0.17) * 4.
//         / m2Q2 / x;

//     double C30 =
//         (1. - f) * beta3 * C2_ps3_highscale_klmv_paper(x, m2Q2, 1., nf, -1)
//         + f * beta3
//               * (-log(eta) / log(x) * C2_ps3_highenergyLL(x, m2Q2, 1.)
//                  + C2_ps30_highenergyNLLB * eta_delta / (D + eta_delta));

//     if (m2mu2 == 1.)
//         return C30;

//     double Lmu = -log(m2mu2);
//     double Lmu2 = Lmu * Lmu;

//     return C30 + C2_ps31(x, m2Q2, nf) * Lmu + C2_ps32(x, m2Q2, nf) * Lmu2;
// }
