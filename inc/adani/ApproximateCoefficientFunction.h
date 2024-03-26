/*
 * =====================================================================================
 *
 *       Filename:  ApproximateCoefficientFunction.h
 *
 *    Description:  Header file for the
 * ApproximateCoefficientFunction.cc file.
 *
 *         Author:  L'artiglio che graffia
 *
 *  In this file there is the approximation for the unknown
 * O(as^3) DIS massive coefficient functions.
 *
 * =====================================================================================
 */

#ifndef Approximate_h
#define Approximate_h

#include "adani/AsymptoticCoefficientFunction.h"
#include "adani/CoefficientFunction.h"
#include "adani/ExactCoefficientFunction.h"
#include "adani/ThresholdCoefficientFunction.h"

//==========================================================================================//
//  struct approximation_parameters: parameters of the damping functions
//------------------------------------------------------------------------------------------//

struct approximation_parameters {
        double A;
        double B;
        double C;
        double D;
};

//==========================================================================================//
//  struct variation_parameters: parameters of the variation of the damping
//  functions
//------------------------------------------------------------------------------------------//

struct variation_parameters {
        double var;
        double fact;
};

//==========================================================================================//
//  struct klmv_params: parameters for klmv approximation
//------------------------------------------------------------------------------------------//

struct klmv_params {
        double gamma;
        double C;
        double log_coeff;
        double log_pow;
        double const_coeff;
};

//==========================================================================================//
//  class AbstractApproximate
//------------------------------------------------------------------------------------------//

class AbstractApproximate : public CoefficientFunction {
    public:
        AbstractApproximate(
            const int &order, const char &kind, const char &channel,
            const double &abserr = 1e-3, const double &relerr = 1e-3,
            const int &dim = 1000, const int &method_flag = 0,
            const int &MCcalls = 25000
        );
        ~AbstractApproximate();

        double MuIndependentTerms(double x, double m2Q2, int nf) const override;

        Value
        fxBand(double x, double m2Q2, double m2mu2, int nf) const override;
        double MuDependentTerms(
            double x, double m2Q2, double m2mu2, int nf
        ) const override;

    private:
        ExactCoefficientFunction *muterms_;
};

//==========================================================================================//
//  class ApproximateCoefficientFunction
//------------------------------------------------------------------------------------------//

class ApproximateCoefficientFunction : public AbstractApproximate {
    public:
        ApproximateCoefficientFunction(
            const int &order, const char &kind, const char &channel,
            const bool &NLL = true, const bool &exact_highscale = false,
            const bool &revised_approx_highscale = true,
            const double &abserr = 1e-3, const double &relerr = 1e-3,
            const int &dim = 1000, const int &method_flag = 0,
            const int &MCcalls = 25000
        );
        ~ApproximateCoefficientFunction() override;

        Value
        MuIndependentTermsBand(double x, double m2Q2, int nf) const override;

    private:
        ThresholdCoefficientFunction *threshold_;
        AsymptoticCoefficientFunction *asymptotic_;

        struct approximation_parameters approximation_;
        struct variation_parameters variation_;

        double Approximation(
            double x, double m2Q2, double asy, double thresh, double A,
            double B, double C, double D
        ) const;
};

//==========================================================================================//
//  class ApproximateCoefficientFunctionKLMV: Approximate coefficient functions
//  from [arXiv:1205.5727] klmv = Kawamura, Lo Presti, Moch, Vogt
//------------------------------------------------------------------------------------------//

class ApproximateCoefficientFunctionKLMV : public AbstractApproximate {
    public:
        ApproximateCoefficientFunctionKLMV(
            const int &order, const char &kind, const char &channel,
            const bool &revised_approx_highscale = true, const bool &lowxi = false,
            const double &abserr = 1e-3, const double &relerr = 1e-3,
            const int &dim = 1000, const int &method_flag = 0,
            const int &MCcalls = 25000
        );
        ~ApproximateCoefficientFunctionKLMV() override;

        Value
        MuIndependentTermsBand(double x, double m2Q2, int nf) const override;

    private:
        ThresholdCoefficientFunction *threshold_;
        HighScaleCoefficientFunction *highscale_;
        HighEnergyCoefficientFunction *highenergy_;

        struct klmv_params params_A_;
        struct klmv_params params_B_;

        double ApproximationA(
            double x, double m2Q2, double he_ll, double he_nll, double hs,
            double thr, double thr_const, double gamma, double C
        ) const;
        double ApproximationB(
            double x, double m2Q2, double he_ll, double he_nll, double hs,
            double thr, double thr_const, double delta, double D
        ) const;
        Value ApproximateNLL(double x, double m2Q2) const;
};

// //==========================================================================================//
// //                      Approximate coefficient functions
// //                      O(as^3)
// //------------------------------------------------------------------------------------------//

// #define default_method 0

// double C2_g3_approximation(
//     double x, double m2Q2, double m2mu2, int nf, int v = 0,
//     int method_flag = default_method
// );
// double C2_g30_approximation(double x, double m2Q2, int nf);
// double C2_g30_approximation_implicit(
//     double x, double m2Q2, int nf, double A, double B, double C, double D,
//     double a, double b, int v1, int v2
// );
// double C2_g30_approximation_BAND(
//     double x, double m2Q2, int nf, int v, double var, double fact
// );
// double
// C2_ps3_approximation(double x, double m2Q2, double m2mu2, int nf, int v = 0);
// double C2_ps30_approximation(double x, double m2Q2, int nf);
// double C2_ps30_approximation_implicit(
//     double x, double m2Q2, int nf, double A, double B, double C, double D,
//     double a, double b, int v
// );
// double C2_ps30_approximation_BAND(
//     double x, double m2Q2, int nf, int v, double var, double fact
// );

// double CL_g3_approximation(
//     double x, double m2Q2, double m2mu2, int nf, int v = 0,
//     int method_flag = default_method
// );
// double CL_g30_approximation(double x, double m2Q2, int nf);
// double CL_g30_approximation_implicit(
//     double x, double m2Q2, int nf, double A, double B, double C, double D,
//     double a, double b, int v
// );
// double CL_g30_approximation_BAND(
//     double x, double m2Q2, int nf, int v, double var, double fact
// );

// double
// CL_ps3_approximation(double x, double m2Q2, double m2mu2, int nf, int v = 0);
// double CL_ps30_approximation(double x, double m2Q2, int nf);
// double CL_ps30_approximation_implicit(
//     double x, double m2Q2, int nf, double A, double B, double C, double D,
//     double a, double b, int v
// );
// double CL_ps30_approximation_BAND(
//     double x, double m2Q2, int nf, int v, double var, double fact
// );

// //==========================================================================================//
// //              Approximate coefficient functions
// //              O(as^2) from [arXiv:1205.5727] klmv =
// //              Kawamura, Lo Presti, Moch, Vogt
// //------------------------------------------------------------------------------------------//

// double C2_g2_approximationA_klmv(double x, double m2Q2, double m2mu2);
// double C2_g2_approximationB_klmv(double x, double m2Q2, double m2mu2);

// double C2_ps2_approximationA_klmv(double x, double m2Q2, double m2mu2);
// double C2_ps2_approximationB_klmv(double x, double m2Q2, double m2mu2);

// //==========================================================================================//
// //              Approximate coefficient functions
// //              O(as^3) from [arXiv:1205.5727] klmv =
// //              Kawamura, Lo Presti, Moch, Vogt
// //------------------------------------------------------------------------------------------//

// //==========================================================================================//
// // the functions labeled with 'paper' use some approximate
// // results for which at the time of the paper
// // [arXiv:1205.5727] the exact result was not known (like
// // aQqPS30) or for which now we have a better approximation
// // (aQg30_B). They are only used as a benchmark against the
// // plots of the paper
// //------------------------------------------------------------------------------------------//

// double C2_g3_approximationA_klmv(
//     double x, double m2Q2, double m2mu2, int nf,
//     int method_flag = default_method
// );
// double C2_g3_approximationB_klmv(
//     double x, double m2Q2, double m2mu2, int nf,
//     int method_flag = default_method
// );
// double C2_g3_approximationB_klmv_paper(
//     double x, double m2Q2, double m2mu2, int nf,
//     int method_flag = default_method
// );
// double C2_g3_approximationBlowxi_klmv(
//     double x, double m2Q2, double m2mu2, int nf,
//     int method_flag = default_method
// );

// double C2_ps3_approximationA_klmv(double x, double m2Q2, double m2mu2, int
// nf); double C2_ps3_approximationB_klmv(double x, double m2Q2, double m2mu2,
// int nf);

// double
// C2_ps3_approximationA_klmv_paper(double x, double m2Q2, double m2mu2, int
// nf); double C2_ps3_approximationB_klmv_paper(double x, double m2Q2, double
// m2mu2, int nf);

#endif
