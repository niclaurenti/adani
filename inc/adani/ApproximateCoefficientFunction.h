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
            const int &dim = 1000, const bool &MCintegral = false,
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
            const bool &NLL = true, const string &highscale_version = "klmv",
            const double &abserr = 1e-3, const double &relerr = 1e-3,
            const int &dim = 1000, const bool &MCintegral = false,
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
            const string &highscale_version = "klmv", const bool &lowxi = false,
            const double &abserr = 1e-3, const double &relerr = 1e-3,
            const int &dim = 1000, const bool &MCintegral = false,
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

#endif
