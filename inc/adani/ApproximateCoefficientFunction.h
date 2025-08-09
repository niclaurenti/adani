/*
 * =====================================================================================
 *
 *       Filename:  ApproximateCoefficientFunction.h
 *
 *         Author:  Daniele Adani
 *
 *    Description:  L'artiglio che graffia
 *
 *  In this file there is the approximation for the unknown
 *  O(as^3) DIS massive coefficient functions.
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
        double var1;
        double var2;
};

//==========================================================================================//
//  struct klmv_params: parameters for klmv approximation
//------------------------------------------------------------------------------------------//

struct klmv_params {
        double eta_exponent;
        double shift;
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
            const int &dim = 1000
        );
        ~AbstractApproximate();

        void SetDoubleIntegralMethod(
            const string &double_int_method, const double &abserr = 1e-3,
            const double &relerr = 1e-3, const int &dim = 1000,
            const int &MCcalls = 25000
        );

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
            const bool &NLL = true, const string &highscale_version = "exact",
            const double &abserr = 1e-3, const double &relerr = 1e-3,
            const int &dim = 1000
        );
        ~ApproximateCoefficientFunction() override;

        void SetLegacyThreshold(const bool &legacy_threshold);
        void SetLegacyPowerTerms(const bool &legacy_pt);
        void SetLegacyApproximation(const bool &legacy_appr);

        Value MuIndependentTermsBand(
            double x, double m2Q2, int nf
        ) const override;

    private:
        ThresholdCoefficientFunction *threshold_;
        AsymptoticCoefficientFunction *asymptotic_;

        Value (ApproximateCoefficientFunction::*fx_)(
            double, double, int
        ) const;

        bool legacy_appr_;

        struct approximation_parameters *approximation_;
        struct variation_parameters *variation_;

        double ApproximationLegacyForm(
            double x, double m2Q2, double asy, double thresh, double A,
            double B, double C, double D
        ) const;

        Value Approximation(double x, double m2Q2, int nf) const;
        Value ApproximationLegacy(double x, double m2Q2, int nf) const;

        void SetLegacyParameters();
};

//==========================================================================================//
//  class ApproximateCoefficientFunctionKLMV: Approximate coefficient functions
//  from [arXiv:1205.5727] klmv = Kawamura, Lo Presti, Moch, Vogt
//------------------------------------------------------------------------------------------//

class ApproximateCoefficientFunctionKLMV : public AbstractApproximate {
    public:
        ApproximateCoefficientFunctionKLMV(
            const int &order, const char &kind, const char &channel,
            const string &highscale_version = "exact",
            const bool &lowxi = false, const double &abserr = 1e-3,
            const double &relerr = 1e-3, const int &dim = 1000
        );
        ~ApproximateCoefficientFunctionKLMV() override;

        bool GetLowXi() const {return lowxi_;}
        void SetLowXi(const bool& lowxi);

        Value MuIndependentTermsBand(
            double x, double m2Q2, int nf
        ) const override;

    private:
        ThresholdCoefficientFunction *threshold_;
        HighScaleCoefficientFunction *highscale_;
        HighEnergyCoefficientFunction *highenergy_;

        Value (ApproximateCoefficientFunctionKLMV::*fx_)(
                    double, double, int
                ) const;

        struct klmv_params *params_A_;
        struct klmv_params *params_B_;

        bool lowxi_;

        Value Order2(double x, double m2Q2, int nf) const;
        Value Order3(double x, double m2Q2, int nf) const;

        Value ApproximateNLL(double x, double m2Q2) const;

        void SetFunctions();
};

#endif
