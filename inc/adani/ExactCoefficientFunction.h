/*
 * =====================================================================================
 *
 *       Filename:  ExactCoefficientFunction.h
 *
 *         Author:  Daniele Adani
 *
 *    Description:  Hanno un cuore differente
 *
 *  In this file there are the exact heavy coefficient
 *  functions (when known)
 *
 * =====================================================================================
 */

#ifndef Exact_h
#define Exact_h

#include "adani/AsymptoticCoefficientFunction.h"
#include "adani/CoefficientFunction.h"
#include "adani/CommonTypes.h"
#include "adani/Convolution.h"
#include "adani/SplittingFunction.h"

#include <memory>
#include <vector>

//==========================================================================================//
//  forward declaration of class ThresholdCoefficientFunction to avoid circular
//  dependencies
//------------------------------------------------------------------------------------------//

class ThresholdCoefficientFunction;

//==========================================================================================//
//  class ExactCoefficientFunction
//------------------------------------------------------------------------------------------//

class ExactCoefficientFunction : public CoefficientFunction {

    public:
        ExactCoefficientFunction(
            const int &order, const char &kind, const char &channel,
            const double &abserr = 1e-3, const double &relerr = 1e-3,
            const int &dim = 1000
        );
        ExactCoefficientFunction(ExactCoefficientFunction &obj);
        ~ExactCoefficientFunction() override = default;

        double GetAbsErr() const;
        double GetRelErr() const;
        int GetDim() const;
        int GetMCcalls() const;
        DoubleIntegralMethod GetDoubleIntegralMethod() const {
            return double_int_meth_;
        };

        void SetDoubleIntegralMethod(
            const DoubleIntegralMethod &double_int_method,
            const double &abserr = 1e-3, const double &relerr = 1e-3,
            const int &dim = 100, const int &MCcalls = 25000
        );

        double fx(double x, double m2Q2, double m2mu2, int nf) const override;
        double MuIndependentTerms(double x, double m2Q2, int nf) const override;
        double MuDependentTerms(
            double x, double m2Q2, double m2mu2, int nf
        ) const override;

        Value
            fxBand(double x, double m2Q2, double m2mu2, int nf) const override;

    private:
        double (ExactCoefficientFunction::*mu_indep_)(
            double, double, int
        ) const;
        double (ExactCoefficientFunction::*mu_dep_)(
            double, double, double, int
        ) const;

        std::unique_ptr<AsymptoticCoefficientFunction> asy_;
        std::unique_ptr<ThresholdCoefficientFunction> thr_;

        std::vector<std::unique_ptr<AbstractConvolution> > convolutions_lmu1_;
        std::vector<std::unique_ptr<AbstractConvolution> > convolutions_lmu2_;

        std::shared_ptr<const ExactCoefficientFunction> gluon_as1_;
        std::shared_ptr<const ExactCoefficientFunction> gluon_as2_;
        std::shared_ptr<const ExactCoefficientFunction> quark_as2_;

        std::shared_ptr<const SplittingFunction> Pgq0_;
        std::shared_ptr<const SplittingFunction> Pgg0_;
        std::shared_ptr<const SplittingFunction> Pgq1_;
        std::shared_ptr<const SplittingFunction> Pqq0_;
        std::shared_ptr<const ConvolutedSplittingFunctions> Pgg0Pgq0_;
        std::shared_ptr<const ConvolutedSplittingFunctions> Pqq0Pgq0_;
        std::shared_ptr<const SplittingFunction> Pgg1_;
        std::shared_ptr<const SplittingFunction> Pqg0_;
        std::shared_ptr<const ConvolutedSplittingFunctions> Pgq0Pqg0_;
        std::shared_ptr<const ConvolutedSplittingFunctions> Pgg0Pgg0_;

        std::shared_ptr<const Delta> delta_;

        DoubleIntegralMethod double_int_meth_;

        void SetFunctions();

        //==========================================================================================//
        //                      Exact massive coefficient functions
        //                      O(as)
        //------------------------------------------------------------------------------------------//

        double C2_g1(double x, double m2Q2, int /*nf*/) const;
        double CL_g1(double x, double m2Q2, int /*nf*/) const;

        //==========================================================================================//
        //  Exact massive coefficient functions O(as^2):
        //  mu-independent terms
        //------------------------------------------------------------------------------------------//

        double C2_g20(double x, double m2Q2, int /*nf*/) const;
        double CL_g20(double x, double m2Q2, int /*nf*/) const;
        double C2_ps20(double x, double m2Q2, int /*nf*/) const;
        double CL_ps20(double x, double m2Q2, int /*nf*/) const;

        //==========================================================================================//
        //  Exact massive coefficient functions O(as^2): terms
        //  proportional to log(mu^2/m^2)
        //------------------------------------------------------------------------------------------//

        double C_g21(double x, double m2Q2) const;
        double C_ps21(double x, double m2Q2) const;

        //==========================================================================================//
        //  Exact massive coefficient functions O(as^2):
        //  mu-dependent terms
        //------------------------------------------------------------------------------------------//

        double
            C_ps2_MuDep(double x, double m2Q2, double m2mu2, int /*nf*/) const;
        double
            C_g2_MuDep(double x, double m2Q2, double m2mu2, int /*nf*/) const;

        //==========================================================================================//
        //  Exact massive coefficient functions O(as^3): terms
        //  proportional to log(mu^2/m^2)
        //------------------------------------------------------------------------------------------//

        double C_ps31(double x, double m2Q2, int nf) const;
        double C_g31(double x, double m2Q2, int nf) const;

        //==========================================================================================//
        //  Exact massive coefficient functions O(as^3): terms
        //  proportional to log(mu^2/m^2)^2
        //------------------------------------------------------------------------------------------//

        double C_ps32(double x, double m2Q2, int nf) const;
        double C_g32(double x, double m2Q2, int nf) const;

        //==========================================================================================//
        //  Exact massive coefficient functions O(as^3):
        //  mu-dependent terms
        //------------------------------------------------------------------------------------------//

        double C_ps3_MuDep(double x, double m2Q2, double m2mu2, int nf) const;
        double C_g3_MuDep(double x, double m2Q2, double m2mu2, int nf) const;

        //==========================================================================================//
        //  Function needed to make muindep_ and mudep_ point to a zero function
        //------------------------------------------------------------------------------------------//

        double ZeroFunction(
            double /*x*/, double /*m2Q2*/, double /*m2mu2*/, int /*nf*/
        ) const {
            return 0.;
        };
        double WarningFunc(double /*x*/, double /*m2Q2*/, int /*nf*/) const;
};

#endif
