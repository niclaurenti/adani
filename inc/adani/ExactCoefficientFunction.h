/*
 * =====================================================================================
 *
 *       Filename:  ExactCoefficientFunction.h
 *
 *    Description:  Header file for the
 * ExactCoefficientFunction.cc file.
 *
 *         Author:  Hanno un cuore differente
 *
 *  In this file there are the exact heavy coefficient
 *  functions (when known)
 *
 * =====================================================================================
 */

#ifndef Exact_h
#define Exact_h

#include "adani/CoefficientFunction.h"
#include "adani/Convolution.h"
#include "adani/SplittingFunction.h"

#include <vector>

//==========================================================================================//
//  class ExactCoefficientFunction
//------------------------------------------------------------------------------------------//

class ExactCoefficientFunction : public CoefficientFunction {

    public:
        ExactCoefficientFunction(
            const int &order, const char &kind, const char &channel,
            const double &abserr = 1e-3, const double &relerr = 1e-3,
            const int &dim = 1000, const bool &MCintegral = false,
            const int &MCcalls = 25000
        );
        ~ExactCoefficientFunction() override;

        double fx(double x, double m2Q2, double m2mu2, int nf) const override;
        double MuIndependentTerms(double x, double m2Q2, int nf) const override;
        double MuDependentTerms(
            double x, double m2Q2, double m2mu2, int nf
        ) const override;

        Value
        fxBand(double x, double m2Q2, double m2mu2, int nf) const override;

        void SetFunctions();

    private:
        double (ExactCoefficientFunction::*mu_indep_)(
            double, double, int
        ) const;
        double (ExactCoefficientFunction::*mu_dep_)(
            double, double, double, int
        ) const;

        std::vector<AbstractConvolution *> convolutions_lmu1_;
        std::vector<AbstractConvolution *> convolutions_lmu2_;

        ExactCoefficientFunction *gluon_as1_;
        ExactCoefficientFunction *gluon_as2_;
        ExactCoefficientFunction *quark_as2_;

        SplittingFunction *Pgq0_;
        SplittingFunction *Pgg0_;
        SplittingFunction *Pgq1_;
        SplittingFunction *Pqq0_;
        ConvolutedSplittingFunctions *Pgg0Pgq0_;
        ConvolutedSplittingFunctions *Pqq0Pgq0_;
        SplittingFunction *Pgg1_;
        SplittingFunction *Pqg0_;
        ConvolutedSplittingFunctions *Pgq0Pqg0_;

        Delta *delta_;

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
