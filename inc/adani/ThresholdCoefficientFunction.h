/*
 * =====================================================================================
 *
 *       Filename:  ThresholdCoefficientFunctions.h
 *
 *         Author:  Daniele Adani
 *
 *    Description:  La garra charruaaaaaaa
 *
 *  In this file there are the coefficient functions in the
 *  threshold limit, i.e. s->4m^2 (where s is the partonic
 *  center of mass energy), or x -> xmax
 *
 * =====================================================================================
 */

#ifndef Threshold_h
#define Threshold_h

#include "adani/CoefficientFunction.h"

//==========================================================================================//
//  forward declaration of class ExactCoefficientFunction to avoid circular
//  dependencies
//------------------------------------------------------------------------------------------//

class ExactCoefficientFunction;

//==========================================================================================//
//  class ThresholdCoefficientFunction
//------------------------------------------------------------------------------------------//

class ThresholdCoefficientFunction : public CoefficientFunction {

    public:
        ThresholdCoefficientFunction(
            const int &order, const char &kind, const char &channel
        );
        ~ThresholdCoefficientFunction() override;

        void SetLegacyThreshold(const bool &legacy_threshold);

        double fx(double x, double m2Q2, double m2mu2, int nf) const override;
        Value
            fxBand(double x, double m2Q2, double m2mu2, int nf) const override;

        double BetaIndependentTerms(double x, double m2Q2, double m2mu2) const;

    private:
        // TODO: fx is the sum of a beta-dependent term and a beta-independent
        // in this way there is some repeated code. Split the pointers into
        // beta-indep and beta-dep
        bool legacy_threshold_;

        double (ThresholdCoefficientFunction::*expansion_beta_)(
            double, double, double, int
        ) const;
        double (ThresholdCoefficientFunction::*expansion_no_beta_)(
            double, double
        ) const;
        double (ThresholdCoefficientFunction::*threshold_as1_)(
            double, double
        ) const;
        Value (ThresholdCoefficientFunction::*fx_)(
            double, double, double, int
        ) const;

        ExactCoefficientFunction *exact_as1_;

        Value Order1(double x, double m2Q2, double m2mu2, int nf) const;
        Value PlainThreshold(double x, double m2Q2, double m2mu2, int nf) const;
        Value ModifiedThreshold(double x, double m2Q2, double m2mu2, int nf) const;

        void SetFunctions();

        //==========================================================================================//
        //                      Threshold (s -> 4m^2) coefficient
        //                      functions O(as)
        //------------------------------------------------------------------------------------------//

        double C2_g1_threshold(double x, double m2Q2) const;
        double CL_g1_threshold(double x, double m2Q2) const;

        //==========================================================================================//
        //                      Threshold (s -> 4m^2) coefficient
        //                      functions O(as^2)
        //------------------------------------------------------------------------------------------//

        double
            C2_g2_threshold_expansion(double x, double m2Q2, double m2mu2, int /*nf*/) const;
        double
            CL_g2_threshold_expansion(double x, double m2Q2, double m2mu2, int /*nf*/) const;
        double C2_g2_threshold_expansion_const(double m2Q2, double m2mu2) const;
        double CL_g2_threshold_expansion_const(double m2Q2, double m2mu2) const;

        //==========================================================================================//
        //                      Threshold (s -> 4m^2) coefficient
        //                      functions O(as^3)
        //------------------------------------------------------------------------------------------//

        double
            C2_g3_threshold_expansion(double x, double m2Q2, double m2mu2, int /*nf*/) const;
        double
            CL_g3_threshold_expansion(double x, double m2Q2, double m2mu2, int /*nf*/) const;
        double C2_g3_threshold_expansion_const(double m2Q2, double m2mu2) const;
        double CL_g3_threshold_expansion_const(double m2Q2, double m2mu2) const;

        //==========================================================================================//
        //  Functions needed for the threshold limit.
        //------------------------------------------------------------------------------------------//

        double c0(double xi) const;
        double c0_bar(double xi) const;

        //==========================================================================================//
        //  Function needed to make fx_ point to a zero function
        //------------------------------------------------------------------------------------------//

        Value ZeroFunction(
            double /*x*/, double /*m2Q2*/, double /*m2mu2*/, int /*nf*/
        ) const {
            return Value(0.);
        };
};

#endif
