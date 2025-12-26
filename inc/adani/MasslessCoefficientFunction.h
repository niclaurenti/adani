/*
 * =====================================================================================
 *
 *       Filename:  MasslessCoefficientFunction.h
 *
 *         Author:  Daniele Adani
 *
 *    Description:  Viva el Futbol
 *
 *  In this file there are the massless coefficient
 *  functions. Observe that for the O(as^2) and
 *  O(as^3) it is used a parameterization of the exact
 *  result, being the latter a very long equation.
 *
 *  OBSERVATION: massless coefficient functions depend on x, Q2/mu2, nf
 *  but in order to achieve polymorfism with the class CoefficientFunction
 *  (the aim is call it inside the class Convolution) for all the methods,
 *  it's always present also a function fx(double, double, double, int)
 *  that returns fx(double, double, int)
 *
 * =====================================================================================
 */

#ifndef Massless_h
#define Massless_h

#include "adani/CoefficientFunction.h"

//==========================================================================================//
//  class MasslessCoefficientFunction
//------------------------------------------------------------------------------------------//

class MasslessCoefficientFunction : public CoefficientFunction {

    public:
        MasslessCoefficientFunction(
            const int &order, const char &kind, const char &channel
        );
        MasslessCoefficientFunction(const MasslessCoefficientFunction &obj);
        ~MasslessCoefficientFunction() override = default;

        double
            fx(double x, double /*m2Q2*/, double Q2mu2, int nf) const override;
        double fx(double x, double Q2mu2, int nf) const;

        double MuIndependentTerms(
            double x, double /*m2mu2*/, int nf
        ) const override;
        double MuIndependentTerms(double x, int nf) const;

        // TODO: to be implemented
        double MuDependentTerms(
            double /*x*/, double /*m2Q2*/, double /*Q2mu2*/, int /*nf*/
        ) const override;
        // TODO: to be implemented
        double MuDependentTerms(
            double /*x*/, double /*Q2mu2*/, int /*nf*/
        ) const;

        Value
            fxBand(double x, double m2Q2, double m2mu2, int nf) const override;
        Value fxBand(double x, double Q2mu2, int nf) const;

    private:
        double (MasslessCoefficientFunction::*mu_indep_)(double, int) const;

        void SetFunctions();

        //==========================================================================================//
        //                      Massless coefficient functions
        //                      O(as)
        //------------------------------------------------------------------------------------------//

        double C2_g1_massless(double x, int nf) const;
        double CL_g1_massless(double x, int nf) const;

        //==========================================================================================//
        //                      Massless coefficient functions
        //                      O(as^2)
        //------------------------------------------------------------------------------------------//

        double C2_g2_massless(double x, int nf) const;
        double C2_ps2_massless(double x, int nf) const;
        double CL_g2_massless(double x, int nf) const;
        double CL_ps2_massless(double x, int nf) const;

        //==========================================================================================//
        //                      Massless coefficient functions
        //                      O(as^3)
        //------------------------------------------------------------------------------------------//

        double C2_g3_massless(double x, int nf) const;
        double C2_ps3_massless(double x, int nf) const;
        double CL_g3_massless(double x, int nf) const;
        double CL_ps3_massless(double x, int nf) const;

        //==========================================================================================//
        //                      Charge factors
        //------------------------------------------------------------------------------------------//

        double fl11g(int nf) const;
        double fl11ps(int nf) const;

        //==========================================================================================//
        //                      Massless coefficient functions (parametrization)
        //                      O(as^2)
        //------------------------------------------------------------------------------------------//

        // double C2_g2_massless_param(double x, int nf) const;
        // double C2_ps2_massless_param(double x, int nf) const;
        // double CL_g2_massless_param(double x, int nf) const;
        // double CL_ps2_massless_param(double x, int nf) const;

        //==========================================================================================//
        //                      Massless coefficient functions (parametrization)
        //                      O(as^3)
        //------------------------------------------------------------------------------------------//

        // double C2_g3_massless_param(double x, int nf) const;
        // double C2_ps3_massless_param(double x, int nf) const;
        // double CL_g3_massless_param(double x, int nf) const;
        // double CL_ps3_massless_param(double x, int nf) const;
};

#endif
