/*
 * =====================================================================================
 *
 *       Filename:  MasslessCoefficientFunction.h
 *
 *    Description:  Header file for the
 * MasslessCoefficientFunction.cc file.
 *
 *         Author:  Viva el Futbol
 *
 *  In this file there are the massless coefficient
 * functions. Observe that for the O(as^2) and
 * O(as^3) it is used a parameterization of the exact
 * result, being the latter a very long equation.
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
        ~MasslessCoefficientFunction() override{};

        double fx(double x, double m2Q2, double m2mu2, int nf) const override;
        double MuIndependentTerms(double /*x*/, double /*m2Q2*/, int /*nf*/)
            const override;
        double MuIndependentTerms(double x, int nf) const;
        double MuDependentTerms(
            double /*x*/, double /*m2Q2*/, double /*m2mu2*/, int /*nf*/
        ) const override;

        Value
        fxBand(double x, double m2Q2, double m2mu2, int nf) const override;

        void SetFunctions();

    private:
        double (MasslessCoefficientFunction::*mu_indep_)(double, int) const;

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
