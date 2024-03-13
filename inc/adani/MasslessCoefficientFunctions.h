/*
 * =====================================================================================
 *
 *       Filename:  MasslessCoefficientFunctions.h
 *
 *    Description:  Header file for the
 * MasslessCoefficientFunctions.cc file.
 *
 *         Author:  Viva el Futbol
 *
 *  In this file there are the massless coefficient
 * functions. Observe that for the O(alpha_s^2) and
 * O(alpha_s^3) it is used a parameterization of the exact
 * result, being the latter a very long equation.
 *
 * =====================================================================================
 */

#ifndef Massless_h
#define Massless_h

#include "adani/CoefficientFunction.h"

class MasslessCoefficientFunction : public CoefficientFunction {

    public:
        MasslessCoefficientFunction(const int& order, const char& kind, const char& channel) : CoefficientFunction(order, kind, channel) {} ;
        ~MasslessCoefficientFunction() {} ;

        double fx(double x, double /*m2Q2*/, double /*m2mu2*/, int nf) const override ;

    private:

        //==========================================================================================//
        //                      Massless coefficient functions
        //                      O(alpha_s)
        //------------------------------------------------------------------------------------------//

        double C2_g1_massless(double x, int nf) const;
        double CL_g1_massless(double x, int nf) const;

        //==========================================================================================//
        //                      Massless coefficient functions
        //                      O(alpha_s^2)
        //------------------------------------------------------------------------------------------//

        double C2_g2_massless(double x, int nf) const;
        double C2_ps2_massless(double x, int nf) const ;
        double CL_g2_massless(double x, int nf) const ;
        double CL_ps2_massless(double x, int nf) const;

        //==========================================================================================//
        //                      Massless coefficient functions
        //                      O(alpha_s^3)
        //------------------------------------------------------------------------------------------//

        double C2_g3_massless(double x, int nf) const;
        double C2_ps3_massless(double x, int nf) const;
        double CL_g3_massless(double x, int nf) const;
        double CL_ps3_massless(double x, int nf) const ;

        //==========================================================================================//
        //                      Charge factors
        //------------------------------------------------------------------------------------------//

        double fl11g(int nf) const;
        double fl11ps(int nf) const;

        //==========================================================================================//
        //                      Massless coefficient functions (parametrization)
        //                      O(alpha_s^2)
        //------------------------------------------------------------------------------------------//

        // double C2_g2_massless_param(double x, int nf) const;
        // double C2_ps2_massless_param(double x, int nf) const;
        // double CL_g2_massless_param(double x, int nf) const;
        // double CL_ps2_massless_param(double x, int nf) const;

        //==========================================================================================//
        //                      Massless coefficient functions (parametrization)
        //                      O(alpha_s^3)
        //------------------------------------------------------------------------------------------//

        // double C2_g3_massless_param(double x, int nf) const;
        // double C2_ps3_massless_param(double x, int nf) const;
        // double CL_g3_massless_param(double x, int nf) const;
        // double CL_ps3_massless_param(double x, int nf) const;

};

#endif
