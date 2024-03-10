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
        // MasslessCoefficientFunction() : CoefficientFunction() {} ;
        ~MasslessCoefficientFunction() {} ;

        double fx(const double x, const double m2Q2, const double m2mu2, const int nf) const ;

    private:

        //==========================================================================================//
        //                      Massless coefficient functions
        //                      O(alpha_s)
        //------------------------------------------------------------------------------------------//

        double C2_g1_massless(const double x, const int nf) const;
        double CL_g1_massless(const double x, const int nf) const;

        //==========================================================================================//
        //                      Massless coefficient functions
        //                      O(alpha_s^2)
        //------------------------------------------------------------------------------------------//

        double C2_g2_massless(const double x, const int nf) const;
        double C2_ps2_massless(const double x, const int nf) const ;
        double CL_g2_massless(const double x, const int nf) const ;
        double CL_ps2_massless(const double x, const int nf) const;

        //==========================================================================================//
        //                      Massless coefficient functions
        //                      O(alpha_s^3)
        //------------------------------------------------------------------------------------------//

        double C2_g3_massless(const double x, const int nf) const;
        double C2_ps3_massless(const double x, const int nf) const;
        double CL_g3_massless(const double x, const int nf) const;
        double CL_ps3_massless(const double x, const int nf) const ;

        //==========================================================================================//
        //                      Charge factors
        //------------------------------------------------------------------------------------------//

        double fl11g(const int nf) const;
        double fl11ps(const int nf) const;

        //==========================================================================================//
        //                      Massless coefficient functions (parametrization)
        //                      O(alpha_s^2)
        //------------------------------------------------------------------------------------------//

        // double C2_g2_massless_param(const double x, const int nf) const;
        // double C2_ps2_massless_param(const double x, const int nf) const;
        // double CL_g2_massless_param(const double x, const int nf) const;
        // double CL_ps2_massless_param(const double x, const int nf) const;

        //==========================================================================================//
        //                      Massless coefficient functions (parametrization)
        //                      O(alpha_s^3)
        //------------------------------------------------------------------------------------------//

        // double C2_g3_massless_param(const double x, const int nf) const;
        // double C2_ps3_massless_param(const double x, const int nf) const;
        // double CL_g3_massless_param(const double x, const int nf) const;
        // double CL_ps3_massless_param(const double x, const int nf) const;

};

#endif
