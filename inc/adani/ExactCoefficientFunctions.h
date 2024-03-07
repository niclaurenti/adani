/*
 * =====================================================================================
 *
 *       Filename:  ExactCoefficientFunctions.h
 *
 *    Description:  Header file for the
 * ExactCoefficientFunctions.cc file.
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

//==========================================================================================//
//                      The notation used is the following:
//                      C^{(k)} = k-th order expansion in
//                      terms of \alpha_s^{[nf]}
//
//------------------------------------------------------------------------------------------//
//==========================================================================================//
//                      Notation:
//                      m2Q2 = m^2/Q^2
//                      m2mu2 = m^2/mu^2
//------------------------------------------------------------------------------------------//

class ExactCoefficientFunction : public CoefficientFunction {

    public:
        ExactCoefficientFunction(const int order, const char kind, const char channel, const int method_flag = 1, const double abserr = 1e-3, const double relerr = 1e-3, const int MCcalls = 25000) ;
        ExactCoefficientFunction() : ExactCoefficientFunction(1, '2', 'g') {} ;
        ~ExactCoefficientFunction() {} ;

        // get methods
        int GetMethodFlag() {return method_flag_;};
        int GetAbserr() {return abserr_;};
        int GetRelerr() {return relerr_;};
        int GetMCcalls() {return MCcalls_;};

        // set methods
        void SetMethodFlag(const int method_flag) ;
        void SetAbserr(const double abserr) ;
        void SetRelerr(const double relerr);
        void SetMCcalls(const int MCcalls);

        double fx(double x, double m2Q2, double m2mu2, int nf) ;

    private:

        int method_flag_ ;
        double abserr_ ;
        double relerr_ ;
        int MCcalls_ ;

        //==========================================================================================//
        //                      Exact massive coefficient functions
        //                      O(alpha_s)
        //------------------------------------------------------------------------------------------//

        double C2_g1(double x, double m2Q2);
        double CL_g1(double x, double m2Q2);

        //==========================================================================================//
        //                      Exact massive coefficient functions
        //                      O(alpha_s^2)
        //------------------------------------------------------------------------------------------//

        double C2_g2(double x, double m2Q2, double m2mu2);
        double C2_ps2(double x, double m2Q2, double m2mu2);

        double CL_g2(double x, double m2Q2, double m2mu2);
        double CL_ps2(double x, double m2Q2, double m2mu2);

        //==========================================================================================//
        //  Exact massive coefficient functions O(alpha_s^2):
        //  mu-independent terms
        //------------------------------------------------------------------------------------------//

        double C2_g20(double x, double m2Q2);
        double CL_g20(double x, double m2Q2);
        double C2_ps20(double x, double m2Q2);
        double CL_ps20(double x, double m2Q2);

        //==========================================================================================//
        //  Exact massive coefficient functions O(alpha_s^2): terms
        //  proportional to log(mu^2/m^2)
        //------------------------------------------------------------------------------------------//

        double C2_g21(double x, double m2Q2);
        double CL_g21(double x, double m2Q2);
        double C2_ps21(double x, double m2Q2);
        double CL_ps21(double x, double m2Q2);

        //==========================================================================================//
        //  Exact massive coefficient functions O(alpha_s^3): terms
        //  proportional to log(mu^2/m^2)
        //------------------------------------------------------------------------------------------//

        double C2_ps31(double x, double m2Q2, int nf);
        double CL_ps31(double x, double m2Q2, int nf);
        double C2_g31(double x, double m2Q2, int nf);
        double CL_g31(double x, double m2Q2, int nf);

        //==========================================================================================//
        //  Exact massive coefficient functions O(alpha_s^3): terms
        //  proportional to log(mu^2/m^2)^2
        //------------------------------------------------------------------------------------------//

        double C2_ps32(double x, double m2Q2, int nf);
        double CL_ps32(double x, double m2Q2, int nf);

        // These two expressions are integrated with montcarlo
        // methods

        double C2_g32(double x, double m2Q2, int nf);
        double CL_g32(double x, double m2Q2, int nf);

};

#endif
