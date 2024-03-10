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
#include "adani/Convolutions.h"

#include <vector>

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
        ExactCoefficientFunction(const int& order, const char& kind, const char& channel, const int& method_flag = 1, const double& abserr = 1e-3, const double& relerr = 1e-3, const int& MCcalls = 25000, const int& dim = 1000) ;
        // ExactCoefficientFunction() : ExactCoefficientFunction(1, '2', 'g') {} ;
        ~ExactCoefficientFunction() {} ;

        // get methods
        int GetMethodFlag() const {return method_flag_;};
        int GetAbserr() const {return abserr_;};
        int GetRelerr() const {return relerr_;};
        int GetMCcalls() const {return MCcalls_;};
        int GetDim() const {return dim_;} ;

        // set methods
        void SetMethodFlag(const int& method_flag) ;
        void SetAbserr(const double& abserr) ;
        void SetRelerr(const double& relerr);
        void SetMCcalls(const int& MCcalls);
        void SetDim(const int& dim);

        double fx(const double x, const double m2Q2, const double m2mu2, const int nf) const ;

        double MuIndependentTerms(const double x, const double m2Q2, const int nf) const override ;
        double MuDependentTerms(const double x, const double m2Q2, const double m2mu2, const int nf) const ;

    private:

        int method_flag_ ;
        double abserr_ ;
        double relerr_ ;
        int MCcalls_ ;
        int dim_ ;

        std::vector<Convolution> convolutions_;
        ExactCoefficientFunction* gluon_leadingorder_ = new ExactCoefficientFunction(1, GetKind(), 'g');

        //==========================================================================================//
        //                      Exact massive coefficient functions
        //                      O(alpha_s)
        //------------------------------------------------------------------------------------------//

        double C2_g1(const double x, const double m2Q2) const;
        double CL_g1(const double x, const double m2Q2) const;

        //==========================================================================================//
        //                      Exact massive coefficient functions
        //                      O(alpha_s^2)
        //------------------------------------------------------------------------------------------//

        double C2_g2(const double x, const double m2Q2, const double m2mu2) const;
        double C2_ps2(const double x, const double m2Q2, const double m2mu2) const;

        double CL_g2(const double x, const double m2Q2, const double m2mu2) const;
        double CL_ps2(const double x, const double m2Q2, const double m2mu2) const;

        //==========================================================================================//
        //  Exact massive coefficient functions O(alpha_s^2):
        //  mu-independent terms
        //------------------------------------------------------------------------------------------//

        double C2_g20(const double x, const double m2Q2) const;
        double CL_g20(const double x, const double m2Q2) const;
        double C2_ps20(const double x, const double m2Q2) const;
        double CL_ps20(const double x, const double m2Q2) const;

        //==========================================================================================//
        //  Exact massive coefficient functions O(alpha_s^2): terms
        //  proportional to log(mu^2/m^2)
        //------------------------------------------------------------------------------------------//

        double C2_g21(const double x, const double m2Q2) const;
        double CL_g21(const double x, const double m2Q2) const;
        double C2_ps21(const double x, const double m2Q2) const;
        double CL_ps21(const double x, const double m2Q2) const;

        //==========================================================================================//
        //  Exact massive coefficient functions O(alpha_s^3): terms
        //  proportional to log(mu^2/m^2)
        //------------------------------------------------------------------------------------------//

        double C2_ps31(const double x, const double m2Q2, const int nf) const;
        double CL_ps31(const double x, const double m2Q2, const int nf) const;
        double C2_g31(const double x, const double m2Q2, const int nf) const;
        double CL_g31(const double x, const double m2Q2, const int nf) const;

        //==========================================================================================//
        //  Exact massive coefficient functions O(alpha_s^3): terms
        //  proportional to log(mu^2/m^2)^2
        //------------------------------------------------------------------------------------------//

        double C2_ps32(const double x, const double m2Q2, const int nf) const;
        double CL_ps32(const double x, const double m2Q2, const int nf) const;

        // These two expressions are integrated with montcarlo
        // methods

        double C2_g32(const double x, const double m2Q2, const int nf) const;
        double CL_g32(const double x, const double m2Q2, const int nf) const;

};

#endif
