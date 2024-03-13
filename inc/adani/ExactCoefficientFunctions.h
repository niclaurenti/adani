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
//                      terms of a_s^{[nf]}
//
//------------------------------------------------------------------------------------------//
//==========================================================================================//
//                      Notation:
//                      m2Q2 = m^2/Q^2
//                      m2mu2 = m^2/mu^2
//------------------------------------------------------------------------------------------//

class ExactCoefficientFunction : public CoefficientFunction {

    public:
        ExactCoefficientFunction(const int& order, const char& kind, const char& channel, const double& abserr = 1e-3, const double& relerr = 1e-3, const int& dim = 1000, const int& method_flag = 1, const int& MCcalls = 25000) ;
        ~ExactCoefficientFunction() override ;

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

        double fx(double x, double m2Q2, double m2mu2, int nf) const override ;

        double MuIndependentTerms(double x, double m2Q2, int nf) const override ;
        double MuDependentTerms(double x, double m2Q2, double m2mu2, int nf) const override ;

        void SetFunctions();

    private:

        int method_flag_ ;
        double abserr_ ;
        double relerr_ ;
        int MCcalls_ ;
        int dim_ ;

        double (ExactCoefficientFunction::*mu_indep_)(double, double, int) const;
        double (ExactCoefficientFunction::*mu_dep_)(double, double, double, int) const;

        std::vector<AbstractConvolution*> convolutions_lmu1_;
        std::vector<AbstractConvolution*> convolutions_lmu2_;

        ExactCoefficientFunction* gluon_lo_ ;
        ExactCoefficientFunction* gluon_nlo_ ;
        ExactCoefficientFunction* quark_nlo_ ;

        SplittingFunction* Pgq0_;
        SplittingFunction* Pgg0_;
        SplittingFunction* Pgq1_;
        SplittingFunction* Pqq0_;
        ConvolutedSplittingFunctions* Pgg0Pgq0_;
        ConvolutedSplittingFunctions* Pqq0Pgq0_;
        SplittingFunction* Pgg1_;
        SplittingFunction* Pqg0_;
        ConvolutedSplittingFunctions* Pgq0Pqg0_;

        Delta* delta_;

        //==========================================================================================//
        //                      Exact massive coefficient functions
        //                      O(alpha_s)
        //------------------------------------------------------------------------------------------//

        double C2_g1(double x, double m2Q2, int /*nf*/) const;
        double CL_g1(double x, double m2Q2, int /*nf*/) const;

        //==========================================================================================//
        //                      Exact massive coefficient functions
        //                      O(alpha_s^2)
        //------------------------------------------------------------------------------------------//

        // double C2_g2(double x, double m2Q2, double m2mu2) const;
        // double C2_ps2(double x, double m2Q2, double m2mu2) const;

        // double CL_g2(double x, double m2Q2, double m2mu2) const;
        // double CL_ps2(double x, double m2Q2, double m2mu2) const;

        //==========================================================================================//
        //  Exact massive coefficient functions O(alpha_s^2):
        //  mu-independent terms
        //------------------------------------------------------------------------------------------//

        double C2_g20(double x, double m2Q2, int /*nf*/) const;
        double CL_g20(double x, double m2Q2, int /*nf*/) const;
        double C2_ps20(double x, double m2Q2, int /*nf*/) const;
        double CL_ps20(double x, double m2Q2, int /*nf*/) const;

        //==========================================================================================//
        //  Exact massive coefficient functions O(alpha_s^2): terms
        //  proportional to log(mu^2/m^2)
        //------------------------------------------------------------------------------------------//

        double C_g21(double x, double m2Q2) const;
        double C_ps21(double x, double m2Q2) const;

        double C_ps2_MuDep(double x, double m2Q2, double m2mu2, int /*nf*/) const ;
        double C_g2_MuDep(double x, double m2Q2, double m2mu2, int /*nf*/) const ;

        //==========================================================================================//
        //  Exact massive coefficient functions O(alpha_s^3): terms
        //  proportional to log(mu^2/m^2)
        //------------------------------------------------------------------------------------------//

        double C_ps31(double x, double m2Q2, int nf) const;
        double C_g31(double x, double m2Q2, int nf) const;

        //==========================================================================================//
        //  Exact massive coefficient functions O(alpha_s^3): terms
        //  proportional to log(mu^2/m^2)^2
        //------------------------------------------------------------------------------------------//

        double C_ps32(double x, double m2Q2, int nf) const;

        // These two expressions are integrated with montcarlo
        // methods

        double C_g32(double x, double m2Q2, int nf) const;

        double C_ps3_MuDep(double x, double m2Q2, double m2mu2, int nf) const ;
        double C_g3_MuDep(double x, double m2Q2, double m2mu2, int nf) const ;

        double ZeroFunction(double /*x*/, double /*m2Q2*/, double /*m2mu2*/, int /*nf*/) const {return 0.;};
        double WarningFunc(double /*x*/, double /*m2Q2*/, int /*nf*/) const ;

};

#endif
