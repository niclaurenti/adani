/*
 * =====================================================================================
 *
 *       Filename:  ThresholdCoefficientFunctions.h
 *
 *    Description:  Header file for the
 * ThresholdCoefficientFunctions.cc file.
 *
 *         Author:  La garra charruaaaaaaa
 *
 *  In this file there are the coefficient functions in the
 * threshold limit, i.e. s->4m^2 (where s is the partonic
 * conter of mass energy), or x -> xmax
 *
 * =====================================================================================
 */

#ifndef Threshold_h
#define Threshold_h

#include "adani/CoefficientFunction.h"
#include "adani/ExactCoefficientFunctions.h"

class ThresholdCoefficientFunction : public CoefficientFunction {

    public:
        ThresholdCoefficientFunction(const int& order, const char& kind, const char& channel) ;
        ~ThresholdCoefficientFunction() ;

        double fx(const double x, const double m2Q2, const double m2mu2, const int nf) const ;

    private:

        ExactCoefficientFunction *exactlo_ ;

        //==========================================================================================//
        //                      Threshold (s -> 4m^2) coefficient
        //                      functions O(alpha_s)
        //------------------------------------------------------------------------------------------//

        double C2_g1_threshold(const double x, const double m2Q2) const;

        //==========================================================================================//
        //                      Threshold (s -> 4m^2) coefficient
        //                      functions O(alpha_s^2)
        //------------------------------------------------------------------------------------------//

        double C2_g2_threshold(const double x, const double m2Q2, const double m2mu2) const;
        double CL_g2_threshold(const double x, const double m2Q2, const double m2mu2) const;

        double threshold_expansion_g2(const double x, const double m2Q2, const double m2mu2) const;
        double threshold_expansion_g2_const(const double m2Q2, const double m2mu2) const;

        double C2_g2_threshold_const(const double x, const double m2Q2, const double m2mu2) const;
        double CL_g2_threshold_const(const double x, const double m2Q2, const double m2mu2) const;

        //==========================================================================================//
        //                      Threshold (s -> 4m^2) coefficient
        //                      functions O(alpha_s^3)
        //------------------------------------------------------------------------------------------//

        double C2_g3_threshold(const double x, const double m2Q2, const double m2mu2, const int nf) const;
        double CL_g3_threshold(const double x, const double m2Q2, const double m2mu2, const int nf) const;

        double threshold_expansion_g3(const double x, const double m2Q2, const double m2mu2, const int nf) const;
        double threshold_expansion_g3_const(const double m2Q2, const double m2mu2) const;

        double C2_g3_threshold_const(const double x, const double m2Q2, const double m2mu2) const;
        double CL_g3_threshold_const(const double x, const double m2Q2, const double m2mu2) const;

        //==========================================================================================//
        //  Functions needed for the threshold limit.
        //------------------------------------------------------------------------------------------//

        double c0(const double xi) const;
        double c0_bar(const double xi) const;

};

#endif
