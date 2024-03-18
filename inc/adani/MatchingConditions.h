/*
 * =====================================================================================
 *
 *       Filename:  MatchingConditions.h
 *
 *    Description:  Header file for the
 * MatchingConditions.cc file.
 *
 *         Author:  Vamo!!
 *
 *  In this file there are the matching conditions.
 *
 * =====================================================================================
 */

#ifndef Match_h
#define Match_h

#include "CoefficientFunction.h"

//==========================================================================================//
//  Matching conditions O(alpha_s)
//------------------------------------------------------------------------------------------//

double K_Qg1(double x, double m2mu2);
double K_gg1_local(double m2mu2);

//==========================================================================================//
//  Matching conditions O(alpha_s^2)
//------------------------------------------------------------------------------------------//

double K_Qg2(double x, double m2mu2);

//==========================================================================================//
//  Matching conditions O(alpha_s^3)
//------------------------------------------------------------------------------------------//

class MatchingCondition {
    public:
        MatchingCondition(const int& order, const char& entry1, const char& entry2, const bool& exact, const bool& revised_approx) ;
        ~MatchingCondition() {};
        
        Value MuIndependentNfIndependentTerm(double x) const ;

    private:
        int order_;
        char entry1_;
        char entry2_;
        
        bool exact_;
        bool revised_approx_;

        double a_Qg_30(double x, int v) const ;
        double a_Qq_PS_30(double x, int v) const ;

};

#endif
