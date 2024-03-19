/*
 * =====================================================================================
 *
 *       Filename:  MatchingCondition.h
 *
 *    Description:  Header file for the
 * MatchingCondition.cc file.
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
//  class MatchingCondition
//------------------------------------------------------------------------------------------//

class MatchingCondition {
  public:
    MatchingCondition(
        const int &order, const char &entry1, const char &entry2,
        const bool &exact, const bool &revised_approx
    );
    ~MatchingCondition(){};

    Value MuIndependentNfIndependentTerm(double x) const;

  private:
    int order_;
    char entry1_;
    char entry2_;

    bool exact_;
    bool revised_approx_;

    //==========================================================================================//
    //  Matching conditions O(as)
    //------------------------------------------------------------------------------------------//

    // double K_Qg1(double x, double m2mu2) const;
    // double K_gg1_local(double m2mu2) const ;

    //==========================================================================================//
    //  Matching conditions O(as^2)
    //------------------------------------------------------------------------------------------//

    // double K_Qg2(double x, double m2mu2) const ;

    //==========================================================================================//
    //  Matching conditions O(as^3)
    //------------------------------------------------------------------------------------------//

    double a_Qg_30(double x, int v) const;
    double a_Qq_PS_30(double x, int v) const;
};

#endif
