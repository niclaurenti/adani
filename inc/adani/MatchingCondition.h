/*
 * =====================================================================================
 *
 *       Filename:  MatchingCondition.h
 *
 *    Description:  Header file for the
 * MatchingCondition.cc file.
 *
 *         Author:  Non credo esista una persona che analizzi il calcio meglio
 * di me.
 *
 *  In this file there are the matching conditions.
 *
 * =====================================================================================
 */

#ifndef Match_h
#define Match_h

#include "CoefficientFunction.h"

#include <string>

using std::string;

//==========================================================================================//
//  class MatchingCondition
//------------------------------------------------------------------------------------------//

class MatchingCondition {
    public:
        MatchingCondition(
            const int &order, const char &entry1, const char &entry2,
            const string &version
        );
        ~MatchingCondition(){};

        int GetOrder() const { return order_; };
        char GetEntry1() const { return entry1_; };
        char GetEntry2() const { return entry2_; };
        string GetVersion() const { return version_; };

        Value MuIndependentNfIndependentTerm(double x) const;
        vector<double> NotOrdered(double x) const;

    private:
        const int order_;
        const char entry1_;
        const char entry2_;
        const string version_;

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
