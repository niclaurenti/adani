/*
 * =====================================================================================
 *
 *       Filename:  MatchingCondition.h
 *
 *         Author:  Daniele Adani
 *
 *    Description:  Non credo esista una persona che analizzi il calcio meglio
 *                  di me.
 *
 *  In this file there are the matching conditions.
 *
 * =====================================================================================
 */

#ifndef Match_h
#define Match_h

#include "CoefficientFunction.h"
#include "adani/CommonTypes.h"

//==========================================================================================//
//  class MatchingCondition
//------------------------------------------------------------------------------------------//

class MatchingCondition {
    public:
        MatchingCondition(
            const int &order, const char &entry1, const char &entry2,
            const HighScaleVersion &version = HighScaleVersion::Exact
        );
        MatchingCondition(const MatchingCondition &obj);
        ~MatchingCondition() = default;

        // copy operator is deleted since the possibility to change
        // the state of an instance of MatchingCondition is not implemented yet
        // For this reason the data members are declared const
        // TODO: implement it!!!
        MatchingCondition &operator=(MatchingCondition &rhs) = delete;

        int GetOrder() const { return order_; };
        char GetEntry1() const { return entry1_; };
        char GetEntry2() const { return entry2_; };
        HighScaleVersion GetHighScaleVersion() const { return version_; };

        Value MuIndependentNfIndependentTerm(double x) const;
        double MuIndependentNfDependentTerm(double x) const;
        Value MuIndependentTerm(double x, int nf) const;
        vector<double> NotOrdered(double x) const;

    private:
        const int order_;
        const char entry1_;
        const char entry2_;
        const HighScaleVersion version_;

        void CheckEntry(char entry) const;

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
        double a_Qg_31(double x) const;
        double a_Qq_PS_30(double x, int v) const;
        double a_Qq_PS_31(double x) const;
};

#endif
