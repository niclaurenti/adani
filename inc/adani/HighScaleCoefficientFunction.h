/*
 * =====================================================================================
 *
 *       Filename:  HighScaleCoefficientFunction.h
 *
 *    Description:  Header file for the
 * HighScaleCoefficientFunction.cc file.
 *
 *         Author:  Da quella posizione trasforma l'acqua in
 * vino
 *
 *  In this file there are the coefficient functions in the
 * high scale limit, i.e. Q^2 >> m^2
 *
 * =====================================================================================
 */

#ifndef HighScale_h
#define HighScale_h

#include "adani/CoefficientFunction.h"
#include "adani/MasslessCoefficientFunction.h"
#include "adani/MatchingCondition.h"

//==========================================================================================//
//                      The notation used is the following:
//                      C^{(k)} = k-th order expansion in
//                      terms of \alpha_s^{[nf]} D^{(k)} =
//                      k-th order expansion in terms of
//                      \alpha_s^{[nf+1]}
//------------------------------------------------------------------------------------------//

//==========================================================================================//
//  class HighScaleCoefficientFunction
//------------------------------------------------------------------------------------------//

class HighScaleCoefficientFunction : public CoefficientFunction {
    public:
        HighScaleCoefficientFunction(
            const int &order, const char &kind, const char &channel,
            const string &version = "klmv"
        );
        ~HighScaleCoefficientFunction() override;

        double fx(double x, double m2Q2, double m2mu2, int nf) const override;

        Value
        fxBand(double x, double m2Q2, double m2mu2, int nf) const override;

        vector<double>
        fxBand_NotOrdered(double x, double m2Q2, double m2mu2, int nf) const;

        void SetFunctions();

    private:
        Value (HighScaleCoefficientFunction::*fx_)(
            double, double, double, int
        ) const;

        MasslessCoefficientFunction *massless_as1_;
        MasslessCoefficientFunction *massless_as2_;
        MasslessCoefficientFunction *massless_as3_;
        MatchingCondition *a_muindep_;

        //==========================================================================================//
        //                      High scale (Q^2 >> m^2) coefficient
        //                      functions O(as)
        //------------------------------------------------------------------------------------------//

        Value
        C2_g1_highscale(double x, double m2Q2, double /*m2mu2*/, int /*nf*/)
            const;
        Value
        CL_g1_highscale(double x, double /*m2Q2*/, double /*m2mu2*/, int /*nf*/)
            const;

        double D2_g1_highscale(double x, double m2Q2) const;
        double DL_g1_highscale(double x) const;

        //==========================================================================================//
        //                      High scale (Q^2 >> m^2) coefficient
        //                      functions O(as^2)
        //------------------------------------------------------------------------------------------//

        Value
        C2_g2_highscale(double x, double m2Q2, double m2mu2, int /*nf*/) const;
        Value
        C2_ps2_highscale(double z, double m2Q2, double m2mu2, int /*nf*/) const;

        Value
        CL_g2_highscale(double x, double m2Q2, double m2mu2, int /*nf*/) const;
        Value
        CL_ps2_highscale(double x, double m2Q2, double m2mu2, int /*nf*/) const;

        double D2_g2_highscale(double x, double m2Q2, double m2mu2) const;
        double D2_ps2_highscale(double x, double m2Q2, double m2mu2) const;

        double DL_g2_highscale(double z, double m2Q2, double m2mu2) const;
        double DL_ps2_highscale(double z, double m2Q2, double m2mu2) const;

        //==========================================================================================//
        //                      High scale (Q^2 >> m^2) coefficient
        //                      functions O(as^3)
        //------------------------------------------------------------------------------------------//

        Value
        C2_g3_highscale(double x, double m2Q2, double m2mu2, int nf) const;
        Value
        C2_ps3_highscale(double x, double m2Q2, double m2mu2, int nf) const;

        Value
        CL_g3_highscale(double x, double m2Q2, double m2mu2, int nf) const;
        Value
        CL_ps3_highscale(double x, double m2Q2, double m2mu2, int nf) const;

        double
        DL_g3_highscale(double z, double m2Q2, double m2mu2, int nf) const;
        double
        DL_ps3_highscale(double z, double m2Q2, double m2mu2, int nf) const;

        Value
        D2_g3_highscale(double x, double m2Q2, double m2mu2, int nf) const;
        Value
        D2_ps3_highscale(double x, double m2Q2, double m2mu2, int nf) const;
};

#endif
