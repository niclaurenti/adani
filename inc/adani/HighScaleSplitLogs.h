/*
 * =====================================================================================
 *
 *       Filename:  HighScaleSplitLogs.h
 *
 *         Author:  Daniele Adani
 *
 *    Description:  I mancini sono veramente l'essenza della qualita' nel calcio
 *
 *  In this file there are the coefficient functions in the
 *  high scale limit, i.e. Q^2 >> m^2, implemented setting mu=Q and splitting
 *  the coefficients of each logarithm
 *
 * =====================================================================================
 */

#ifndef HighScaleLogs_h
#define HighScaleLogs_h

#include "adani/CoefficientFunction.h"
#include "adani/MasslessCoefficientFunction.h"
#include "adani/MatchingCondition.h"

//==========================================================================================//
//           Legend:
//           at order O(alpha_s^n)
//           N^(k)LL is the coefficient of log(m^2/Q^2)^(n-k)
//           (The expressions do not contain the log)
//------------------------------------------------------------------------------------------//

//==========================================================================================//
//  class HighScaleSplitLogs
//------------------------------------------------------------------------------------------//

class HighScaleSplitLogs : public CoefficientFunction {
    public:
        HighScaleSplitLogs(
            const int &order, const char &kind, const char &channel,
            const string &version = "klmv"
        );
        ~HighScaleSplitLogs() override;

        double fx(double x, double m2Q2, int nf) const;
        Value fxBand(double x, double m2Q2, int nf) const;

        [[deprecated(
            "This function is deprecated and should not be used since "
            "HighScaleSplitLogs is implemented only for Q=mu. It will throw a "
            "NotImplementedException. Use HighScaleSplitLogs::fx(double x, "
            "double m2Q2, int nf) instead"
        )]] double
            fx(double /*x*/, double /*m2Q2*/, double /*m2mu2*/,
               int /*nf*/) const override;
        [[deprecated(
            "This function is deprecated and should not be used since "
            "HighScaleSplitLogs is implemented only for Q=mu. It will throw a "
            "NotImplementedException. Use HighScaleSplitLogs::fxBand(double x, "
            "double m2Q2, int nf) instead"
        )]] Value
            fxBand(
                double x, double m2Q2, double /*m2mu2*/, int nf
            ) const override;

        // division of the total result in the log terms
        double LL(double x, int nf) const { return (this->*LL_)(x, nf); }
        double NLL(double x, int nf) const { return (this->*NLL_)(x, nf); }
        double N2LL(double x, int nf) const { return (this->*N2LL_)(x, nf); }
        Value N3LL(double x, int nf) const { return (this->*N3LL_)(x, nf); }

    private:
        MasslessCoefficientFunction *massless_as1_;
        MasslessCoefficientFunction *massless_;
        MatchingCondition *a_muindep_;

        double (HighScaleSplitLogs::*LL_)(double, int) const;
        double (HighScaleSplitLogs::*NLL_)(double, int) const;
        double (HighScaleSplitLogs::*N2LL_)(double, int) const;
        Value (HighScaleSplitLogs::*N3LL_)(double, int) const;

        void SetFunctions();

        //==========================================================================================//
        //           Gluon channel, F2
        //------------------------------------------------------------------------------------------//

        double C2_g3_highscale_LL(double x, int nf) const;
        double C2_g3_highscale_NLL(double x, int nf) const;
        double C2_g3_highscale_N2LL(double x, int nf) const;
        Value C2_g3_highscale_N3LL(double x, int nf) const;

        //==========================================================================================//
        //           Pure singlet channel, F2
        //------------------------------------------------------------------------------------------//

        double C2_ps3_highscale_LL(double x, int nf) const;
        double C2_ps3_highscale_NLL(double x, int nf) const;
        double C2_ps3_highscale_N2LL(double x, int nf) const;
        Value C2_ps3_highscale_N3LL(double x, int nf) const;

        //==========================================================================================//
        //           Gluon channel, FL
        //------------------------------------------------------------------------------------------//

        double CL_g3_highscale_NLL(double x, int /*nf*/) const;
        double CL_g3_highscale_N2LL(double x, int nf) const;
        Value CL_g3_highscale_N3LL(double x, int nf) const;

        //==========================================================================================//
        //           Pure singlet channel, FL
        //------------------------------------------------------------------------------------------//

        double CL_ps3_highscale_NLL(double x, int /*nf*/) const;
        double CL_ps3_highscale_N2LL(double x, int /*nf*/) const;
        Value CL_ps3_highscale_N3LL(double x, int nf) const;

        double ZeroFunction(double /*x*/, int /*nf*/) const { return 0.; };
};

#endif
