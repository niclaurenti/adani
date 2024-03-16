/*
 * =====================================================================================
 *
 *       Filename:  AsymptoticCoefficientFunctions.h
 *
 *    Description:  Header file for the
 * AsymptoticCoefficientFunctions.cc file.
 *
 *         Author:  Anche tu sei un teorico Max
 *
 *  In this file there are the coefficient functions in the
 * asymptotic limit, i.e. x -> 0 and/or Q^2 >> m^2
 *
 * =====================================================================================
 */

#ifndef Asymptotic_h
#define Asymptotic_h

#include "adani/CoefficientFunction.h"
#include "adani/HighEnergyCoefficientFunctions.h"
#include "adani/HighScaleCoefficientFunctions.h"

class AsymptoticCoefficientFunction : public AbstractHighEnergyCoefficientFunction {
    public:
        AsymptoticCoefficientFunction(const int& order, const char& kind, const char& channel, const bool& NLL = true, const bool& exact_highscale = false, const bool& revised_approx_highscale = true);
        ~AsymptoticCoefficientFunction() ;

        double fx(double x, double m2Q2, double m2mu2, int nf) const override;
        // double MuIndependentTerms(double x, double m2Q2, int nf) const override ;
        // double MuDependentTerms(double x, double m2Q2, double m2mu2, int nf) const override ;

        Value fxBand(double x, double m2Q2, double m2mu2, int nf) const override;

    private:
        HighScaleCoefficientFunction* highscale_ ;
        PowerTermsCoefficientFunction* powerterms_ ;

};

#endif
