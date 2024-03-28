/*
 * =====================================================================================
 *
 *       Filename:  AsymptoticCoefficientFunction.h
 *
 *    Description:  Header file for the
 * AsymptoticCoefficientFunction.cc file.
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
#include "adani/HighEnergyCoefficientFunction.h"
#include "adani/HighScaleCoefficientFunction.h"

//==========================================================================================//
//  class AsymptoticCoefficientFunction
//------------------------------------------------------------------------------------------//

class AsymptoticCoefficientFunction
    : public AbstractHighEnergyCoefficientFunction {
    public:
        AsymptoticCoefficientFunction(
            const int &order, const char &kind, const char &channel,
            const bool &NLL = true, const string &highscale_version = "klmv"
        );
        ~AsymptoticCoefficientFunction();

        double fx(double x, double m2Q2, double m2mu2, int nf) const override;

        Value
        fxBand(double x, double m2Q2, double m2mu2, int nf) const override;

        vector<double>
        AllVariations(double x, double m2Q2, double m2mu2, int nf) const;

    private:
        HighScaleCoefficientFunction *highscale_;
        PowerTermsCoefficientFunction *powerterms_;
};

#endif
