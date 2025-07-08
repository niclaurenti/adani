/*
 * =====================================================================================
 *
 *       Filename:  AsymptoticCoefficientFunction.h
 *
 *         Author:  Daniele Adani
 *
 *    Description:  Anche tu sei un teorico Max
 *
 *  In this file there are the coefficient functions in the
 *  asymptotic limit, i.e. x -> 0 and/or Q^2 >> m^2
 *
 * =====================================================================================
 */

#ifndef Asymptotic_h
#define Asymptotic_h

#include "adani/CoefficientFunction.h"
#include "adani/HighEnergyCoefficientFunction.h"
#include "adani/HighScaleCoefficientFunction.h"
#include "adani/PowerTermsCoefficientFunction.h"

//==========================================================================================//
//  class AsymptoticCoefficientFunction
//------------------------------------------------------------------------------------------//

class AsymptoticCoefficientFunction : public CoefficientFunction {
    public:
        AsymptoticCoefficientFunction(
            const int &order, const char &kind, const char &channel,
            const bool &NLL = true, const string &highscale_version = "exact"
        );
        ~AsymptoticCoefficientFunction();

        void SetLegacyPowerTerms(const bool &legacy_pt);

        Value
            fxBand(double x, double m2Q2, double m2mu2, int nf) const override;

        vector<double>
            AllVariations(double x, double m2Q2, double m2mu2, int nf) const;

    private:
        const bool NLL_;
        HighScaleCoefficientFunction *highscale_;
        AbstractPowerTerms *powerterms_;

        Value (AsymptoticCoefficientFunction::*fx_)(
            double, double, double, int
        ) const;

        Value
            AdditiveMatching(double x, double m2Q2, double m2mu2, int nf) const;
        Value MultiplicativeMatching(
            double x, double m2Q2, double m2mu2, int nf
        ) const;
};

#endif
