/*
 * =====================================================================================
 *
 *       Filename:  PowerTermsCoefficientFunction.h
 *
 *         Author:  Daniele Adani
 *
 *    Description:  L'allenamento si fa!
 *
 *  In this file there are the coefficient functions in the
 *  high energy limit, i.e. x -> 0
 *
 * =====================================================================================
 */

#ifndef PowerTerms_h
#define PowerTerms_h

#include "adani/HighEnergyCoefficientFunction.h"

//==========================================================================================//
//  class AbstractPowerTerms
//------------------------------------------------------------------------------------------//

class AbstractPowerTerms : public CoefficientFunction {

    public:
        AbstractPowerTerms(
            const int &order, const char &kind, const char &channel,
            const bool &NLL = true
        );
        ~AbstractPowerTerms() override;

        HighEnergyCoefficientFunction* GetHighEnergy() const {return highenergy_;};
        HighEnergyHighScaleCoefficientFunction* GetHighEnergyHighScale() const {return highenergyhighscale_;};

    private:
        HighEnergyCoefficientFunction *highenergy_;
        HighEnergyHighScaleCoefficientFunction *highenergyhighscale_;
};


//==========================================================================================//
//  class PowerTermsCoefficientFunction
//------------------------------------------------------------------------------------------//

class PowerTermsCoefficientFunction
    : public AbstractPowerTerms {

    public:
        PowerTermsCoefficientFunction(
            const int &order, const char &kind, const char &channel,
            const bool &NLL = true
        ) : AbstractPowerTerms(order, kind, channel, NLL) {} ;
        ~PowerTermsCoefficientFunction() override {};

        Value
            fxBand(double x, double m2Q2, double m2mu2, int nf) const override;
};

//==========================================================================================//
//  class MultiplicativeAsymptotic
//------------------------------------------------------------------------------------------//

class MultiplicativeAsymptotic
    : public AbstractPowerTerms {

    public:
        MultiplicativeAsymptotic(
            const int &order, const char &kind, const char &channel,
            const bool &NLL = true
        );
        ~MultiplicativeAsymptotic() override {};

        Value
            fxBand(double x, double m2Q2, double m2mu2, int nf) const override;

    private:
        Value (MultiplicativeAsymptotic::*fx_)(
            double, double, double, int
        ) const;

        Value PlainRatio(double x, double m2Q2, double m2mu2, int nf) const;
        Value RegoularizedRatio(double x, double m2Q2, double m2mu2, int nf) const;

        Value OneFunction(double /*x*/, double /*m2Q2*/, double /*m2mu2*/, int /*nf*/) const {
            return Value(1.);
        }

};


#endif
