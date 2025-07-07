/*
 * =====================================================================================
 *
 *       Filename:  PowerTermsCoefficientFunction.h
 *
 *         Author:  Daniele Adani
 *
 *    Description:  Tutti in piedi! Tutti in piedi per il Brasile! Tutti a danzare
 *    con la Selecao, per onorare il calcio. Non c'è posto per i deboli di spirito:
 *    chi ama il calcio deve godere, gioire e, se può, ringraziare.
 *
 *  In this file there are the terms of the coefficient functions needed to match the
 *  high scale coefficient function to the high energy coefficient functions
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
