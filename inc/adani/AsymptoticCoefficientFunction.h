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
#include "adani/CommonTypes.h"
#include "adani/HighEnergyCoefficientFunction.h"
#include "adani/HighScaleCoefficientFunction.h"

#include <memory>

//==========================================================================================//
//  class AsymptoticCoefficientFunction
//------------------------------------------------------------------------------------------//

class AsymptoticCoefficientFunction : public CoefficientFunction {
    public:
        AsymptoticCoefficientFunction(
            const int &order, const char &kind, const char &channel,
            const bool &NLL = true,
            const HighScaleVersion &highscale_version = HighScaleVersion::Exact
        );
        AsymptoticCoefficientFunction(const AsymptoticCoefficientFunction &obj);
        ~AsymptoticCoefficientFunction() override = default;

        bool GetNLL() const { return highenergy_->GetNLL(); };
        HighScaleVersion GetHighScaleVersion() const {
            return highscale_->GetHighScaleVersion();
        };

        bool IsLegacyPowerTerms() const { return legacy_pt_; };
        void SetLegacyPowerTerms(const bool &legacy_pt);

        Value
            fxBand(double x, double m2Q2, double m2mu2, int nf) const override;

    private:
        bool legacy_pt_;
        double a_fact_;
        std::unique_ptr<HighScaleCoefficientFunction> highscale_;
        std::unique_ptr<HighEnergyCoefficientFunction> highenergy_;
        std::unique_ptr<HighEnergyHighScaleCoefficientFunction>
            highenergyhighscale_;

        Value (AsymptoticCoefficientFunction::*fx_)(
            double, double, double, int
        ) const;

        void SetFunctions();

        Value PlainAdditiveMatching(
            double x, double m2Q2, double m2mu2, int nf
        ) const;
        Value PlainMultiplicativeMatching(
            double x, double m2Q2, double m2mu2, int nf
        ) const;
        Value ModifiedMultiplicativeMatching1(
            double x, double m2Q2, double m2mu2, int nf
        ) const;
        Value ModifiedMultiplicativeMatching2(
            double x, double m2Q2, double m2mu2, int nf
        ) const;
        Value
            C2_2_asymptotic(double x, double m2Q2, double m2mu2, int nf) const;
        Value
            CL_2_asymptotic(double x, double m2Q2, double m2mu2, int nf) const;
        Value
            C2_3_asymptotic(double x, double m2Q2, double m2mu2, int nf) const;
        Value
            CL_3_asymptotic(double x, double m2Q2, double m2mu2, int nf) const;
        Value Delta2(
            Value central, Value variation
        ) const;
        Value Delta3(
            Value central, Value variation1, Value variation2
        ) const;
        double C_highenergy_lim(
            double highenergy_ll, double highscalehighenergy_ll, double a_fact
        ) const;
        double ComputeDampDelta(double m2Q2, double m2mu2) const;
};

#endif
