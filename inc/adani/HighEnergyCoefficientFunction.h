/*
 * =====================================================================================
 *
 *       Filename:  HighEnergyCoefficientFunction.h
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

#ifndef HighEnergy_h
#define HighEnergy_h

#include "adani/CoefficientFunction.h"

//==========================================================================================//
//                      Notation:
//      High energy: small x limit
//      High energy high scale: Q^2 >> m^2 limit of the
//      small x limit (or the opposite)
//------------------------------------------------------------------------------------------//

//==========================================================================================//
//  class AbstractHighEnergyCoefficientFunction
//------------------------------------------------------------------------------------------//

class AbstractHighEnergyCoefficientFunction : public CoefficientFunction {
    public:
        AbstractHighEnergyCoefficientFunction(
            const int &order, const char &kind, const char &channel,
            const bool &NLL = true
        );
        ~AbstractHighEnergyCoefficientFunction() override = default;

        // get methods
        bool GetNLL() const { return NLL_; };

        virtual double LL(double m2Q2, double m2mu2) const = 0;
        virtual Value NLL(double m2Q2, double m2mu2, int nf) const = 0;

        Value
            fxBand(double x, double m2Q2, double m2mu2, int nf) const override;

    private:
        const bool NLL_;

        Value (AbstractHighEnergyCoefficientFunction::*fx_)(
            double, double, double, int
        ) const;

        Value Order2(double x, double m2Q2, double m2mu2, int nf) const;
        Value Order3(double x, double m2Q2, double m2mu2, int nf) const;
        Value Order3LL(double x, double m2Q2, double m2mu2, int nf) const;

        Value ZeroFunctionBand(
            double /*x*/, double /*m2Q2*/, double /*m2mu2*/, int /*nf*/
        ) const {
            return Value(0.);
        };

        //==========================================================================================//
        //                  Color factors O(as^3)
        //------------------------------------------------------------------------------------------//

    protected:
        double a_10(int nf) const;
        double a_11() const;
        double a_21(int nf) const;
        double a_21_new(int nf) const;
};

//==========================================================================================//
//  class HighEnergyCoefficientFunction
//------------------------------------------------------------------------------------------//

class HighEnergyCoefficientFunction
    : public AbstractHighEnergyCoefficientFunction {

    public:
        HighEnergyCoefficientFunction(
            const int &order, const char &kind, const char &channel,
            const bool &NLL = true
        );
        HighEnergyCoefficientFunction(const HighEnergyCoefficientFunction &obj);
        ~HighEnergyCoefficientFunction() override = default;

        double LL(double m2Q2, double m2mu2) const override;
        Value NLL(double m2Q2, double m2mu2, int nf) const override;

    private:
        double (HighEnergyCoefficientFunction::*LL_)(double, double) const;
        Value (HighEnergyCoefficientFunction::*NLL_)(double, double, int) const;

        void SetFunctions();

        //==========================================================================================//
        //                      High energy coefficient functions
        //                      O(as^2)
        //------------------------------------------------------------------------------------------//

        double C2_g2_highenergyLL(double m2Q2, double m2mu2) const;
        double C2_ps2_highenergyLL(double m2Q2, double m2mu2) const;
        double CL_g2_highenergyLL(double m2Q2, double m2mu2) const;
        double CL_ps2_highenergyLL(double m2Q2, double m2mu2) const;

        //==========================================================================================//
        //                      High energy coefficient functions
        //                      O(as^3) at leading log
        //------------------------------------------------------------------------------------------//

        double C2_g3_highenergyLL(double m2Q2, double m2mu2) const;
        double CL_g3_highenergyLL(double m2Q2, double m2mu2) const;
        double C2_ps3_highenergyLL(double m2Q2, double m2mu2) const;
        double CL_ps3_highenergyLL(double m2Q2, double m2mu2) const;

        //==========================================================================================//
        //                      High energy coefficient functions
        //                      O(as^3) at next-to-leading log
        //------------------------------------------------------------------------------------------//

        Value C2_g3_highenergyNLL(double m2Q2, double m2mu2, int nf) const;
        Value CL_g3_highenergyNLL(double m2Q2, double m2mu2, int nf) const;
        Value C2_ps3_highenergyNLL(double m2Q2, double m2mu2, int nf) const;
        Value CL_ps3_highenergyNLL(double m2Q2, double m2mu2, int nf) const;

        double C2_g3_highenergyNLL(
            double m2Q2, double m2mu2, double a11, double a10, double a21,
            double beta0
        ) const;
        double CL_g3_highenergyNLL(
            double m2Q2, double m2mu2, double a11, double a10, double a21,
            double beta0
        ) const;

        Value ThrowException(double /*m2Q2*/, double /*m2mu2*/, int /*nf*/) const;
};

//==========================================================================================//
//  class HighEnergyHighScaleCoefficientFunction
//------------------------------------------------------------------------------------------//

class HighEnergyHighScaleCoefficientFunction
    : public AbstractHighEnergyCoefficientFunction {

    public:
        HighEnergyHighScaleCoefficientFunction(
            const int &order, const char &kind, const char &channel,
            const bool &NLL = true
        );
        HighEnergyHighScaleCoefficientFunction(
            const HighEnergyHighScaleCoefficientFunction &obj
        );
        ~HighEnergyHighScaleCoefficientFunction() override = default;

        double LL(double m2Q2, double m2mu2) const override;
        Value NLL(double m2Q2, double m2mu2, int nf) const override;

    private:
        double (HighEnergyHighScaleCoefficientFunction::*LL_)(
            double, double
        ) const;
        Value (HighEnergyHighScaleCoefficientFunction::*NLL_)(
            double, double, int
        ) const;

        void SetFunctions();

        //==========================================================================================//
        //                      Q>>m limit of the high energy
        //                      coefficient functions O(as^2)
        //------------------------------------------------------------------------------------------//

        double C2_g2_highenergy_highscaleLL(double m2Q2, double m2mu2) const;
        double C2_ps2_highenergy_highscaleLL(double m2Q2, double m2mu2) const;
        double CL_g2_highenergy_highscaleLL(double m2Q2, double m2mu2) const;
        double CL_ps2_highenergy_highscaleLL(double m2Q2, double m2mu2) const;

        //==========================================================================================//
        //  Q>>m limit of the high energy coefficient functions
        //  O(as^3) at leading log
        //------------------------------------------------------------------------------------------//

        double C2_g3_highenergy_highscaleLL(double m2Q2, double m2mu2) const;
        double CL_g3_highenergy_highscaleLL(double m2Q2, double m2mu2) const;
        double C2_ps3_highenergy_highscaleLL(double m2Q2, double m2mu2) const;
        double CL_ps3_highenergy_highscaleLL(double m2Q2, double m2mu2) const;

        //==========================================================================================//
        //  Q>>m limit of the high energy coefficient functions
        //  O(as^3) at next-to-leading log
        //------------------------------------------------------------------------------------------//

        Value C2_g3_highenergy_highscaleNLL(
            double m2Q2, double m2mu2, int nf
        ) const;
        Value C2_ps3_highenergy_highscaleNLL(
            double m2Q2, double m2mu2, int nf
        ) const;
        Value CL_g3_highenergy_highscaleNLL(
            double m2Q2, double m2mu2, int nf
        ) const;
        Value CL_ps3_highenergy_highscaleNLL(
            double m2Q2, double m2mu2, int nf
        ) const;

        double C2_g3_highenergy_highscaleNLL(
            double m2Q2, double m2mu2, double a11, double a10, double a21,
            double beta0
        ) const;
        double CL_g3_highenergy_highscaleNLL(
            double m2Q2, double m2mu2, double a11, double a10, double a21,
            double beta0
        ) const;

        Value ThrowException(double /*m2Q2*/, double /*m2mu2*/, int /*nf*/) const;
};

#endif
