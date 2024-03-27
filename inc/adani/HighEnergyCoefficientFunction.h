/*
 * =====================================================================================
 *
 *       Filename:  HighEnergyCoefficientFunction.h
 *
 *    Description:  Header file for the
 * HighEnergyCoefficientFunction.cc file.
 *
 *         Author:  L'allenamento si fa!
 *
 *  In this file there are the coefficient functions in the
 * high energy limit, i.e. x -> 0
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
//      small x limit (or the opposite) Power terms: power
//      terms in the small x liimit, obtained via
//      C_powerterms = C_highenergy - C_highenergy_highscale
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
        ~AbstractHighEnergyCoefficientFunction() override{};

        // get methods
        bool GetNLL() const { return NLL_; };

        // set methods
        void SetNLL(const bool &NLL) { NLL_ = NLL; };

    private:
        bool NLL_;

        //==========================================================================================//
        //                  Color factors O(as^3)
        //------------------------------------------------------------------------------------------//

    protected:
        double a_10(int nf) const;
        double a_11() const;
        double a_21(int nf) const;
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
        ~HighEnergyCoefficientFunction() override{};

        Value
        fxBand(double x, double m2Q2, double m2mu2, int nf) const override;

        void SetFunctions();

    private:
        double (HighEnergyCoefficientFunction::*LL_)(
            double, double, double
        ) const;
        Value (HighEnergyCoefficientFunction::*NLL_)(
            double, double, double, int
        ) const;

        //==========================================================================================//
        //                      High energy coefficient functions
        //                      O(as^2)
        //------------------------------------------------------------------------------------------//

        double C2_g2_highenergy(double x, double m2Q2, double m2mu2) const;
        double C2_ps2_highenergy(double x, double m2Q2, double m2mu2) const;
        double CL_g2_highenergy(double x, double m2Q2, double m2mu2) const;
        double CL_ps2_highenergy(double x, double m2Q2, double m2mu2) const;

        //==========================================================================================//
        //                      High energy coefficient functions
        //                      O(as^3) at leading log
        //------------------------------------------------------------------------------------------//

        double C2_g3_highenergyLL(double x, double m2Q2, double m2mu2) const;
        double CL_g3_highenergyLL(double x, double m2Q2, double m2mu2) const;
        double C2_ps3_highenergyLL(double x, double m2Q2, double m2mu2) const;
        double CL_ps3_highenergyLL(double x, double m2Q2, double m2mu2) const;

        Value
        C2_g3_highenergyNLL(double x, double m2Q2, double m2mu2, int nf) const;
        Value
        CL_g3_highenergyNLL(double x, double m2Q2, double m2mu2, int nf) const;
        Value
        C2_ps3_highenergyNLL(double x, double m2Q2, double m2mu2, int nf) const;
        Value
        CL_ps3_highenergyNLL(double x, double m2Q2, double m2mu2, int nf) const;

        //==========================================================================================//
        //                      High energy coefficient functions
        //                      O(as^3)
        //------------------------------------------------------------------------------------------//

        Value
        C2_g3_highenergy(double x, double m2Q2, double m2mu2, int nf) const;
        Value
        C2_ps3_highenergy(double x, double m2Q2, double m2mu2, int nf) const;
        Value
        CL_g3_highenergy(double x, double m2Q2, double m2mu2, int nf) const;
        Value
        CL_ps3_highenergy(double x, double m2Q2, double m2mu2, int nf) const;

        //==========================================================================================//
        //  Function needed to make LL_ and NLL_ point to a zero function
        //------------------------------------------------------------------------------------------//

        double
        ZeroFunction(double /*x*/, double /*m2Q2*/, double /*m2mu2*/) const {
            return 0.;
        };
        Value ZeroFunctionBand(
            double /*x*/, double /*m2Q2*/, double /*m2mu2*/, int /*nf*/
        ) const {
            return Value(0.);
        };
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
        ~HighEnergyHighScaleCoefficientFunction() override{};

        Value
        fxBand(double x, double m2Q2, double m2mu2, int nf) const override;

        void SetFunctions();

    private:
        double (HighEnergyHighScaleCoefficientFunction::*LL_)(
            double, double, double
        ) const;
        Value (HighEnergyHighScaleCoefficientFunction::*NLL_)(
            double, double, double, int
        ) const;

        //==========================================================================================//
        //                      Q>>m limit of the high energy
        //                      coefficient functions O(as^2)
        //------------------------------------------------------------------------------------------//

        double
        C2_g2_highenergy_highscale(double x, double m2Q2, double m2mu2) const;
        double
        C2_ps2_highenergy_highscale(double x, double m2Q2, double m2mu2) const;
        double
        CL_g2_highenergy_highscale(double x, double m2Q2, double m2mu2) const;
        double
        CL_ps2_highenergy_highscale(double x, double m2Q2, double m2mu2) const;

        //==========================================================================================//
        //  Q>>m limit of the high energy coefficient functions
        //  O(as^3) at leading log
        //------------------------------------------------------------------------------------------//

        double
        C2_g3_highenergy_highscaleLL(double x, double m2Q2, double m2mu2) const;
        double
        CL_g3_highenergy_highscaleLL(double x, double m2Q2, double m2mu2) const;
        double C2_ps3_highenergy_highscaleLL(
            double x, double m2Q2, double m2mu2
        ) const;
        double CL_ps3_highenergy_highscaleLL(
            double x, double m2Q2, double m2mu2
        ) const;

        Value C2_g3_highenergy_highscaleNLL(
            double x, double m2Q2, double m2mu2, int nf
        ) const;
        Value C2_ps3_highenergy_highscaleNLL(
            double x, double m2Q2, double m2mu2, int nf
        ) const;
        Value CL_g3_highenergy_highscaleNLL(
            double x, double m2Q2, double m2mu2, int nf
        ) const;
        Value CL_ps3_highenergy_highscaleNLL(
            double x, double m2Q2, double m2mu2, int nf
        ) const;

        //==========================================================================================//
        //                  Q^2>>m^2 limit of the high energy
        //                  coefficient functions O(as^3)
        //------------------------------------------------------------------------------------------//

        Value C2_g3_highenergy_highscale(
            double x, double m2Q2, double m2mu2, int nf
        ) const;
        Value C2_ps3_highenergy_highscale(
            double x, double m2Q2, double m2mu2, int nf
        ) const;
        Value CL_g3_highenergy_highscale(
            double x, double m2Q2, double m2mu2, int nf
        ) const;
        Value CL_ps3_highenergy_highscale(
            double x, double m2Q2, double m2mu2, int nf
        ) const;

        //==========================================================================================//
        //  Function needed to make LL_ and NLL_ point to a zero function
        //------------------------------------------------------------------------------------------//

        double
        ZeroFunction(double /*x*/, double /*m2Q2*/, double /*m2mu2*/) const {
            return 0.;
        };
        Value ZeroFunctionBand(
            double /*x*/, double /*m2Q2*/, double /*m2mu2*/, int /*nf*/
        ) const {
            return Value(0.);
        };
};

//==========================================================================================//
//  class PowerTermsCoefficientFunction
//------------------------------------------------------------------------------------------//

class PowerTermsCoefficientFunction
    : public AbstractHighEnergyCoefficientFunction {

    public:
        PowerTermsCoefficientFunction(
            const int &order, const char &kind, const char &channel,
            const bool &NLL = true
        );
        ~PowerTermsCoefficientFunction() override;

        Value
        fxBand(double x, double m2Q2, double m2mu2, int nf) const override;

    private:
        HighEnergyCoefficientFunction *highenergy_;
        HighEnergyHighScaleCoefficientFunction *highenergyhighscale_;
};

#endif
