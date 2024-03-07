/*
 * =====================================================================================
 *
 *       Filename:  HighEnergyCoefficientFunctions.h
 *
 *    Description:  Header file for the
 * HighEnergyCoefficientFunctions.cc file.
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
//          C_powerterms = C_highenergy -
//          C_highenergy_highscale
//------------------------------------------------------------------------------------------//

class HighTmpCoefficientFunction : public CoefficientFunction {
    public:
        HighTmpCoefficientFunction(const int order, const char kind, const char channel, const bool NLL = true) ;
        HighTmpCoefficientFunction() : HighTmpCoefficientFunction(1, '2', 'g', true) {} ;
        ~HighTmpCoefficientFunction() {};

        // get methods
        bool GetNLL() {return NLL_;};

        // set methods
        void SetNLL(const bool NLL) {NLL_ = NLL ; } ;

    private:
        bool NLL_ ;

        //==========================================================================================//
        //                  Color factors O(alpha_s^3)
        //------------------------------------------------------------------------------------------//

    protected:
        double a_10(int nf);
        double a_11();
        double a_21(int nf);
};

class HighEnergyCoefficientFunction : public HighTmpCoefficientFunction {

    public:
        HighEnergyCoefficientFunction(const int order, const char kind, const char channel, const bool NLL = true) : HighTmpCoefficientFunction(order, kind, channel, NLL) {};
        HighEnergyCoefficientFunction() : HighTmpCoefficientFunction() {} ;
        ~HighEnergyCoefficientFunction() {} ;

        double fx(double x, double m2Q2, double m2mu2, int nf) ;

    private:

        //==========================================================================================//
        //                      High energy coefficient functions
        //                      O(alpha_s^2)
        //------------------------------------------------------------------------------------------//

        double C2_g2_highenergy(double x, double m2Q2, double m2mu2);
        double C2_ps2_highenergy(double x, double m2Q2, double m2mu2);
        double CL_g2_highenergy(double x, double m2Q2, double m2mu2);
        double CL_ps2_highenergy(double x, double m2Q2, double m2mu2);

        //==========================================================================================//
        //                      High energy coefficient functions
        //                      O(alpha_s^3) at leading log
        //------------------------------------------------------------------------------------------//

        double C2_g3_highenergyLL(double x, double m2Q2, double m2mu2);
        double CL_g3_highenergyLL(double x, double m2Q2, double m2mu2);
        double C2_ps3_highenergyLL(double x, double m2Q2, double m2mu2);
        double CL_ps3_highenergyLL(double x, double m2Q2, double m2mu2);

        double C2_g3_highenergyNLL(double x, double m2Q2, double m2mu2, int nf, int v);
        double CL_g3_highenergyNLL(double x, double m2Q2, double m2mu2, int nf, int v);

        //==========================================================================================//
        //                      High energy coefficient functions
        //                      O(alpha_s^3)
        //------------------------------------------------------------------------------------------//

        double C2_g3_highenergy(double x, double m2Q2, double m2mu2, int nf, int v);
        double C2_ps3_highenergy(double x, double m2Q2, double m2mu2, int nf, int v);
        double CL_g3_highenergy(double x, double m2Q2, double m2mu2, int nf, int v);
        double CL_ps3_highenergy(double x, double m2Q2, double m2mu2, int nf, int v);

};

class HighEnergyHighScaleCoefficientFunction : public HighTmpCoefficientFunction {

    public:
        HighEnergyHighScaleCoefficientFunction(const int order, const char kind, const char channel, const bool NLL = true) : HighTmpCoefficientFunction(order, kind, channel, NLL) {};
        HighEnergyHighScaleCoefficientFunction() : HighTmpCoefficientFunction() {} ;
        ~HighEnergyHighScaleCoefficientFunction() {} ;

        double fx(double x, double m2Q2, double m2mu2, int nf) override {return 0.;};

    private:

        //==========================================================================================//
        //                      Q>>m limit of the high energy
        //                      coefficient functions O(alpha_s^2)
        //------------------------------------------------------------------------------------------//

        double C2_g2_highenergy_highscale(double x, double m2Q2, double m2mu2);
        double C2_ps2_highenergy_highscale(double x, double m2Q2, double m2mu2);
        double CL_g2_highenergy_highscale(double x, double m2Q2, double m2mu2);
        double CL_ps2_highenergy_highscale(double x, double m2Q2, double m2mu2);

        //==========================================================================================//
        //  Q>>m limit of the high energy coefficient functions
        //  O(alpha_s^3) at leading log
        //------------------------------------------------------------------------------------------//

        double C2_g3_highenergy_highscaleLL(double x, double m2Q2, double m2mu2);
        double CL_g3_highenergy_highscaleLL(double x, double m2Q2, double m2mu2);
        double C2_ps3_highenergy_highscaleLL(double x, double m2Q2, double m2mu2);
        double CL_ps3_highenergy_highscaleLL(double x, double m2Q2, double m2mu2);

        double C2_g3_highenergy_highscaleNLL(
            double x, double m2Q2, double m2mu2, int nf, int v
        );
        // double C2_ps3_highenergy_highscaleNLL(double x, double
        // m2Q2, double m2mu2, int nf, int v);
        double CL_g3_highenergy_highscaleNLL(
            double x, double m2Q2, double m2mu2, int nf, int v
        );
        // double CL_ps3_highenergy_highscaleNLL(double x, double
        // m2Q2, double m2mu2, int nf, int v);

        //==========================================================================================//
        //                  Q^2>>m^2 limit of the high energy
        //                  coefficient functions O(alpha_s^3)
        //------------------------------------------------------------------------------------------//

        double
        C2_g3_highenergy_highscale(double x, double m2Q2, double m2mu2, int nf, int v);
        double
        C2_ps3_highenergy_highscale(double x, double m2Q2, double m2mu2, int nf, int v);
        double
        CL_g3_highenergy_highscale(double x, double m2Q2, double m2mu2, int nf, int v);
        double
        CL_ps3_highenergy_highscale(double x, double m2Q2, double m2mu2, int nf, int v);

};

class PowerTermsCoefficientFunction : public HighTmpCoefficientFunction {

    public:
        PowerTermsCoefficientFunction(const int order, const char kind, const char channel, const bool NLL = true);
        PowerTermsCoefficientFunction() : HighTmpCoefficientFunction() {} ;
        ~PowerTermsCoefficientFunction() {} ;

        double fx(double x, double m2Q2, double m2mu2, int nf);

    private:

        HighEnergyCoefficientFunction highenergy_ ;
        HighEnergyHighScaleCoefficientFunction highenergyhighscale_ ;

        // //==========================================================================================//
        // //                  Power terms of the coefficient function
        // //                  O(alpha_s^2)
        // //------------------------------------------------------------------------------------------//

        // double C2_g2_power_terms(double x, double m2Q2, double m2mu2);
        // double C2_ps2_power_terms(double x, double m2Q2, double m2mu2);
        // double CL_g2_power_terms(double x, double m2Q2, double m2mu2);
        // double CL_ps2_power_terms(double x, double m2Q2, double m2mu2);

        // //==========================================================================================//
        // //                  Power terms of the coefficient function
        // //                  O(alpha_s^3) at leading log
        // //------------------------------------------------------------------------------------------//

        // double C2_g3_power_termsLL(double x, double m2Q2 , double m2mu2);
        // double C2_ps3_power_termsLL(double x, double m2Q2, double m2mu2);

        // //==========================================================================================//
        // //                  Power terms of the coefficient function
        // //                  O(alpha_s^3)
        // //------------------------------------------------------------------------------------------//

        // double C2_g3_power_terms(double x, double m2Q2, double m2mu2, int nf, int v);
        // double C2_ps3_power_terms(double x, double m2Q2, double m2mu2, int nf, int v);
        // double CL_g3_power_terms(double x, double m2Q2, double m2mu2, int nf, int v);
        // double CL_ps3_power_terms(double x, double m2Q2, double m2mu2, int nf, int v);

};

#endif
