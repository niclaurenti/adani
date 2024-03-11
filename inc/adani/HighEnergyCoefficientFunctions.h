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

class AbstractHighEnergyCoefficientFunction : public CoefficientFunction {
    public:
        AbstractHighEnergyCoefficientFunction(const int& order, const char& kind, const char& channel, const bool& NLL = true) ;
        // AbstractHighEnergyCoefficientFunction() : AbstractHighEnergyCoefficientFunction(1, '2', 'g', true) {} ;
        ~AbstractHighEnergyCoefficientFunction() {};

        // get methods
        bool GetNLL() const {return NLL_;};

        // set methods
        void SetNLL(const bool& NLL) {NLL_ = NLL ; } ;

    private:
        bool NLL_ ;

    //==========================================================================================//
    //                  Color factors O(alpha_s^3)
    //------------------------------------------------------------------------------------------//

    protected:
        double a_10(const int nf) const;
        double a_11() const;
        double a_21(const int nf) const;
};

class HighEnergyCoefficientFunction : public AbstractHighEnergyCoefficientFunction {

    public:
        HighEnergyCoefficientFunction(const int& order, const char& kind, const char& channel, const bool& NLL = true) : AbstractHighEnergyCoefficientFunction(order, kind, channel, NLL) {};
        // HighEnergyCoefficientFunction() : AbstractHighEnergyCoefficientFunction() {} ;
        ~HighEnergyCoefficientFunction() {} ;

        double fx(const double x, const double m2Q2, const double m2mu2, const int nf) const ;

    private:

        //==========================================================================================//
        //                      High energy coefficient functions
        //                      O(alpha_s^2)
        //------------------------------------------------------------------------------------------//

        double C2_g2_highenergy(const double x, const double m2Q2, const double m2mu2) const ;
        double C2_ps2_highenergy(const double x, const double m2Q2, const double m2mu2) const ;
        double CL_g2_highenergy(const double x, const double m2Q2, const double m2mu2) const ;
        double CL_ps2_highenergy(const double x, const double m2Q2, const double m2mu2) const ;

        //==========================================================================================//
        //                      High energy coefficient functions
        //                      O(alpha_s^3) at leading log
        //------------------------------------------------------------------------------------------//

        double C2_g3_highenergyLL(const double x, const double m2Q2, const double m2mu2) const ;
        double CL_g3_highenergyLL(const double x, const double m2Q2, const double m2mu2) const ;
        double C2_ps3_highenergyLL(const double x, const double m2Q2, const double m2mu2) const ;
        double CL_ps3_highenergyLL(const double x, const double m2Q2, const double m2mu2) const ;

        double C2_g3_highenergyNLL(const double x, const double m2Q2, const double m2mu2, const int nf, const int v) const ;
        double CL_g3_highenergyNLL(const double x, const double m2Q2, const double m2mu2, const int nf, const int v) const ;

        //==========================================================================================//
        //                      High energy coefficient functions
        //                      O(alpha_s^3)
        //------------------------------------------------------------------------------------------//

        double C2_g3_highenergy(const double x, const double m2Q2, const double m2mu2, const int nf, const int v) const ;
        double C2_ps3_highenergy(const double x, const double m2Q2, const double m2mu2, const int nf, const int v) const ;
        double CL_g3_highenergy(const double x, const double m2Q2, const double m2mu2, const int nf, const int v) const ;
        double CL_ps3_highenergy(const double x, const double m2Q2, const double m2mu2, const int nf, const int v) const ;

};

class HighEnergyHighScaleCoefficientFunction : public AbstractHighEnergyCoefficientFunction {

    public:
        HighEnergyHighScaleCoefficientFunction(const int& order, const char& kind, const char& channel, const bool& NLL = true) : AbstractHighEnergyCoefficientFunction(order, kind, channel, NLL) {};
        // HighEnergyHighScaleCoefficientFunction() : AbstractHighEnergyCoefficientFunction() {} ;
        ~HighEnergyHighScaleCoefficientFunction() {} ;

        double fx(const double x, const double m2Q2, const double m2mu2, const int nf) const ;

    private:

        //==========================================================================================//
        //                      Q>>m limit of the high energy
        //                      coefficient functions O(alpha_s^2)
        //------------------------------------------------------------------------------------------//

        double C2_g2_highenergy_highscale(const double x, const double m2Q2, const double m2mu2) const ;
        double C2_ps2_highenergy_highscale(const double x, const double m2Q2, const double m2mu2) const ;
        double CL_g2_highenergy_highscale(const double x, const double m2Q2, const double m2mu2) const ;
        double CL_ps2_highenergy_highscale(const double x, const double m2Q2, const double m2mu2) const ;

        //==========================================================================================//
        //  Q>>m limit of the high energy coefficient functions
        //  O(alpha_s^3) at leading log
        //------------------------------------------------------------------------------------------//

        double C2_g3_highenergy_highscaleLL(const double x, const double m2Q2, const double m2mu2) const;
        double CL_g3_highenergy_highscaleLL(const double x, const double m2Q2, const double m2mu2) const;
        double C2_ps3_highenergy_highscaleLL(const double x, const double m2Q2, const double m2mu2) const;
        double CL_ps3_highenergy_highscaleLL(const double x, const double m2Q2, const double m2mu2) const;

        double C2_g3_highenergy_highscaleNLL(
            const double x, const double m2Q2, const double m2mu2, const int nf, const int v
        ) const;
        // double C2_ps3_highenergy_highscaleNLL(double x, double
        // m2Q2, double m2mu2, int nf, int v);
        double CL_g3_highenergy_highscaleNLL(
            const double x, const double m2Q2, const double m2mu2, const int nf, const int v
        ) const;
        // double CL_ps3_highenergy_highscaleNLL(double x, double
        // m2Q2, double m2mu2, int nf, int v);

        //==========================================================================================//
        //                  Q^2>>m^2 limit of the high energy
        //                  coefficient functions O(alpha_s^3)
        //------------------------------------------------------------------------------------------//

        double
        C2_g3_highenergy_highscale(const double x, const double m2Q2, const double m2mu2, const int nf, const int v) const;
        double
        C2_ps3_highenergy_highscale(const double x, const double m2Q2, const double m2mu2, const int nf, const int v) const;
        double
        CL_g3_highenergy_highscale(const double x, const double m2Q2, const double m2mu2, const int nf, const int v) const;
        double
        CL_ps3_highenergy_highscale(const double x, const double m2Q2, const double m2mu2, const int nf, const int v) const;

};

class PowerTermsCoefficientFunction : public AbstractHighEnergyCoefficientFunction {

    public:
        PowerTermsCoefficientFunction(const int& order, const char& kind, const char& channel, const bool& NLL = true);
        // PowerTermsCoefficientFunction() : AbstractHighEnergyCoefficientFunction() {} ;
        ~PowerTermsCoefficientFunction() ;

        double fx(const double x, const double m2Q2, const double m2mu2, const int nf) const;

    private:

        HighEnergyCoefficientFunction *highenergy_ ;
        HighEnergyHighScaleCoefficientFunction *highenergyhighscale_ ;

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
