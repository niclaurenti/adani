/*
 * =====================================================================================
 *
 *       Filename:  HighEnergyCoefficientFunctions.h
 *
 *    Description:  Header file for the HighEnergyCoefficientFunctions.cc file.
 *
 *         Author:  L'allenamento si fa!
 *
 *  In this file there are the coefficient functions in the high energy limit, i.e. x -> 0
 *
 * =====================================================================================
 */

#ifndef HighEnergy_h
#define HighEnergy_h

//==========================================================================================//
//                      Notation:
//      High energy: small x limit
//      High energy high scale: Q^2 >> m^2 limit of the small x limit (or the opposite)
//      Power terms: power terms in the small x liimit, obtained via
//          C_powerterms = C_highenergy - C_highenergy_highscale
//------------------------------------------------------------------------------------------//

//==========================================================================================//
//                      High energy coefficient functions O(alpha_s^2)
//------------------------------------------------------------------------------------------//

double C2_g2_highenergy(double x, double mQ, double mMu);
double C2_ps2_highenergy(double x, double mQ, double mMu);
double CL_g2_highenergy(double x, double mQ, double mMu);
double CL_ps2_highenergy(double x, double mQ, double mMu);

//==========================================================================================//
//                      Q>>m limit of the high energy coefficient functions O(alpha_s^2)
//------------------------------------------------------------------------------------------//

double C2_g2_highenergy_highscale(double x, double mQ, double mMu);
double C2_ps2_highenergy_highscale(double x, double mQ, double mMu);
double CL_g2_highenergy_highscale(double x, double mQ, double mMu);
double CL_ps2_highenergy_highscale(double x, double mQ, double mMu);

//==========================================================================================//
//                  Power terms of the coefficient function O(alpha_s^2)
//------------------------------------------------------------------------------------------//

double C2_g2_power_terms(double x, double mQ, double mMu);
double C2_ps2_power_terms(double x, double mQ, double mMu);
double CL_g2_power_terms(double x, double mQ, double mMu);
double CL_ps2_power_terms(double x, double mQ, double mMu);

//==========================================================================================//
//                      High energy coefficient functions O(alpha_s^3) at leading log
//------------------------------------------------------------------------------------------//

double C2_g3_highenergyLL(double x, double mQ, double mMu);
double CL_g3_highenergyLL(double x, double mQ, double mMu);
double C2_ps3_highenergyLL(double x, double mQ, double mMu);
// double CL_ps3_highenergyLL(double x, double mQ, double mMu);

double C2_g3_highenergyNLL(double x, double mQ, double mMu, int nf, int v);
// double C2_ps3_highenergyNLL(double x, double mQ, double mMu, int nf, int v);
double CL_g3_highenergyNLL(double x, double mQ, double mMu, int nf, int v);
// double CL_ps3_highenergyNLL(double x, double mQ, double mMu, int nf, int v);

//==========================================================================================//
//                      High energy coefficient functions O(alpha_s^3)
//------------------------------------------------------------------------------------------//

double C2_g3_highenergy(double x, double mQ, double mMu, int nf, int v);
double C2_ps3_highenergy(double x, double mQ, double mMu, int nf, int v);
double CL_g3_highenergy(double x, double mQ, double mMu, int nf, int v);
double CL_ps3_highenergy(double x, double mQ, double mMu, int nf, int v);

//==========================================================================================//
//  Q>>m limit of the high energy coefficient functions O(alpha_s^3) at leading log
//------------------------------------------------------------------------------------------//

double C2_g3_highenergy_highscaleLL(double x, double mQ, double mMu);
double CL_g3_highenergy_highscaleLL(double x, double mQ, double mMu);
// double C2_ps3_highenergy_highscaleLL(double x, double mQ, double mMu);
// double CL_ps3_highenergy_highscaleLL(double x, double mQ, double mMu);

double C2_g3_highenergy_highscaleNLL(double x, double mQ, double mMu, int nf, int v);
// double C2_ps3_highenergy_highscaleNLL(double x, double mQ, double mMu, int nf, int v);
double CL_g3_highenergy_highscaleNLL(double x, double mQ, double mMu, int nf, int v);
// double CL_ps3_highenergy_highscaleNLL(double x, double mQ, double mMu, int nf, int v);

//==========================================================================================//
//                  Q^2>>m^2 limit of the high energy coefficient functions O(alpha_s^3)
//------------------------------------------------------------------------------------------//

double C2_g3_highenergy_highscale(double x, double mQ, double mMu, int nf, int v);
double C2_ps3_highenergy_highscale(double x, double mQ, double mMu, int nf, int v);
double CL_g3_highenergy_highscale(double x, double mQ, double mMu, int nf, int v);
double CL_ps3_highenergy_highscale(double x, double mQ, double mMu, int nf, int v);

//==========================================================================================//
//                  Power terms of the coefficient function O(alpha_s^3) at leading log
//------------------------------------------------------------------------------------------//

// double C2_g3_power_termsLL(double x, double mQ , double mMu);
// double C2_ps3_power_termsLL(double x, double mQ , double mMu);

//==========================================================================================//
//                  Power terms of the coefficient function O(alpha_s^3)
//------------------------------------------------------------------------------------------//

double C2_g3_power_terms(double x, double mQ , double mMu, int nf, int v);
double C2_ps3_power_terms(double x, double mQ , double mMu, int nf, int v);
double CL_g3_power_terms(double x, double mQ , double mMu, int nf, int v);
double CL_ps3_power_terms(double x, double mQ , double mMu, int nf, int v);

#endif
