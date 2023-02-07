/*
 * =====================================================================================
 *
 *       Filename:  AsymptoticCoefficientFunctions.h
 *
 *    Description:  Header file for the AsymptoticCoefficientFunctions.cc file.
 *
 *         Author:  Anche tu sei un teorico Max
 *
 *  In this file there are the coefficient functions in the asymptotic limit, i.e. x -> 0
 *  and/or Q^2 >> m^2
 *
 * =====================================================================================
 */

#ifndef Asymptotic_h
#define Asymptotic_h

//==========================================================================================//
//                      Asymptotic coefficient functions O(alpha_s)
//------------------------------------------------------------------------------------------//

double C2_g1_asymptotic(double x, double mQ);
double CL_g1_asymptotic(double x, double mQ);

//==========================================================================================//
//                      Asymptotic coefficient functions O(alpha_s^2)
//------------------------------------------------------------------------------------------//

double C2_g2_asymptotic(double x, double mQ, double mMu);
double C2_ps2_asymptotic(double x, double mQ, double mMu);

double CL_g2_asymptotic(double x, double mQ, double mMu);
double CL_ps2_asymptotic(double x, double mQ, double mMu);

//==========================================================================================//
//                      Asymptotic coefficient functions O(alpha_s^3) at leading log
//------------------------------------------------------------------------------------------//

double C2_g3_asymptoticLL(double x, double mQ, double mMu, int nf, int v);

//==========================================================================================//
//                      Asymptotic coefficient functions O(alpha_s^3) at leading log
//------------------------------------------------------------------------------------------//

double C2_g3_asymptotic(double x, double mQ, double mMu, int nf, int v1=0, int v2=0);
double C2_ps3_asymptotic(double x, double mQ, double mMu, int nf);
double CL_g3_asymptotic(double x, double mQ, double mMu, int nf);
double CL_ps3_asymptotic(double x, double mQ, double mMu, int nf);

#endif
