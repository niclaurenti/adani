/*
 * =====================================================================================
 *
 *       Filename:  AsymptoticCoefficientFunctions.h
 *
 *    Description:  Header file for the AsymptoticCoefficientFunctions.cc file.
 *
 *         Author:  Sei un teorico
 *
 *  In this file there are the coefficient functions in the asymptotic limit, i.e. x -> 0
 *  and/or Q^2 >> m^2
 *
 * =====================================================================================
 */

#ifndef Asymptotic_h
#define Asymptotic_h

//==========================================================================================//
//                      Asymptotic coefficient funtions O(alpha_s)
//------------------------------------------------------------------------------------------//

double C2m_g1_asymptotic(double x, double mQ);
double CLm_g1_asymptotic(double x, double mQ);

//==========================================================================================//
//                      Asymptotic coefficient funtions O(alpha_s^2)
//------------------------------------------------------------------------------------------//

double C2m_g2_asymptotic(double x, double mQ, double mMu);
double C2m_ps2_asymptotic(double x, double mQ, double mMu);

double CLm_g2_asymptotic(double x, double mQ, double mMu);
double CLm_ps2_asymptotic(double x, double mQ, double mMu);

//==========================================================================================//
//                      Asymptotic coefficient funtions O(alpha_s^3) at leading log
//------------------------------------------------------------------------------------------//

double C2m_g3_asymptoticLL(double x, double mQ, double mMu, int nf, int v);

//==========================================================================================//
//                      Asymptotic coefficient funtions O(alpha_s^3) at leading log
//------------------------------------------------------------------------------------------//

double C2m_g3_asymptotic(double x, double mQ, double mMu, int nf, int v1=0, int v2=0);
double C2m_ps3_asymptotic(double x, double mQ, double mMu, int nf);
double CLm_g3_asymptotic(double x, double mQ, double mMu, int nf);
double CLm_ps3_asymptotic(double x, double mQ, double mMu, int nf);

#endif
