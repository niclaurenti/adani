/*
 * =====================================================================================
 *
 *       Filename:  ThresholdCoefficientFunctions.h
 *
 *    Description:  Header file for the ThresholdCoefficientFunctions.cc file.
 *
 *         Author:  La garra charruaaaaaaa
 *
 *  In this file there are the coefficient functions in the threshold limit, i.e. s->4m^2
 *  (where s is the partonic conter of mass energy), or x -> xmax
 *
 * =====================================================================================
 */

#ifndef Threshold_h
#define Threshold_h

//==========================================================================================//
//                      Threshold (s -> 4m^2) coefficient functions O(alpha_s)
//------------------------------------------------------------------------------------------//

double C2m_g1_threshold(double x, double mQ);

//==========================================================================================//
//                      Threshold (s -> 4m^2) coefficient functions O(alpha_s^2)
//------------------------------------------------------------------------------------------//

double C2m_g2_threshold(double x, double mQ, double mMu);
double CLm_g2_threshold(double x, double mQ, double mMu);

double threshold_expansion_g2(double x, double mQ, double mMu);
double threshold_expansion_g2_const(double x, double mQ, double mMu);

double C2m_g2_threshold_const(double x, double mQ, double mMu);
double CLm_g2_threshold_const(double x, double mQ, double mMu);

//==========================================================================================//
//                      Threshold (s -> 4m^2) coefficient functions O(alpha_s^3)
//------------------------------------------------------------------------------------------//

double C2m_g3_threshold(double x, double mQ, double mMu, int nf);
double CLm_g3_threshold(double x, double mQ, double mMu, int nf);

double threshold_expansion_g3(double mQ, double mMu, int nf);
double threshold_expansion_g3_const(double mQ, double mMu);

double C2m_g3_threshold_const(double x, double mQ, double mMu);

//==========================================================================================//
//  Function needed for the threshold limit.
//------------------------------------------------------------------------------------------//

double c0(double xi);
double c0_bar(double xi);

#endif
