/*
 * =====================================================================================
 *
 *       Filename:  ThresholdCoefficientFunctions.h
 *
 *    Description:  Header file for the ThresholdCoefficientFunctions.cc file.
 *
 *         Author:  LeBron James
 *   Organization:  Los Angeles Lakers (at the time of writing)
 *
 *  In this file there are the exact heavy coefficient functions
 *
 * =====================================================================================
 */
#ifndef Threshold_h
#define Threshold_h

//==========================================================================================//
//                      Threshold (s -> 4m^2) coefficient funtions O(alpha_s)
//------------------------------------------------------------------------------------------//
double C2m_g1_threshold(double x, double mQ);

//==========================================================================================//
//                      Threshold (s -> 4m^2) coefficient funtions O(alpha_s^2)
//------------------------------------------------------------------------------------------//

double C2m_g2_threshold(double x, double mQ, double mMu);
double CLm_g2_threshold(double x, double mQ, double mMu);

//==========================================================================================//
//                      Threshold (s -> 4m^2) coefficient funtions O(alpha_s^3)
//------------------------------------------------------------------------------------------//
double C2m_g3_threshold(double x, double mQ, double mMu, int nf);
double CLm_g3_threshold(double x, double mQ, double mMu, int nf);

#endif
