/*
 * =====================================================================================
 *
 *       Filename:  ThresholdCoefficientFunctions.h
 *
 *    Description:  Header file for the
 * ThresholdCoefficientFunctions.cc file.
 *
 *         Author:  La garra charruaaaaaaa
 *
 *  In this file there are the coefficient functions in the
 * threshold limit, i.e. s->4m^2 (where s is the partonic
 * conter of mass energy), or x -> xmax
 *
 * =====================================================================================
 */

#ifndef Threshold_h
#define Threshold_h

//==========================================================================================//
//                      Threshold (s -> 4m^2) coefficient
//                      functions O(alpha_s)
//------------------------------------------------------------------------------------------//

double C2_g1_threshold(double x, double m2Q2);

//==========================================================================================//
//                      Threshold (s -> 4m^2) coefficient
//                      functions O(alpha_s^2)
//------------------------------------------------------------------------------------------//

double C2_g2_threshold(double x, double m2Q2, double m2mu2);
double CL_g2_threshold(double x, double m2Q2, double m2mu2);

double threshold_expansion_g2(double x, double m2Q2, double m2mu2);
double threshold_expansion_g2_const(double x, double m2Q2, double m2mu2);

double C2_g2_threshold_const(double x, double m2Q2, double m2mu2);
double CL_g2_threshold_const(double x, double m2Q2, double m2mu2);

//==========================================================================================//
//                      Threshold (s -> 4m^2) coefficient
//                      functions O(alpha_s^3)
//------------------------------------------------------------------------------------------//

double C2_g3_threshold(double x, double m2Q2, double m2mu2, int nf);
double CL_g3_threshold(double x, double m2Q2, double m2mu2, int nf);

double threshold_expansion_g3(double m2Q2, double m2mu2, int nf);
double threshold_expansion_g3_const(double m2Q2, double m2mu2);

double C2_g3_threshold_const(double x, double m2Q2, double m2mu2);

//==========================================================================================//
//  Function needed for the threshold limit.
//------------------------------------------------------------------------------------------//

double c0(double xi);
double c0_bar(double xi);

#endif
