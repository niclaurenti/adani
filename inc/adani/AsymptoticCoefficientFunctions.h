/*
 * =====================================================================================
 *
 *       Filename:  AsymptoticCoefficientFunctions.h
 *
 *    Description:  Header file for the
 * AsymptoticCoefficientFunctions.cc file.
 *
 *         Author:  Anche tu sei un teorico Max
 *
 *  In this file there are the coefficient functions in the
 * asymptotic limit, i.e. x -> 0 and/or Q^2 >> m^2
 *
 * =====================================================================================
 */

#ifndef Asymptotic_h
#define Asymptotic_h

//==========================================================================================//
//                      Asymptotic coefficient functions
//                      O(alpha_s)
//------------------------------------------------------------------------------------------//

double C2_g1_asymptotic(double x, double m2Q2);
double CL_g1_asymptotic(double x, double m2Q2);

//==========================================================================================//
//                      Asymptotic coefficient functions
//                      O(alpha_s^2)
//------------------------------------------------------------------------------------------//

double C2_g2_asymptotic(double x, double m2Q2, double m2mu2);
double C2_ps2_asymptotic(double x, double m2Q2, double m2mu2);

double CL_g2_asymptotic(double x, double m2Q2, double m2mu2);
double CL_ps2_asymptotic(double x, double m2Q2, double m2mu2);

//==========================================================================================//
//                      Asymptotic coefficient functions
//                      O(alpha_s^3) at leading log
//------------------------------------------------------------------------------------------//

// double C2_g3_asymptoticLL(double x, double m2Q2, double
// m2mu2, int nf, int v);

//==========================================================================================//
//                      Asymptotic coefficient functions
//                      O(alpha_s^3) at leading log
//------------------------------------------------------------------------------------------//

double
C2_g3_asymptotic(double x, double m2Q2, double m2mu2, int nf, int v1, int v2);
double C2_ps3_asymptotic(double x, double m2Q2, double m2mu2, int nf, int v);
double CL_g3_asymptotic(double x, double m2Q2, double m2mu2, int nf, int v);
double CL_ps3_asymptotic(double x, double m2Q2, double m2mu2, int nf, int v);

#endif
