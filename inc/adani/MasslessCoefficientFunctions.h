/*
 * =====================================================================================
 *
 *       Filename:  MasslessCoefficientFunctions.h
 *
 *    Description:  Header file for the
 * MasslessCoefficientFunctions.cc file.
 *
 *         Author:  Viva el Futbol
 *
 *  In this file there are the massless coefficient
 * functions. Observe that for the O(alpha_s^2) and
 * O(alpha_s^3) it is used a parameterization of the exact
 * result, being the latter a very long equation.
 *
 * =====================================================================================
 */

#ifndef Massless_h
#define Massless_h

//==========================================================================================//
//                      Massless coefficient functions
//                      O(alpha_s)
//------------------------------------------------------------------------------------------//

double C2_g1_massless(double x, int nf);
double CL_g1_massless(double x, int nf);

//==========================================================================================//
//                      Massless coefficient functions
//                      O(alpha_s^2)
//------------------------------------------------------------------------------------------//

double C2_g2_massless(double x, int nf);
double C2_ps2_massless(double x, int nf);
double CL_g2_massless(double x, int nf);
double CL_ps2_massless(double x, int nf);

//==========================================================================================//
//                      Massless coefficient functions
//                      O(alpha_s^3)
//------------------------------------------------------------------------------------------//

double C2_g3_massless(double x, int nf);
double C2_ps3_massless(double x, int nf);
double CL_g3_massless(double x, int nf);
double CL_ps3_massless(double x, int nf);

//==========================================================================================//
//                      Charge factors
//------------------------------------------------------------------------------------------//

double fl11g(int nf);
double fl11ps(int nf);

//==========================================================================================//
//                      Massless coefficient functions (parametrization)
//                      O(alpha_s^2)
//------------------------------------------------------------------------------------------//

// double C2_g2_massless_param(double x, int nf);
// double C2_ps2_massless_param(double x, int nf);
// double CL_g2_massless_param(double x, int nf);
// double CL_ps2_massless_param(double x, int nf);

//==========================================================================================//
//                      Massless coefficient functions (parametrization)
//                      O(alpha_s^3)
//------------------------------------------------------------------------------------------//

// double C2_g3_massless_param(double x, int nf);
// double C2_ps3_massless_param(double x, int nf);
// double CL_g3_massless_param(double x, int nf);
// double CL_ps3_massless_param(double x, int nf);

#endif
