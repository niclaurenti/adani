/*
 * =====================================================================================
 *
 *       Filename:  MasslessCoefficientFunctions.h
 *
 *    Description:  Header file for the MasslessCoefficientFunctions.cc file.
 *
 *         Author:  
 *
 *  In this file there are the massless coefficient functions. Observe that for the
 *  O(alpha_s^2) and O(alpha_s^3) it is used a parameterization of the exact result,
 *  being the latter a very long equation.
 *
 * =====================================================================================
 */

#ifndef Massless_h
#define Massless_h

//==========================================================================================//
//                      Massless coefficient funtions O(alpha_s)
//------------------------------------------------------------------------------------------//

double C2_g1(double x, int nf);
double CL_g1(double x, int nf);

//==========================================================================================//
//                      Massless coefficient funtions O(alpha_s^2)
//------------------------------------------------------------------------------------------//

double C2_g2(double x, int nf);
double C2_ps2(double x, int nf);

double CL_g2(double x, int nf);
double CL_ps2(double x, int nf);

//==========================================================================================//
//                      Massless coefficient funtions O(alpha_s^3)
//------------------------------------------------------------------------------------------//

double C2_g3(double x, int nf);
double C2_ps3(double x, int nf);

double CL_g3(double x, int nf);
double CL_ps3(double x, int nf);


#endif
