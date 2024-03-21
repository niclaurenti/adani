/*
 * =====================================================================================
 *
 *       Filename:  ExactCoefficientFunctions.h
 *
 *    Description:  Header file for the
 * ExactCoefficientFunctions.cc file.
 *
 *         Author:  Hanno un cuore differente
 *
 *  In this file there are the exact heavy coefficient
 *  functions (when known)
 *
 * =====================================================================================
 */

#ifndef Exact_h
#define Exact_h

//==========================================================================================//
//                      The notation used is the following:
//                      C^{(k)} = k-th order expansion in
//                      terms of \alpha_s^{[nf]}
//
//------------------------------------------------------------------------------------------//
//==========================================================================================//
//                      Notation:
//                      m2Q2 = m^2/Q^2
//                      m2mu2 = m^2/mu^2
//------------------------------------------------------------------------------------------//

//==========================================================================================//
//                      Exact massive coefficient functions
//                      O(alpha_s)
//------------------------------------------------------------------------------------------//

double C2_g1(double x, double m2Q2);
double CL_g1(double x, double m2Q2);

//==========================================================================================//
//                      Exact massive coefficient functions
//                      O(alpha_s^2)
//------------------------------------------------------------------------------------------//

/// @cond UNNECESSARY
/**
 * @name Fortran massive coefficient functions
 * Fortran functions for the O(alpha_s^2)
 * coefficient functions from 'src/hqcoef.f'.
 */
///@{
extern "C" {
    // double c2log_(double *wr,double *xi);
    double c2nlog_(double *wr, double *xi);
    double clnlog_(double *wr, double *xi);
    double c2nloq_(double *wr, double *xi);
    double clnloq_(double *wr, double *xi);
    // double c2nlobarg_(double *wr,double *xi);
    // double clnlobarg_(double *wr,double *xi);
    // double c2nlobarq_(double *wr,double *xi);
    // double clnlobarq_(double *wr,double *xi);
}
///@}
/// \endcond

double C2_g2(double x, double m2Q2, double m2mu2);
double C2_ps2(double x, double m2Q2, double m2mu2);

double CL_g2(double x, double m2Q2, double m2mu2);
double CL_ps2(double x, double m2Q2, double m2mu2);

//==========================================================================================//
//  Exact massive coefficient functions O(alpha_s^2):
//  mu-independent terms
//------------------------------------------------------------------------------------------//

double C2_g20(double x, double m2Q2);
double CL_g20(double x, double m2Q2);
double C2_ps20(double x, double m2Q2);
double CL_ps20(double x, double m2Q2);

//==========================================================================================//
//  Exact massive coefficient functions O(alpha_s^2): terms
//  proportional to log(mu^2/m^2)
//------------------------------------------------------------------------------------------//

double C2_g21(double x, double m2Q2);
double CL_g21(double x, double m2Q2);
double C2_ps21(double x, double m2Q2);
double CL_ps21(double x, double m2Q2);

//==========================================================================================//
//  Exact massive coefficient functions O(alpha_s^3): terms
//  proportional to log(mu^2/m^2)
//------------------------------------------------------------------------------------------//

double C2_ps31(double x, double m2Q2, int nf);
double CL_ps31(double x, double m2Q2, int nf);
double C2_g31(double x, double m2Q2, int nf);
double CL_g31(double x, double m2Q2, int nf);

//==========================================================================================//
//  Exact massive coefficient functions O(alpha_s^3): terms
//  proportional to log(mu^2/m^2)^2
//------------------------------------------------------------------------------------------//

double C2_ps32(double x, double m2Q2, int nf);
double CL_ps32(double x, double m2Q2, int nf);

// These two expressions are integrated with montcarlo
// methods

double C2_g32(double x, double m2Q2, int nf, int method_flag);
double CL_g32(double x, double m2Q2, int nf, int method_flag);

#endif
