/*
 * =====================================================================================
 *
 *       Filename:  MassiveCoefficientFunctions.h
 *
 *    Description:  Header file for the MassiveCoefficientFunctions.cc file.
 *
 *         Author:  LeBron James
 *   Organization:  Los Angeles Lakers (at the time of writing)
 *
 *  In this file there are the exact heavy coefficient functions
 *
 * =====================================================================================
 */
#ifndef Massive_h
#define Massive_h

//==========================================================================================//
//                      The notation used is the following:
//                      C^{(k)} = k-th order expansion in terms of \alpha_s^{[nf]}
//  D^{(k)}_4 = k-th order expansion in terms of \alpha_s^{[nf+1]} in thhe nf flavor scheme
//  D^{(k)}_5 = k-th order expansion in terms of \alpha_s^{[nf+1]} in thhe nf+1 flavor scheme
//------------------------------------------------------------------------------------------//
//==========================================================================================//
//                      Notation:
//                      mQ = m^2/Q^2
//                      mMu = m^2/mu^2
//------------------------------------------------------------------------------------------//

//==========================================================================================//
//                      Exact massive coefficient funtions O(alpha_s)
//------------------------------------------------------------------------------------------//

double C2m_g1(double x, double mQ);
double CLm_g1(double x, double mQ);

double D2m_g1_4(double x, double mQ);
double DLm_g1_4(double x, double mQ);
double D2m_g1_5(double x, double mQ);
double DLm_g1_5(double x, double mQ);

//==========================================================================================//
//                      Exact massive coefficient funtions O(alpha_s^2)
//------------------------------------------------------------------------------------------//

double C2m_g2(double x, double mQ, double mMu);
double C2m_ps2(double x, double mQ, double mMu);

double CLm_g2(double x, double mQ, double mMu);
double CLm_ps2(double x, double mQ, double mMu);

double D2m_g2_4(double x, double mQ, double mMu);
double DLm_g2_4(double x, double mQ, double mMu);

double D2m_g2_5(double x, double mQ, double mMu);
double DLm_g2_5(double x, double mQ, double mMu);

//==========================================================================================//
//  Exact massive coefficient funtions O(alpha_s^2): terms proportional to log(mu^2/m^2)
//------------------------------------------------------------------------------------------//

double C2m_g21(double x, double mQ);
double CLm_g21(double x, double mQ);
double C2m_ps21(double x, double mQ);
double CLm_ps21(double x, double mQ);

//==========================================================================================//
//  Exact massive coefficient funtions O(alpha_s^3): terms proportional to log(mu^2/m^2)
//------------------------------------------------------------------------------------------//

double C2m_ps31(double x, double mQ, int nf);
double CLm_ps31(double x, double mQ, int nf);
double C2m_g31(double x, double mQ, int nf);
double CLm_g31(double x, double mQ, int nf);

//==========================================================================================//
//  Exact massive coefficient funtions O(alpha_s^3): terms proportional to log(mu^2/m^2)^2
//------------------------------------------------------------------------------------------//

double C2m_ps32(double x, double mQ, int nf) ;
double CLm_ps32(double x, double mQ, int nf) ;

double C2m_g31(double x, double mQ, int nf) ;
double CLm_g31(double x, double mQ, int nf) ;

// These two expressions are integrated with montcarlo techniques
double C2m_g32(double x, double mQ, int nf, int method_flag = 0, int calls = 500000);
double CLm_g32(double x, double mQ, int nf, int method_flag = 0, int calls = 500000) ;

#endif
