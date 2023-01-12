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
 *  In this file there are the matching conditions.
 *
 * =====================================================================================
 */

#ifndef Match_h
#define Match_h

double K_Qg1(double x, double mMu);
double K_gg1_local(double mMu);

double K_Qg2(double x, double mMu);
double K_Qg2_apfel(double x, double mMu);

double K_Qg3(double x, double mQ, int nf);
double a_Qg_30(double x, int v);

#endif
