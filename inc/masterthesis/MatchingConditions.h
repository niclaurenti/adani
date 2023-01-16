/*
 * =====================================================================================
 *
 *       Filename:  MassiveCoefficientFunctions.h
 *
 *    Description:  Header file for the MassiveCoefficientFunctions.cc file.
 *
 *         Author:  Vamo!!
 *
 *  In this file there are the matching conditions.
 *
 * =====================================================================================
 */

#ifndef Match_h
#define Match_h

//==========================================================================================//
//  Matching conditions O(alpha_s)
//------------------------------------------------------------------------------------------//

double K_Qg1(double x, double mMu);
double K_gg1_local(double mMu);

//==========================================================================================//
//  Matching conditions O(alpha_s^2)
//------------------------------------------------------------------------------------------//

double K_Qg2(double x, double mMu);
double K_Qg2_apfel(double x, double mMu);

//==========================================================================================//
//  Matching conditions O(alpha_s^3)
//------------------------------------------------------------------------------------------//

double a_Qg_30(double x, int v);

#endif
