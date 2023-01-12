/*
 * =====================================================================================
 *
 *       Filename:  SplittingFunctions.h
 *
 *    Description:  Header file for the SplittingFunctions.cc file.
 *
 *         Author:  LeBron James
 *   Organization:  Los Angeles Lakers (at the time of writing)
 *
 *  In this file there are the exact heavy coefficient functions
 *
 * =====================================================================================
 */
#ifndef Split
#define Split

//==========================================================================================//
//                      Splitting funtions O(alpha_s) without color factors
//------------------------------------------------------------------------------------------//

double pgq(double x);
double pqg(double x);
double pggreg(double x);
double pggsing(double x);

//==========================================================================================//
//                      Splitting funtions O(alpha_s)
//------------------------------------------------------------------------------------------//

double Pgq0(double x);
double Pqg0(double x, int nf);

double Pgg0reg(double x);
double Pgg0loc(int nf);
double Pgg0sing(double x);

double Pqq0reg(double x);
double Pqq0loc();
double Pqq0sing(double x);

//==========================================================================================//
//                      Splitting funtions O(alpha_s^2)
//------------------------------------------------------------------------------------------//

double Pgq1(double x, int nf);

double Pgg1reg(double x, int nf);
double Pgg1sing (double x, int nf);
double Pgg1loc(int nf);

#endif
