/*
 * =====================================================================================
 *
 *       Filename:  SplittingFunctions.h
 *
 *    Description:  Header file for the
 * SplittingFunctions.cc file.
 *
 *         Author:  Al livello di servilismo come siamo
 * messi?
 *
 *  In this file there are the splitting functions.
 *
 * =====================================================================================
 */

#ifndef Split
#define Split

//==========================================================================================//
//                      Splitting functions O(alpha_s)
//                      without color factors
//------------------------------------------------------------------------------------------//

double pgq(double x);
double pqg(double x);
double pggreg(double x);
double pggsing(double x);

//==========================================================================================//
//                      Splitting functions O(alpha_s)
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
//                      Splitting functions O(alpha_s^2)
//------------------------------------------------------------------------------------------//

double Pgq1(double x, int nf);

double Pgg1reg(double x, int nf);
double Pgg1sing(double x, int nf);
double Pgg1loc(int nf);

#endif
