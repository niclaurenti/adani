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
 *  In this file there are some useful functions.
 *
 * =====================================================================================
 */

#ifndef Special_h
#define Special_h

double zeta(int i);
double beta(int ord, int nf);

double theta(double x) ;

double Li2(double x);
double Li3(double x);
double Li4(double x);

double S12(double x);

double c0(double xi);
double c0_bar(double xi);

double H(double x, int i);
double H(double x, int i, int j);
double H(double x, int i, int j, int k);

#endif
