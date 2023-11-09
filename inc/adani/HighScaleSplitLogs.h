/*
 * =====================================================================================
 *
 *       Filename:  HighScaleSplitLogs.h
 *
 *    Description:  Header file for the
 * HighScaleSplitLogs.cc file.
 *
 *         Author:  Find Quote
 *
 *  In this file there are the coefficient functions in the
 *  high scale limit, i.e. Q^2 >> m^2, implemented setting mu=Q and splitting
 *  the coefficients of each logarithm
 *
 * =====================================================================================
 */

#ifndef HighScaleLogs_h
#define HighScaleLogs_h

//==========================================================================================//
//           Legend:
//           at order O(alpha_s^n)
//           N^(k)LL is the coefficient of log(m^2/Q^2)^(n-k)
//           (The expressions do not contain the log)
//------------------------------------------------------------------------------------------//

//==========================================================================================//
//           Gluon channel, F2
//------------------------------------------------------------------------------------------//

double C2_g3_highscale_LL(double x, int nf);
double C2_g3_highscale_NLL(double x, int nf);
double C2_g3_highscale_N2LL(double x, int nf);
double C2_g3_highscale_N3LL(double x, int nf, int v);

//==========================================================================================//
//           Pure singlet channel, F2
//------------------------------------------------------------------------------------------//

double C2_ps3_highscale_LL(double x, int nf);
double C2_ps3_highscale_NLL(double x, int nf);
double C2_ps3_highscale_N2LL(double x, int nf);
double C2_ps3_highscale_N3LL(double x, int nf);

//==========================================================================================//
//           Gluon channel, FL
//------------------------------------------------------------------------------------------//

double CL_g3_highscale_NLL(double x);
double CL_g3_highscale_N2LL(double x, int nf);
double CL_g3_highscale_N3LL(double x, int nf);

//==========================================================================================//
//           Pure singlet channel, FL
//------------------------------------------------------------------------------------------//

double CL_ps3_highscale_NLL(double x);
double CL_ps3_highscale_N2LL(double x);
double CL_ps3_highscale_N3LL(double x, int nf);

#endif
