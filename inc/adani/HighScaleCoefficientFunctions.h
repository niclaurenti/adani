/*
 * =====================================================================================
 *
 *       Filename:  HighScaleCoefficientFunctions.h
 *
 *    Description:  Header file for the HighScaleCoefficientFunctions.cc file.
 *
 *         Author:  Da quella posizione trasforma l'acqua in vino
 *
 *  In this file there are the coefficient functions in the high scale limit, i.e. Q^2 >> m^2
 *
 * =====================================================================================
 */

#ifndef HighScale_h
#define HighScale_h

//==========================================================================================//
//                      The notation used is the following:
//                      C^{(k)} = k-th order expansion in terms of \alpha_s^{[nf]}
//                      D^{(k)} = k-th order expansion in terms of \alpha_s^{[nf+1]}
//------------------------------------------------------------------------------------------//

//==========================================================================================//
//                      High scale (Q^2 >> m^2) coefficient functions O(alpha_s)
//------------------------------------------------------------------------------------------//

double C2_g1_highscale(double x, double mQ);
double CL_g1_highscale(double x);

double D2_g1_highscale(double x, double mQ);
double DL_g1_highscale(double x, double mQ);

//==========================================================================================//
//                      High scale (Q^2 >> m^2) coefficient functions O(alpha_s^2)
//------------------------------------------------------------------------------------------//

double C2_g2_highscale(double x, double mQ, double mMu);
double C2_ps2_highscale(double z, double mQ, double mMu);

double CL_g2_highscale(double x, double mQ, double mMu);
double CL_ps2_highscale(double x, double mQ, double mMu);

double D2_g2_highscale(double x, double mQ, double mMu);
double D2_ps2_highscale(double x, double mQ, double mMu);

double DL_g2_highscale(double z, double mQ, double mMu);
double DL_ps2_highscale(double z, double mQ, double mMu);

//==========================================================================================//
//                      High scale (Q^2 >> m^2) coefficient functions O(alpha_s^3)
//------------------------------------------------------------------------------------------//

double C2_g3_highscale(double x, double mQ, double mMu, int nf, int v=0);
double C2_ps3_highscale(double x, double mQ, double mMu, int nf);

double CL_g3_highscale(double x, double mQ, double mMu, int nf);
double CL_ps3_highscale(double x, double mQ, double mMu, int nf);

double DL_g3_highscale(double z, double mQ, double mMu, int nf);
double DL_ps3_highscale(double z, double mQ, double mMu, int nf);

double D2_g3_highscale(double x, double mQ, double mMu, int nf, int v);
double D2_ps3_highscale(double x, double mQ, double mMu, int nf, int v=0);


//==========================================================================================//
//  High scale (Q^2 >> m^2) coefficient functions at O(alpha_s^3) used in [arXiv:1205.5727]
//  klmv = Kawamura, Lo Presti, Moch, Vogt
//  This functions uses the approximation for aQqPS30 (that now is exactly known) from Eq. ??
//  of [arXiv:1205.5727]. It is only used for benchmark against the plots on the paper
//------------------------------------------------------------------------------------------//

double C2_ps3_highscale_klmv_paper(double x, double mQ, double mMu, int nf, int v);
double D2_ps3_highscale_klmv_paper(double x, double mQ, double mMu, int nf, int v);

#endif
