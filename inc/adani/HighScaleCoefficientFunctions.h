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

double C2m_g1_highscale(double x, double mQ);
double CLm_g1_highscale(double x);

double D2m_g1_highscale(double x, double mQ);
double DLm_g1_highscale(double x, double mQ);

//==========================================================================================//
//                      High scale (Q^2 >> m^2) coefficient functions O(alpha_s^2)
//------------------------------------------------------------------------------------------//

double C2m_g2_highscale(double x, double mQ, double mMu);
double C2m_ps2_highscale(double x, double mQ, double mMu);

double CLm_g2_highscale(double x, double mQ, double mMu);
double CLm_ps2_highscale(double x, double mQ, double mMu);

double D2m_g2_highscale(double x, double mQ, double mMu);
double D2m_ps2_highscale(double x, double mQ, double mMu);

double DLm_g2_highscale(double x, double mQ, double mMu);
double DLm_ps2_highscale(double x, double mQ, double mMu);

//==========================================================================================//
//                      High scale (Q^2 >> m^2) coefficient functions O(alpha_s^3)
//------------------------------------------------------------------------------------------//

double C2m_g3_highscale(double x, double mQ, double mMu, int nf, int v=0);
double C2m_ps3_highscale(double x, double mQ, double mMu, int nf);

double CLm_g3_highscale(double x, double mQ, double mMu, int nf);
double CLm_ps3_highscale(double x, double mQ, double mMu, int nf);

double DLm_g3_highscale(double x, double mQ, double mMu, int nf);
double DLm_ps3_highscale(double x, double mQ, double mMu, int nf);

double D2m_g3_highscale(double x, double mQ, double mMu, int nf, int v);
double D2m_ps3_highscale(double x, double mQ, double mMu, int nf);


//==========================================================================================//
//  High scale (Q^2 >> m^2) coefficient functions at O(alpha_s^3) used in [arXiv:1205.5727]
//  klmv = Kawamura, Lo Presti, Moch, Vogt
//  This functions uses the approximation for aQqPS30 (that now is exactly known) from Eq. ??
//  of [arXiv:1205.5727]. It is only used for benchmark against the plots on the paper
//------------------------------------------------------------------------------------------//

double C2m_ps3_highscale_klmv_paper(double x, double mQ, double mMu, int nf, int v);
double D2m_ps3_highscale_klmv_paper(double x, double mQ, double mMu, int nf, int v);


/*
* @name Fortran harmonic polylogarithms
* @brief Harmonic polylogarithms up to weight five
* @param x: real input argument
* @param nw: maximum number of weights requested
* @param Hr1: weight 1 harmonic polylogs (1D array)
* @param Hr2: weight 2 harmonic polylogs (2D array)
* @param Hr3: weight 3 harmonic polylogs (3D array)
* @param Hr4: weight 4 harmonic polylogs (4D array)
* @param Hr5: weight 5 harmonic polylogs (5D array)
* @param n1: lower bound of the weight index requested
* @param n2: upper bound of the weight index requested
* @note This is just a suitably formatted wrapper of the original
* fortran function (see src/hplog.f) to facilitate the call
* of the harmonic logarithms from a C++ code.
*/
extern "C"
{
    double apf_hplog_(double *wx, int *wnw, double *Hr1, double *Hr2, double *Hr3, double *Hr4, double *Hr5, int *wn1, int *wn2);
}

#endif
