/*
 * =====================================================================================
 *
 *       Filename:  ApproximateCoefficientFunctions.h
 *
 *    Description:  Header file for the ApproximateCoefficientFunctions.cc file.
 *
 *         Author:  Dwyane Wade
 *   Organization:  retired (after a glorious carrer with the Miami Heat)
 *
 *  In this file there is the approximation for the unknown O(alpha_s^3) DIS massive
 *  coefficient functions.
 *
 * =====================================================================================
 */

#ifndef Approximate_h
#define Approximate_h

//==========================================================================================//
//                      Approximate coefficient funtions O(alpha_s)
//------------------------------------------------------------------------------------------//

double C2m_g1_approximation(double x, double mQ, double k, double h);
double C2m_g1_approximation(double x, double mQ);
double CLm_g1_approximation(double x, double mQ);

//==========================================================================================//
//                      Approximate coefficient funtions O(alpha_s^2)
//------------------------------------------------------------------------------------------//

double C2m_g2_approximation(double x, double mQ, double mMu);
double C2m_g2_approximation(double x, double mQ, double mMu, double A, double B, double C, double D, double a, double b);
double C2m_g2_approximation_BAND(double x, double mQ, double mMu, double var, double fact, int v);

double C2m_ps2_approximation(double x, double mQ, double mMu);
double C2m_ps2_approximation(double x, double mQ, double mMu, double A, double B, double C, double D, double a, double b);

double CLm_g2_approximation(double x, double mQ, double mMu, double A, double B, double C, double D, double a, double b);
double CLm_g2_approximation(double x, double mQ, double mMu);

double CLm_ps2_approximation(double x, double mQ, double mMu);
double CLm_ps2_approximation(double x, double mQ, double mMu, double A, double B, double C, double D, double a, double b);

//==========================================================================================//
//                      Approximate coefficient funtions O(alpha_s^3)
//------------------------------------------------------------------------------------------//

double C2m_g3_approximation(double x, double mQ, double mMu, int nf, int method_flag = 0, int calls = 500000);
double C2m_g3_approximation(double x, double mQ, double mMu, int nf, double A, double B, double C, double D, double a, double b, int v1,int v2, int method_flag = 0, int calls = 500000);
double C2m_g3_approximation_BAND(double x, double mQ, double mMu, int nf, double var, double fact, int v, int method_flag = 0, int calls = 500000);

double C2m_ps3_approximation(double x, double mQ, double mMu, int nf);

double CLm_g3_approximation(double x, double mQ, double mMu, int nf, int method_flag = 0, int calls = 500000);

double CLm_ps3_approximation(double x, double mQ, double mMu, int nf);

//==========================================================================================//
//              Approximate coefficient funtions O(alpha_s^2) from [arXiv:1205.5727]
//              klmv = Kawamura, Lo Presti, Moch, Vogt
//------------------------------------------------------------------------------------------//

double C2m_g2_approximationA_klmv(double x, double mQ, double mMu);
double C2m_g2_approximationB_klmv(double x, double mQ, double mMu);

double C2m_ps2_approximationA_klmv(double x, double mQ, double mMu);
double C2m_ps2_approximationB_klmv(double x, double mQ, double mMu);

//==========================================================================================//
//              Approximate coefficient funtions O(alpha_s^3) from [arXiv:1205.5727]
//              klmv = Kawamura, Lo Presti, Moch, Vogt
//------------------------------------------------------------------------------------------//

// the functions labeled with 'paper' use some approximate results for which at
// the time of the paper [arXiv:1205.5727] the exact result was not known.
// They are only used as a benchmark against the plots of the paper

double C2m_g3_approximationA_klmv(double x, double mQ, double mMu, int nf);
double C2m_g3_approximationB_klmv(double x, double mQ, double mMu, int nf);
double C2m_g3_approximationB_klmv_paper(double x, double mQ, double mMu, int nf) ;
double C2m_g3_approximationBlowxi_klmv(double x, double mQ, double mMu, int nf);

double C2m_ps3_approximationA_klmv(double x, double mQ, double mMu, int nf);
double C2m_ps3_approximationB_klmv(double x, double mQ, double mMu, int nf);

double C2m_ps3_approximationA_klmv_paper(double x, double mQ, double mMu, int nf);
double C2m_ps3_approximationB_klmv_paper(double x, double mQ, double mMu, int nf);

#endif
