/*
 * =====================================================================================
 *
 *       Filename:  ApproximateCoefficientFunctions.h
 *
 *    Description:  Header file for the
 * ApproximateCoefficientFunctions.cc file.
 *
 *         Author:  L'artiglio che graffia
 *
 *  In this file there is the approximation for the unknown
 * O(alpha_s^3) DIS massive coefficient functions.
 *
 * =====================================================================================
 */

#ifndef Approximate_h
#define Approximate_h

//==========================================================================================//
//                      Approximate coefficient functions
//                      O(alpha_s)
//------------------------------------------------------------------------------------------//

double C2_g1_approximation_implicit(double x, double m2Q2, double k, double h);
double C2_g1_approximation(double x, double m2Q2);
double CL_g1_approximation(double x, double m2Q2);

//==========================================================================================//
//                      Approximate coefficient functions
//                      O(alpha_s^2)
//------------------------------------------------------------------------------------------//

struct approximation_parameters {
    double A;
    double B;
    double C;
    double D;
    double a;
    double b;
};

double C2_g2_approximation(double x, double m2Q2, double m2mu2, int v = 0);
double C2_g20_approximation(double x, double m2Q2);
double C2_g20_approximation_implicit(
    double x, double m2Q2, double A, double B, double C, double D, double a,
    double b
);
double C2_g20_approximation_BAND(
    double x, double m2Q2, int v, double var, double fact
);

double C2_ps2_approximation(double x, double m2Q2, double m2mu2, int v = 0);
double C2_ps20_approximation(double x, double m2Q2);
double C2_ps20_approximation_implicit(
    double x, double m2Q2, double A, double B, double C, double D, double a,
    double b
);
double C2_ps20_approximation_BAND(
    double x, double m2Q2, int v, double var, double fact
);

double CL_g2_approximation(double x, double m2Q2, double m2mu2, int v = 0);
double CL_g20_approximation(double x, double m2Q2);
double CL_g20_approximation_implicit(
    double x, double m2Q2, double A, double B, double C, double D, double a,
    double b
);
double CL_g20_approximation_BAND(
    double x, double m2Q2, int v, double var, double fact
);

double CL_ps2_approximation(double x, double m2Q2, double m2mu2, int v = 0);
double CL_ps20_approximation(double x, double m2Q2);
double CL_ps20_approximation_implicit(
    double x, double m2Q2, double A, double B, double C, double D, double a,
    double b
);
double CL_ps20_approximation_BAND(
    double x, double m2Q2, int v, double var, double fact
);

//==========================================================================================//
//                      Approximate coefficient functions
//                      O(alpha_s^3)
//------------------------------------------------------------------------------------------//

#define default_method 0

double C2_g3_approximation(
    double x, double m2Q2, double m2mu2, int nf, int v = 0,
    int method_flag = default_method
);
double C2_g30_approximation(double x, double m2Q2, int nf);
double C2_g30_approximation_implicit(
    double x, double m2Q2, int nf, double A, double B, double C, double D,
    double a, double b, int v1, int v2
);
double C2_g30_approximation_BAND(
    double x, double m2Q2, int nf, int v, double var, double fact
);
double
C2_ps3_approximation(double x, double m2Q2, double m2mu2, int nf, int v = 0);
double C2_ps30_approximation(double x, double m2Q2, int nf);
double C2_ps30_approximation_implicit(
    double x, double m2Q2, int nf, double A, double B, double C, double D,
    double a, double b, int v
);
double C2_ps30_approximation_BAND(
    double x, double m2Q2, int nf, int v, double var, double fact
);

double CL_g3_approximation(
    double x, double m2Q2, double m2mu2, int nf, int v = 0,
    int method_flag = default_method
);
double CL_g30_approximation(double x, double m2Q2, int nf);
double CL_g30_approximation_implicit(
    double x, double m2Q2, int nf, double A, double B, double C, double D,
    double a, double b, int v
);
double CL_g30_approximation_BAND(
    double x, double m2Q2, int nf, int v, double var, double fact
);

double
CL_ps3_approximation(double x, double m2Q2, double m2mu2, int nf, int v = 0);
double CL_ps30_approximation(double x, double m2Q2, int nf);
double CL_ps30_approximation_implicit(
    double x, double m2Q2, int nf, double A, double B, double C, double D,
    double a, double b, int v
);
double CL_ps30_approximation_BAND(
    double x, double m2Q2, int nf, int v, double var, double fact
);

//==========================================================================================//
//              Approximate coefficient functions
//              O(alpha_s^2) from [arXiv:1205.5727] klmv =
//              Kawamura, Lo Presti, Moch, Vogt
//------------------------------------------------------------------------------------------//

double C2_g2_approximationA_klmv(double x, double m2Q2, double m2mu2);
double C2_g2_approximationB_klmv(double x, double m2Q2, double m2mu2);

double C2_ps2_approximationA_klmv(double x, double m2Q2, double m2mu2);
double C2_ps2_approximationB_klmv(double x, double m2Q2, double m2mu2);

//==========================================================================================//
//              Approximate coefficient functions
//              O(alpha_s^3) from [arXiv:1205.5727] klmv =
//              Kawamura, Lo Presti, Moch, Vogt
//------------------------------------------------------------------------------------------//

//==========================================================================================//
// the functions labeled with 'paper' use some approximate
// results for which at the time of the paper
// [arXiv:1205.5727] the exact result was not known (like
// aQqPS30) or for which now we have a better approximation
// (aQg30_B). They are only used as a benchmark against the
// plots of the paper
//------------------------------------------------------------------------------------------//

double C2_g3_approximationA_klmv(
    double x, double m2Q2, double m2mu2, int nf,
    int method_flag = default_method
);
double C2_g3_approximationB_klmv(
    double x, double m2Q2, double m2mu2, int nf,
    int method_flag = default_method
);
double C2_g3_approximationB_klmv_paper(
    double x, double m2Q2, double m2mu2, int nf,
    int method_flag = default_method
);
double C2_g3_approximationBlowxi_klmv(
    double x, double m2Q2, double m2mu2, int nf,
    int method_flag = default_method
);

double C2_ps3_approximationA_klmv(double x, double m2Q2, double m2mu2, int nf);
double C2_ps3_approximationB_klmv(double x, double m2Q2, double m2mu2, int nf);

double
C2_ps3_approximationA_klmv_paper(double x, double m2Q2, double m2mu2, int nf);
double
C2_ps3_approximationB_klmv_paper(double x, double m2Q2, double m2mu2, int nf);

#endif
