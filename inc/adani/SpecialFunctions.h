/*
 * =====================================================================================
 *
 *       Filename:  SpecialFunctions.h
 *
 *    Description:  Header file for the SpecialFunctions.cc
 * file.
 *
 *         Author:  Allison? No, Franco Armani
 *
 *  In this file there are some useful functions.
 *
 * =====================================================================================
 */

#ifndef Special_h
#define Special_h

//==========================================================================================//
//  Beta function.
//------------------------------------------------------------------------------------------//

double beta(int ord, int nf);

//==========================================================================================//
//  Theta function. (no longer used)
//------------------------------------------------------------------------------------------//

// double theta(double x);

//==========================================================================================//
//  Polylogarithms.
//------------------------------------------------------------------------------------------//

double Li2(double x);
double Li3(double x);

//==========================================================================================//
//  Nielsen function.
//------------------------------------------------------------------------------------------//

double S12(double x);

//==========================================================================================//
//  Harmonic polylogarithms up to weight 3.
//------------------------------------------------------------------------------------------//

double H0(double x);
double Hm1(double x);
double H1(double x);

double Hm1m1(double x);
double Hm10(double x);
double Hm11(double x);
double H0m1(double x);
double H00(double x);
double H01(double x);
double H1m1(double x);
double H10(double x);
double H11(double x);

double Hm1m1m1(double x);
double Hm1m10(double x);
double Hm1m11(double x);
double Hm10m1(double x);
double Hm100(double x);
double Hm101(double x);
double Hm11m1(double x);
double Hm110(double x);
double Hm111(double x);
double H0m1m1(double x);
double H0m10(double x);
double H0m11(double x);
double H00m1(double x);
double H000(double x);
double H001(double x);
double H01m1(double x);
double H010(double x);
double H011(double x);
double H1m1m1(double x);
double H1m10(double x);
double H1m11(double x);
double H10m1(double x);
double H100(double x);
double H101(double x);
double H11m1(double x);
double H110(double x);
double H111(double x);

//==========================================================================================//
//  Harmonic polylogarithms up to weight 5.
//------------------------------------------------------------------------------------------//

/**
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
 * @note This is just a suitably formatted wrapper of the
 * original fortran function (see src/hplog.f) to facilitate
 * the call of the harmonic logarithms from a C++ code.
 */
extern "C" {
    double apf_hplog_(
        double *wx, int *wnw, double *Hr1, double *Hr2, double *Hr3,
        double *Hr4, double *Hr5, int *wn1, int *wn2
    );
}

#endif
