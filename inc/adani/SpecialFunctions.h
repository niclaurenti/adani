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

double H_0(double x);
double H_m1(double x);
double H_1(double x);

double H_m1m1(double x);
double H_m10(double x);
double H_m11(double x);
double H_0m1(double x);
double H_00(double x);
double H_01(double x);
double H_1m1(double x);
double H_10(double x);
double H_11(double x);

double H_m1m1m1(double x);
double H_m1m10(double x);
double H_m1m11(double x);
double H_m10m1(double x);
double H_m100(double x);
double H_m101(double x);
double H_m11m1(double x);
double H_m110(double x);
double H_m111(double x);
double H_0m1m1(double x);
double H_0m10(double x);
double H_0m11(double x);
double H_00m1(double x);
double H_000(double x);
double H_001(double x);
double H_01m1(double x);
double H_010(double x);
double H_011(double x);
double H_1m1m1(double x);
double H_1m10(double x);
double H_1m11(double x);
double H_10m1(double x);
double H_100(double x);
double H_101(double x);
double H_11m1(double x);
double H_110(double x);
double H_111(double x);

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
