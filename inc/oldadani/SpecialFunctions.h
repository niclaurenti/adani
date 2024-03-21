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

double H(double x, int i);
double H(double x, int i, int j);
double H(double x, int i, int j, int k);

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
