/*
 * =====================================================================================
 *
 *       Filename:  Convolutions.h
 *
 *    Description:  Header file for the Convolutions.cc
 * file.
 *
 *         Author:  Dribbla tutti, anche i cammelli del
 * deserto
 *
 *  In this file there are the convolutions between the
 * functions defined in this library.
 *
 *  For the convolution between a regular function and a
 * function containing regular, singular (plus distribution)
 * and local (delta) parts see Eq. (4.7) of Ref.
 * [arXiv:hep-ph/0504242v1]
 *
 * =====================================================================================
 */

#ifndef Convolutions_h
#define Convolutions_h

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

struct function_params {
    double x;
    double m2Q2;
    int nf;
};

//==========================================================================================//
//  Convolution between first order massless quark
//  coefficient functions and first order matchings
//------------------------------------------------------------------------------------------//

double C2_Q1_x_K_Qg1(double x, double m2mu2);
double CL_Q1_x_K_Qg1(double x, double m2mu2);

//==========================================================================================//
//  WARNING: from now on all the convolutions are numerical.
//  In all the integrals we have switched off the ERROR
//  raised by gsl when the numerical precision is not
//  achieved. This is not ideal but as a wise man once said:
//  "If one has no better options, then he has sex even with
//  his own wife" (the identity of this enlightened person
//  will not be revealed)
//------------------------------------------------------------------------------------------//

//==========================================================================================//
//  Convolution between first order massive gluon
//  coefficient functions and splitting function Pgq0
//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgq0_integrand(double z, void *p);
double C2_g1_x_Pgq0(double x, double m2Q2);
double CL_g1_x_Pgq0_integrand(double z, void *p);
double CL_g1_x_Pgq0(double x, double m2Q2);

//==========================================================================================//
//  Convolution between first order massive gluon
//  coefficient functions and splitting function Pgg0
//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_reg_integrand(double z, void *p);
double C2_g1_x_Pgg0_sing_integrand(double z, void *p);
double CL_g1_x_Pgg0_reg_integrand(double z, void *p);
double CL_g1_x_Pgg0_sing_integrand(double z, void *p);
double Pgg0sing_integrated(double x);
double C2_g1_x_Pgg0(double x, double m2Q2, int nf);
double CL_g1_x_Pgg0(double x, double m2Q2, int nf);

//==========================================================================================//
//  Convolution between first order massive gluon
//  coefficient functions and splitting function Pgq1
//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgq1_integrand(double z, void *p);
double C2_g1_x_Pgq1(double x, double m2Q2, int nf);
double CL_g1_x_Pgq1_integrand(double z, void *p);
double CL_g1_x_Pgq1(double x, double m2Q2, int nf);

//==========================================================================================//
//  Convolution between the mu-independent term of the
//  second order massive gluon coefficient functions and
//  splitting function Pgq0
//------------------------------------------------------------------------------------------//

double C2_g20_x_Pgq0_integrand(double z, void *p);
double C2_g20_x_Pgq0(double x, double m2Q2);
double CL_g20_x_Pgq0_integrand(double z, void *p);
double CL_g20_x_Pgq0(double x, double m2Q2);

//==========================================================================================//
//  Convolution between the mu-independent term of the
//  second order massive quark coefficient functions and
//  splitting function Pqq0
//------------------------------------------------------------------------------------------//

double C2_ps20_x_Pqq0_reg_integrand(double z, void *p);
double C2_ps20_x_Pqq0_sing_integrand(double z, void *p);
double CL_ps20_x_Pqq0_reg_integrand(double z, void *p);
double CL_ps20_x_Pqq0_sing_integrand(double z, void *p);
double Pqq0sing_integrated(double x);
double C2_ps20_x_Pqq0(double x, double m2Q2, int nf);
double CL_ps20_x_Pqq0(double x, double m2Q2, int nf);

//==========================================================================================//
//  Convolution between the splitting functions Pgg0/Pqq0
//  and Pgq0
//------------------------------------------------------------------------------------------//

double Pgg0_x_Pgq0(double x, int nf);
double Pqq0_x_Pgq0(double x);

//==========================================================================================//
//  Convolution between the first order massive gluon
//  coefficient functions and the convolution between the
//  splitting functions Pgg0/Pqq0 and Pgq0
//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgq0_integrand(double z, void *p);
double C2_g1_x_Pgg0_x_Pgq0(double x, double m2Q2, int nf);
double C2_g1_x_Pqq0_x_Pgq0_integrand(double z, void *p);
double C2_g1_x_Pqq0_x_Pgq0(double x, double m2Q2, int nf);
double CL_g1_x_Pgg0_x_Pgq0_integrand(double z, void *p);
double CL_g1_x_Pgg0_x_Pgq0(double x, double m2Q2, int nf);
double CL_g1_x_Pqq0_x_Pgq0_integrand(double z, void *p);
double CL_g1_x_Pqq0_x_Pgq0(double x, double m2Q2, int nf);

//==========================================================================================//
//  Convolution between first order massive gluon
//  coefficient functions and splitting function Pgg1
//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg1_reg_integrand(double z, void *p);
double C2_g1_x_Pgg1_sing_integrand(double z, void *p);
double Pgg1sing_integrated(double x, int nf);
double C2_g1_x_Pgg1(double x, double m2Q2, int nf);
double CL_g1_x_Pgg1_reg_integrand(double z, void *p);
double CL_g1_x_Pgg1_sing_integrand(double z, void *p);
double CL_g1_x_Pgg1(double x, double m2Q2, int nf);

//==========================================================================================//
//  Convolution between the mu-independent term of the
//  second order massive quark coefficient functions and
//  splitting function Pqg0
//------------------------------------------------------------------------------------------//

double C2_ps20_x_Pqg0_integrand(double z, void *p);
double C2_ps20_x_Pqg0(double x, double m2Q2, int nf);
double CL_ps20_x_Pqg0_integrand(double z, void *p);
double CL_ps20_x_Pqg0(double x, double m2Q2, int nf);

//==========================================================================================//
//  Convolution between the mu-independent term of the
//  second order massive gluon coefficient functions and
//  splitting function Pgg0
//------------------------------------------------------------------------------------------//

double C2_g20_x_Pgg0_reg_integrand(double z, void *p);
double C2_g20_x_Pgg0_sing_integrand(double z, void *p);
double C2_g20_x_Pgg0(double x, double m2Q2, int nf);
double CL_g20_x_Pgg0_reg_integrand(double z, void *p);
double CL_g20_x_Pgg0_sing_integrand(double z, void *p);
double CL_g20_x_Pgg0(double x, double m2Q2, int nf);

//==========================================================================================//
//  Convolution between the first order massive gluon
//  coefficient functions and the convolution between the
//  splitting functions Pqg0 and Pgq0
//------------------------------------------------------------------------------------------//

double Pqg0_x_Pgq0(double x, int nf);
double C2_g1_x_Pqg0_x_Pgq0_integrand(double z, void *p);
double C2_g1_x_Pqg0_x_Pgq0(double x, double m2Q2, int nf);
double CL_g1_x_Pqg0_x_Pgq0_integrand(double z, void *p);
double CL_g1_x_Pqg0_x_Pgq0(double x, double m2Q2, int nf);

//==========================================================================================//
//  Convolution between the first order massive gluon
//  coefficient functions and the convolution between the
//  splitting functions Pgg0 and Pgg0
//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgg0_reg_integrand(double z, void *p);
double C2_g1_x_Pgg0_x_Pgg0_sing_integrand(double z, void *p);
double C2_g1_x_Pgg0_x_Pgg0(double x, double m2Q2, int nf);
double CL_g1_x_Pgg0_x_Pgg0_reg_integrand(double z, void *p);
double CL_g1_x_Pgg0_x_Pgg0_sing_integrand(double z, void *p);
double CL_g1_x_Pgg0_x_Pgg0(double x, double m2Q2, int nf);

//==========================================================================================//
//  Convolution between the first order massive gluon
//  coefficient functions and the convolution between the
//  splitting functions Pgg0 and Pgg0 using monte carlo
//  mathods
//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgg0_reg1_integrand(double z[], size_t dim, void *p);
double C2_g1_x_Pgg0_x_Pgg0_reg2_integrand(double z[], size_t dim, void *p);
double C2_g1_x_Pgg0_x_Pgg0_reg3_integrand(double z[], size_t dim, void *p);

double C2_g1_x_Pgg0_x_Pgg0_reg(double x, double m2Q2, int nf);

double C2_g1_x_Pgg0_x_Pgg0_sing1_integrand(double z[], size_t dim, void *p);
double C2_g1_x_Pgg0_x_Pgg0_sing2_integrand(double z[], size_t dim, void *p);
double C2_g1_x_Pgg0_x_Pgg0_sing3_integrand(double z[], size_t dim, void *p);

double C2_g1_x_Pgg0_x_Pgg0_sing(double x, double m2Q2, int nf);

double C2_g1_x_Pgg0_x_Pgg0_MC(double x, double m2Q2, int nf);

double CL_g1_x_Pgg0_x_Pgg0_reg1_integrand(double z[], size_t dim, void *p);
double CL_g1_x_Pgg0_x_Pgg0_reg2_integrand(double z[], size_t dim, void *p);
double CL_g1_x_Pgg0_x_Pgg0_reg3_integrand(double z[], size_t dim, void *p);

double CL_g1_x_Pgg0_x_Pgg0_reg(double x, double m2Q2, int nf);

double CL_g1_x_Pgg0_x_Pgg0_sing1_integrand(double z[], size_t dim, void *p);
double CL_g1_x_Pgg0_x_Pgg0_sing2_integrand(double z[], size_t dim, void *p);
double CL_g1_x_Pgg0_x_Pgg0_sing3_integrand(double z[], size_t dim, void *p);

double CL_g1_x_Pgg0_x_Pgg0_sing(double x, double m2Q2, int nf);

double CL_g1_x_Pgg0_x_Pgg0_MC(double x, double m2Q2, int nf);

#endif
