/*
 * =====================================================================================
 *
 *       Filename:  MassiveCoefficientFunctions.h
 *
 *    Description:  Header file for the MassiveCoefficientFunctions.cc file.
 *
 *         Author:  Dribbla tutti, anche i cammelli del deserto
 *
 *  In this file there are the convolutions between the functions defined in this library.
 *
 * =====================================================================================
 */

#ifndef Convolutions_h
#define Convolutions_h

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

struct function_params {double x; double mQ; int nf;};

//==========================================================================================//
//  Convolution between first order massless quark coefficient functions and first order matchings
//------------------------------------------------------------------------------------------//

double C2_b1_x_K_Qg1(double x, double mMu);
double CL_b1_x_K_Qg1(double x, double mMu);

double K_Qg1_x_K_gg2(double x, double mMu);

double C2_b2_x_K_bg1(double x, double mMu, int nf);

//==========================================================================================//
//  Convolution between first order massive gluon coefficient functions and splitting 
//  function Pgq0
//------------------------------------------------------------------------------------------//

double C2m_g1_x_Pgq0_integrand(double z, void * p) ;
double C2m_g1_x_Pgq0(double x, double mQ) ;
double CLm_g1_x_Pgq0_integrand(double z, void * p) ;
double CLm_g1_x_Pgq0(double x, double mQ) ;

//==========================================================================================//
//  Convolution between first order massive gluon coefficient functions and splitting
//  function Pgg0
//------------------------------------------------------------------------------------------//

double C2m_g1_x_Pgg0_reg_integrand(double z, void * p);
double C2m_g1_x_Pgg0_sing_integrand(double z, void * p) ;
double CLm_g1_x_Pgg0_reg_integrand(double z, void * p) ;
double CLm_g1_x_Pgg0_sing_integrand(double z, void * p);
double Pgg0sing_integrand(double z, void * p);
double C2m_g1_x_Pgg0(double x, double mQ, int nf);
double CLm_g1_x_Pgg0(double x, double mQ, int nf);

//==========================================================================================//
//  Convolution between first order massive gluon coefficient functions and splitting
//  function Pgq1
//------------------------------------------------------------------------------------------//

double C2m_g1_x_Pgq1_integrand(double z, void * p);
double C2m_g1_x_Pgq1(double x, double mQ, int nf) ;
double CLm_g1_x_Pgq1_integrand(double z, void * p);
double CLm_g1_x_Pgq1(double x, double mQ, int nf) ;

//==========================================================================================//
//  Convolution between the mu-independent term of the second order massive gluon coefficient
//  functions and splitting function Pgq0
//------------------------------------------------------------------------------------------//

double C2m_g20_x_Pgq0_integrand(double z, void * p);
double C2m_g20_x_Pgq0(double x, double mQ);
double CLm_g20_x_Pgq0_integrand(double z, void * p);
double CLm_g20_x_Pgq0(double x, double mQ);

//==========================================================================================//
//  Convolution between the mu-independent term of the second order massive quark coefficient
//  functions and splitting function Pqq0
//------------------------------------------------------------------------------------------//

double C2m_ps20_x_Pqq0_reg_integrand(double z, void * p);
double C2m_ps20_x_Pqq0_sing_integrand(double z, void * p);
double CLm_ps20_x_Pqq0_reg_integrand(double z, void * p);
double CLm_ps20_x_Pqq0_sing_integrand(double z, void * p);
double Pqq0sing_integrand(double z, void * p);
double C2m_ps20_x_Pqq0(double x, double mQ, int nf);
double CLm_ps20_x_Pqq0(double x, double mQ, int nf) ;

//==========================================================================================//
//  Convolution between the splitting functions Pgg0/Pqq0 and Pgq0
//------------------------------------------------------------------------------------------//

double Pgg0_x_Pgq0(double x, int nf);
double Pqq0_x_Pgq0(double x) ;

//==========================================================================================//
//  Convolution between thefirst order massive gluon coefficient functions and the convolution
//  between the splitting functions Pgg0/Pqq0 and Pgq0
//------------------------------------------------------------------------------------------//

double C2m_g1_x_Pgg0_x_Pgq0_integrand(double z, void * p);
double C2m_g1_x_Pgg0_x_Pgq0(double x, double mQ, int nf) ;
double C2m_g1_x_Pqq0_x_Pgq0_integrand(double z, void * p);
double C2m_g1_x_Pqq0_x_Pgq0(double x, double mQ, int nf);
double CLm_g1_x_Pgg0_x_Pgq0_integrand(double z, void * p);
double CLm_g1_x_Pgg0_x_Pgq0(double x, double mQ, int nf);
double CLm_g1_x_Pqq0_x_Pgq0_integrand(double z, void * p);
double CLm_g1_x_Pqq0_x_Pgq0(double x, double mQ, int nf);

//==========================================================================================//
//  Convolution between first order massive gluon coefficient functions and splitting
//  function Pgg1
//------------------------------------------------------------------------------------------//

double C2m_g1_x_Pgg1_reg_integrand(double z, void * p);
double C2m_g1_x_Pgg1_sing_integrand(double z, void * p);
double Pgg1sing_integrand(double z, void * p);
double C2m_g1_x_Pgg1(double x, double mQ, int nf) ;
double CLm_g1_x_Pgg1_reg_integrand(double z, void * p);
double CLm_g1_x_Pgg1_sing_integrand(double z, void * p);
double CLm_g1_x_Pgg1(double x, double mQ, int nf);

//==========================================================================================//
//  Convolution between the mu-independent term of the second order massive quark coefficient
//  functions and splitting function Pqg0
//------------------------------------------------------------------------------------------//

double C2m_ps20_x_Pqg0_integrand(double z, void * p);
double C2m_ps20_x_Pqg0(double x, double mQ, int nf);
double CLm_ps20_x_Pqg0_integrand(double z, void * p);
double CLm_ps20_x_Pqg0(double x, double mQ, int nf) ;

//==========================================================================================//
//  Convolution between the mu-independent term of the second order massive gluon coefficient
//  functions and splitting function Pgg0
//------------------------------------------------------------------------------------------//

double C2m_g20_x_Pgg0_reg_integrand(double z, void * p);
double C2m_g20_x_Pgg0_sing_integrand(double z, void * p);
double C2m_g20_x_Pgg0(double x, double mQ, int nf);
double CLm_g20_x_Pgg0_reg_integrand(double z, void * p);
double CLm_g20_x_Pgg0_sing_integrand(double z, void * p);
double CLm_g20_x_Pgg0(double x, double mQ, int nf) ;

//==========================================================================================//
//  Convolution between the first order massive gluon coefficient functions and the convolution
//  between the splitting functions Pqg0 and Pgq0
//------------------------------------------------------------------------------------------//

double Pqg0_x_Pgq0(double x, int nf) ;
double C2m_g1_x_Pqg0_x_Pgq0_integrand(double z, void * p);
double C2m_g1_x_Pqg0_x_Pgq0(double x, double mQ, int nf);
double CLm_g1_x_Pqg0_x_Pgq0_integrand(double z, void * p);
double CLm_g1_x_Pqg0_x_Pgq0(double x, double mQ, int nf);

//==========================================================================================//
//  Convolution between thefirst order massive gluon coefficient functions and the convolution
//  between the splitting functions Pgg0 and Pgg0
//------------------------------------------------------------------------------------------//

double C2m_g1_x_Pgg0_x_Pgg0_reg_integrand(double z, void * p);
double C2m_g1_x_Pgg0_x_Pgg0_sing_integrand(double z, void * p);
double C2m_g1_x_Pgg0_x_Pgg0(double x, double mQ, int nf);
double CLm_g1_x_Pgg0_x_Pgg0_reg_integrand(double z, void * p);
double CLm_g1_x_Pgg0_x_Pgg0_sing_integrand(double z, void * p);
double CLm_g1_x_Pgg0_x_Pgg0(double x, double mQ, int nf);

//==========================================================================================//
//  Convolution between thefirst order massive gluon coefficient functions and the convolution
//  between the splitting functions Pgg0 and Pgg0 using monte carlo mathods
//------------------------------------------------------------------------------------------//

double C2m_g1_x_Pgg0_x_Pgg0_reg1_integrand(double z[], size_t dim, void * p) ;
double C2m_g1_x_Pgg0_x_Pgg0_reg2_integrand(double z[], size_t dim, void * p) ;
double C2m_g1_x_Pgg0_x_Pgg0_reg3_integrand(double z[], size_t dim, void * p);

double C2m_g1_x_Pgg0_x_Pgg0_reg(double x, double mQ, int nf, size_t calls);

double C2m_g1_x_Pgg0_x_Pgg0_sing1_integrand(double z[], size_t dim, void * p) ;
double C2m_g1_x_Pgg0_x_Pgg0_sing2_integrand(double z[], size_t dim, void * p) ;
double C2m_g1_x_Pgg0_x_Pgg0_sing3_integrand(double z[], size_t dim, void * p);

double C2m_g1_x_Pgg0_x_Pgg0_sing(double x, double mQ, int nf, size_t calls);

double C2m_g1_x_Pgg0_x_Pgg0_MC(double x, double mQ, int nf, size_t calls) ;

double CLm_g1_x_Pgg0_x_Pgg0_reg1_integrand(double z[], size_t dim, void * p) ;
double CLm_g1_x_Pgg0_x_Pgg0_reg2_integrand(double z[], size_t dim, void * p) ;
double CLm_g1_x_Pgg0_x_Pgg0_reg3_integrand(double z[], size_t dim, void * p);

double CLm_g1_x_Pgg0_x_Pgg0_reg(double x, double mQ, int nf, size_t calls);

double CLm_g1_x_Pgg0_x_Pgg0_sing1_integrand(double z[], size_t dim, void * p) ;
double CLm_g1_x_Pgg0_x_Pgg0_sing2_integrand(double z[], size_t dim, void * p) ;
double CLm_g1_x_Pgg0_x_Pgg0_sing3_integrand(double z[], size_t dim, void * p);

double CLm_g1_x_Pgg0_x_Pgg0_sing(double x, double mQ, int nf, size_t calls);

double CLm_g1_x_Pgg0_x_Pgg0_MC(double x, double mQ, int nf, size_t calls) ;

#endif