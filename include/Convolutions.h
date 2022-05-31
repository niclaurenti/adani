#ifndef Convolutions_h
#define Convolutions_h

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

struct function_params {double x; double mQ; int nf;};

double C2_b1_x_K_bg1(double x, double mMu);
double CL_b1_x_K_bg1(double x, double mMu);

double K_bg1_x_K_gg2(double x, double mMu);

double C2_b2_x_K_bg1(double x, double mMu, int nf);

double C2m_g1_x_Pgq0_integrand(double z, void * p) ;
double C2m_g1_x_Pgq0(double x, double mQ) ;
double CLm_g1_x_Pgq0_integrand(double z, void * p) ;
double CLm_g1_x_Pgq0(double x, double mQ) ;

double C2m_g1_x_Pgg0_reg_integrand(double z, void * p);
double C2m_g1_x_Pgg0_sing_integrand(double z, void * p) ;
double CLm_g1_x_Pgg0_reg_integrand(double z, void * p) ;
double CLm_g1_x_Pgg0_sing_integrand(double z, void * p);
double Pgg0sing_integrand(double z, void * p);
double C2m_g1_x_Pgg0(double x, double mQ, int nf);
double CLm_g1_x_Pgg0(double x, double mQ, int nf);


double C2m_g1_x_Pgq1_integrand(double z, void * p);
double C2m_g1_x_Pgq1(double x, double mQ, int nf) ;
double CLm_g1_x_Pgq1_integrand(double z, void * p);
double CLm_g1_x_Pgq1(double x, double mQ, int nf) ;

double C2m_g20_x_Pgq0_integrand(double z, void * p);
double C2m_g20_x_Pgq0(double x, double mQ);
double CLm_g20_x_Pgq0_integrand(double z, void * p);
double CLm_g20_x_Pgq0(double x, double mQ);

double C2m_ps20_x_Pqq0_reg_integrand(double z, void * p);
double C2m_ps20_x_Pqq0_sing_integrand(double z, void * p);
double CLm_ps20_x_Pqq0_reg_integrand(double z, void * p);
double CLm_ps20_x_Pqq0_sing_integrand(double z, void * p);
double Pqq0sing_integrand(double z, void * p);
double C2m_ps20_x_Pqq0(double x, double mQ, int nf);
double CLm_ps20_x_Pqq0(double x, double mQ, int nf) ;

double Pgg0_x_Pgq0(double x, int nf);
double Pqq0_x_Pgq0(double x) ;

double C2m_g1_x_Pgg0_x_Pgq0_integrand(double z, void * p);
double C2m_g1_x_Pgg0_x_Pgq0(double x, double mQ, int nf) ;
double C2m_g1_x_Pqq0_x_Pgq0_integrand(double z, void * p);
double C2m_g1_x_Pqq0_x_Pgq0(double x, double mQ, int nf);
double CLm_g1_x_Pgg0_x_Pgq0_integrand(double z, void * p);
double CLm_g1_x_Pgg0_x_Pgq0(double x, double mQ, int nf);
double CLm_g1_x_Pqq0_x_Pgq0_integrand(double z, void * p);
double CLm_g1_x_Pqq0_x_Pgq0(double x, double mQ, int nf);

double C2m_g1_x_Pgg1_reg_integrand(double z, void * p);
double C2m_g1_x_Pgg1_sing_integrand(double z, void * p);
double CLm_g1_x_Pgg1_reg_integrand(double z, void * p);
double CLm_g1_x_Pgg1_sing_integrand(double z, void * p);
double Pgg1sing_integrand(double z, void * p);
double C2m_g1_x_Pgg1(double x, double mQ, int nf) ;
double CLm_g1_x_Pgg1(double x, double mQ, int nf);

double C2m_ps20_x_Pqg0_integrand(double z, void * p);
double C2m_ps20_x_Pqg0(double x, double mQ, int nf);

double C2m_g20_x_Pgg0_reg_integrand(double z, void * p);
double C2m_g20_x_Pgg0_sing_integrand(double z, void * p);
double C2m_g20_x_Pgg0(double x, double mQ, int nf);

double CLm_ps20_x_Pqg0_integrand(double z, void * p);
double CLm_ps20_x_Pqg0(double x, double mQ, int nf) ;

double CLm_g20_x_Pgg0_reg_integrand(double z, void * p);
double CLm_g20_x_Pgg0_sing_integrand(double z, void * p);
double CLm_g20_x_Pgg0(double x, double mQ, int nf) ;

double Pqg0_x_Pgq0(double x, int nf) ;
double C2m_g1_x_Pqg0_x_Pgq0_integrand(double z, void * p);
double C2m_g1_x_Pqg0_x_Pgq0(double x, double mQ, int nf);
double CLm_g1_x_Pqg0_x_Pgq0_integrand(double z, void * p);
double CLm_g1_x_Pqg0_x_Pgq0(double x, double mQ, int nf);

double C2m_g1_x_Pgg0_x_Pgg0_reg_integrand(double z, void * p);
double C2m_g1_x_Pgg0_x_Pgg0_sing_integrand(double z, void * p);
double C2m_g1_x_Pgg0_x_Pgg0(double x, double mQ, int nf);

double CLm_g1_x_Pgg0_x_Pgg0_reg_integrand(double z, void * p);
double CLm_g1_x_Pgg0_x_Pgg0_sing_integrand(double z, void * p);
double CLm_g1_x_Pgg0_x_Pgg0(double x, double mQ, int nf);

double C2m_g1_x_Pgg0_x_Pgg0_reg1_integrand(double z[], size_t dim, void * p) ;
double C2m_g1_x_Pgg0_x_Pgg0_reg2_integrand(double z[], size_t dim, void * p) ;
double C2m_g1_x_Pgg0_x_Pgg0_reg3_integrand(double z[], size_t dim, void * p);

double C2m_g1_x_Pgg0_x_Pgg0_reg(double x, double mQ, int nf, size_t calls);

double C2m_g1_x_Pgg0_x_Pgg0_sing1_integrand(double z[], size_t dim, void * p) ;
double C2m_g1_x_Pgg0_x_Pgg0_sing2_integrand(double z[], size_t dim, void * p) ;
double C2m_g1_x_Pgg0_x_Pgg0_sing3_integrand(double z[], size_t dim, void * p);

double C2m_g1_x_Pgg0_x_Pgg0_sing(double x, double mQ, int nf, size_t calls);

double C2m_g1_x_Pgg0_x_Pgg0_MC(double x, double mQ, int nf, size_t calls) ;

#endif
