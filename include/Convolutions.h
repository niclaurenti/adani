#ifndef Convolutions_h
#define Convolutions_h

struct function_params {double x; double mQ; int nf;};

double C2_b1_x_K_bg1(double x, double mMu);
double CL_b1_x_K_bg1(double x, double mMu);

double K_bg1_x_K_gg2(double x, double mMu);

double C2_b2_x_K_bg1(double x, double mMu, int nf);

double C2m_g1_x_Pgq0(double z, void * p) ;
double CLm_g1_x_Pgq0(double z, void * p) ;
double C2m_g1_x_Pgg0_reg(double z, void * p);
double C2m_g1_x_Pgg0_sing(double z, void * p) ;
double CLm_g1_x_Pgg0_reg(double z, void * p) ;
double CLm_g1_x_Pgg0_sing(double z, void * p);

double Pgg0sing_int(double z, void * p);

double C2m_g1_x_Pgq1(double z, void * p);
double CLm_g1_x_Pgq1(double z, void * p);
double C2m_g20_x_Pgq0(double z, void * p);
double CLm_g20_x_Pgq0(double z, void * p);

double C2m_ps20_x_Pqq0_reg(double z, void * p);
double C2m_ps20_x_Pqq0_sing(double z, void * p);
double CLm_ps20_x_Pqq0_reg(double z, void * p);
double CLm_ps20_x_Pqq0_sing(double z, void * p);

double Pqq0sing_int(double z, void * p);

double Pgg0_x_Pgq0(double x, int nf);
double Pqq0_x_Pgq0(double x) ;

double C2m_g1_x_Pgg0_x_Pgq0(double z, void * p);
double C2m_g1_x_Pqq0_x_Pgq0(double z, void * p);
double CLm_g1_x_Pgg0_x_Pgq0(double z, void * p);
double CLm_g1_x_Pqq0_x_Pgq0(double z, void * p);

double C2m_g1_x_Pgg1_reg(double z, void * p);
double C2m_g1_x_Pgg1_sing(double z, void * p);
double CLm_g1_x_Pgg1_reg(double z, void * p);
double CLm_g1_x_Pgg1_sing(double z, void * p);

double Pgg1sing_int(double z, void * p);

double C2m_ps21_x_Pqg0(double z, void * p);

double C2m_g21_x_Pgg0_reg(double z, void * p);
double C2m_g21_x_Pgg0_sing(double z, void * p);

double CLm_ps21_x_Pqg0(double z, void * p);

double CLm_g21_x_Pgg0_reg(double z, void * p);
double CLm_g21_x_Pgg0_sing(double z, void * p);

double Pqg0_x_Pgq0(double x, int nf) ;

double C2m_g1_x_Pgg0(double x, double mQ, int nf);

double C2m_g1_x_Pgg0_Pgg0_reg(double z, void * p);
double C2m_g1_x_Pgg0_Pgg0_sing(double z, void * p);

double CLm_g1_x_Pgg0(double x, double mQ, int nf);

double CLm_g1_x_Pgg0_Pgg0_reg(double z, void * p);
double CLm_g1_x_Pgg0_Pgg0_sing(double z, void * p);

#endif
