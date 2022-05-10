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

#endif
