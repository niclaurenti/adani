#ifndef Massive_h
#define Massive_h


//NLO => O(\alpha_s^1)
//massive coeff func of the gluon at first order in alpha_s^[nf] in the nf flavor scheme
double C2m_g1(double x, double mQ);
double CLm_g1(double x, double mQ);

//massive coeff func of the gluon at first order in alpha_s^[nf+1] in the nf flavor scheme
double D2m_g1_4(double x, double mQ);
double DLm_g1_4(double x, double mQ);
//massive coeff func of the gluon at first order in alpha_s^[nf+1] in the nf+1 flavor scheme
double D2m_g1_5(double x, double mQ);
double DLm_g1_5(double x, double mQ);

//NNLO => O(\alpha_s^2)
double C2m_g2(double x, double mQ, double mMu);//mMu=m^2/mu^2
double C2m_ps2(double x, double mQ, double mMu);

double CLm_g2(double x, double mQ, double mMu);
double CLm_ps2(double x, double mQ, double mMu);

double D2m_g2_4(double x, double mQ, double mMu);
double DLm_g2_4(double x, double mQ, double mMu);

double D2m_g2_5(double x, double mQ, double mMu);
double DLm_g2_5(double x, double mQ, double mMu);

double C2m_g21(double x, double mQ);
double CLm_g21(double x, double mQ);
double C2m_ps21(double x, double mQ);
double CLm_ps21(double x, double mQ);

double C2m_ps31(double x, double mQ, int nf);
double CLm_ps31(double x, double mQ, int nf);
double C2m_g31(double x, double mQ, int nf);
double CLm_g31(double x, double mQ, int nf);


double C2m_ps32(double x, double mQ, int nf) ;
double CLm_ps32(double x, double mQ, int nf) ;

double C2m_g31(double x, double mQ, int nf) ;
double CLm_g31(double x, double mQ, int nf) ;

double C2m_g32(double x, double mQ, int nf, int method_flag = 0, int calls = 500000);
double CLm_g32(double x, double mQ, int nf, int method_flag = 0, int calls = 500000) ;

#endif
