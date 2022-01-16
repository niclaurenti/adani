#ifndef Asymptotic_h
#define Asymptotic_h

double C2m_g1_asymptotic(double x, double mQ);
double C2m_g2_asymptotic(double x, double mQ, double mMu);
double C2m_ps2_asymptotic(double x, double mQ, double mMu);
//double C2m_g2_asymptoticV2(double x, double mQ, double mMu);

double CLm_g1_asymptotic(double x, double mQ);
double CLm_g2_asymptotic(double x, double mQ, double mMu);
//double CLm_g2_asymptoticV2(double x, double mQ, double mMu);

double C2m_g3_asymptotic(double x, double mQ, double mMu, int nf, int v);
double C2m_g3_asymptoticNLL(double x, double mQ, double mMu, int nf, int v1=0, int v2=0);
double C2m_ps3_asymptoticNLL(double x, double mQ, double mMu, int nf);

#endif
