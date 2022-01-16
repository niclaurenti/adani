#ifndef Combination_h
#define Combination_h

//approximated coefficient funtions O(alpha_s)

double C2m_g1_approximation(double x, double mQ, double k, double h);
double C2m_g1_approximation(double x, double mQ);
double CLm_g1_approximation(double x, double mQ);


//approximated coefficient funtions O(alpha_s^2)

double C2m_g2_approximation(double x, double mQ, double mMu);
double C2m_g2_approximation(double x, double mQ, double mMu, double A, double B, double C, double D, double a, double b);
double C2m_g2_approximation_BAND(double x, double mQ, double mMu, double var, double fact, int v);

double C2m_ps2_approximation(double x, double mQ, double mMu);
double C2m_ps2_approximation(double x, double mQ, double mMu, double A, double B, double C, double D, double a, double b);

//approximated coefficient funtions O(alpha_s^3)

double C2m_g30_approximation(double x, double mQ, double mMu, int nf);
double C2m_ps30_approximation(double x, double mQ, double mMu, int nf);
double C2m_g30_approximation(double x, double mQ, double mMu, int nf, double A, double B, double C, double D, double a, double b, int v1,int v2);
double C2m_g30_approximation_BAND(double x, double mQ, double mMu, int nf, double var, double fact, int v);

//approximated coefficient funtions O(alpha_s^2) from [arXiv:1205.5727] 

double C2m_g2_approximationA_vogt(double x, double mQ, double mMu);
double C2m_g2_approximationB_vogt(double x, double mQ, double mMu);

//approximated coefficient funtions O(alpha_s^3) from [arXiv:1205.5727] 

double C2m_g30_approximationA_vogt(double x, double mQ, double mMu, int nf);
double C2m_g30_approximationB_vogt(double x, double mQ, double mMu, int nf);
double C2m_g30_approximationBlowxi_vogt(double x, double mQ, double mMu, int nf);

#endif

