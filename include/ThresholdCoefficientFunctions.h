#ifndef Threshold_h
#define Threshold_h

//NLO => O(\alpha_s^1)
double C2m_g1_threshold(double x, double mQ);

//NNLO => O(\alpha_s^2)
double C2m_g2_threshold(double x, double mQ, double mMu);
double CLm_g2_threshold(double x, double mQ, double mMu);

//N3LO => O(\alpha_s^3)
double C2m_g3_threshold(double x, double mQ, double mMu, int nf);

#endif
