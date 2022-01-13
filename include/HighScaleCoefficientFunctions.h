#ifndef HighScale_h
#define HighScale_h

//TODO: change name to HighScaleCoefficientFunctions

// C^{(k)} = k-th order expansion in terms of \alpha_s^{[nf]}
// D^{(k)} = k-th order expansion in terms of \alpha_s^{[nf+1]}

//NLO => O(\alpha_s^1)
double C2m_g1_highscale(double x, double mQ);
double CLm_g1_highscale(double x, double mQ);

double D2m_g1_highscale(double x, double mQ);
double DLm_g1_highscale(double x, double mQ);

//NNLO => O(\alpha_s^2)
double C2m_g2_highscale(double x, double mQ, double mMu);
double CLm_g2_highscale(double x, double mQ, double mMu);

double C2m_ps2_highscale(double x, double mQ, double mMu);
double CLm_ps2_highscale(double x, double mQ, double mMu);

double D2m_g2_highscale(double x, double mQ, double mMu);
double DLm_g2_highscale(double x, double mQ, double mMu);

double D2m_ps2_highscale(double x, double mQ, double mMu);
double DLm_ps2_highscale(double x, double mQ, double mMu);


//double C2m0_g2(double x, double mQ);
//double CLm0_g2(double x, double mQ);
//double D2m0_g2(double x, double mQ);
//double DLm0_g2(double x, double mQ);


//N3LO => O(\alpha_s^3)
double C2m_g3_highscale(double x, double mQ, double mMu, int nf, int v);
double D2m_g3_highscale(double x, double mQ, double mMu, int nf, int v);

double D2m_ps3_highscale(double x, double mQ, double mMu, int nf);
double C2m_ps3_highscale(double x, double mQ, double mMu, int nf);

double CLm_g3_highscale(double x, double mQ, double mMu, int nf); //TODO: check if it is really C and not D

double CLm_ps3_highscale(double x, double mQ, double mMu, int nf); //TODO: check if it is really C and not D




//double D2m0_g3(double x, double mQ, int nf);
//double C2m0_g3(double x, double mQ, int nf);

#endif
