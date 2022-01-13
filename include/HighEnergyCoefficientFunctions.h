#ifndef HighEnergy_h
#define HighEnergy_h

//NNLO => O(\alpha_s^2)
double C2m_g2_highenergy(double x, double mQ, double mMu);
double CLm_g2_highenergy(double x, double mQ, double mMu);
double C2m_ps2_highenergy(double x, double mQ, double mMu);
double CLm_ps2_highenergy(double x, double mQ, double mMu);

double C2m_g2_highenergy_highscale(double x, double mQ , double mMu); //cerca un altro nome
double CLm_g2_highenergy_highscale(double x, double mQ , double mMu);
double C2m_ps2_highenergy_highscale(double x, double mQ, double mMu);
double CLm_ps2_highenergy_highscale(double x, double mQ, double mMu);

double C2m_g2_power_terms(double x, double mQ , double mMu);
double CLm_g2_power_terms(double x, double mQ , double mMu);
double C2m_ps2_power_terms(double x, double mQ , double mMu);
double CLm_ps2_power_terms(double x, double mQ , double mMu);

//N3LO => O(\alpha_s^3)

double C2m_g3_highenergyNLL(double x, double mQ, double mMu,int nf);
double C2m_ps3_highenergyNLL(double x, double mQ, double mMu,int nf);

double C2m_g3_highenergy_highscaleNLL(double x, double mQ, double mMu, int nf);
double C2m_ps3_highenergy_highscaleNLL(double x, double mQ, double mMu, int nf);

double C2m_g3_highenergyNLL_ERR(double x, double mQ, double mMu, int nf);
double C2m_g3_highenergy_highscaleNLL_ERR(double x, double mQ, double mMu, int nf);

double C2m_g3_power_termsNLL(double x, double mQ , double mMu, int nf,int v=0);
double C2m_ps3_power_termsNLL(double x, double mQ , double mMu, int nf);

/////////////////////
/////////////////////

double C2m_g3_highenergy(double x, double mQ, double mMu);//term without the NLL approximation
double C2m_g3_highenergy_highscale(double x, double mQ , double mMu);
double C2m_g3_power_terms(double x, double mQ , double mMu);

double C2m_g3_highenergyNLL(double x, double mQ, double mMu, int nf, int v);
double C2m_g3_highenergy_highscaleNLL(double x, double mQ, double mMu, int nf, int v);

double C2m_g3_highenergyNLL_NEW(double x, double mQ, double mMu, int nf);
double C2m_g3_highenergyNLL_ERR_NEW(double x, double mQ, double mMu, int nf);

double C2m_g3_highenergy_highscaleNLL_NEW(double x, double mQ, double mMu, int nf);
double C2m_g3_highenergy_highscaleNLL_ERR_NEW(double x, double mQ, double mMu, int nf);

#endif
