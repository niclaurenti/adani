#ifndef HighEnergy_h
#define HighEnergy_h

//==========================================================================================//
//                      High energy coefficient funtions O(alpha_s^2)
//------------------------------------------------------------------------------------------//

double C2m_g2_highenergy(double x, double mQ, double mMu);
double C2m_ps2_highenergy(double x, double mQ, double mMu);
double CLm_g2_highenergy(double x, double mQ, double mMu);
double CLm_ps2_highenergy(double x, double mQ, double mMu);

//==========================================================================================//
//                      Q>>m limit of the high energy coefficient funtions O(alpha_s^2)
//------------------------------------------------------------------------------------------//

double C2m_g2_highenergy_highscale(double x, double mQ, double mMu);
double C2m_ps2_highenergy_highscale(double x, double mQ, double mMu);
double CLm_g2_highenergy_highscale(double x, double mQ, double mMu);
double CLm_ps2_highenergy_highscale(double x, double mQ, double mMu);

double C2m_g2_power_terms(double x, double mQ, double mMu);
double C2m_ps2_power_terms(double x, double mQ, double mMu);
double CLm_g2_power_terms(double x, double mQ, double mMu);
double CLm_ps2_power_terms(double x, double mQ, double mMu);

//==========================================================================================//
//                      High energy coefficient funtions O(alpha_s^3)
//------------------------------------------------------------------------------------------//

double C2m_g3_highenergyLL(double x, double mQ, double mMu);
double C2m_g3_highenergyNLL(double x, double mQ, double mMu,int nf);//NLL=LL+NLL TODO: change it
double C2m_ps3_highenergyLL(double x, double mQ, double mMu);
double C2m_ps3_highenergyNLL(double x, double mQ, double mMu,int nf);
double CLm_g3_highenergyNLL(double x, double mQ, double mMu, int nf);
double CLm_ps3_highenergyNLL(double x, double mQ, double mMu, int nf);

//==========================================================================================//
//                      Q>>m limit of the high energy coefficient funtions O(alpha_s^3)
//------------------------------------------------------------------------------------------//

double C2m_g3_highenergy_highscaleLL(double x, double mQ, double mMu);
double C2m_g3_highenergy_highscaleNLL(double x, double mQ, double mMu, int nf);
double C2m_ps3_highenergy_highscaleLL(double x, double mQ, double mMu);
double C2m_ps3_highenergy_highscaleNLL(double x, double mQ, double mMu, int nf);
double CLm_g3_highenergy_highscaleNLL(double x, double mQ, double mMu, int nf);
double CLm_ps3_highenergy_highscaleNLL(double x, double mQ, double mMu, int nf);

double C2m_g3_highenergyNLL_ERR(double x, double mQ, double mMu, int nf);
double C2m_g3_highenergy_highscaleNLL_ERR(double x, double mQ, double mMu, int nf);

double C2m_g3_power_termsLL(double x, double mQ , double mMu);
double C2m_g3_power_termsNLL(double x, double mQ , double mMu, int nf,int v=0);
double C2m_ps3_power_termsLL(double x, double mQ , double mMu);
double C2m_ps3_power_termsNLL(double x, double mQ , double mMu, int nf);
double CLm_g3_power_termsNLL(double x, double mQ , double mMu, int nf);
double CLm_ps3_power_termsNLL(double x, double mQ , double mMu, int nf);

/////////////////////
/////////////////////

//term with band
double C2m_g3_highenergyNLL(double x, double mQ, double mMu, int nf, int v);
double C2m_g3_highenergy_highscaleNLL(double x, double mQ, double mMu, int nf, int v);



#endif
