#include "../include/ThresholdCoefficientFunctions.h"
#include "../include/MassiveCoefficientFunctions.h"
#include "../include/SpecialFunctions.h"
#include "../include/ColorFactors.h"
#include <cmath>
//#include<gsl/gsl_sf.h>
#include "../include/SpecialFunctions.h"
#include "apfel/massivecoefficientfunctionsunp_sl.h"

//#include <iostream>

//using namespace std;

using namespace apfel;

//_____________________________________________________________

double C2m_g1_threshold(double x, double mQ) {

  double xmax=1/(1+4*mQ);

  if(x>xmax || x<0) return 0;
  
  double beta=sqrt(1-4*mQ*x/(1-x));
  double xi=1./mQ;
  double overall=xi/(4*M_PI*M_PI);
  
  return overall*M_PI*TR*beta/(1+xi/4)/x;
  
}//in order to pass to Vogt normalization multiply mQ*4*M_PI*M_PI*x 

//__________________________________________________________

double C2m_g2_threshold(double x, double mQ, double mMu) {

  double xmax=1/(1+4*mQ);

  if(x>xmax || x<0) return 0;
  	
  double beta=sqrt(1-4*mQ*x/(1-x));
  double xi=1./mQ;
    //double overall=xi/(4*M_PI*M_PI);

  double logb=log(beta);
  double log2b=logb*logb;
    
  double C_log2b=16*CA;
    
  double C_logb=48*CA*log(2) - 40*CA 
    + 8*CA*log(mMu);
    
  double C_fracb=(2*CF-CA)*M_PI*M_PI;

  double C_const= c0(xi)+36*CA*log(2)*log(2)-60*CA*log(2)
    +log(mMu)*(8*CA*log(2)-c0_bar(xi));

  return C2m_g1(x,mQ)/4/M_PI*(C_log2b*log2b + C_logb*logb + C_fracb/beta + C_const);

}//in order to pass to Vogt normalization multiply mQ*M_PI*x and put mu^2=Q^2+4m^2

//_________________________________________________________

double CLm_g2_threshold(double x, double mQ, double mMu) {

  double xmax=1/(1+4*mQ);

  if(x>xmax || x<0) return 0;
  	
  double beta=sqrt(1-4*mQ*x/(1-x));
  double xi=1./mQ;
    //double overall=xi/(4*M_PI*M_PI);

  double logb=log(beta);
  double log2b=logb*logb;
    
  double C_log2b=16*CA;
    
  double C_logb=48*CA*log(2) - 40*CA 
    + 8*CA*log(mMu);
    
  double C_fracb=(2*CF-CA)*M_PI*M_PI;

  double C_const= c0(xi)+36*CA*log(2)*log(2)-60*CA*log(2)
    +log(mMu)*(8*CA*log(2)-c0_bar(xi));

  return CLm_g1(x,mQ)/4/M_PI*(C_log2b*log2b + C_logb*logb + C_fracb/beta + C_const);

}//in order to pass to Vogt normalization multiply mQ*M_PI*x and put mu^2=Q^2+4m^2

//_________________________________________________________

double C2m_g3_threshold(double x, double mQ, double mMu, int nf) {

  double x_max=1./(1+4*mQ);

  if(x>x_max || x<0) return 0;

  double xi=1/mQ;  
  double beta=sqrt(1-4*mQ*x/(1-x));
  //double eta=1./4.*xi*(1-x)/x - 1;
  //double rho=4*mQ*x/(1-x);

  // double L=log(mQ);
  double Lm=log(mMu);
  //double L2=L*L;    
  double Lm2=Lm*Lm; 
  double l=log(beta);
  double l2=l*l;
  double l3=l2*l;
  double l4=l3*l;
  
  double log2=log(2);
  
  double z3=zeta(3);
  double pi2=M_PI*M_PI;
  double pi4=pi2*pi2;

  //double overall=xi/(4*pi2);

  double c_log4= 128.*CA*CA;

  double c_log3=(768.*log2 - 6464./9.)*CA*CA + 128./9.*CA*nf + 128.*CA*CA*Lm;

  double c_log2=
  	(1728.*log2*log2 - 3232.*log2 - 208./3.*pi2 + 15520./9.)*CA*CA + (64.*log2 - 640./9.)*CA*nf + 16.*CA*c0(xi) + 
  	32.*CA*(CF - CA/2)*pi2/beta - 
  	((-512*log2 + 1136./3.)*CA*CA - 32./3.*CA*nf + 16*CA*c0_bar(xi))*Lm+ 32*CA*CA*Lm2;

  double c_log_const=
  	(1728.*log2*log2*log2 - 4848*log2*log2 + 15520./3.*log2 - 208*pi2*log2 + 936*z3 + 608./3.*pi2 - 88856./27.)*CA*CA + 
  	(96*log2*log2 - 640./3.*log2 - 16./3.*pi2 + 4592./27.)*CA*nf - 
  	32*CF*(CF - CA/2)*pi2 + (48*log2 - 40)*CA*c0(xi);

  double c_log_fracbeta=((-92./3. + 32*log2)*CA + 8./3.*nf)*(CF - CA/2)*pi2;

  double c_log_Lm= - ((-672*log2*log2 + 976*log2 + 104./3.*pi2 - 4160./9.)*CA*CA + (-32*log2 + 320./9.)*CA*nf + (48*log2 - 40)*CA*c0_bar(xi) - 8*CA*c0(xi) - 16*CA*(CF - CA/2)*pi2/beta);

  double c_log_Lm2=((64*log2 - 44./3.)*CA*CA + 8./3.*CA*nf - 8*CA*c0_bar(xi));

  double c_log = c_log_const + c_log_fracbeta/beta + c_log_Lm*Lm + c_log_Lm2*Lm2;

  double c_fracbeta=((8*log2*log2 - 68./3.*log2 + 8./3.*pi2 - 658./9.)*CA + (8./3.*log2 - 20./9.)*nf + 2*c0(xi) + (26./3.*CA + 4./3.*nf - 2*c0_bar(xi))*Lm)*(CF - CA/2)*pi2;

  double c_fracbeta2= 4./3.*(CF - CA/2)*(CF - CA/2)*pi4;

  double c_const_sqrt=c0(xi) + 36*CA*log2*log2 - 60*CA*log2 + Lm*(8*CA*log2 - c0_bar(xi));

  double c_const= c_const_sqrt*c_const_sqrt;

  return C2m_g1(x,mQ)/pi2/16.*(c_log4*l4 + c_log3*l3 + c_log2*l2 + c_log*l + c_fracbeta/beta + c_fracbeta2/beta/beta + c_const);

}

//_________________________________________________________

double CLm_g3_threshold(double x, double mQ, double mMu, int nf) {

  double x_max=1./(1+4*mQ);

  if(x>x_max || x<0) return 0;

  double xi=1/mQ;  
  double beta=sqrt(1-4*mQ*x/(1-x));
  //double eta=1./4.*xi*(1-x)/x - 1;
  //double rho=4*mQ*x/(1-x);

  // double L=log(mQ);
  double Lm=log(mMu);
  //double L2=L*L;    
  double Lm2=Lm*Lm; 
  double l=log(beta);
  double l2=l*l;
  double l3=l2*l;
  double l4=l3*l;
  
  double log2=log(2);
  
  double z3=zeta(3);
  double pi2=M_PI*M_PI;
  double pi4=pi2*pi2;

  //double overall=xi/(4*pi2);

  double c_log4= 128.*CA*CA;

  double c_log3=(768.*log2 - 6464./9.)*CA*CA + 128./9.*CA*nf + 128.*CA*CA*Lm;

  double c_log2=
  	(1728.*log2*log2 - 3232.*log2 - 208./3.*pi2 + 15520./9.)*CA*CA + (64.*log2 - 640./9.)*CA*nf + 16.*CA*c0(xi) + 
  	32.*CA*(CF - CA/2)*pi2/beta - 
  	((-512*log2 + 1136./3.)*CA*CA - 32./3.*CA*nf + 16*CA*c0_bar(xi))*Lm+ 32*CA*CA*Lm2;

  double c_log_const=
  	(1728.*log2*log2*log2 - 4848*log2*log2 + 15520./3.*log2 - 208*pi2*log2 + 936*z3 + 608./3.*pi2 - 88856./27.)*CA*CA + 
  	(96*log2*log2 - 640./3.*log2 - 16./3.*pi2 + 4592./27.)*CA*nf - 
  	32*CF*(CF - CA/2)*pi2 + (48*log2 - 40)*CA*c0(xi);

  double c_log_fracbeta=((-92./3. + 32*log2)*CA + 8./3.*nf)*(CF - CA/2)*pi2;

  double c_log_Lm= - ((-672*log2*log2 + 976*log2 + 104./3.*pi2 - 4160./9.)*CA*CA + (-32*log2 + 320./9.)*CA*nf + (48*log2 - 40)*CA*c0_bar(xi) - 8*CA*c0(xi) - 16*CA*(CF - CA/2)*pi2/beta);

  double c_log_Lm2=((64*log2 - 44./3.)*CA*CA + 8./3.*CA*nf - 8*CA*c0_bar(xi));

  double c_log = c_log_const + c_log_fracbeta/beta + c_log_Lm*Lm + c_log_Lm2*Lm2;

  double c_fracbeta=((8*log2*log2 - 68./3.*log2 + 8./3.*pi2 - 658./9.)*CA + (8./3.*log2 - 20./9.)*nf + 2*c0(xi) + (26./3.*CA + 4./3.*nf - 2*c0_bar(xi))*Lm)*(CF - CA/2)*pi2;

  double c_fracbeta2= 4./3.*(CF - CA/2)*(CF - CA/2)*pi4;

  double c_const_sqrt=c0(xi) + 36*CA*log2*log2 - 60*CA*log2 + Lm*(8*CA*log2 - c0_bar(xi));

  double c_const= c_const_sqrt*c_const_sqrt;

  return CLm_g1(x,mQ)/pi2/16.*(c_log4*l4 + c_log3*l3 + c_log2*l2 + c_log*l + c_fracbeta/beta + c_fracbeta2/beta/beta + c_const);

}

