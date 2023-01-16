#include "adani/SpecialFunctions.h"
#include "adani/ColorFactors.h"

//#include<gsl/gsl_sf.h>

#include <iostream>
#include <cmath>

using namespace std;

//______________________________________________________

//const double ZETA2  =  1.64493406684822643647;
/*
const double ZETA2 = M_PI*M_PI/6;
const double ZETA3  =  1.2020569031595942855;

//const double ZETA4  =  1.0823232337111381915;
const double ZETA4  = M_PI*M_PI*M_PI*M_PI/90;
const double ZETA5  =  1.03692775514336992633136548646;
*/
//______________________________________________________________

double beta(int ord, int nf) {
    if(ord == 0) return (11. /3 * CA - 2. /3 * nf) / 4. / M_PI ;
    if(ord == 1) {
        double TF = TR * nf ;
        double b_ca2 = 34. / 3. * CA * CA ;
        double b_ca = -20. / 3. * CA * TF ;
        double b_cf = -4. * CF * TF ;
        double beta_1 = b_ca2 + b_ca + b_cf ;
        return beta_1 / (16. * M_PI * M_PI);
    }
    else {
        cout << "beta("<<ord<<") is not implemented"<<endl;
        exit(-1);
    }
}

//_____________________________________________________________

double zeta(int i) {
    
    double pi2 = M_PI*M_PI;
    double pi4=pi2*pi2;
    double pi6=pi4*pi2;
    
    if(i==2) return pi2/6;
    if(i==3) return 1.2020569031595942855;
    if(i==4) return pi4/90;
    if(i==5) return 1.03692775514336992633136548646;
    if(i==6) return pi6/945;
    
    else {
        cout    << "zeta(<<"<<i<<") is not implemented! \n Exiting..." 
                << endl;
        exit(-1);
    }

}

//___________________________________________________________

double theta(double x) {

    return 0.5 * (1. + x / fabs(x)) ;
    
}

//___________________________________________________________

double Li2(double x) {
    
    double ZETA2=zeta(2);
    double x_0 = -0.30;
    double x_1 = 0.25;
    double x_2 = 0.51;
    
    if (x == 1.) return ZETA2;
    
    if (x <= x_0) { 
        double temp = log(fabs(1.0-x));
        return -Li2(-x/(1.0-x)) - temp*temp/2 ; 
    }
  
    else if (x < x_1) {
        double z = - log(1.0-x);
        double temp = (
            z*(1.0-z/4.0*(1.0-z/9.0*(1.0-z*z/100.0
            *(1.0-5.0*z*z/294.0*(1.0-7.0*z*z/360.0
            *(1.0-5.0*z*z/242.0*(1.0-7601.0*z*z/354900.0
            *(1.0-91.0*z*z/4146.0*(1.0-3617.0*z*z/161840.0)
            ))))))))
        ) ;
        return temp; 
    }
    
    else if (x < x_2) return - Li2(-x) + Li2(x*x)/2.0 ;
  
    else {return ZETA2 - Li2(1.0-x) - log(fabs(x))*log(fabs(1.0-x)) ;} ; 
}

//___________________________________________________________

double Li3(double x) {

    double ZETA2=zeta(2);
    double ZETA3=zeta(3);
    double x_0 = -1.0;
    double x_1 = -0.85;
    double x_2 = 0.25;
    double x_3 = 0.63;
    double x_4 =  1.0;
    
    if (x == 1.) return ZETA3;
    
    if (x == -1.) return - 0.75 * ZETA3;
    
    if (x <= x_0){ 
        double lnx = log(-x);
        return Li3(1.0/x) - ZETA2*lnx - lnx*lnx*lnx/6.0; 
    }
  
    else if (x < x_1){
    return Li3(x*x)/4.0 - Li3(-x); }
  
    else if (x < x_2){
        double z = - log(1.0-x);
        double temp = (
            z*(1.0-3.0*z/8.0*(1.0-17.0*z/81.0*(1.0-15*z/136.0
            *(1.0-28.0*z/1875.0*(1.0+5.0*z/8.0*(1.0-304.0*z/7203.0
            *(1.0+945.0*z/2432.0*(1.0-44.0*z/675.0*(1.0+7.0*z/24.0
            *(1.0-26104.0*z/307461.0*(1.0+1925.0*z/8023.0
            *(1.0-53598548.0*z/524808375.0
            *(1.0+22232925.0*z/107197096.0
            )))))))))))))
        ) ;
        return temp;
    }
    
    else if (x < x_3) {
        return Li3(x*x)/4.0 - Li3(-x);
    }
  
    else if (x < x_4) {
        double ln1x = log(1.0-x); 
        return (
            -Li3(1.0-x) - Li3(-x/(1.0-x)) + ZETA3 + ZETA2*ln1x
            - log(x)*ln1x*ln1x/2.0 + ln1x*ln1x*ln1x/6.0
        ) ;
    }
       
    else { 
        double lnx = log(x);
        return Li3(1./x) + 2.0*ZETA2*lnx - lnx*lnx*lnx/6.0;
    }
  
}

//______________________________________________________________

double Li4(double x) {
    
    cout<< "Li4(x) is not implemented yet!!\n Exiting..."<<endl;
    exit(-1);

}

//_______________________________________________________________

double S12(double x){
    
    double ZETA3=zeta(3);
    double c;
    if (x > 1) { 
        return (
            - Li3(1.-x) + ZETA3 + log(x-1.) * Li2(1.-x) 
            + 0.5 * log(x) * ( log(x-1.)*log(x-1.) - M_PI*M_PI )
        ) ;
    }
           
    else if (x==1){ return ZETA3; }
  
    else if ((0 < x) && (x < 1)){ 
        return (
            - Li3(1.-x) + ZETA3 + log(1.-x) * Li2(1.-x) 
            + 0.5 * log(x) * log(1.-x)*log(1.-x)
        ); 
    }
           
    else if (x==0){ return 0.; }
  
    else if (x<0){ 
        c = 1./(1.-x);
        return (
            - Li3(c) + ZETA3 + log(c) * Li2(c) 
            + 0.5 * log(c)*log(c) * log(1.-c) 
            - 1./6. * log(c)*log(c)*log(c)
        ); 
    }
           
    else { return 0.; }
  
}

//_______________________________________________________

double c0(double xi) {

    double y=sqrt(1+4/xi);

    double L1=log(1+xi/2);
    double L2=log(2+xi/2);  
    double L3=log(sqrt(xi)*(y-1)/2);

    double xp2=2+xi;
    double xp4=4+xi;

    double Li_2=Li2(-2/xp2);
    double z2=zeta(2);
    double pi2=M_PI*M_PI;

    double c_CA= (
        50 - pi2 + 12*L3/y + 4*L3*L3 + L1*L1 + 6*L2 - 
        4*L2*L2 + 2*Li_2 + 48/xp2 - 4*L2/xp2 + 64*L2/xp2/xp2 - 
        128*L2/(xp2*xp2*xp4) - 160/xp2/xp4 - 64*L2/xp2/xp4 + 
        128/(xp2*xp4*xp4) - 12*(4+z2)/xp4 - 8*L3*L3/xp4 + 64/xp4/xp4
    );

    double c_CF= (
        -18 - 2./3.*pi2 - 24*L3/y - 8*L3*L3 + 2*L1*L1 - 6*L2 + 
        4*Li_2 - 48/xp2 + 8*L2/xp2 + 360/xp2/xp4 + 128*L2/xp2/xp4 - 
        544/(xp2*xp4*xp4) + 48*L3*L3/xp4 - 8*L1*L1/xp4 + 
        (44+40*z2)/xp4 - 120*L2/xp2/xp2 + 256*L2/(xp2*xp2*xp4) - 
        16*Li_2/xp4 - 272/xp4/xp4
    ) ;

  return CA*c_CA + CF*c_CF ;

}

//________________________________________________________

double c0_bar(double xi) {

    return 4*CA*(2+log(1+xi/4))-4./3.*TR;
  
}

//________________________________________________________

//Harmonic Polylogarithms

double H(double x, int i) {
    
    if(i==1) return -log(1-x);
    if(i==0) return log(x);
    if(i==-1) return log(1+x);
    
    else {
        
        cout << "H(x,i) is defined only for i = -1,0,1! \n Exiting..." << endl;
        exit(-1);
    }
    
}

//_______________________________________________________

double H(double x, int i, int j) {
    
    double xp=1+x;
    double xm=1-x;
    
    double L=log(x);
    double L2=L*L;
    double Lp=log(xp);
    double Lp2=Lp*Lp;
    double Lm=log(xm);
    double Lm2=Lm*Lm;
    
    if(i==0 && j==0) return 0.5*L2;
    if(i==0 && j==1) return Li2(x);
    if(i==0 && j==-1) return -Li2(-x);
    if(i==1 && j==0) return -L*Lm - Li2(x);
    if(i==1 && j==1) return 0.5*Lm2;
    if(i==1 && j==-1) return Li2(xm/2) - log(2)*Lm - Li2(0.5);
    if(i==-1 && j==0) return L*Lp + Li2(-x);
    if(i==-1 && j==1) return Li2(xp/2) - log(2)*Lp - Li2(0.5);
    if(i==-1 && j==-1) return 0.5*Lp2;
    
    else  {
        cout    << "H(x,"<<i<<","<<j<<") is not defined."<<endl
                << "Exiting..."<<endl;
        exit(-1);
    }
    
}
    
//____________________________________________________________


double H(double x, int i, int j, int k) {
    
    double xp=1+x;
    double xm=1-x;
    
    double L=log(x);
    double L2=L*L;
    double L3=L2*L;
    
    double Lp=log(xp);
    double Lp2=Lp*Lp;
    double Lp3=Lp2*Lp;
    
    double Lm=log(xm);
    double Lm2=Lm*Lm;
    double Lm3=Lm2*Lm;
    
    double pi2=M_PI*M_PI;
    double z3=zeta(3);
    
    //double l=log(2);
    //double l2=l*l;
    //double l3=l2*l;
    
    
    if(i==1 && j==1 && k==1) return -1./6.*Lm3;
    
    if(i==1 && j==1 && k==0) return 1./6*pi2*Lm - Li3(xm) + z3;
    
    if(i==1 && j==1 && k==-1) return 0.5*H(x,1)*H(x,1,-1) - 0.5*H(x,1,-1,1);
    
    if(i==1 && j==0 && k==1) return -2*S12(x) - Lm*Li2(x);
    
    if(i==1 && j==0 && k==0) return -1./2*L2*Lm - L*Li2(x) + Li3(x);
    
    if(i==1 && j==0 && k==-1) return H(x,1)*H(x,0,-1) - H(x,0,1,-1) - H(x,0,-1,1);
    
    if(i==1 && j==-1 && k==1) return H(x,1)*H(x,-1,1) - 2*H(x,-1,1,1);
    
    if(i==1 && j==-1 && k==0) return H(x,0)*H(x,1,-1) - H(x,0,1,-1) - H(x,1,0,-1);
    
    if(i==1 && j==-1 && k==-1) return -0.5*log(xm/2)*Lp2 - Li3(0.5) - Lp*Li2(xp/2)+ Li3(xp/2);
    
    
    if(i==0 && j==1 && k==1) return S12(x);
    
    if(i==0 && j==1 && k==0) return L*Li2(x) - 2*Li3(x);
    
    if(i==0 && j==1 && k==-1) return Li3(2*x/xp) - Li3(x/xp) - Li3(xp/2) - Li3(x) + Lp*Li2(0.5) + Lp*Li2(x) + 0.5*log(2)*Lp2+ Li3(0.5);
    
    if(i==0 && j==0 && k==1) return Li3(x);
    
    if(i==0 && j==0 && k==0) return 1./6.*L3;
    
    if(i==0 && j==0 && k==-1) return -Li3(-x);
    
    if(i==0 && j==-1 && k==1) return -S12(x) + Li3(-2*x/xm) - Li3(xm/2) - Li3(-x) + Li3(0.5) + Li3(x) + Lm*Li2(-x) + Lm*Li2(0.5) - Lm*Li2(x)    + 0.5*log(2)*Lm2 - 1./6.*Lm3;
    
    if(i==0 && j==-1 && k==0) return -L*Li2(-x) + 2*Li3(-x);
    
    if(i==0 && j==-1 && k==-1) return S12(-x);
    
    
    if(i==-1 && j==1 && k==1) return 0.5*log(xp/2)*Lm2 + Li3(0.5) + Lm*Li2(xm/2) - Li3(xm/2);
    
    if(i==-1 && j==1 && k==0) return H(x,0)*H(x,-1,1) - H(x,0,-1,1) - H(x,-1,0,1);
    
    if(i==-1 && j==1 && k==-1) return H(x,-1)*H(x,1,-1) - 2*H(x,1,-1,-1);
    
    if(i==-1 && j==0 && k==1) return H(x,-1)*H(x,0,1) - H(x,0,-1,1) - H(x,0,1,-1);
    
    if(i==-1 && j==0 && k==0) return 0.5*L2*Lp + L*Li2(-x) - Li3(-x);
    
    if(i==-1 && j==0 && k==-1) return H(x,0,-1)*H(x,-1) - 2*H(x,0,-1,-1);
    
    if(i==-1 && j==-1 && k==1) return 0.5*H(x,-1)*H(x,-1,1) - 0.5*H(x,-1,1,-1);
        
    if(i==-1 && j==-1 && k==0) return H(x,0)*H(x,-1,-1) - H(x,0,-1,-1) - H(x,-1,0,-1);
    
    if(i==-1 && j==-1 && k==-1) return 1./6.*Lp3;
        
    
    else  {
        cout    << "H(x,"<<i<<","<<j<<","<<k<<") is not defined."<<endl
                << "Exiting..."<<endl;
        exit(-1);
    }
    
        
}    

//__________________________________________________________