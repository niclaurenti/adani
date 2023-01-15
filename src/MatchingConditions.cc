#include "masterthesis/MatchingConditions.h"
#include "masterthesis/ColorFactors.h"
#include "masterthesis/SpecialFunctions.h"
#include "apfel/matchingconditions_sl.h"
#include<cmath>
#include <iostream>

using namespace apfel;
using namespace std;

double K_Qg1(double x, double mMu) {

    return 2 * TR * (
        x * x + (x - 1) * (x - 1)
    ) * log(1./mMu) / 4. / M_PI ;
}


//___________________________________________________________

double K_Qg2(double x, double mMu) {
    
    double x2=x*x;
    
    double L=log(x);
    double L2=L*L;
    double L3=L2*L;
    
    double Lm=log(1-x);
    double Lm2=Lm*Lm;
    double Lm3=Lm2*Lm;
    
    double Lp=log(1+x);
    double Lp2=Lp*Lp;
    
    double Li2xm=Li2(1-x);
    double Li3xm=Li3(1-x);
    double Li2minus=Li2(-x);
    double Li3minus=Li3(-x);
    double S12xm=S12(1-x);
    double S12minus=S12(-x);
    
    double z2=zeta(2);
    double z3=zeta(3);
    
    double pi2=M_PI*M_PI;
    
    double Lmu=log(mMu);
    double Lmu2=Lmu*Lmu;
    
    double logmu2_CFTR = Lm*(8 - 16*x + 16*x2) - L*(4 - 8*x + 16*x2) - (2 - 8*x);
    
    double logmu2_CATR = -Lm*(8 - 16*x + 16*x2) - (8 + 32*x)*L - 16./3./x - 4 - 32*x + 124./3.*x2;
    
    double logmu2_TR2 = -16./3.*(2*x2 - 2*x + 1);
    
    double logmu2 = CF*TR*logmu2_CFTR + CA*TR*logmu2_CATR + TR*TR*logmu2_TR2;
    
    double logmu_CFTR = (
        (8 - 16*x + 16*x2)*(2*L*Lm - Lm2 + 2*z2) - 
        (4 - 8*x + 16*x2)*L2 - 32*x*(1-x)*Lm - 
        (12 - 16*x + 32*x2)*L - 56 + 116*x - 80*x2
    );
    
    double logmu_CATR = (
        (16 + 32*x + 32*x2)*(Li2minus + L*Lp) + 
        (8 - 16*x + 16*x2)*Lm2 + (8 + 16*x)*L2 + 32*x*z2 + 
        32*x*(1-x)*Lm - (8 + 64*x + 352./3.*x2)*L - 
        160./9./x + 16 -200*x + 1744./9.*x2
    );
    
    double logmu = CF*TR*logmu_CFTR + CA*TR*logmu_CATR;
    
    double const_CFTR = (
        (1 - 2*x + 2*x2)*(8*z3 + 4./3.*Lm3 - 8*Lm*Li2xm + 8*z2*L - 4*L*Lm2 + 2./3.*L3 - 8*L*Li2xm + 8*Li3xm - 24*S12xm) + 
        x2*(-16*z2*L + 4./3.*L3 + 16*L*Li2xm + 32*S12xm) -
        (4 + 96*x - 64*x2)*Li2xm - (4 - 48*x + 40*x2)*z2 -
        (8 + 48*x - 24*x2)*L*Lm + (4 + 8*x - 12*x2)*Lm2 -
        (1 + 12*x - 20*x2)*L2 - (52*x - 48*x2)*Lm -
        (16 + 18*x + 48*x2)*L + 26 - 82*x + 80*x2
    );
    
    double const_CATR = (
        (1 - 2*x + 2*x2)*(-4./3.*Lm3 + 8*Lm*Li2xm- 8*Li3xm) +
        (1 + 2*x + 2*x2)*(-8*z2*Lp - 16*Lp*Li2minus - 8*L*Lp2 + 4*L2*Lp + 8*L*Li2minus - 8*Li3minus - 16*S12minus) +
        (16 + 64*x)*(2*S12xm + L*Li2xm) -
        (4./3. + 8./3.*x)*L3 + (8 - 32*x + 16*x2)*z3 -
        (16 + 64*x)*z2*L + (16*x + 16*x2)*(Li2minus + L*Lp) +//there is a typo here in the e-Print ( 16+16*x2 -> 16*x+16*x2 )
        (32./3./x + 12 + 64*x - 272./3.*x2)*Li2xm -
        (12 + 48*x - 260./3.*x2 + 32./3./x)*z2 - 4*x2*L*Lm -
        (2 + 8*x - 10*x2)*Lm2 + (2 + 8*x + 46./3.*x2)*L2 +
        (4 + 16*x - 16*x2)*Lm - (56./3. + 172./3.*x + 1600./9.*x2)*L -
        448./27./x -4./3. - 628./3.*x + 6352./27.*x2
    );
    
    double const_tot = CF*TR*const_CFTR + CA*TR*const_CATR;
    
        
    double tmp = const_tot + logmu*Lmu + logmu2*Lmu2;
    
    return 0.5*tmp/16/pi2;
    
}

//_________________________________________________________


double K_gg1_local(double mMu) {
    return -4./3.*TR*log(1./mMu)/4./M_PI;
}

//___________________________________________________________

double a_Qg_30(double x, int v) {
    
    double L=log(x);
    double L2=L*L;
    
    double x1=1-x;
    
    double L1=log(x1);
    double L12=L1*L1;
    double L13=L12*L1;
    
    if(v==0) {
        return 0.5*( a_Qg_30(x,1) + a_Qg_30(x,2));
    }
    
    if(v==1) {
        return 354.1002*L13 + 479.3838*L12 - 7856.784*(2-x) - 6233.530*L2 + 9416.621/x + 1548.891/x*L;
    }
    
    if(v==2) {
        return 226.3840*L13 - 652.2045*L12 - 2686.387*L1 - 7714.786*(2-x) - 2841.851*L2 + 7721.120/x + 1548.891/x*L ;
    }
    
    if(v==3) {//Updated version w.r.t v==4
        return L/x*CA*CA*(41984./243 + 160./9*zeta(2) - 224./9*zeta(3));
    }
    
    if(v==4) { //Version of the paper (used only for benchamrk)
        return -2658.323*L12 - 7449.948*L1 - 7460.002*(2-x) + 3178.819*L2 + 4710.725/x + 1548.891/x*L;
    }
    
    else {
        cout<<"Choose either v=0, v=1, v=2, v=3 or v=4 !!\nExiting!!\n"<<endl;
        exit(-1);
    }

}











