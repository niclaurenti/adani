#include "../include/MassiveCoefficientFunctions.h"
#include "../include/MatchingConditions.h"
#include "../include/Convolutions.h"
#include "../include/ColorFactors.h"
#include "../include/Convolutions.h"
#include "../include/SplittingFunctions.h"
#include "../include/SpecialFunctions.h"
#include "apfel/massivecoefficientfunctionsunp_sl.h"
#include <gsl/gsl_integration.h>
#include<cmath>

using namespace apfel;

//NLO => O(\alpha_s^1)
//____________________________________________________


double C2m_g1(double x, double mQ) { //mQ=m^2/Q^2
  
  double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0;

  double beta=sqrt(1-4*mQ*x/(1-x));
  double x2=x*x;
  double mQ2=mQ*mQ;
  double L=log((1+beta)/(1-beta ));
    
  return 4*TR*(L*(-8*x2*mQ2 - 4*x*mQ*(3*x-1) + 2*x2-2*x+1) + 
  			 beta*(4*mQ*x*(x-1)-(8*x2-8*x+1)))/4./M_PI;
 
}

//______________________________________________________

double CLm_g1(double x, double mQ) {
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0;

  double beta=sqrt(1-4*mQ*x/(1-x));
  double x2=x*x;
  double L=log((1+beta)/(1-beta));
  
  return 16*TR*(x*(1-x)*beta - 2*x2*mQ*L)/4/M_PI;

}

//__________________________________________________________

double D2m_g1_4(double x, double mQ) {
	
	return C2m_g1(x,mQ);
	
}

//__________________________________________________________

double DLm_g1_4(double x, double mQ) {
	
	return CLm_g1(x,mQ);
	
}


//__________________________________________________________

double D2m_g1_5(double x, double mQ) {

  return D2m_g1_4(x,mQ)-2*K_Qg1(x,mQ);

}

//__________________________________________________________

double DLm_g1_5(double x, double mQ) {

  return DLm_g1_4(x,mQ);

}


//______________________________________________________
//NNLO => O(\alpha_s^2)

double C2m_g2(double x, double mQ, double mMu) {
	
	double xi=1/mQ;
	double eta= 0.25*xi*(1-x)/x - 1;
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0; 
	
	//if(eta > 1e6 || eta < 1e-6 || xi<1e-3 || xi>1e5) return __builtin_nan("");
  if(eta > 1e6 || eta < 1e-6 || xi<1e-3 || xi>1e5) return 0.;
	

  double pi2=M_PI*M_PI;
  
  Cm22gNC cm(1./(1+4*mQ));
  Cm22bargNC cm_log(1./(1+4*mQ));
  
  return (cm.Regular(x*(1+4*mQ)) + cm_log.Regular(x*(1+4*mQ))*log(1./mMu) )/16./pi2;
	
}

//________________________________________________________

double C2m_ps2(double x, double mQ, double mMu) {
	
	double xi=1/mQ;
	double eta= 0.25*xi*(1-x)/x - 1;
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0; 
	
	//if(eta > 1e6 || eta < 1e-6 || xi<1e-3 || xi>1e5) return __builtin_nan("");
  if(eta > 1e6 || eta < 1e-6 || xi<1e-3 || xi>1e5) return 0.;
	

  double pi2=M_PI*M_PI;
  
  Cm22psNC cm(1./(1+4*mQ));
  Cm22barpsNC cm_log(1./(1+4*mQ));
  
  return (cm.Regular(x*(1+4*mQ)) + cm_log.Regular(x*(1+4*mQ))*log(1./mMu) )/16./pi2;
	
}

//______________________________________________________

double CLm_g2(double x, double mQ, double mMu) {
	
	double xi=1/mQ;
	double eta= 0.25*xi*(1-x)/x - 1;
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0; 
	
	if(eta > 1e6 || eta < 1e-6 || xi<1e-3 || xi>1e5) return __builtin_nan("");
	

  double pi2=M_PI*M_PI;
  
  CmL2gNC cm(1./(1+4*mQ));
  CmL2bargNC cm_log(1./(1+4*mQ));
  
  return ( cm.Regular(x*(1+4*mQ)) + cm_log.Regular(x*(1+4*mQ))*log(1./mMu) )/16./pi2;
	
}

//________________________________________________________

double CLm_ps2(double x, double mQ, double mMu) {
	
	double xi=1/mQ;
	double eta= 0.25*xi*(1-x)/x - 1;
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0; 
	
	if(eta > 1e6 || eta < 1e-6 || xi<1e-3 || xi>1e5) return __builtin_nan("");
	

  double pi2=M_PI*M_PI;
  
  CmL2psNC cm(1./(1+4*mQ));
  CmL2barpsNC cm_log(1./(1+4*mQ));
  
  return (cm.Regular(x*(1+4*mQ)) + cm_log.Regular(x*(1+4*mQ))*log(1./mMu) )/16./pi2;
	
}

//________________________________________________________

double D2m_g2_4(double x, double mQ, double mMu) {
	
	double Lmu=log(1./mMu);
	
	return C2m_g2(x,mQ,mMu) - Lmu/6/M_PI*C2m_g1(x,mQ);
	
}

//_______________________________________________________

double DLm_g2_4(double x, double mQ, double mMu) {
	
	double Lmu=log(1./mMu);
	
	return CLm_g2(x,mQ,mMu) - Lmu/6/M_PI*CLm_g1(x,mQ);
	
}

//_______________________________________________________

double D2m_g2_5(double x, double mQ, double mMu) {

  return 
  	 D2m_g2_4(x, mQ, mMu) 
  	-D2m_g1_4(x,mQ)*K_gg1_local(mMu)
    -2*(K_Qg2(x,mMu)-K_Qg1(x,mMu)*K_gg1_local(mMu))
    -2*C2_b1_x_K_bg1(x,mQ);
  
}


//_____________________________________________________________


double DLm_g2_5(double x, double mQ, double mMu) {

  return 
   DLm_g2_4(x, mQ, mMu) 
  -DLm_g1_4(x,mQ)*K_gg1_local(mMu)
  -2*CL_b1_x_K_bg1(x,mQ);
  
}

//___________________________________________________________

double C2m_ps21(double x, double mQ) {
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0; 	

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  double result, error, relerr = 0.0001;
  struct function_params params ={x, mQ, 1};

  gsl_function F;
  F.function = &C2m_g1_x_Pgq0;
  F.params = &params;

  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &result, &error);

  gsl_integration_workspace_free (w);
  
  return -result;
	
}

//_________________________________________________________

double CLm_ps21(double x, double mQ) {
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0; 	

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  double result, error, relerr = 0.0001;
  struct function_params params ={x, mQ, 1};

  gsl_function F;
  F.function = &CLm_g1_x_Pgq0;
  F.params = &params;

  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &result, &error);

  gsl_integration_workspace_free (w);
  
  return -result;
	
}

//__________________________________________________________

double C2m_g21(double x, double mQ) {
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0; 	

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  double regular, singular1, singular2, error, relerr = 0.0001;
  struct function_params params ={x, mQ, 1};

  gsl_function F;
  F.function = &C2m_g1_x_Pgg0_reg;
  F.params = &params;

  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &regular, &error);

  F.function = &C2m_g1_x_Pgg0_sing;
  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &singular1, &error);

  F.function = &Pgg0sing_int;
  gsl_integration_qag(&F, 0, x, 0, relerr, 1000, 4, w, &singular2, &error);

  singular2 *= - C2m_g1(x, mQ) ;

  gsl_integration_workspace_free (w);
  
  return -(regular + singular1 + singular2);
	
}

//__________________________________________________________

double CLm_g21(double x, double mQ) {
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0; 	

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  double regular, singular1, singular2, error, relerr = 0.0001;
  struct function_params params ={x, mQ, 1};

  gsl_function F;
  F.function = &CLm_g1_x_Pgg0_reg;
  F.params = &params;

  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &regular, &error);

  F.function = &CLm_g1_x_Pgg0_sing;
  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &singular1, &error);

  F.function = &Pgg0sing_int;
  gsl_integration_qag(&F, 0, x, 0, relerr, 1000, 4, w, &singular2, &error);

  singular2 *= - CLm_g1(x, mQ) ;

  gsl_integration_workspace_free (w);
  
  return -(regular + singular1 + singular2);
	
}

//______________________________________________________________

double C2m_ps31(double x, double mQ, int nf) {
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0; 	

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  double regular1, regular2, regular3, singular1, singular2, local, error, relerr = 0.0001;
  struct function_params params ={x, mQ, nf};

  gsl_function F;
  F.function = &C2m_g1_x_Pgq1;
  F.params = &params;

  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &regular1, &error);

  F.function = &C2m_g20_x_Pgq0;
  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &regular2, &error);

  F.function = &C2m_ps20_x_Pqq0_reg;
  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &regular3, &error);

  F.function = &C2m_ps20_x_Pqq0_sing;
  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &singular1, &error);

  F.function = &Pqq0sing_int;
  gsl_integration_qag(&F, 0, x, 0, relerr, 1000, 4, w, &singular2, &error);

  singular2 *= - C2m_ps2(x, mQ, 1) ;

  local = C2m_ps2(x, mQ, 1) * (Pqq0loc() - 2 * beta(0, nf));

  gsl_integration_workspace_free (w);
  
  return -(regular1 + regular2 + regular3 + singular1 + singular2 + local);
	
}

//__________________________________________________________


double CLm_ps31(double x, double mQ, int nf) {
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0; 	

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  double regular1, regular2, regular3, singular1, singular2, local, error, relerr = 0.0001;
  struct function_params params ={x, mQ, nf};

  gsl_function F;
  F.function = &CLm_g1_x_Pgq1;
  F.params = &params;

  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &regular1, &error);

  F.function = &CLm_g20_x_Pgq0;
  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &regular2, &error);

  F.function = &CLm_ps20_x_Pqq0_reg;
  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &regular3, &error);

  F.function = &CLm_ps20_x_Pqq0_sing;
  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &singular1, &error);

  F.function = &Pqq0sing_int;
  gsl_integration_qag(&F, 0, x, 0, relerr, 1000, 4, w, &singular2, &error);

  singular2 *= - CLm_ps2(x, mQ, 1) ;

  local = CLm_ps2(x, mQ, 1) * (Pqq0loc() - 2 * beta(0, nf));

  gsl_integration_workspace_free (w);
  
  return -(regular1 + regular2 + regular3 + singular1 + singular2 + local);
	
}

//__________________________________________________________

double C2m_ps32(double x, double mQ, int nf) {
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0; 	

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  double Pgg, Pqq, Pgq, error, relerr = 0.0001;
  struct function_params params ={x, mQ, nf};

  gsl_function F;
  F.function = &C2m_g1_x_Pgg0_x_Pgq0;
  F.params = &params;

  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &Pgg, &error);

  F.function = &C2m_g1_x_Pqq0_x_Pgq0;
  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &Pqq, &error);

  F.function = &C2m_g1_x_Pgq0;
  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &Pgq, &error);

  gsl_integration_workspace_free (w);
  
  return 0.5 * (Pgg + Pqq) - 3. / 2 * beta(0, nf) * Pgq ;
	
}

//__________________________________________________________

double CLm_ps32(double x, double mQ, int nf) {
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0; 	

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);

  double Pgg, Pqq, Pgq, error, relerr = 0.0001;
  struct function_params params ={x, mQ, nf};

  gsl_function F;
  F.function = &CLm_g1_x_Pgg0_x_Pgq0;
  F.params = &params;

  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &Pgg, &error);

  F.function = &CLm_g1_x_Pqq0_x_Pgq0;
  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &Pqq, &error);

  F.function = &CLm_g1_x_Pgq0;
  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &Pgq, &error);

  gsl_integration_workspace_free (w);
  
  return 0.5 * (Pgg + Pqq) - 3. / 2 * beta(0, nf) * Pgq ;
	
}

//__________________________________________________________

double C2m_g31(double x, double mQ, int nf) {
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0; 	

  gsl_integration_workspace * w = gsl_integration_workspace_alloc(100000);

  double regular1, regular2, regular3, singular1, singular2, singular3, singular4, local1, local2, local3, error, relerr = 0.0001;
  struct function_params params ={x, mQ, nf};

  gsl_function F;
  F.function = &C2m_g1_x_Pgg1_reg;
  F.params = &params;

  gsl_integration_qag(&F, x, 1, 0, relerr, 100000, 4, w, &regular1, &error);

  F.function = &C2m_g1_x_Pgg1_sing;
  gsl_integration_qag(&F, x, 1, 0, relerr, 100000, 4, w, &singular1, &error);

  F.function = &Pgg1sing_int;
  gsl_integration_qag(&F, 0, x, 0, relerr, 100000, 4, w, &singular2, &error);

  singular2 *= - C2m_g1(x, mQ) ;

  local1 = C2m_g1(x, mQ) * Pgg1loc(nf) ;

  local2 = C2m_g1(x, mQ) * ( - beta(1,nf) ) ;

  F.function = &C2m_ps20_x_Pqg0;
  gsl_integration_qag(&F, x, 1, 0, relerr, 100000, 4, w, &regular2, &error);

  F.function = &C2m_g20_x_Pgg0_reg;
  gsl_integration_qag(&F, x, 1, 0, relerr, 100000, 4, w, &regular3, &error);

  F.function = &C2m_g20_x_Pgg0_sing;
  gsl_integration_qag(&F, x, 1, 0, relerr, 100000, 4, w, &singular3, &error);

  F.function = &Pgg0sing_int;
  gsl_integration_qag(&F, 0, x, 0, relerr, 100000, 4, w, &singular4, &error);

  singular4 *= - C2m_g2(x, mQ, 1) ;

  local3 = C2m_g2(x, mQ, 1) * ( - 2 * beta(0, nf)) ;

  gsl_integration_workspace_free (w);
  
  //return -(regular1 + singular1 + singular2 + local1 + local2 + regular2 + regular3 + singular3 + singular4 + local3 );
  return -(regular1 + singular1 + singular2 + local1 + local2 + singular4 + local3 );
}

//__________________________________________________________

double CLm_g31(double x, double mQ, int nf) {
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0; 	

  gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

  double regular1, regular2, regular3, singular1, singular2, singular3, singular4, local1, local2, local3, error, relerr = 0.0001;
  struct function_params params ={x, mQ, nf};

  gsl_function F;
  F.function = &CLm_g1_x_Pgg1_reg;
  F.params = &params;

  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &regular1, &error);

  F.function = &CLm_g1_x_Pgg1_sing;
  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &singular1, &error);

  F.function = &Pgg1sing_int;
  gsl_integration_qag(&F, 0, x, 0, relerr, 1000, 4, w, &singular2, &error);

  singular2 *= - CLm_g1(x, mQ) ;

  local1 = CLm_g1(x, mQ) * Pgg1loc(nf) ;

  local2 = CLm_g1(x, mQ) * ( - beta(1,nf) ) ;

  F.function = &CLm_ps20_x_Pqg0;
  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &regular2, &error);

  F.function = &CLm_g20_x_Pgg0_reg;
  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &regular3, &error);

  F.function = &CLm_g20_x_Pgg0_sing;
  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &singular3, &error);

  F.function = &Pgg0sing_int;
  gsl_integration_qag(&F, 0, x, 0, relerr, 1000, 4, w, &singular4, &error);

  singular4 *= - CLm_g2(x, mQ, 1) ;

  local3 = CLm_g2(x, mQ, 1) * ( - 2 * beta(0, nf)) ;

  gsl_integration_workspace_free (w);
  
  return -(regular1 + singular1 + singular2 + local1 + local2 + regular2 + regular3 + singular3 + singular4 + local3 );
	
}

//__________________________________________________________


double C2m_g32(double x, double mQ, int nf) {
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0; 	

  gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

  double regular1, regular2, singular1, singular2, local1, error, relerr = 0.0001;
  struct function_params params ={x, mQ, nf};

  double C2m_g1xPgg0 = C2m_g1_x_Pgg0(x, mQ, nf) ;

  gsl_function F;
  F.function = &C2m_g1_x_Pgg0_Pgg0_reg;
  F.params = &params;

  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &regular1, &error);

  F.function = &C2m_g1_x_Pgg0_Pgg0_sing;
  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &singular1, &error);

  F.function = &Pgg0sing_int;
  gsl_integration_qag(&F, 0, x, 0, relerr, 1000, 4, w, &singular2, &error);

  singular2 *= - C2m_g1xPgg0 ;

  local1 = C2m_g1xPgg0 * Pgg0loc(nf) ;

  regular2 = 0.5 * Pqg0_x_Pgq0(x, nf);  

  gsl_integration_workspace_free (w);
  
  return (regular1 + singular1 + singular2 + local1 + regular2 - 3. / 2 * beta(0, nf) * C2m_g1xPgg0 + beta(0,nf) * beta(0,nf)*C2m_g1(x,mQ) );
	
}

//__________________________________________________________

double CLm_g32(double x, double mQ, int nf) {
	
	double x_max=1./(1+4*mQ);
  
  if (x>x_max || x<0) return 0; 	

  gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);

  double regular1, regular2, singular1, singular2, local1, error, relerr = 0.0001;
  struct function_params params ={x, mQ, nf};

  double CLm_g1xPgg0 = CLm_g1_x_Pgg0(x, mQ, nf) ;

  gsl_function F;
  F.function = &CLm_g1_x_Pgg0_Pgg0_reg;
  F.params = &params;

  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &regular1, &error);

  F.function = &CLm_g1_x_Pgg0_Pgg0_sing;
  gsl_integration_qag(&F, x, 1, 0, relerr, 1000, 4, w, &singular1, &error);

  F.function = &Pgg0sing_int;
  gsl_integration_qag(&F, 0, x, 0, relerr, 1000, 4, w, &singular2, &error);

  singular2 *= - CLm_g1xPgg0 ;

  local1 = CLm_g1xPgg0 * Pgg0loc(nf) ;

  regular2 = 0.5 * Pqg0_x_Pgq0(x, nf);  

  gsl_integration_workspace_free (w);
  
  return (regular1 + singular1 + singular2 + local1 + regular2 - 3. / 2 * beta(0, nf) * CLm_g1xPgg0 + beta(0,nf) * beta(0,nf)*CLm_g1(x,mQ) );
	
}

//__________________________________________________________