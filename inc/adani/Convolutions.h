/*
 * =====================================================================================
 *
 *       Filename:  Convolutions.h
 *
 *    Description:  Header file for the Convolutions.cc
 * file.
 *
 *         Author:  Dribbla tutti, anche i cammelli del
 * deserto
 *
 *  In this file there are the convolutions between the
 * functions defined in this library.
 *
 *  For the convolution between a regular function and a
 * function containing regular, singular (plus distribution)
 * and local (delta) parts see Eq. (4.7) of Ref.
 * [arXiv:hep-ph/0504242v1]
 *
 * =====================================================================================
 */

#ifndef Convolutions_h
#define Convolutions_h

#include "adani/CoefficientFunction.h"
#include "adani/SplittingFunctions.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

struct function_params {
    double x;
    double m2Q2;
    int nf;
};

class Convolution {
    public:
        Convolution(const CoefficientFunction& coefffunc, const SplittingFunction& splitfunc, const double& abserr = 1e-3, const double& relerr = 1e-3, const int& dim = 1000);
        Convolution(const CoefficientFunction* coefffunc, const SplittingFunction* splitfunc, const double& abserr = 1e-3, const double& relerr = 1e-3, const int& dim = 1000);
        // Convolution() : Convolution(ExactCoefficientFunction(), SplittingFunction()) {} ;
        ~Convolution();

        // get methods
        int GetAbserr() const {return abserr_;};
        int GetRelerr() const {return relerr_;};
        int GetDim() const {return dim_;} ;

        // set methods
        void SetAbserr(const double& abserr) ;
        void SetRelerr(const double& relerr);
        void SetDim(const int& dim);

        double convolute(double x, double m2Q2, int nf) const ;

    private:
        double abserr_;
        double relerr_;
        int dim_;

        CoefficientFunction *coefffunc_ ;
        SplittingFunction *splitfunc_ ;

        double regular_integrand(double z, void *p) ;
        double singular_integrand(double z, void *p) ;
};

class MonteCarloDoubleConvolution {
    public:
        MonteCarloDoubleConvolution(const CoefficientFunction& coefffunc, const SplittingFunction& splitfunc, const int& MCcalls = 25000, const int& dim = 1000) ;
        MonteCarloDoubleConvolution(const CoefficientFunction* coefffunc, const SplittingFunction* splitfunc, const int& MCcalls = 25000, const int& dim = 1000) ;
        ~MonteCarloDoubleConvolution() ;

        double convolute(double x, double m2Q2, int nf) const ;

    private:
        int dim_;
        int MCcalls_;

        CoefficientFunction *coefffunc_ ;
        SplittingFunction *splitfunc_ ;

        double regular1_integrand(double z[], size_t dim, void *p) const ;
        double regular2_integrand(double z[], size_t dim, void *p) const ;
        double regular3_integrand(double z[], size_t dim, void *p) const ;

        double singular1_integrand(double z[], size_t dim, void *p) const ;
        double singular2_integrand(double z[], size_t dim, void *p) const ;
        double singular3_integrand(double z[], size_t dim, void *p) const ;
};

//==========================================================================================//
//  Convolution between the first order massive gluon
//  coefficient functions and the convolution between the
//  splitting functions Pgg0 and Pgg0 using monte carlo
//  mathods
//------------------------------------------------------------------------------------------//

double C2_g1_x_Pgg0_x_Pgg0_reg1_integrand(double z[], size_t dim, void *p);
double C2_g1_x_Pgg0_x_Pgg0_reg2_integrand(double z[], size_t dim, void *p);
double C2_g1_x_Pgg0_x_Pgg0_reg3_integrand(double z[], size_t dim, void *p);

double C2_g1_x_Pgg0_x_Pgg0_reg(double x, double m2Q2, int nf);

double C2_g1_x_Pgg0_x_Pgg0_sing1_integrand(double z[], size_t dim, void *p);
double C2_g1_x_Pgg0_x_Pgg0_sing2_integrand(double z[], size_t dim, void *p);
double C2_g1_x_Pgg0_x_Pgg0_sing3_integrand(double z[], size_t dim, void *p);

double C2_g1_x_Pgg0_x_Pgg0_sing(double x, double m2Q2, int nf);

double C2_g1_x_Pgg0_x_Pgg0_MC(double x, double m2Q2, int nf);

double CL_g1_x_Pgg0_x_Pgg0_reg1_integrand(double z[], size_t dim, void *p);
double CL_g1_x_Pgg0_x_Pgg0_reg2_integrand(double z[], size_t dim, void *p);
double CL_g1_x_Pgg0_x_Pgg0_reg3_integrand(double z[], size_t dim, void *p);

double CL_g1_x_Pgg0_x_Pgg0_reg(double x, double m2Q2, int nf);

double CL_g1_x_Pgg0_x_Pgg0_sing1_integrand(double z[], size_t dim, void *p);
double CL_g1_x_Pgg0_x_Pgg0_sing2_integrand(double z[], size_t dim, void *p);
double CL_g1_x_Pgg0_x_Pgg0_sing3_integrand(double z[], size_t dim, void *p);

double CL_g1_x_Pgg0_x_Pgg0_sing(double x, double m2Q2, int nf);

double CL_g1_x_Pgg0_x_Pgg0_MC(double x, double m2Q2, int nf);

#endif
