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

class AbstractConvolution {
    public:
        // AbstractConvolution(const CoefficientFunction& coefffunc, const SplittingFunction& splitfunc, const double& abserr = 1e-3, const double& relerr = 1e-3, const int& dim = 1000);
        AbstractConvolution(CoefficientFunction* coefffunc, SplittingFunction* splitfunc, const double& abserr = 1e-3, const double& relerr = 1e-3, const int& dim = 1000);
        ~AbstractConvolution() {} ;

        double Convolute(double x, double m2Q2, int nf) const ;

        virtual double RegularPart(double x, double m2Q2, int nf) const = 0;
        virtual double SingularPart(double x, double m2Q2, int nf) const = 0;
        virtual double LocalPart(double x, double m2Q2, int nf) const = 0;

        // get methods
        int GetAbserr() const {return abserr_;};
        int GetRelerr() const {return relerr_;};
        int GetDim() const {return dim_;} ;

        // CoefficientFunction* GetCoeffFunc() const {return coefffunc_;};
        // SplittingFunction* GetSplitFunc() const {return splitfunc_;};

        // set methods
        void SetAbserr(const double& abserr) ;
        void SetRelerr(const double& relerr);
        void SetDim(const int& dim);

    private:
        double abserr_;
        double relerr_;
        int dim_;

    protected:
        CoefficientFunction *coefffunc_ ;
        SplittingFunction *splitfunc_ ;

};

class Convolution : public AbstractConvolution {
    public:
        Convolution(CoefficientFunction* coefffunc, SplittingFunction* splitfunc, const double& abserr = 1e-3, const double& relerr = 1e-3, const int& dim = 1000) : AbstractConvolution(coefffunc, splitfunc, abserr, relerr, dim) {} ;
        ~Convolution() {} ;

        // double Convolute(double x, double m2Q2, int nf) const ;

        double RegularPart(double x, double m2Q2, int nf) const override;
        double SingularPart(double x, double m2Q2, int nf) const override;
        double LocalPart(double x, double m2Q2, int nf) const override;

        // friend double regular_integrand(double z, void *p) ;
        // friend double singular_integrand(double z, void *p) ;

    private:
        double regular_integrand(double z, void *p) const ;
        double singular_integrand(double z, void *p) const ;

        // static double static_regular_integrand(double z, void *params) {
        //     // Call the member function through the params pointer
        //     return static_cast<Convolution*>(params)->regular_integrand(z, nullptr);
        // }

};

class MonteCarloDoubleConvolution : public AbstractConvolution {
    public:
        MonteCarloDoubleConvolution(CoefficientFunction* coefffunc, SplittingFunction* splitfunc, const double& abserr = 1e-3, const double& relerr = 1e-3, const int& dim = 1000, const int& MCcalls = 25000) ;
        ~MonteCarloDoubleConvolution() ;

        double RegularPart(double x, double m2Q2, int nf) const override;
        double SingularPart(double x, double m2Q2, int nf) const override;
        double LocalPart(double x, double m2Q2, int nf) const override;

        // get methods
        int GetMCcalls() const {return MCcalls_;};

        // set methods
        void SetMCcalls(const int& MCcalls) ;

    private:
        int MCcalls_;

        Convolution* convolution_;

        double regular1_integrand(double z[], size_t dim, void *p) const ;
        double regular2_integrand(double z[], size_t dim, void *p) const ;
        double regular3_integrand(double z, void *p) const ;

        double singular1_integrand(double z[], size_t dim, void *p) const ;
        double singular2_integrand(double z[], size_t dim, void *p) const ;
        double singular3_integrand(double z, void *p) const ;
};

//==========================================================================================//
//  Convolution between the first order massive gluon
//  coefficient functions and the convolution between the
//  splitting functions Pgg0 and Pgg0 using monte carlo
//  mathods
//------------------------------------------------------------------------------------------//

#endif
