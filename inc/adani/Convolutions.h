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

class AbstractConvolution {
    public:
        AbstractConvolution(CoefficientFunction* coefffunc, AbstractSplittingFunction* splitfunc, const double& abserr = 1e-3, const double& relerr = 1e-3, const int& dim = 1000);
        virtual ~AbstractConvolution() = 0;

        double Convolute(double x, double m2Q2, int nf) const ;

        virtual double RegularPart(double x, double m2Q2, int nf) const = 0;
        virtual double SingularPart(double x, double m2Q2, int nf) const = 0;
        virtual double LocalPart(double x, double m2Q2, int nf) const = 0;

        // get methods
        int GetAbserr() const {return abserr_;};
        int GetRelerr() const {return relerr_;};
        int GetDim() const {return dim_;} ;

        CoefficientFunction* GetCoeffFunc() const {return coefffunc_;};
        AbstractSplittingFunction* GetSplitFunc() const {return splitfunc_;};

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
        AbstractSplittingFunction *splitfunc_ ;

};

class Convolution : public AbstractConvolution {
    public:
        Convolution(CoefficientFunction* coefffunc, AbstractSplittingFunction* splitfunc, const double& abserr = 1e-3, const double& relerr = 1e-3, const int& dim = 1000) : AbstractConvolution(coefffunc, splitfunc, abserr, relerr, dim) {} ;
        ~Convolution() override {} ;

        double RegularPart(double x, double m2Q2, int nf) const override;
        double SingularPart(double x, double m2Q2, int nf) const override;
        double LocalPart(double x, double m2Q2, int nf) const override;

        static double regular_integrand(double z, void *p) ;
        static double singular_integrand(double z, void *p) ;

};


class ConvolutedCoefficientFunction : public CoefficientFunction {
    public:
        ConvolutedCoefficientFunction(Convolution* conv) : CoefficientFunction(conv -> GetCoeffFunc()) {conv_ = conv;};
        ~ConvolutedCoefficientFunction() override {} ;

        double MuIndependentTerms(double x, double m2Q2, int nf) const override {return conv_ -> Convolute(x, m2Q2, nf);};
        // double MuDependentTerms(double x, double m2Q2, double m2mu2, int nf) const {return 0.;};
        double fx(double x, double m2Q2, double /*m2mu2*/, int nf) const override {return MuIndependentTerms(x, m2Q2, nf);};
    private:
        Convolution* conv_;

};

class MonteCarloDoubleConvolution : public AbstractConvolution {
    public:
        MonteCarloDoubleConvolution(CoefficientFunction* coefffunc, AbstractSplittingFunction* splitfunc, const double& abserr = 1e-3, const double& relerr = 1e-3, const int& dim = 1000, const int& method_flag = 1, const int& MCcalls = 25000) ;
        ~MonteCarloDoubleConvolution() override ;

        double RegularPart(double x, double m2Q2, int nf) const override;
        double SingularPart(double x, double m2Q2, int nf) const override;
        double LocalPart(double x, double m2Q2, int nf) const override;

        // get methods
        int GetMethodFlag() const {return method_flag_;};
        int GetMCcalls() const {return MCcalls_;};

        // set methods
        void SetMethodFlag(const int& method_flag) ;
        void SetMCcalls(const int& MCcalls) ;

        static double regular1_integrand(double z[], size_t /*dim*/, void *p) ;
        static double regular2_integrand(double z[], size_t /*dim*/, void *p) ;
        static double regular3_integrand(double z, void *p) ;

        static double singular1_integrand(double z[], size_t /*dim*/, void *p) ;
        static double singular2_integrand(double z[], size_t /*dim*/, void *p) ;
        static double singular3_integrand(double z, void *p) ;

    private:
        int method_flag_;
        int MCcalls_;

        Convolution* convolution_;
        ConvolutedCoefficientFunction* conv_coeff_;

};

struct function_params {
    double x;
    double m2Q2;
    int nf;
    const AbstractConvolution* conv;
};

//==========================================================================================//
//  Convolution between the first order massive gluon
//  coefficient functions and the convolution between the
//  splitting functions Pgg0 and Pgg0 using monte carlo
//  mathods
//------------------------------------------------------------------------------------------//

#endif
