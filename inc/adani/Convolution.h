/*
 * =====================================================================================
 *
 *       Filename:  Convolution.h
 *
 *    Description:  Header file for the Convolution.cc
 * file.
 *
 *         Author:  Dribbla tutti, anche i cammelli del
 * deserto
 *
 *  In this file there are the convolutions between the coefficeint functions
 * and splitting functions.
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
#include "adani/SplittingFunction.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

//==========================================================================================//
//  class AbstractConvolution
//------------------------------------------------------------------------------------------//

class AbstractConvolution {
    public:
        AbstractConvolution(
            CoefficientFunction *coefffunc,
            AbstractSplittingFunction *splitfunc, const double &abserr = 1e-3,
            const double &relerr = 1e-3, const int &dim = 1000
        );
        virtual ~AbstractConvolution() = 0;

        // result of the convolution
        double Convolute(double x, double m2Q2, int nf) const;

        // integrals of the reular, singular and local parts of the splittings
        virtual double RegularPart(double x, double m2Q2, int nf) const = 0;
        virtual double SingularPart(double x, double m2Q2, int nf) const = 0;
        virtual double LocalPart(double x, double m2Q2, int nf) const = 0;

        // get methods
        double GetAbserr() const { return abserr_; };
        double GetRelerr() const { return relerr_; };
        int GetDim() const { return dim_; };
        gsl_integration_workspace *GetWorkspace() const { return w_; };
        CoefficientFunction *GetCoeffFunc() const { return coefffunc_; };
        AbstractSplittingFunction *GetSplitFunc() const { return splitfunc_; };

        // set methods
        void SetAbserr(const double &abserr);
        void SetRelerr(const double &relerr);
        void SetDim(const int &dim);

    private:
        double abserr_;
        double relerr_;
        int dim_;
        gsl_integration_workspace *w_;

    protected:
        CoefficientFunction *coefffunc_;
        AbstractSplittingFunction *splitfunc_;
};

//==========================================================================================//
//  class Convolution
//------------------------------------------------------------------------------------------//

class Convolution : public AbstractConvolution {
    public:
        Convolution(
            CoefficientFunction *coefffunc,
            AbstractSplittingFunction *splitfunc, const double &abserr = 1e-3,
            const double &relerr = 1e-3, const int &dim = 1000
        )
            : AbstractConvolution(coefffunc, splitfunc, abserr, relerr, dim){};
        ~Convolution() override{};

        double RegularPart(double x, double m2Q2, int nf) const override;
        double SingularPart(double x, double m2Q2, int nf) const override;
        double LocalPart(double x, double m2Q2, int nf) const override;

        // support function for the integral. It is static in order to be passed
        // to gsl
        static double regular_integrand(double z, void *p);
        static double singular_integrand(double z, void *p);
};

//==========================================================================================//
//  class ConvolutedCoefficientFunction: convolution between a
//  CoefficientFunction and a SplittingFunction
//------------------------------------------------------------------------------------------//

class ConvolutedCoefficientFunction : public CoefficientFunction {
    public:
        ConvolutedCoefficientFunction(
            CoefficientFunction *coefffunc,
            AbstractSplittingFunction *splitfunc, const double &abserr = 1e-3,
            const double &relerr = 1e-3, const int &dim = 1000
        );
        ~ConvolutedCoefficientFunction() override;

        // get method
        Convolution *GetConv() const { return conv_; };

        double MuIndependentTerms(double x, double m2Q2, int nf) const override;
        double MuDependentTerms(
            double /*x*/, double /*m2Q2*/, double /*m2mu2*/, int /*nf*/
        ) const override {
            return 0.;
        };
        double fx(double x, double m2Q2, double m2mu2, int nf) const override;
        Value
        fxBand(double x, double m2Q2, double m2mu2, int nf) const override;

    private:
        Convolution *conv_;
};

//==========================================================================================//
//  class DoubleConvolution
//------------------------------------------------------------------------------------------//

class DoubleConvolution : public AbstractConvolution {
    public:
        DoubleConvolution(
            CoefficientFunction *coefffunc,
            AbstractSplittingFunction *splitfunc, const double &abserr = 1e-3,
            const double &relerr = 1e-3, const int &dim = 1000,
            const bool &MCintegral = false, const int &MCcalls = 25000
        );
        ~DoubleConvolution() override;

        double RegularPart(double x, double m2Q2, int nf) const override;
        double SingularPart(double x, double m2Q2, int nf) const override;
        double LocalPart(double x, double m2Q2, int nf) const override;

        // get methods
        bool GetMCintegral() const { return MCintegral_; };
        int GetMCcalls() const { return MCcalls_; };

        // set methods
        void SetMCintegral(const bool &MCintegral);
        void SetMCcalls(const int &MCcalls);

        // support function for the integral. it is static in order to be passed
        // to gsl
        static double regular1_integrand(double z[], size_t /*dim*/, void *p);
        static double regular2_integrand(double z[], size_t /*dim*/, void *p);
        static double regular3_integrand(double z, void *p);

        static double singular1_integrand(double z[], size_t /*dim*/, void *p);
        static double singular2_integrand(double z[], size_t /*dim*/, void *p);
        static double singular3_integrand(double z, void *p);

    private:
        bool MCintegral_;
        int MCcalls_;
        gsl_monte_vegas_state *s_;
        gsl_rng *r_;

        Convolution *convolution_;
        ConvolutedCoefficientFunction *conv_coeff_;
};

//==========================================================================================//
//  struct function_params to be passed to gsl
//------------------------------------------------------------------------------------------//

struct function_params {
        double x;
        double m2Q2;
        int nf;
        const AbstractConvolution *conv;
};

#endif
