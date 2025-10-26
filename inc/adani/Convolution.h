/*
 * =====================================================================================
 *
 *       Filename:  Convolution.h
 *
 *         Author:  Daniele Adani
 *
 *    Description:  Dribbla tutti, anche i cammelli del deserto
 *
 *  In this file there are the convolutions between the coefficeint functions
 *  and splitting functions.
 *
 *  For the convolution between a regular function and a
 *  function containing regular, singular (plus distribution)
 *  and local (delta) parts see Eq. (4.7) of Ref.
 *  [arXiv:hep-ph/0504242v1]
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

#include <memory>

//==========================================================================================//
//  class AbstractConvolution
//------------------------------------------------------------------------------------------//

class AbstractConvolution {
    public:
        AbstractConvolution(
            std::shared_ptr<const CoefficientFunction> coefffunc,
            std::shared_ptr<const AbstractSplittingFunction> splitfunc,
            const double &abserr = 1e-3, const double &relerr = 1e-3,
            const int &dim = 1000
        );
        virtual ~AbstractConvolution() = 0;

        AbstractConvolution &operator=(const AbstractConvolution &obj);

        // get methods
        double GetAbsErr() const { return abserr_; };
        double GetRelErr() const { return relerr_; };
        int GetDim() const { return dim_; };
        std::shared_ptr<const CoefficientFunction> GetCoeffFunc() const {
            return coefffunc_;
        };
        std::shared_ptr<const AbstractSplittingFunction> GetSplitFunc() const {
            return splitfunc_;
        };

        // set methods
        void SetAbserr(const double &abserr);
        void SetRelerr(const double &relerr);
        void SetDim(const int &dim);

        // result of the convolution
        double Convolute(double x, double m2Q2, int nf) const;

        // integrals of the reular, singular and local parts of the splittings
        virtual double RegularPart(double x, double m2Q2, int nf) const = 0;
        virtual double SingularPart(double x, double m2Q2, int nf) const = 0;
        virtual double LocalPart(double x, double m2Q2, int nf) const = 0;

    private:
        double abserr_;
        double relerr_;
        int dim_;

        std::shared_ptr<const CoefficientFunction> coefffunc_;
        std::shared_ptr<const AbstractSplittingFunction> splitfunc_;
};

//==========================================================================================//
//  class Convolution
//------------------------------------------------------------------------------------------//

class Convolution : public AbstractConvolution {
    public:
        Convolution(
            std::shared_ptr<const CoefficientFunction> coefffunc,
            std::shared_ptr<const AbstractSplittingFunction> splitfunc,
            const double &abserr = 1e-3, const double &relerr = 1e-3,
            const int &dim = 1000
        );
        Convolution(const Convolution &obj);
        ~Convolution() override;

        Convolution &operator=(const Convolution &obj);

        double RegularPart(double x, double m2Q2, int nf) const override;
        double SingularPart(double x, double m2Q2, int nf) const override;
        double LocalPart(double x, double m2Q2, int nf) const override;

    private:
        static int number_of_instances;
        static gsl_error_handler_t *old_handler_;

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
            std::shared_ptr<const CoefficientFunction> coefffunc,
            std::shared_ptr<const AbstractSplittingFunction> splitfunc,
            const double &abserr = 1e-3, const double &relerr = 1e-3,
            const int &dim = 1000
        );
        ConvolutedCoefficientFunction(const ConvolutedCoefficientFunction &obj);
        ~ConvolutedCoefficientFunction() override = default;

        // get method
        std::shared_ptr<const Convolution> GetConv() const { return conv_; };

        double MuIndependentTerms(double x, double m2Q2, int nf) const override;

        // TODO: these three last functions should be marked as deprecated since
        // not implemented
        double MuDependentTerms(
            double /*x*/, double /*m2Q2*/, double /*m2mu2*/, int /*nf*/
        ) const override {
            return 0.;
        };
        double fx(double x, double m2Q2, double m2mu2, int nf) const override;
        Value
            fxBand(double x, double m2Q2, double m2mu2, int nf) const override;

    private:
        std::shared_ptr<const Convolution> conv_;
};

//==========================================================================================//
//  class DoubleConvolution
//------------------------------------------------------------------------------------------//

class DoubleConvolution : public AbstractConvolution {
    public:
        DoubleConvolution(
            std::shared_ptr<const CoefficientFunction> coefffunc,
            std::shared_ptr<const AbstractSplittingFunction> splitfunc,
            const double &abserr = 1e-3, const double &relerr = 1e-3,
            const int &dim = 1000, const bool &MCintegral = false,
            const int &MCcalls = 25000
        );
        DoubleConvolution(const DoubleConvolution &obj);
        ~DoubleConvolution() override = default;

        DoubleConvolution &operator=(const DoubleConvolution &obj);

        // get methods
        bool GetMCintegral() const { return MCintegral_; };
        int GetMCcalls() const { return MCcalls_; };

        // set methods
        void SetMCcalls(const int &MCcalls);
        void SetMCintegral(const bool &MCintegral);

        double RegularPart(double x, double m2Q2, int nf) const override;
        double SingularPart(double x, double m2Q2, int nf) const override;
        double LocalPart(double x, double m2Q2, int nf) const override;

    private:
        bool MCintegral_;
        int MCcalls_;

        std::unique_ptr<const Convolution> convolution_;
        std::shared_ptr<const ConvolutedCoefficientFunction> conv_coeff_;

        // support function for the integral. it is static in order to be passed
        // to gsl
        static double regular1_integrand(double z[], size_t /*dim*/, void *p);
        static double regular2_integrand(double z[], size_t /*dim*/, void *p);
        static double regular3_integrand(double z, void *p);

        static double singular1_integrand(double z[], size_t /*dim*/, void *p);
        static double singular2_integrand(double z[], size_t /*dim*/, void *p);
        static double singular3_integrand(double z, void *p);
};

#endif
