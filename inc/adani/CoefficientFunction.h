#ifndef CoeffFunc
#define CoeffFunc

#include "adani/SplittingFunctions.h"

class CoefficientFunction {

    public:
        CoefficientFunction(const int order, const char kind, const char channel) ;
        CoefficientFunction() : CoefficientFunction(1, '2', 'g') {} ;
        virtual ~CoefficientFunction() = 0 ;

        virtual double fx(double x, double m2Q2, double m2mu2, int nf) = 0;

        // get methods
        int GetOrder() { return order_; } ;
        char GetKind() { return kind_; } ;
        char GetChannel() { return channel_; } ;

        // set methods
        void SetOrder(const int order) ;
        void SetKind(const char kind) ;
        void SetChannel(const char channel) ;

        // virtual void Set_fx() = 0;

        double ConvoluteWithSplittingFunction(const SplittingFunction split) ;

    private:
        int order_ ; // order = 1, 2, or 3
        char kind_ ; // kind_ = '2' for F2 and 'L' for FL
        char channel_ ; // channel_ = 'g' for Cg and 'q' for Cq

    // TODO: IMPLEMENT POINTER TO THE RIGHT FUNCTION
    // protected:
    //     double (*fx_)(double, double, double, int);

};

#endif
