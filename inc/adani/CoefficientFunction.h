/*
 * =====================================================================================
 *
 *       Filename:  CoefficientFunctions.h
 *
 *    Description:  Header file for the
 *                  CoefficientFunctions.cc file.
 *
 *         Author:  La tigre e il leone saranno anche più forti,
 *                  ma non vedrete mai un lupo al circo
 *
 *  In this file there is the abstract class CoefficientFunction
 *
 * =====================================================================================
 */

#ifndef CoeffFunc
#define CoeffFunc

#include "adani/SplittingFunctions.h"

class Value {
    public:
        Value(const double& central, const double& higher, const double& lower) ;
        Value(const double& central);
        Value(const Value& value) ;

        double GetCentral() const {return central_;};
        double GetHigher() const {return higher_;};
        double GetLower() const {return lower_;};

        double* ToArray() const ;

        Value operator+(const Value& rhs) const;
        // Value operator-(const Value& rhs) const;

        Value operator+(const double& rhs) const;
        friend Value operator+(const double& lhs, const Value& rhs);

        Value operator-(const double& rhs) const;
        friend Value operator-(const double& lhs, const Value& rhs);

        Value operator*(const double& rhs) const;
        friend Value operator*(const double& lhs, const Value& rhs);

        Value operator/(const double& rhs) const;
        friend Value operator/(const double& lhs, const Value& rhs);

        const Value& operator=(const Value& rhs);

        const Value& operator*=(const double& rhs);
        const Value& operator/=(const double& rhs);

    private:
        double central_;
        double higher_;
        double lower_;
};

class CoefficientFunction {

    public:
        CoefficientFunction(const int& order, const char& kind, const char& channel) ;
        CoefficientFunction(CoefficientFunction* coeff) : CoefficientFunction(coeff -> GetOrder(), coeff -> GetKind(), coeff -> GetChannel()) {};

        virtual ~CoefficientFunction() = 0 ;

        // virtual Value fx(double x, double m2Q2, double m2mu2, int nf) const = 0 ;
        // virtual Value MuIndependentTerms(double x, double m2Q2, int nf) const {return fx(x, m2Q2, 1., nf);} ;
        // virtual Value MuDependentTerms(double x, double m2Q2, double m2mu2, int nf) const {
        //     return fx(x, m2Q2, m2mu2, nf) - fx(x, m2Q2, 1., nf);
        // }

        // get methods
        int GetOrder() const { return order_; } ;
        char GetKind() const { return kind_; } ;
        char GetChannel() const { return channel_; } ;

        // set methods
        void SetOrder(const int& order) ;
        void SetKind(const char& kind) ;
        void SetChannel(const char& channel) ;

    private:
        int order_ ; // order = 1, 2, or 3
        char kind_ ; // kind_ = '2' for F2 and 'L' for FL
        char channel_ ; // channel_ = 'g' for Cg and 'q' for Cq

};

#endif
