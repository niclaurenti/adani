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

#include "adani/Value.h"

#include <iostream>
#include <vector>

//==========================================================================================//
//                      The notation used is the following:
//                      C^{(k)} = k-th order expansion in
//                      terms of a_s^{[nf]}
//
//------------------------------------------------------------------------------------------//
//==========================================================================================//
//                      Notation:
//                      m2Q2 = m^2/Q^2
//                      m2mu2 = m^2/mu^2
//------------------------------------------------------------------------------------------//

//==========================================================================================//
//  class CoefficientFunction: abstract base class for all the
//  CoefficeintFunction
//------------------------------------------------------------------------------------------//

class CoefficientFunction {

    public:
        CoefficientFunction(
            const int &order, const char &kind, const char &channel
        );
        CoefficientFunction(CoefficientFunction *coeff)
            : CoefficientFunction(
                coeff->GetOrder(), coeff->GetKind(), coeff->GetChannel()
            ){};

        virtual ~CoefficientFunction() = 0;

        virtual double MuIndependentTerms(double x, double m2Q2, int nf) const;
        virtual double fx(double x, double m2Q2, double m2mu2, int nf) const;
        virtual double
        MuDependentTerms(double x, double m2Q2, double m2mu2, int nf) const;
        virtual Value
        fxBand(double x, double m2Q2, double m2mu2, int nf) const = 0;
        virtual Value
        MuIndependentTermsBand(double x, double m2Q2, int nf) const;
        virtual Value
        MuDependentTermsBand(double x, double m2Q2, double m2mu2, int nf) const;

        // get methods
        int GetOrder() const { return order_; };
        char GetKind() const { return kind_; };
        char GetChannel() const { return channel_; };

        // set methods
        void SetOrder(const int &order);
        void SetKind(const char &kind);
        void SetChannel(const char &channel);

    private:
        int order_;    // order = 1, 2, or 3
        char kind_;    // kind_ = '2' for F2 and 'L' for FL
        char channel_; // channel_ = 'g' for Cg and 'q' for Cq
};

#endif
