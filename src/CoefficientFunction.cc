#include "adani/CoefficientFunction.h"
#include "adani/Exceptions.h"

//==========================================================================================//
//  CoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

CoefficientFunction::CoefficientFunction(
    const int &order, const char &kind, const char &channel
)
    : order_(order), kind_(kind), channel_(channel) {

    try {
        // check order
        if (order_ < 1) {
            throw NotPresentException(
                "massive coefficient functions below order 1 are not present! "
                "Got order="
                    + to_string(order_),
                __PRETTY_FUNCTION__, __LINE__
            );
        }
        if (order_ > 3) {
            throw NotKnownException(
                "massive coefficient functions below order 3 are not known! "
                "Got order="
                    + to_string(order_),
                __PRETTY_FUNCTION__, __LINE__
            );
        }

        // check kind
        if (kind_ != '2' && kind_ != 'L') {
            throw NotValidException(
                "kind must be '2' or 'L'! Got '" + string(1, kind_) + "'",
                __PRETTY_FUNCTION__, __LINE__
            );
        }

        // check channel
        if (channel_ != 'g' && channel_ != 'q') {
            throw NotValidException(
                "channel must be 'g' or 'q'! Got '" + string(1, channel_) + "'",
                __PRETTY_FUNCTION__, __LINE__
            );
        }
        if (channel_ == 'q' && order_ == 1) {
            throw NotPresentException(
                "quark massive coefficient functions at order 1 are not "
                "present! Got channel="
                    + string(1, channel_) + ",order=" + to_string(order_),
                __PRETTY_FUNCTION__, __LINE__
            );
        }
    } catch (const NotPresentException &e) {
        e.runtime_error();
    } catch (const NotKnownException &e) {
        e.runtime_error();
    } catch (const NotValidException &e) {
        e.runtime_error();
    }
}

//==========================================================================================//
//  CoefficientFunction: destructor
//------------------------------------------------------------------------------------------//

CoefficientFunction::~CoefficientFunction(){};

//==========================================================================================//
//  CoefficientFunction: central value of fx
//------------------------------------------------------------------------------------------//

double
    CoefficientFunction::fx(double x, double m2Q2, double m2mu2, int nf) const {
    return fxBand(x, m2Q2, m2mu2, nf).GetCentral();
}

//==========================================================================================//
//  CoefficientFunction: central value of mu independent terms
//------------------------------------------------------------------------------------------//

double CoefficientFunction::MuIndependentTerms(
    double x, double m2Q2, int nf
) const {
    return fx(x, m2Q2, 1., nf);
}

//==========================================================================================//
//  CoefficientFunction: central value mu dependent terms
//------------------------------------------------------------------------------------------//

double CoefficientFunction::MuDependentTerms(
    double x, double m2Q2, double m2mu2, int nf
) const {
    return fx(x, m2Q2, m2mu2, nf) - MuIndependentTerms(x, m2Q2, nf);
}

//==========================================================================================//
//  CoefficientFunction: band for mu independent terms
//------------------------------------------------------------------------------------------//

Value CoefficientFunction::MuIndependentTermsBand(
    double x, double m2Q2, int nf
) const {
    return fxBand(x, m2Q2, 1., nf);
    ;
}

//==========================================================================================//
//  CoefficientFunction: band for mu dependent terms
//------------------------------------------------------------------------------------------//

Value CoefficientFunction::MuDependentTermsBand(
    double x, double m2Q2, double m2mu2, int nf
) const {

    double central = fxBand(x, m2Q2, m2mu2, nf).GetCentral()
                     - MuIndependentTermsBand(x, m2Q2, nf).GetCentral();
    double higher = fxBand(x, m2Q2, m2mu2, nf).GetHigher()
                    - MuIndependentTermsBand(x, m2Q2, nf).GetHigher();
    double lower = fxBand(x, m2Q2, m2mu2, nf).GetLower()
                   - MuIndependentTermsBand(x, m2Q2, nf).GetLower();

    return Value(central, higher, lower);
}
