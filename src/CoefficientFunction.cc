#include "adani/CoefficientFunction.h"

using std::cout;
using std::endl;

//==========================================================================================//
//  CoefficientFunction: constructor
//------------------------------------------------------------------------------------------//

CoefficientFunction::CoefficientFunction(
    const int &order, const char &kind, const char &channel
) {

    SetOrder(order);
    SetKind(kind);
    SetChannel(channel);
}

//==========================================================================================//
//  CoefficientFunction: destructor
//------------------------------------------------------------------------------------------//

CoefficientFunction::~CoefficientFunction(){};

//==========================================================================================//
//  CoefficientFunction: set method for order: order = 1, 2, 3
//------------------------------------------------------------------------------------------//

void CoefficientFunction::SetOrder(const int &order) {
    // check order
    if (order < 1 || order > 3) {
        cout << "Error: order must be 1,2 or 3. Got: " << order << endl;
        exit(-1);
    }
    order_ = order;
}

//==========================================================================================//
//  CoefficientFunction: set method for kind: kind = '2', 'L'
//------------------------------------------------------------------------------------------//

void CoefficientFunction::SetKind(const char &kind) {
    // check kind
    if (kind != '2' && kind != 'L') {
        cout << "Error: kind must be 2 or L. Got: " << kind << endl;
        exit(-1);
    }
    kind_ = kind;
}

//==========================================================================================//
//  CoefficientFunction: set method for channel: channel = 'g', 'q'
//------------------------------------------------------------------------------------------//

void CoefficientFunction::SetChannel(const char &channel) {
    // check channel
    if (channel != 'g' && channel != 'q') {
        cout << "Error: channel must be g or q. Got: " << channel << endl;
        exit(-1);
    }
    if (channel_ == 'q' && order_ == 1) {
        cout << "Error: quark coefficeint function at O(as) doesn't exist!"
             << endl;
        exit(-1);
    }
    channel_ = channel;
}

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

double
CoefficientFunction::MuIndependentTerms(double x, double m2Q2, int nf) const {
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
