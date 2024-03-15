#include "adani/CoefficientFunction.h"
#include <iostream>

using std::cout ;
using std::endl ;

Value::Value(const double& central, const double& higher, const double& lower) {
    central_ = central;

    if (higher < central) {
        cout << "Error: class Value initialized with higher < central!" << endl;
        exit(-1);
    }
    higher_ = higher;

    if (lower > central) {
        cout << "Error: class Value initialized with lower > central!" << endl;
        exit(-1);
    }
    lower_ = lower;

}

Value Value::operator*(const double& rhs) const {
    return Value(rhs * central_, rhs * higher_, rhs * lower_);
}

Value operator*(const double& lhs, const Value& rhs){
    return Value(lhs * rhs.central_, lhs * rhs.higher_, lhs * rhs.lower_);
}

Value Value::operator/(const double& rhs) const {
    return Value(central_ / rhs, higher_ / rhs, lower_ / rhs);
}

Value operator/(const double& lhs, const Value& rhs){
    return Value(rhs.central_ / lhs, rhs.higher_ / lhs, rhs.lower_ / lhs);
}

CoefficientFunction::CoefficientFunction(const int& order, const char& kind, const char& channel) {

    SetOrder(order);
    SetKind(kind);
    SetChannel(channel);

}

CoefficientFunction::~CoefficientFunction() {};

void CoefficientFunction::SetOrder(const int& order) {
    // check order
    if (order < 1 || order > 3) {
        cout << "Error: order must be 1,2 or 3. Got: " << order << endl ;
        exit(-1);
    }
    order_ = order ;
}

void CoefficientFunction::SetKind(const char& kind) {
    // check kind
    if (kind != '2' && kind !='L') {
        cout << "Error: kind must be 2 or L. Got: " << kind << endl ;
        exit(-1);
    }
    kind_ = kind ;
}

void CoefficientFunction::SetChannel(const char& channel) {
    //check channel
    if (channel != 'g' && channel != 'q') {
        cout << "Error: channel must be g or q. Got: " << channel << endl ;
        exit(-1);
    }
    if (channel_ == 'q' && order_ == 1) {
        cout << "Error: quark coefficeint function at O(as) doesn't exist!" << endl ;
        exit(-1);
    }
    channel_ = channel ;
}
