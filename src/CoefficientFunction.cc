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

Value::Value(const double& central) {
    central_ = central;
    higher_ = central;
    lower_ = central;
}

Value::Value(const Value& value) {
    central_ = value.central_;
    higher_ = value.higher_;
    lower_ = value.lower_;
}

double* Value::ToArray() const {
    double res[] = {central_, higher_, lower_};
    return res;
}

Value Value::operator+(const Value& rhs) const {
    double res_central = central_ + rhs.central_;
    double res_higher = res_central, res_lower = res_central;
    
    double vec_lhs[3] = {central_, higher_, lower_};
    double vec_rhs[3] = {rhs.central_, rhs.higher_, rhs.lower_};

    double tmp;
    for (int i=0; i<3; i++) {
        for(int j=0; j<3; j++) {
            tmp = vec_lhs[i] + vec_rhs[j];
            if (tmp > res_higher) res_higher = tmp;
            if (tmp < res_lower) res_higher = tmp;
        }
    }
    return Value(res_central, res_higher, res_lower);
}

// Value Value::operator-(const Value& rhs) const {
    
// }

Value Value::operator+(const double& rhs) const {
    return Value(rhs + central_, rhs + higher_, rhs + lower_);
}

Value operator+(const double& lhs, const Value& rhs){
    return Value(lhs + rhs.central_, lhs + rhs.higher_, lhs + rhs.lower_);
}

Value Value::operator-(const double& rhs) const {
    return Value(rhs - central_, rhs - higher_, rhs - lower_);
}

Value operator-(const double& lhs, const Value& rhs){
    return Value(lhs - rhs.central_, lhs - rhs.higher_, lhs - rhs.lower_);
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

const Value& Value::operator=(const Value& rhs) {
    central_ = rhs.central_;
    higher_ = rhs.higher_;
    lower_ = rhs.lower_;
    
    return *this;
}

const Value& Value::operator*=(const double& rhs) {
    central_ *= rhs;
    higher_ *= rhs;
    lower_ *= rhs;

    return *this;
}

const Value& Value::operator/=(const double& rhs) {
    central_ /= rhs;
    higher_ /= rhs;
    lower_ /= rhs;

    return *this;
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
