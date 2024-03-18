#include "adani/CoefficientFunction.h"

using std::cout ;
using std::endl ;

Value::Value(const double& central, const double& higher, const double& lower) {
    central_ = central;

    if (higher < central) {
        cout << "Error: class Value initialized with higher < central!" << endl;
        cout << "Got: central=" << central << ", higher=" << higher << endl;
        exit(-1);
    }
    higher_ = higher;

    if (lower > central) {
        cout << "Error: class Value initialized with lower > central!" << endl;
        cout << "Got: central=" << central << ", lower=" << lower << endl;
        exit(-1);
    }
    lower_ = lower;
}

Value::Value(const double& central) {
    central_ = central;
    higher_ = central;
    lower_ = central;
}

Value::Value(const double& higher, const double& lower) {

    if (higher < lower) {
        cout << "Error: class Value initialized with higher < lower!" << endl;
        cout << "Got: higher=" << higher << ", lower=" << lower << endl;
        exit(-1);
    }
    higher_ = higher;
    lower_ = lower;
    central_ = (higher + lower) / 2;
}

Value::Value(const Value& value) {
    central_ = value.central_;
    higher_ = value.higher_;
    lower_ = value.lower_;
}

vector<double> Value::ToVect() const {
    return {central_, higher_, lower_};
}

Value Value::operator+(const Value& rhs) const {
    return Value(central_ + rhs.central_, higher_ + rhs.higher_, lower_ + rhs.lower_);
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
    return Value(central_ - rhs, higher_ - rhs, lower_ - rhs);
}

// Value operator-(const double& lhs, const Value& rhs){
//     return Value(lhs - rhs.central_, lhs - rhs.higher_, lhs - rhs.lower_);
// }

Value Value::operator*(const double& rhs) const {
    if (rhs > 0) return Value(rhs * central_, rhs * higher_, rhs * lower_);
    else return Value(rhs * central_, rhs * lower_, rhs * higher_);
}

Value operator*(const double& lhs, const Value& rhs){
    return Value(rhs.central_, rhs.higher_, rhs.lower_) * lhs;
}

Value Value::operator/(const double& rhs) const {
    if (rhs > 0) return Value(central_ / rhs, higher_ / rhs, lower_ / rhs);
    else return Value(central_ / rhs, lower_ / rhs, higher_ / rhs);
}

Value operator/(const double& lhs, const Value& rhs){
    return Value(rhs.central_, rhs.higher_, rhs.lower_) / lhs ;
}

const Value& Value::operator=(const Value& rhs) {
    central_ = rhs.central_;
    higher_ = rhs.higher_;
    lower_ = rhs.lower_;
    
    return *this;
}

const Value& Value::operator*=(const double& rhs) {
    if (rhs > 0) {
        central_ *= rhs;
        higher_ *= rhs;
        lower_ *= rhs;
    } else {
        central_ *= rhs;
        higher_ = lower_ * rhs;
        lower_ = higher_ * rhs;  
    }

    return *this;
}

const Value& Value::operator/=(const double& rhs) {
    if (rhs > 0) {
        central_ /= rhs;
        higher_ /= rhs;
        lower_ /= rhs;
    } else {
        central_ /= rhs;
        higher_ = lower_ / rhs;
        lower_ = higher_ / rhs;  
    }

    return *this;
}

ostream& operator<<(ostream& os, const Value& rhs){
    os << "(" << rhs.central_ << ", " << rhs.higher_ << ", " << rhs.lower_ << ")" ;
    return os;
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

double CoefficientFunction::fx(double x, double m2Q2, double m2mu2, int nf) const {
    return fxBand(x, m2Q2, m2mu2, nf).GetCentral();
}

double CoefficientFunction::MuIndependentTerms(double x, double m2Q2, int nf) const {
    return fx(x, m2Q2, 1., nf);
}

double CoefficientFunction::MuDependentTerms(double x, double m2Q2, double m2mu2, int nf) const {
    return fx(x, m2Q2, m2mu2, nf) - MuIndependentTerms(x, m2Q2, nf);
}

Value CoefficientFunction::MuIndependentTermsBand(double x, double m2Q2, int nf) const {
    return fxBand(x, m2Q2, 1., nf);;
}

Value CoefficientFunction::MuDependentTermsBand(double x, double m2Q2, double m2mu2, int nf) const {
    
    double central = fxBand(x, m2Q2, m2mu2, nf).GetCentral() -  MuIndependentTermsBand(x, m2Q2, nf).GetCentral();
    double higher = fxBand(x, m2Q2, m2mu2, nf).GetHigher() -  MuIndependentTermsBand(x, m2Q2, nf).GetHigher();
    double lower = fxBand(x, m2Q2, m2mu2, nf).GetLower() -  MuIndependentTermsBand(x, m2Q2, nf).GetLower();

    return Value(central, higher, lower);
}