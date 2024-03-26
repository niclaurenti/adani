#include "adani/Value.h"

using std::cout;
using std::endl;

//==========================================================================================//
//  Value: constructor
//------------------------------------------------------------------------------------------//

Value::Value(const double &central, const double &higher, const double &lower) {
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

//==========================================================================================//
//  Value: constructor
//------------------------------------------------------------------------------------------//

Value::Value(const double &central) {
    central_ = central;
    higher_ = central;
    lower_ = central;
}

//==========================================================================================//
//  Value: constructor
//------------------------------------------------------------------------------------------//

Value::Value(const double &higher, const double &lower) {

    if (higher < lower) {
        cout << "Error: class Value initialized with higher < lower!" << endl;
        cout << "Got: higher=" << higher << ", lower=" << lower << endl;
        exit(-1);
    }
    higher_ = higher;
    lower_ = lower;
    central_ = (higher + lower) / 2;
}

//==========================================================================================//
//  Value: constructor
//------------------------------------------------------------------------------------------//

Value::Value(const Value &value) {
    central_ = value.central_;
    higher_ = value.higher_;
    lower_ = value.lower_;
}

//==========================================================================================//
//  Value: export to vector<double>
//------------------------------------------------------------------------------------------//

vector<double> Value::ToVect() const { return { central_, higher_, lower_ }; }

//==========================================================================================//
//  Value: overload of operator + double
//------------------------------------------------------------------------------------------//

Value Value::operator+(const Value &rhs) const {
    return Value(
        central_ + rhs.central_, higher_ + rhs.higher_, lower_ + rhs.lower_
    );
}

//==========================================================================================//
//  Value: overload of operator - Value
//------------------------------------------------------------------------------------------//

// Value Value::operator-(const Value& rhs) const {

// }

//==========================================================================================//
//  Value: overload of operator double +
//------------------------------------------------------------------------------------------//

Value Value::operator+(const double &rhs) const {
    return Value(rhs + central_, rhs + higher_, rhs + lower_);
}

//==========================================================================================//
//  Value: overload of operator + Value
//------------------------------------------------------------------------------------------//

Value operator+(const double &lhs, const Value &rhs) {
    return Value(lhs + rhs.central_, lhs + rhs.higher_, lhs + rhs.lower_);
}

//==========================================================================================//
//  Value: overload of operator - double
//------------------------------------------------------------------------------------------//

Value Value::operator-(const double &rhs) const {
    return Value(central_ - rhs, higher_ - rhs, lower_ - rhs);
}

//==========================================================================================//
//  Value: overload of operator * double
//------------------------------------------------------------------------------------------//

Value Value::operator*(const double &rhs) const {
    if (rhs > 0)
        return Value(rhs * central_, rhs * higher_, rhs * lower_);
    else if (rhs < 0)
        return Value(rhs * central_, rhs * lower_, rhs * higher_);
    else
        return Value(0.);
}

//==========================================================================================//
//  Value: overload of operator double *
//------------------------------------------------------------------------------------------//

Value operator*(const double &lhs, const Value &rhs) {
    return Value(rhs.central_, rhs.higher_, rhs.lower_) * lhs;
}

//==========================================================================================//
//  Value: overload of operator / double
//------------------------------------------------------------------------------------------//

Value Value::operator/(const double &rhs) const {
    if (rhs > 0)
        return Value(central_ / rhs, higher_ / rhs, lower_ / rhs);
    else
        return Value(central_ / rhs, lower_ / rhs, higher_ / rhs);
}

//==========================================================================================//
//  Value: overload of operator =
//------------------------------------------------------------------------------------------//

const Value &Value::operator=(const Value &rhs) {
    central_ = rhs.central_;
    higher_ = rhs.higher_;
    lower_ = rhs.lower_;

    return *this;
}

//==========================================================================================//
//  Value: overload of operator *= double
//------------------------------------------------------------------------------------------//

const Value &Value::operator*=(const double &rhs) {
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

//==========================================================================================//
//  Value: overload of operator /= double
//------------------------------------------------------------------------------------------//

const Value &Value::operator/=(const double &rhs) {
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

//==========================================================================================//
//  Value: overload of operator <<
//------------------------------------------------------------------------------------------//

ostream &operator<<(ostream &os, const Value &rhs) {
    os << "(" << rhs.central_ << ", " << rhs.higher_ << ", " << rhs.lower_
       << ")";
    return os;
}
