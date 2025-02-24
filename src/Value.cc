#include "adani/Value.h"
#include "adani/Exceptions.h"

//==========================================================================================//
//  Value: constructor
//------------------------------------------------------------------------------------------//

Value::Value(const double &central, const double &higher, const double &lower) {
    central_ = central;
    try {
        higher_ = higher;
        if (higher < central) {
            throw NotValidException(
                "class Value initialized with higher < central! Got: higher="
                    + to_string(higher) + ", central=" + to_string(central),
                __PRETTY_FUNCTION__, __LINE__
            );
        }

        lower_ = lower;
        if (lower > central) {
            throw NotValidException(
                "class Value initialized with lower > central! Got: lower="
                    + to_string(lower) + ", central=" + to_string(central),
                __PRETTY_FUNCTION__, __LINE__
            );
        }
    } catch (NotValidException &e) {
        e.runtime_error();
    }
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
    try {
        higher_ = higher;
        lower_ = lower;
        if (higher < lower) {
            throw NotValidException(
                "class Value initialized with lower > higher! Got: lower="
                    + to_string(lower) + ", higher=" + to_string(higher),
                __PRETTY_FUNCTION__, __LINE__
            );
        }
    } catch (NotValidException &e) {
        e.runtime_error();
    }

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
