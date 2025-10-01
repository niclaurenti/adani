#include "adani/Value.h"
#include "adani/Exceptions.h"

#include <cmath>

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
//  Value: get the delta of higher band
//------------------------------------------------------------------------------------------//

double Value::GetHigherDelta() const { return higher_ - central_; }

//==========================================================================================//
//  Value: get the delta of lower band
//------------------------------------------------------------------------------------------//

double Value::GetLowerDelta() const { return central_ - lower_; }

//==========================================================================================//
//  Value: get max between higher and lower delta
//------------------------------------------------------------------------------------------//

double Value::GetMaxDelta() const {
    double high = GetHigherDelta();
    double low = GetLowerDelta();
    return high > low ? high : low;
}

//==========================================================================================//
//  Value: get min between higher and lower delta
//------------------------------------------------------------------------------------------//

double Value::GetMinDelta() const {
    double high = GetHigherDelta();
    double low = GetLowerDelta();
    return high < low ? high : low;
}

//==========================================================================================//
//  Value: get average between higher and lower delta
//------------------------------------------------------------------------------------------//

double Value::GetAvgDelta() const {
    return 0.5 * (GetHigherDelta() + GetLowerDelta());
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
//  the error on the result is the larger between the errors of rhs and lhs
//------------------------------------------------------------------------------------------//

Value Value::operator-(const Value &rhs) const {
    double higher = higher_ - rhs.higher_;
    double lower = lower_ - rhs.lower_;
    if (higher > lower) {
        return Value(central_ - rhs.central_, higher, lower);
    } else {
        return Value(central_ - rhs.central_, lower, higher);
    }
}

//==========================================================================================//
//  Value: overload of operator * Value
//  since when some elements of the Value class are negative, it gives probles
//  in the ordering of central, higher and lower, we compute the error as sum of
//  the errors
//------------------------------------------------------------------------------------------//

// Value Value::operator*(const Value &rhs) const {
//     double central = central_ * rhs.central_;

//     double delta_low_lhs = GetLowerDelta();
//     double delta_up_lhs = GetHigherDelta();

//     double delta_low_rhs = rhs.GetLowerDelta();
//     double delta_up_rhs = rhs.GetHigherDelta();

//     double delta_low = delta_low_lhs + delta_low_rhs;
//     double delta_up = delta_up_lhs + delta_up_rhs;

//     return Value(central, central + delta_up, central - delta_low);
// }

//==========================================================================================//
//  Value: overload of operator / Value
//  since when some elements of the Value class are negative, it gives probles
//  in the ordering of central, higher and lower, we compute the error as the
//  larger error between lhs and rhs
//------------------------------------------------------------------------------------------//

// Value Value::operator/(const Value &rhs) const {
//     double central = central_ / rhs.central_;

//     double delta_low_lhs = GetLowerDelta();
//     double delta_up_lhs = GetHigherDelta();

//     double delta_low_rhs = rhs.GetLowerDelta();
//     double delta_up_rhs = rhs.GetHigherDelta();

//     double delta_low =
//         delta_low_lhs > delta_low_rhs ? delta_low_lhs : delta_low_rhs;
//     double delta_up = delta_up_lhs > delta_up_rhs ? delta_up_lhs :
//     delta_up_rhs;

//     return Value(central, central + delta_up, central - delta_low);
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
