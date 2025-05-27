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

Value Value::operator-(const Value& rhs) const {
    double higher = higher_ - rhs.higher_;
    double lower = lower_ - rhs.lower_;
    if (higher > lower) {
        return Value(
            central_ - rhs.central_,
            higher,
            lower
        );
    } else {
        return Value(
            central_ - rhs.central_,
            lower,
            higher
        );
    }
    // TODO: in this way the error is very small: should take all the
    // combinations?
}

//==========================================================================================//
//  Value: overload of operator * Value
//  since when some elements of the Value class are negative, it gives probles in the
//  ordering of central, higher and lower, we compute the error as sum of the relative errors
//------------------------------------------------------------------------------------------//

Value Value::operator*(const Value& rhs) const {
    double central = central_ * rhs.central_;
    double res = fabs(central);
    double delta_lhs, delta_rhs;

    double d_lhs_h=fabs(central_ - higher_);
    double d_lhs_l = fabs(central_ - lower_);
    double d_rhs_h = fabs(rhs.central_ - rhs.higher_);
    double d_rhs_l = fabs(rhs.central_ - rhs.lower_);

    if (d_lhs_h > d_lhs_l) delta_lhs = d_lhs_h;
    else delta_lhs = d_lhs_l;

    if (d_rhs_h > d_rhs_l) delta_rhs = d_rhs_h;
    else delta_rhs = d_rhs_l;

    // TODO: should I return the average of the errors?

    double delta_res = res * (delta_lhs / fabs(central_) + delta_rhs / fabs(rhs.central_));

    return Value(central, central + delta_res, central - delta_res);

}

//==========================================================================================//
//  Value: overload of operator / Value
//  since when some elements of the Value class are negative, it gives probles in the
//  ordering of central, higher and lower, we compute the error as sum of the relative errors
//------------------------------------------------------------------------------------------//

Value Value::operator/(const Value& rhs) const {
    double central = central_ / rhs.central_;
    double res = fabs(central);
    double delta_lhs, delta_rhs;

    double d_lhs_h=fabs(central_ - higher_);
    double d_lhs_l = fabs(central_ - lower_);
    double d_rhs_h = fabs(rhs.central_ - rhs.higher_);
    double d_rhs_l = fabs(rhs.central_ - rhs.lower_);

    if (d_lhs_h > d_lhs_l) delta_lhs = d_lhs_h;
    else delta_lhs = d_lhs_l;

    if (d_rhs_h > d_rhs_l) delta_rhs = d_rhs_h;
    else delta_rhs = d_rhs_l;

    // TODO: should I return the average of the errors?

    double delta_res = res * (delta_lhs / fabs(central_) + delta_rhs / fabs(rhs.central_));

    return Value(central, central + delta_res, central - delta_res);
}

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
