/*
 * =====================================================================================
 *
 *       Filename:  Value.h
 *
 *         Author:  Daniele Adani
 *
 *    Description: La fatica chiama e noi fedeli dobbiamo abbracciarla
 *
 *  In this file there is the class Value
 *
 * =====================================================================================
 */

#ifndef Value_h
#define Value_h

#include <iostream>
#include <array>

using std::ostream;
using std::array;

//==========================================================================================//
//  class Value
//------------------------------------------------------------------------------------------//

class Value {
    public:
        Value(const double &central, const double &higher, const double &lower);
        explicit Value(const double &central);
        Value(const double &higher, const double &lower);
        Value(const Value &value);

        // get methods
        double GetCentral() const { return central_; };
        double GetHigher() const { return higher_; };
        double GetLower() const { return lower_; };

        double GetHigherDelta() const;
        double GetLowerDelta() const;
        double GetMaxDelta() const;
        double GetMinDelta() const;
        double GetAvgDelta() const;

        array<double, 3> ToVect() const;

        // overload of operators

        Value operator+(const Value &rhs) const;
        Value operator-(const Value &rhs) const;
        // Value operator*(const Value &rhs) const;
        //  Value operator/(const Value &rhs) const;

        Value operator+(const double &rhs) const;
        friend Value operator+(const double &lhs, const Value &rhs);

        Value operator-(const double &rhs) const;

        Value operator*(const double &rhs) const;
        friend Value operator*(const double &lhs, const Value &rhs);

        Value operator/(const double &rhs) const;

        Value &operator=(const Value &rhs);

        Value &operator*=(const double &rhs);
        Value &operator/=(const double &rhs);
        Value &operator+=(const double &rhs);
        Value &operator-=(const double &rhs);

        friend ostream &operator<<(ostream &os, const Value &rhs);

    private:
        double central_;
        double higher_;
        double lower_;
};

#endif
