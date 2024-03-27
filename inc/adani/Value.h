/*
 * =====================================================================================
 *
 *       Filename:  Value.h
 *
 *    Description:  Header file for the
 *                  Value.cc file.
 *
 *         Author: Non credo esista una persona che analizzi il calcio meglio di
 * me
 *
 *  In this file there is the class Value
 *
 * =====================================================================================
 */

#ifndef Value_h
#define Value_h

#include <iostream>
#include <vector>

using std::ostream;
using std::vector;

//==========================================================================================//
//  class Value
//------------------------------------------------------------------------------------------//

class Value {
    public:
        Value(const double &central, const double &higher, const double &lower);
        Value(const double &central);
        Value(const double &higher, const double &lower);
        Value(const Value &value);

        // get methods
        double GetCentral() const { return central_; };
        double GetHigher() const { return higher_; };
        double GetLower() const { return lower_; };

        vector<double> ToVect() const;

        // overload of operators

        Value operator+(const Value &rhs) const;
        // Value operator-(const Value& rhs) const;

        Value operator+(const double &rhs) const;
        friend Value operator+(const double &lhs, const Value &rhs);

        Value operator-(const double &rhs) const;

        Value operator*(const double &rhs) const;
        friend Value operator*(const double &lhs, const Value &rhs);

        Value operator/(const double &rhs) const;

        const Value &operator=(const Value &rhs);

        const Value &operator*=(const double &rhs);
        const Value &operator/=(const double &rhs);

        friend ostream &operator<<(ostream &os, const Value &rhs);

    private:
        double central_;
        double higher_;
        double lower_;
};

#endif
