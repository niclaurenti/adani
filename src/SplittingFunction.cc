#include "adani/SplittingFunction.h"
#include "adani/Constants.h"
#include "adani/SpecialFunctions.h"

#include <cmath>
#include <iostream>

using std::cout;
using std::endl;

//==========================================================================================//
//  AbstractSplittingFunction: destructor
//------------------------------------------------------------------------------------------//

AbstractSplittingFunction::~AbstractSplittingFunction(){};

//==========================================================================================//
//  SplittingFunction: constructor
//------------------------------------------------------------------------------------------//

SplittingFunction::SplittingFunction(
    const int &order, const char &entry1, const char &entry2
)
    : AbstractSplittingFunction(),
      order_(order),
      entry1_(entry1),
      entry2_(entry2) {

    // check order
    if (order != 0 && order != 1) {
        cout << "Error: order must be 0 or 1. Got " << order << endl;
        exit(-1);
    }

    // check entry1
    if (entry1 != 'g' && entry1 != 'q') {
        cout << "Error: entry1 must be g or q. Got " << entry1 << endl;
        exit(-1);
    }

    // check entry2
    if (entry2 != 'g' && entry2 != 'q') {
        cout << "Error: entry2 must be g or q. Got " << entry2 << endl;
        exit(-1);
    }

    SetFunctions();
}

//==========================================================================================//
//  SplittingFunction: Regular part
//------------------------------------------------------------------------------------------//

double SplittingFunction::Regular(double x, int nf) const {
    return GetMultFact() * (this->*reg_)(x, nf);
}

//==========================================================================================//
//  SplittingFunction: Singular part
//------------------------------------------------------------------------------------------//

double SplittingFunction::Singular(double x, int nf) const {
    return GetMultFact() * (this->*sing_)(x, nf);
}

//==========================================================================================//
//  SplittingFunction: Integral from 0 to x of the singular part
//------------------------------------------------------------------------------------------//

double SplittingFunction::SingularIntegrated(double x, int nf) const {
    return GetMultFact() * (this->*sing_int_)(x, nf);
}

//==========================================================================================//
//  SplittingFunction: Local part
//------------------------------------------------------------------------------------------//

double SplittingFunction::Local(int nf) const {
    return GetMultFact() * (this->*loc_)(nf);
}

//==========================================================================================//
//  SplittingFunction: overload of operator * double
//-------------------------------------------------------------constructor-----------------------------//

SplittingFunction SplittingFunction::operator*(const double &rhs) const {
    SplittingFunction res(order_, entry1_, entry2_);
    res.SetMultFact(GetMultFact() * rhs);
    return res;
}

//==========================================================================================//
//  SplittingFunction: overload of operator double *
//------------------------------------------------------------------------------------------//

SplittingFunction operator*(const double &lhs, const SplittingFunction &rhs) {
    SplittingFunction res(rhs.order_, res.entry1_, res.entry2_);
    res.SetMultFact(rhs.GetMultFact() * lhs);
    return res;
}

//==========================================================================================//
//  SplittingFunction: overload of operator / double
//------------------------------------------------------------------------------------------//

SplittingFunction SplittingFunction::operator/(const double &rhs) const {
    SplittingFunction res(order_, entry1_, entry2_);
    res.SetMultFact(GetMultFact() / rhs);
    return res;
}

//==========================================================================================//
//  SplittingFunction: function that sets all pointers to the right function
//------------------------------------------------------------------------------------------//

void SplittingFunction::SetFunctions() {

    if (order_ == 0) {
        if (entry1_ == 'g' && entry2_ == 'q') {

            reg_ = &SplittingFunction::Pgq0;
            sing_ = &SplittingFunction::ZeroFunction_x_nf;
            loc_ = &SplittingFunction::ZeroFunction_nf;
            sing_int_ = &SplittingFunction::ZeroFunction_x_nf;

        } else if (entry1_ == 'q' && entry2_ == 'g') {

            reg_ = &SplittingFunction::Pqg0;

            sing_ = &SplittingFunction::ZeroFunction_x_nf;
            loc_ = &SplittingFunction::ZeroFunction_nf;
            sing_int_ = &SplittingFunction::ZeroFunction_x_nf;

        } else if (entry1_ == 'g' && entry2_ == 'g') {

            reg_ = &SplittingFunction::Pgg0reg;
            sing_ = &SplittingFunction::Pgg0sing;
            loc_ = &SplittingFunction::Pgg0loc;
            sing_int_ = &SplittingFunction::Pgg0sing_integrated;

        } else if (entry1_ == 'q' && entry2_ == 'q') {

            reg_ = &SplittingFunction::Pqq0reg;
            sing_ = &SplittingFunction::Pqq0sing;
            loc_ = &SplittingFunction::Pqq0loc;
            sing_int_ = &SplittingFunction::Pqq0sing_integrated;

        } else {
            cout << "Error: something has gone wrong in "
                    "SplittingFunction::SetFunctions!"
                 << endl;
            exit(-1);
        }
    } else if (order_ == 1) {
        if (entry1_ == 'g' && entry2_ == 'q') {

            reg_ = &SplittingFunction::Pgq1;
            sing_ = &SplittingFunction::ZeroFunction_x_nf;
            loc_ = &SplittingFunction::ZeroFunction_nf;
            sing_int_ = &SplittingFunction::ZeroFunction_x_nf;

        } else if (entry1_ == 'g' && entry2_ == 'g') {

            reg_ = &SplittingFunction::Pgg1reg;
            sing_ = &SplittingFunction::Pgg1sing;
            loc_ = &SplittingFunction::Pgg1loc;
            sing_int_ = &SplittingFunction::Pgg1sing_integrated;

        } else {
            cout << "Error: something has gone wrong in "
                    "SplittingFunction::SetFunctions!"
                 << endl;
            exit(-1);
        }
    } else {
        cout << "Error: P" << entry1_ << entry2_ << order_
             << " is not implemented!" << endl;
        exit(-1);
    }
}

//==========================================================================================//
//  ConvolutedSplittingFunctions: constructor
//------------------------------------------------------------------------------------------//

ConvolutedSplittingFunctions::ConvolutedSplittingFunctions(
    const int &order1, const char &entry1, const char &entry2,
    const int &order2, const char &entry3, const char &entry4
)
    : AbstractSplittingFunction(),
      order1_(order1),
      entry1_(entry1),
      entry2_(entry2),
      order2_(order2),
      entry3_(entry3),
      entry4_(entry4) {

    // check order
    if (order1 != 0 && order1 != 1) {
        cout << "Error: order1 must be 0 or 1. Got " << order1 << endl;
        exit(-1);
    }

    // check order
    if (order2 != 0 && order2 != 1) {
        cout << "Error: order2 must be 0 or 1. Got " << order2 << endl;
        exit(-1);
    }

    // check entry1
    if (entry1 != 'g' && entry1 != 'q') {
        cout << "Error: entry1 must be g or q. Got " << entry1 << endl;
        exit(-1);
    }

    // check entry2
    if (entry2 != 'g' && entry2 != 'q') {
        cout << "Error: entry2 must be g or q. Got " << entry2 << endl;
        exit(-1);
    }

    // check entry3
    if (entry3 != 'g' && entry3 != 'q') {
        cout << "Error: entry3 must be g or q. Got " << entry3 << endl;
        exit(-1);
    }

    // check entry4
    if (entry4 != 'g' && entry4 != 'q') {
        cout << "Error: entry3 must be g or q. Got " << entry4 << endl;
        exit(-1);
    }

    if (order1_ == 0 && entry1_ == 'g' && entry2_ == 'g' && order2_ == 0
        && entry3_ == 'g' && entry4_ == 'g') {
        Pgg0_ = new SplittingFunction(0, 'g', 'g');
    } else {
        Pgg0_ = nullptr;
    }

    SetFunctions();
}

ConvolutedSplittingFunctions::~ConvolutedSplittingFunctions() { delete Pgg0_; }

//==========================================================================================//
//  ConvolutedSplittingFunctions: regular part
//------------------------------------------------------------------------------------------//

double ConvolutedSplittingFunctions::Regular(double x, int nf) const {
    return GetMultFact() * (this->*reg_)(x, nf);
}

//==========================================================================================//
//  ConvolutedSplittingFunctions: singular part
//------------------------------------------------------------------------------------------//

double ConvolutedSplittingFunctions::Singular(double x, int nf) const {
    return GetMultFact() * (this->*sing_)(x, nf);
}

//==========================================================================================//
//  ConvolutedSplittingFunctions: integral from 0 to x of the singular part
//------------------------------------------------------------------------------------------//

double
    ConvolutedSplittingFunctions::SingularIntegrated(double x, int nf) const {
    return GetMultFact() * (this->*sing_int_)(x, nf);
}

//==========================================================================================//
//  ConvolutedSplittingFunctions: local part
//------------------------------------------------------------------------------------------//

double ConvolutedSplittingFunctions::Local(int nf) const {
    return GetMultFact() * (this->*loc_)(nf);
}

//==========================================================================================//
//  ConvolutedSplittingFunctions: overload of operator * double
//------------------------------------------------------------------------------------------//

ConvolutedSplittingFunctions
    ConvolutedSplittingFunctions::operator*(const double &rhs) const {
    ConvolutedSplittingFunctions res(
        order1_, entry1_, entry2_, order2_, entry3_, entry4_
    );
    res.SetMultFact(GetMultFact() * rhs);
    return res;
}

//==========================================================================================//
//  ConvolutedSplittingFunctions: overload of operator double *
//------------------------------------------------------------------------------------------//

ConvolutedSplittingFunctions
    operator*(const double &lhs, const ConvolutedSplittingFunctions &rhs) {
    ConvolutedSplittingFunctions res(
        rhs.order1_, rhs.entry1_, rhs.entry2_, rhs.order2_, rhs.entry3_,
        rhs.entry4_
    );
    res.SetMultFact(rhs.GetMultFact() * lhs);
    return res;
}

//==========================================================================================//
//  ConvolutedSplittingFunctions: overload of operator / double
//------------------------------------------------------------------------------------------//

ConvolutedSplittingFunctions
    ConvolutedSplittingFunctions::operator/(const double &rhs) const {
    ConvolutedSplittingFunctions res(
        order1_, entry1_, entry2_, order2_, entry3_, entry4_
    );
    res.SetMultFact(GetMultFact() / rhs);
    return res;
}

//==========================================================================================//
//  ConvolutedSplittingFunctions: function that sets all pointers to the right
//  function
//------------------------------------------------------------------------------------------//

void ConvolutedSplittingFunctions::SetFunctions() {

    if (order1_ == 0 && order2_ == 0) {
        if (entry1_ == 'g' && entry2_ == 'q' && entry3_ == 'q'
            && entry4_ == 'g') {
            reg_ = &ConvolutedSplittingFunctions::Pgq0_x_Pqg0;
            sing_ = &ConvolutedSplittingFunctions::ZeroFunction_x_nf;
            sing_int_ = &ConvolutedSplittingFunctions::ZeroFunction_x_nf;
            loc_ = &ConvolutedSplittingFunctions::ZeroFunction_nf;
        } else if (entry1_ == 'g' && entry2_ == 'g' && entry3_ == 'g'
                   && entry4_ == 'q') {
            reg_ = &ConvolutedSplittingFunctions::Pgg0_x_Pgq0;
            sing_ = &ConvolutedSplittingFunctions::ZeroFunction_x_nf;
            sing_int_ = &ConvolutedSplittingFunctions::ZeroFunction_x_nf;
            loc_ = &ConvolutedSplittingFunctions::ZeroFunction_nf;
        } else if (entry1_ == 'q' && entry2_ == 'q' && entry3_ == 'g'
                   && entry4_ == 'q') {
            reg_ = &ConvolutedSplittingFunctions::Pqq0_x_Pgq0;
            sing_ = &ConvolutedSplittingFunctions::ZeroFunction_x_nf;
            sing_int_ = &ConvolutedSplittingFunctions::ZeroFunction_x_nf;
            loc_ = &ConvolutedSplittingFunctions::ZeroFunction_nf;
        } else if (entry1_ == 'g' && entry2_ == 'g' && entry3_ == 'g'
                   && entry4_ == 'g') {
            reg_ = &ConvolutedSplittingFunctions::Pgg0_x_Pgg0_reg;
            sing_ = &ConvolutedSplittingFunctions::Pgg0_x_Pgg0_sing;
            sing_int_ =
                &ConvolutedSplittingFunctions::Pgg0_x_Pgg0_sing_integrated;
            loc_ = &ConvolutedSplittingFunctions::Pgg0_x_Pgg0_loc;
        } else {
            cout << "Error: P" << entry1_ << entry2_ << order1_ << " x P"
                 << entry3_ << entry4_ << order2_ << " is not implemented!"
                 << endl;
            exit(-1);
        }
    } else {
        cout << "Error: P" << entry1_ << entry2_ << order1_ << " x P" << entry3_
             << entry4_ << order2_ << " is not implemented!" << endl;
        exit(-1);
    }
}

//==========================================================================================//
//  Delta: overload of operator * double
//------------------------------------------------------------------------------------------//

Delta Delta::operator*(const double &rhs) const {
    Delta res = Delta();
    res.SetMultFact(GetMultFact() * rhs);
    return res;
}

//==========================================================================================//
//  Delta: overload of operator double *
//------------------------------------------------------------------------------------------//

Delta operator*(const double &lhs, const Delta &rhs) {
    Delta res = Delta();
    res.SetMultFact(rhs.GetMultFact() * lhs);
    return res;
}

//==========================================================================================//
//  Delta: overload of operator double /
//------------------------------------------------------------------------------------------//

Delta Delta::operator/(const double &rhs) const {
    Delta res = Delta();
    res.SetMultFact(GetMultFact() / rhs);
    return res;
}

//==========================================================================================//
//  Gluon-quark splitting functions O(as) without color factors
//
//  Eq. (4.11) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double SplittingFunction::pgq(double x) const { return 2. / x - 2. + x; }

//==========================================================================================//
//  Quark-gluon splitting functions O(as) without color factors
//
//  Eq. (4.11) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double SplittingFunction::pqg(double x) const {
    return 1. - 2. * x + 2. * x * x;
}

//==========================================================================================//
//  Regular part of the gluon-gluon splitting functions O(as) without color
//  factors
//
//  Eq. (4.11) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double SplittingFunction::pggreg(double x) const {
    return 1. / x - 2. + x - x * x;
}

//==========================================================================================//
//  Singular part of the gluon-gluon splitting functions O(as) without
//  color factors
//
//  Eq. (4.11) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double SplittingFunction::pggsing(double x) const { return 4. / (1. - x); }

//==========================================================================================//
//  Gluon-quark splitting functions O(as)
//
//  Eq. (4.6) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double SplittingFunction::Pgq0(double x, int /*nf*/) const {
    return 2. * CF * pgq(x);
}

//==========================================================================================//
//  Quark-gluon splitting functions O(as)
//
//  Eq. (4.6) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double SplittingFunction::Pqg0(double x, int nf) const {
    return 2. * nf * pqg(x);
}

//==========================================================================================//
//  Regular part of the gluon-gluon splitting functions O(as)
//
//  Eq. (4.6) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double SplittingFunction::Pgg0reg(double x, int /*nf*/) const {
    return CA * 4. * pggreg(x);
}

//==========================================================================================//
//  Local part of the gluon-gluon splitting functions O(as)
//
//  Eq. (4.6) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double SplittingFunction::Pgg0loc(int nf) const {
    // return 11. / 3 * CA - 2. / 3 * nf;
    return 3.6666666666666665 * CA - 0.6666666666666666 * nf;
}

//==========================================================================================//
//  Singular part of the gluon-gluon splitting functions O(as)
//
//  Eq. (4.6) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double SplittingFunction::Pgg0sing(double x, int /*nf*/) const {
    return CA * pggsing(x);
}

//==========================================================================================//
//  Integral from 0 to x of the Singular part of the gluon-gluon splitting
//  functions O(as)
//------------------------------------------------------------------------------------------//

double SplittingFunction::Pgg0sing_integrated(double x, int /*nf*/) const {
    return -4 * CA * log(1. - x);
}

//==========================================================================================//
//  Regular part of the quark-quark splitting functions O(as)
//
//  Eq. (4.6) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double SplittingFunction::Pqq0reg(double x, int /*nf*/) const {
    return CF * 2. * (-1. - x);
}

//==========================================================================================//
//  Local part of the quark-quark splitting functions O(as)
//------------------------------------------------------------------------------------------//

double SplittingFunction::Pqq0loc(int /*nf*/) const { return CF * 3.; }

//==========================================================================================//
//  Singular part of the quark-quark splitting functions O(as)
//
//  Eq. (4.6) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double SplittingFunction::Pqq0sing(double x, int /*nf*/) const {
    return CF * 4. / (1. - x);
}

//==========================================================================================//
//  Integral from 0 to x of the Singular part of the quark-quark splitting
//  functions O(as)
//------------------------------------------------------------------------------------------//

double SplittingFunction::Pqq0sing_integrated(double x, int /*nf*/) const {
    return -4. * CF * log(1. - x);
}

//==========================================================================================//
//  Gluon-quark splitting functions O(as^2)
//
//  Eq. (4.9) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double SplittingFunction::Pgq1(double x, int nf) const {

    double H0 = H_0(x);
    double H1 = H_1(x);

    double H00 = H_00(x);
    double H10 = H_10(x);
    double Hm10 = H_m10(x);
    double H01 = H_01(x);
    double H11 = H_11(x);

    return CF * nf
               * (8.88888888888889 - 8.88888888888889 / x
                  - 7.111111111111111 * x
                  + H1
                        * (-5.333333333333333 + 5.333333333333333 / x
                           + 2.6666666666666665 * x))
           + CF * CF
                 * (-10. - 24. * H1 + 16. * H11 + (24. * H1 - 16. * H11) / x
                    - 14. * x + 20. * H1 * x - 8. * H11 * x
                    + H00 * (-8. + 4. * x) + H0 * (8. + 14. * x))
           + (CA * CF
              * (4. - 29.333333333333332 * H1 + 16. * H10 + 16. * H11
                 + 16. * Hm10
                 + (8.444444444444445 - 47.99999999999999 * H0 + 16. * H00) * x
                 + H01 * (16. + x * (-16. + 8. * x))
                 + x
                       * (-16. * H10 - 16. * H11 + 16. * Hm10
                          + H1 * (29.333333333333332 - 22.666666666666664 * x)
                          + (16.444444444444443 - 20. * H0 + 8. * H00) * x
                          + x
                                * (8. * H10 + 8. * H11 + 8. * Hm10
                                   + 19.555555555555554 * x
                                   - 10.666666666666666 * H0 * x)
                          + 16. * zeta2)))
                 / x;
}

//==========================================================================================//
//  Regular part of the gluon-gluon splitting functions O(as^2)
//
//  Eq. (4.10) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double SplittingFunction::Pgg1reg(double x, int nf) const {

    double x2 = x * x;

    double H0 = H_0(x);

    double H00 = H_00(x);
    double H10 = H_10(x);
    double Hm10 = H_m10(x);
    double H01 = H_01(x);

    return CA * nf
               * (12.888888888888888
                  + H0 * (-2.6666666666666665 - 2.6666666666666665 * x)
                  - 10.222222222222221 / x - 8.444444444444445 * x
                  + 10.222222222222221 * x2)
           + CF * nf
                 * (-32. + H0 * (-12. - 20. * x) + H00 * (-8. - 8. * x)
                    + 2.6666666666666665 / x + 16. * x
                    + 13.333333333333332 * x2)
           + CA * CA
                 * (-5.555555555555555 - 33.33333333333333 * H0
                    + 2. * (H00 + 2 * H01 + 2. * H10) * pggsing(x)
                    - 24.22222222222222 * x
                    + H0 * (14.666666666666666 - 58.666666666666664 * x) * x
                    + H00 * (32. - 16. * x) * x + (8. * H00) / (1. + x)
                    + H01 * (-32. + 16. / x + 16. * x - 16. * x2)
                    + H10 * (-32. + 16. / x + 16. * x - 16. * x2)
                    + (16. * Hm10 * (1. + x * (1. + 1. * x))
                       * (1. + x * (1. + 1. * x)))
                          / (x * (1. + x))
                    + (32. + 16. * x2 - 8. / (1. + x)) * zeta2);
}

//==========================================================================================//
//  Local part of the gluon-gluon splitting functions O(as^2)
//
//  Eq. (4.10) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double SplittingFunction::Pgg1loc(int nf) const {

    return -2.6666666666666665 * CA * nf - 2. * CF * nf
           + CA * CA * (10.666666666666666 + 12. * zeta3);
}

//==========================================================================================//
//  Singular part of the gluon-gluon splitting functions O(as^2)
//
//  Eq. (4.10) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double SplittingFunction::Pgg1sing(double x, int nf) const {

    return (-1.1111111111111112 * CA * nf
            + CA * CA * (7.444444444444444 - 2. * zeta2))
           * pggsing(x);
}

//==========================================================================================//
//  Integral from o to x of the Singular part of the gluon-gluon splitting
//  functions O(as^2)
//
//  Eq. (4.10) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double SplittingFunction::Pgg1sing_integrated(double x, int nf) const {

    return -(-4.444444444444445 * CA * nf
             + CA * CA * (29.777777777777775 - 8. * zeta2))
           * log(1. - x);
}

//==========================================================================================//
//  Analytical convolution between the splitting functions Pgg0 and Pgq0
//------------------------------------------------------------------------------------------//

double ConvolutedSplittingFunctions::Pgg0_x_Pgq0(double x, int nf) const {

    return (-4. * CF * nf * (2. + (-2. + x) * x)
            + 2. * CA * CF
                  * (-40. + x * (26. + x * (17. + 8. * x))
                     + 12. * (2. + (-2. + x) * x) * log(1. - x)
                     - 24. * (1. + x + x * x) * log(x)))
           / 3. / x;
}

//==========================================================================================//
//  Analytical onvolution between the splitting functions Pqq0 and Pgq0
//------------------------------------------------------------------------------------------//

double ConvolutedSplittingFunctions::Pqq0_x_Pgq0(double x, int /*nf*/) const {

    return (2. * CF * CF
            * (4. * (2. + (-2. + x) * x) * log(1. - x)
               - x * (-4. + x + 2. * (-2. + x) * log(x))))
           / x;
}

//==========================================================================================//
//  Analytical convolution between the splitting functions Pgq0 and Pqg0
//------------------------------------------------------------------------------------------//

double ConvolutedSplittingFunctions::Pgq0_x_Pqg0(double x, int nf) const {

    return 4. * CF * nf
           * (1. + 4. / 3 / x - x - 4. * x * x / 3 + 2. * (1 + x) * log(x));
}

//==========================================================================================//
//  Analytical convolution between the splitting functions Pgg0reg and Pgg0reg
//------------------------------------------------------------------------------------------//

double ConvolutedSplittingFunctions::Pgg0reg_x_Pgg0reg(double x) const {

    return 16 * CA * CA
           * ((-1. + x) * (11. + x * (5. + 2. * x))
              - 3. * (1. + x * (4. + x + x * x)) * log(x))
           / (3. * x);
}

//==========================================================================================//
//  Analytical convolution between the splitting functions Pgg0reg and Pgg0sing
//------------------------------------------------------------------------------------------//

double ConvolutedSplittingFunctions::Pgg0reg_x_Pgg0sing(double x) const {

    return 16. * CA * CA
           * (0.5 * (1. - 4. * x + 3. * x * x)
              + (-2. + 1. / x + x - x * x) * log(1. - x)
              + (2. - x + x * x) * log(x));
}

//==========================================================================================//
//  Regular part of the analytical convolution between the splitting functions
//  Pgg0sing and Pgg0sing
//------------------------------------------------------------------------------------------//

double ConvolutedSplittingFunctions::Pgg0sing_x_Pgg0sing_reg(double x) const {

    return -16 * CA * CA * log(x) / (1. - x);
}

//==========================================================================================//
//  Singular part of the analytical convolution between the splitting functions
//  Pgg0sing and Pgg0sing
//------------------------------------------------------------------------------------------//

double ConvolutedSplittingFunctions::Pgg0sing_x_Pgg0sing_sing(double x) const {

    return 16 * CA * CA * 2. * log(1. - x) / (1. - x);
}

//==========================================================================================//
//  Integral from 0 to x of the singular part of the analytical convolution
//  between the splitting functions Pgg0sing and Pgg0sing
//------------------------------------------------------------------------------------------//

double ConvolutedSplittingFunctions::Pgg0sing_x_Pgg0sing_sing_integrated(
    double x
) const {

    double L1 = log(1. - x);
    return -16. * CA * CA * L1 * L1;
}

//==========================================================================================//
//  Local part of the analytical convolution between the splitting functions
//  Pgg0sing and Pgg0sing
//------------------------------------------------------------------------------------------//

double ConvolutedSplittingFunctions::Pgg0sing_x_Pgg0sing_loc() const {

    return -16 * CA * CA * zeta2;
}

//==========================================================================================//
//  Regular part of the analytical convolution between the splitting functions
//  Pgg0 and Pgg0
//------------------------------------------------------------------------------------------//

double ConvolutedSplittingFunctions::Pgg0_x_Pgg0_reg(double x, int nf) const {

    return Pgg0reg_x_Pgg0reg(x) + 2 * Pgg0reg_x_Pgg0sing(x)
           + 2 * Pgg0_->Regular(x, nf) * Pgg0_->Local(nf)
           + Pgg0sing_x_Pgg0sing_reg(x);
}

//==========================================================================================//
//  Singular part of the analytical convolution between the splitting functions
//  Pgg0 and Pgg0
//------------------------------------------------------------------------------------------//

double ConvolutedSplittingFunctions::Pgg0_x_Pgg0_sing(double x, int nf) const {

    return Pgg0sing_x_Pgg0sing_sing(x)
           + 2 * Pgg0_->Singular(x, nf) * Pgg0_->Local(nf);
}

//==========================================================================================//
//  Integral from 0 to x of the singular part of the analytical convolution
//  between the splitting functions Pgg0 and Pgg0
//------------------------------------------------------------------------------------------//

double ConvolutedSplittingFunctions::Pgg0_x_Pgg0_sing_integrated(
    double x, int nf
) const {

    return Pgg0sing_x_Pgg0sing_sing_integrated(x)
           + 2 * Pgg0_->SingularIntegrated(x, nf) * Pgg0_->Local(nf);
}

//==========================================================================================//
//  Local part of the analytical convolution between the splitting functions
//  Pgg0 and Pgg0
//------------------------------------------------------------------------------------------//

double ConvolutedSplittingFunctions::Pgg0_x_Pgg0_loc(int nf) const {

    double loc = Pgg0_->Local(nf);
    return Pgg0sing_x_Pgg0sing_loc() + loc * loc;
}
