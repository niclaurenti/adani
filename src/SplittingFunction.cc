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
    : AbstractSplittingFunction() {

    // check order
    if (order != 0 && order != 1) {
        cout << "Error: order must be 0 or 1. Got " << order << endl;
        exit(-1);
    }
    order_ = order;

    // check entry1
    if (entry1 != 'g' && entry1 != 'q') {
        cout << "Error: entry1 must be g or q. Got " << entry1 << endl;
        exit(-1);
    }
    entry1_ = entry1;

    // check entry2
    if (entry2 != 'g' && entry2 != 'q') {
        cout << "Error: entry2 must be g or q. Got " << entry2 << endl;
        exit(-1);
    }
    entry2_ = entry2;

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
    : AbstractSplittingFunction() {

    // check order
    if (order1 != 0 && order1 != 1) {
        cout << "Error: order1 must be 0 or 1. Got " << order1 << endl;
        exit(-1);
    }
    order1_ = order1;

    // check order
    if (order2 != 0 && order2 != 1) {
        cout << "Error: order2 must be 0 or 1. Got " << order2 << endl;
        exit(-1);
    }
    order2_ = order2;

    // check entry1
    if (entry1 != 'g' && entry1 != 'q') {
        cout << "Error: entry1 must be g or q. Got " << entry1 << endl;
        exit(-1);
    }
    entry1_ = entry1;

    // check entry2
    if (entry2 != 'g' && entry2 != 'q') {
        cout << "Error: entry2 must be g or q. Got " << entry2 << endl;
        exit(-1);
    }
    entry2_ = entry2;

    // check entry3
    if (entry3 != 'g' && entry3 != 'q') {
        cout << "Error: entry3 must be g or q. Got " << entry3 << endl;
        exit(-1);
    }
    entry3_ = entry3;

    // check entry4
    if (entry4 != 'g' && entry4 != 'q') {
        cout << "Error: entry3 must be g or q. Got " << entry4 << endl;
        exit(-1);
    }
    entry4_ = entry4;

    SetFunctions();
}

//==========================================================================================//
//  ConvolutedSplittingFunctions: regular part
//------------------------------------------------------------------------------------------//

double ConvolutedSplittingFunctions::Regular(double x, int nf) const {
    return GetMultFact() * (this->*reg_)(x, nf);
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
            && entry4_ == 'g')
            reg_ = &ConvolutedSplittingFunctions::Pgq0_x_Pqg0;
        else if (entry1_ == 'g' && entry2_ == 'g' && entry3_ == 'g' && entry4_ == 'q')
            reg_ = &ConvolutedSplittingFunctions::Pgg0_x_Pgq0;
        else if (entry1_ == 'q' && entry2_ == 'q' && entry3_ == 'g' && entry4_ == 'q')
            reg_ = &ConvolutedSplittingFunctions::Pqq0_x_Pgq0;
        else {
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

double SplittingFunction::pggsing(double x) const { return 1. / (1. - x); }

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
    return 11. / 3 * CA - 2. / 3 * nf;
}

//==========================================================================================//
//  Singular part of the gluon-gluon splitting functions O(as)
//
//  Eq. (4.6) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double SplittingFunction::Pgg0sing(double x, int /*nf*/) const {
    return 4. * CA * pggsing(x);
}

//==========================================================================================//
//  Integral from 0 to x of the Singular part of the gluon-gluon splitting
//  functions O(as)
//------------------------------------------------------------------------------------------//

double SplittingFunction::Pgg0sing_integrated(double x, int /*nf*/) const {
    return -Pgg0sing(0., 0) * log(1. - x);
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
    return CF * 2. * 2. / (1. - x);
}

//==========================================================================================//
//  Integral from 0 to x of the Singular part of the quark-quark splitting
//  functions O(as)
//------------------------------------------------------------------------------------------//

double SplittingFunction::Pqq0sing_integrated(double x, int /*nf*/) const {
    return -Pqq0sing(0., 0) * log(1. - x);
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

    double tmp_CACF =
        zeta2 * 16. + 76. / 9 + 4. / x + 148. / 9 * x + 176. / 9 * x * x
        + Hm10 * (+16. + 16. / x + 8. * x)
        + H0 * (-48. - 20. * x - 32. / 3 * x * x) + H00 * (+16. + 8. * x)
        + H1 * (+88. / 3 - 88. / 3 / x - 68. / 3 * x)
        + H10 * (-16. + 16. / x + 8. * x) + H11 * (-16. + 16. / x + 8. * x)
        + H01 * (-16. + 16. / x + 8. * x);

    double tmp_CFnf = +80. / 9 - 80. / 9 / x - 64. / 9 * x
                      + H1 * (-16. / 3 + 16. / 3 / x + 8. / 3 * x);

    double tmp_CFCF = -10. - 14. * x + H0 * (+8. + 14. * x)
                      + H00 * (-8. + 4. * x) + H1 * (-24. + 24. / x + 20. * x)
                      + H11 * (+16. - 16. / x - 8. * x);

    return CA * CF * tmp_CACF + CF * nf * tmp_CFnf + CF * CF * tmp_CFCF;
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

    double gx = (67. / 18 - zeta2 + H00 + 2. * H10 + 2 * H01);
    double g1 = 67. / 18 - zeta2;

    double tmp_CAnf = 116. / 9 - 92. / 9 / x - 76. / 9 * x + 92. / 9 * x2
                      + H0 * (-8. / 3 - 8. / 3 * x);

    double tmp_CACA =
        zeta2 * (32. - 8. / (1. + x) + 16. * x2) - 50. / 9 - 218. / 9 * x
        + Hm10 * (+32. - 16. / (1. + x) + 16. / x + 16. * x + 16. * x2)
        + H0 * (-100. / 3 + 44. / 3 * x - 176. / 3 * x2)
        + H00 * (8. / (1. + x) + 32. * x - 16. * x2)
        + H10 * (-32. + 16. / x + 16. * x - 16. * x2)
        + H01 * (-32. + 16. / x + 16. * x - 16. * x2)
        + 8. * (gx - g1) * pggsing(x);
    // the last term comes from expanding g(z)[f(z)]_+ = g(1)[f(z)]_+ +
    // (g(z)-g(1))f(z) where (g(z)-g(1))f(z) is regular

    double tmp_CFnf = -32. + 8. / 3 * 1. / x + 16. * x + 40. / 3 * x2
                      + H0 * (-12. - 20. * x) + H00 * (-8. - 8. * x);

    return tmp_CAnf * CA * nf + tmp_CACA * CA * CA + tmp_CFnf * CF * nf;
}

//==========================================================================================//
//  Local part of the gluon-gluon splitting functions O(as^2)
//
//  Eq. (4.10) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double SplittingFunction::Pgg1loc(int nf) const {

    double tmp_CAnf = -2. / 3;
    double tmp_CACA = 8. / 3 + 3. * zeta3;
    double tmp_CFnf = -1. / 2;

    return 4. * (tmp_CAnf * CA * nf + tmp_CACA * CA * CA + tmp_CFnf * CF * nf);
}

//==========================================================================================//
//  Singular part of the gluon-gluon splitting functions O(as^2)
//
//  Eq. (4.10) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double SplittingFunction::Pgg1sing(double x, int nf) const {

    double g1 = 67. / 18 - zeta2;

    double tmp_CAnf = -10. / 9 * pggsing(x);
    double tmp_CACA = 2. * g1 * pggsing(x);

    return 4. * (tmp_CAnf * CA * nf + tmp_CACA * CA * CA);
}

//==========================================================================================//
//  Integral from o to x of the Singular part of the gluon-gluon splitting
//  functions O(as^2)
//
//  Eq. (4.10) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double SplittingFunction::Pgg1sing_integrated(double x, int nf) const {

    return -Pgg1sing(0., nf) * log(1. - x);
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
