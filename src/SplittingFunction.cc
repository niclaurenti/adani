#include "adani/SplittingFunction.h"
#include "adani/Constants.h"
#include "adani/Exceptions.h"
#include "adani/SpecialFunctions.h"

#include <cmath>

//==========================================================================================//
//  AbstractSplittingFunction: destructor
//------------------------------------------------------------------------------------------//

AbstractSplittingFunction::~AbstractSplittingFunction(){};

//==========================================================================================//
//  AbstractSplittingFunction: CheckOrder
//------------------------------------------------------------------------------------------//

void AbstractSplittingFunction::CheckOrder(int order) const {
    if (order < 0) {
        throw NotValidException(
            "order must be non negative. Got order=" + to_string(order),
            __PRETTY_FUNCTION__, __LINE__
        );
    }
    if (order > 1) {
        throw NotImplementedException(
            "order greater than 1 is not implemented. Got order="
                + to_string(order),
            __PRETTY_FUNCTION__, __LINE__
        );
    }
}

//==========================================================================================//
//  AbstractSplittingFunction: CheckEntry
//------------------------------------------------------------------------------------------//

void AbstractSplittingFunction::CheckEntry(char entry) const {
    if (entry != 'g' && entry != 'q') {
        throw NotValidException(
            "entry must be non 'q' or 'g'. Got '" + string(1, entry) + "'",
            __PRETTY_FUNCTION__, __LINE__
        );
    }
}

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

    try {
        // check order
        CheckOrder(order);
        // check entry1
        CheckEntry(entry1);
        // check entry2
        CheckEntry(entry2);

        SetFunctions();

    } catch (NotImplementedException &e) {
        e.runtime_error();
    } catch (NotValidException &e) {
        e.runtime_error();
    } catch (UnexpectedException &e) {
        e.runtime_error();
    }
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
            throw UnexpectedException(
                "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
            );
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
            throw UnexpectedException(
                "Unexpected exception!", __PRETTY_FUNCTION__, __LINE__
            );
        }
    } else {
        throw UnexpectedException("Unexpected exception!", __PRETTY_FUNCTION__, __LINE__);
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

    try {
        // check order 1
        CheckOrder(order1);
        // check order
        CheckOrder(order2);
        // check entry1
        CheckEntry(entry1);
        // check entry2
        CheckEntry(entry2);
        // check entry3
        CheckEntry(entry3);
        // check entry4
        CheckEntry(entry4);

        if (order1_ == 0 && entry1_ == 'g' && entry2_ == 'g' && order2_ == 0
            && entry3_ == 'g' && entry4_ == 'g') {
            Pgg0_ = new SplittingFunction(0, 'g', 'g');
        } else {
            Pgg0_ = nullptr;
        }

        SetFunctions();

    } catch (NotImplementedException &e) {
        e.runtime_error();
    } catch (NotValidException &e) {
        e.runtime_error();
    } catch (UnexpectedException &e) {
        e.runtime_error();
    }
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
            throw NotImplementedException(
                "P" + string(1, entry1_) + string(1, entry2_)
                    + to_string(order1_) + " x " + "P" + string(1, entry3_)
                    + string(1, entry4_) + to_string(order2_)
                    + " is not implemented",
                __PRETTY_FUNCTION__, __LINE__
            );
        }
    } else {
        throw NotImplementedException(
            "P" + string(1, entry1_) + string(1, entry2_) + to_string(order1_)
                + " x " + "P" + string(1, entry3_) + string(1, entry4_)
                + to_string(order2_) + " is not implemented",
            __PRETTY_FUNCTION__, __LINE__
        );
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
