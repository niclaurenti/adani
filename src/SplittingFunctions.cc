#include "adani/SplittingFunctions.h"
#include "adani/Constants.h"
#include "adani/SpecialFunctions.h"
#include <cmath>

using namespace std;

//==========================================================================================//
//  Gluon-quark splitting functions O(alpha_s) without color factors
//
//  Eq. (4.11) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double pgq(double x) { return 2. / x - 2. + x; }

//==========================================================================================//
//  Quark-gluon splitting functions O(alpha_s) without color factors
//
//  Eq. (4.11) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double pqg(double x) { return 1. - 2. * x + 2. * x * x; }

//==========================================================================================//
//  Regular part of the gluon-gluon splitting functions O(alpha_s) without color
//  factors
//
//  Eq. (4.11) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double pggreg(double x) { return 1. / x - 2. + x - x * x; }

//==========================================================================================//
//  Singular part of the gluon-gluon splitting functions O(alpha_s) without
//  color factors
//
//  Eq. (4.11) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double pggsing(double x) { return 1. / (1. - x); }

//==========================================================================================//
//  Gluon-quark splitting functions O(alpha_s)
//
//  Eq. (4.6) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double Pgq0(double x) { return 2. * CF * pgq(x); }

//==========================================================================================//
//  Quark-gluon splitting functions O(alpha_s)
//
//  Eq. (4.6) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double Pqg0(double x, int nf) { return 2. * nf * pqg(x); }

//==========================================================================================//
//  Regular part of the gluon-gluon splitting functions O(alpha_s)
//
//  Eq. (4.6) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double Pgg0reg(double x) { return CA * 4. * pggreg(x); }

//==========================================================================================//
//  Local part of the gluon-gluon splitting functions O(alpha_s)
//
//  Eq. (4.6) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double Pgg0loc(int nf) { return 11. / 3 * CA - 2. / 3 * nf; }

//==========================================================================================//
//  Singular part of the gluon-gluon splitting functions O(alpha_s)
//
//  Eq. (4.6) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double Pgg0sing(double x) { return 4. * CA * pggsing(x); }

//==========================================================================================//
//  Regular part of the quark-quark splitting functions O(alpha_s)
//
//  Eq. (4.6) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double Pqq0reg(double x) { return CF * 2. * (-1. - x); }

//==========================================================================================//
//  Local part of the quark-quark splitting functions O(alpha_s)
//------------------------------------------------------------------------------------------//

double Pqq0loc() { return CF * 3.; }

//==========================================================================================//
//  Singular part of the quark-quark splitting functions O(alpha_s)
//
//  Eq. (4.6) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double Pqq0sing(double x) { return CF * 2. * 2. / (1. - x); }

//==========================================================================================//
//  Gluon-quark splitting functions O(alpha_s^2)
//
//  Eq. (4.9) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double Pgq1(double x, int nf) {

    double H0 = H(x, 0);
    double H1 = H(x, 1);

    double H00 = H(x, 0, 0);
    double H10 = H(x, 1, 0);
    double Hm10 = H(x, -1, 0);
    double H01 = H(x, 0, 1);
    double H11 = H(x, 1, 1);

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
//  Regular part of the gluon-gluon splitting functions O(alpha_s^2)
//
//  Eq. (4.10) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double Pgg1reg(double x, int nf) {

    double x2 = x * x;

    double H0 = H(x, 0);

    double H00 = H(x, 0, 0);
    double H10 = H(x, 1, 0);
    double Hm10 = H(x, -1, 0);
    double H01 = H(x, 0, 1);

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
//  Local part of the gluon-gluon splitting functions O(alpha_s^2)
//
//  Eq. (4.10) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double Pgg1loc(int nf) {

    double tmp_CAnf = -2. / 3;
    double tmp_CACA = 8. / 3 + 3. * zeta3;
    double tmp_CFnf = -1. / 2;

    return 4. * (tmp_CAnf * CA * nf + tmp_CACA * CA * CA + tmp_CFnf * CF * nf);
}

//==========================================================================================//
//  Singular part of the gluon-gluon splitting functions O(alpha_s^2)
//
//  Eq. (4.10) from Ref. [arXiv:hep-ph/0404111]
//------------------------------------------------------------------------------------------//

double Pgg1sing(double x, int nf) {

    double g1 = 67. / 18 - zeta2;

    double tmp_CAnf = -10. / 9 * pggsing(x);
    double tmp_CACA = 2. * g1 * pggsing(x);

    return 4. * (tmp_CAnf * CA * nf + tmp_CACA * CA * CA);
}
