import adani as ad
import numpy as np
from scipy.special import spence

CF = 4./3
CA = 3
ln2 = 0.6931471805599453
Li2_1_2 = 0.5822405265
zeta2 = 1.6449340668482264
zeta3 = 1.2020569031595942
zeta4 = 1.0823232337111381
zeta5 = 1.03692775514337

def pgq(x):
    return 2. / x - 2. + x

def pqg(x):
    return 1. - 2. * x + 2. * x * x

def pggreg(x):
    return 1. / x - 2. + x - x * x

def pggsing(x):
    return 1. / (1. - x)

def Li2(x):
    return spence(1-x)

def H2(x, k):

    (i,j) = k

    if (i == 0 and j == 0):
        L = np.log(x)
        return 0.5 * L * L
    if (i == 0 and j == 1) :
        return Li2(x)
    if (i == 0 and j == -1):
        return -Li2(-x)
    if (i == 1 and j == 0):
        return -np.log(x) * np.log(1. - x) - Li2(x)
    # there's a typo in the paper: +Li2(x)
    if (i == 1 and j == 1):
        Lm = np.log(1. - x)
        return 0.5 * Lm * Lm
    if (i == 1 and j == -1):
        xm = 1. - x
        return Li2(xm / 2.) - ln2 * np.log(xm) - Li2_1_2
    if (i == -1 and j == 0):
        return np.log(x) * np.log(1. + x) + Li2(-x)
    if (i == -1 and j == 1) :
        xp = 1. + x
        return Li2(xp / 2) - ln2 * np.log(xp) - Li2_1_2
    if (i == -1 and j == -1):
        Lp = np.log(1. + x)
        return 0.5 * Lp * Lp

def test_Pgq0():
    Pgq0 = lambda x: 2. * CF * pgq(x)

    sf = ad.SplittingFunction(0, 'g', 'q')

    for x in np.geomspace(1e-5, 1., 10, endpoint=True):
        for nf in range(1, 6 + 1):
            np.testing.assert_allclose(Pgq0(x), sf.Regular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.Singular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.SingularIntegrated(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.Local(nf), rtol = 1e-7)

def test_Pqg0():
    Pqg0 = lambda x, nf: 2. * nf * pqg(x)

    sf = ad.SplittingFunction(0, 'q', 'g')

    for x in np.geomspace(1e-5, 1., 10, endpoint=True):
        for nf in range(1, 6 + 1):
            np.testing.assert_allclose(Pqg0(x, nf), sf.Regular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.Singular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.SingularIntegrated(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.Local(nf), rtol = 1e-7)

def test_Pgg0():
    Pgg0reg = lambda x, nf: CA * 4. * pggreg(x)
    Pgg0sing = lambda x, nf: 4. * CA * pggsing(x)
    Pgg0loc = lambda nf: 11. / 3 * CA - 2. / 3 * nf
    Pgg0sing_int = lambda x, nf: -Pgg0sing(0., 0) * np.log(1. - x)

    sf = ad.SplittingFunction(0, 'g', 'g')

    for x in np.geomspace(1e-5, 1., 10, endpoint=False):
        for nf in range(1, 6 + 1):
            np.testing.assert_allclose(Pgg0reg(x, nf), sf.Regular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(Pgg0sing(x, nf), sf.Singular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(Pgg0sing_int(x, nf), sf.SingularIntegrated(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(Pgg0loc(nf), sf.Local(nf), rtol = 1e-7)

def test_Pqq0():
    Pqq0reg = lambda x, nf:  CF * 2. * (-1. - x)
    Pqq0sing = lambda x, nf: CF * 2. * 2. / (1. - x)
    Pqq0loc = lambda nf: CF * 3.
    Pqq0sing_int = lambda x, nf: -Pqq0sing(0., 0) * np.log(1. - x)

    sf = ad.SplittingFunction(0, 'q', 'q')

    for x in np.geomspace(1e-5, 1., 10, endpoint=False):
        for nf in range(1, 6 + 1):
            np.testing.assert_allclose(Pqq0reg(x, nf), sf.Regular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(Pqq0sing(x, nf), sf.Singular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(Pqq0sing_int(x, nf), sf.SingularIntegrated(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(Pqq0loc(nf), sf.Local(nf), rtol = 1e-7)

def Pgq1(x, nf):

    H0 = np.log(x)
    H1 = - np.log(1 - x)

    H00 = H2(x, (0, 0))
    H10 = H2(x, (1, 0))
    Hm10 = H2(x, (-1, 0))
    H01 = H2(x, (0, 1))
    H11 = H2(x, (1, 1))

    tmp_CACF = (
        zeta2 * 16. + 76. / 9 + 4. / x + 148. / 9 * x + 176. / 9 * x * x
        + Hm10 * (+16. + 16. / x + 8. * x)
        + H0 * (-48. - 20. * x - 32. / 3 * x * x) + H00 * (+16. + 8. * x)
        + H1 * (+88. / 3 - 88. / 3 / x - 68. / 3 * x)
        + H10 * (-16. + 16. / x + 8. * x) + H11 * (-16. + 16. / x + 8. * x)
        + H01 * (-16. + 16. / x + 8. * x))

    tmp_CFnf = (+80. / 9 - 80. / 9 / x - 64. / 9 * x
                      + H1 * (-16. / 3 + 16. / 3 / x + 8. / 3 * x))

    tmp_CFCF = (-10. - 14. * x + H0 * (+8. + 14. * x)
                      + H00 * (-8. + 4. * x) + H1 * (-24. + 24. / x + 20. * x)
                      + H11 * (+16. - 16. / x - 8. * x))

    return CA * CF * tmp_CACF + CF * nf * tmp_CFnf + CF * CF * tmp_CFCF

def test_Pgq1():

    sf = ad.SplittingFunction(1, 'g', 'q')

    for x in np.geomspace(1e-5, 1., 10, endpoint=False):
        for nf in range(1, 6 + 1):
            np.testing.assert_allclose(Pgq1(x, nf), sf.Regular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.Singular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.SingularIntegrated(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.Local(nf), rtol = 1e-7)

def Pgg1reg(x, nf):
    x2 = x * x

    H0 = np.log(x)

    H00 = H2(x, (0, 0))
    H10 = H2(x, (1, 0))
    Hm10 = H2(x, (-1, 0))
    H01 = H2(x, (0, 1))

    gx = (67. / 18 - zeta2 + H00 + 2. * H10 + 2 * H01)
    g1 = 67. / 18 - zeta2

    tmp_CAnf = (116. / 9 - 92. / 9 / x - 76. / 9 * x + 92. / 9 * x2
                      + H0 * (-8. / 3 - 8. / 3 * x))

    tmp_CACA = (
        zeta2 * (32. - 8. / (1. + x) + 16. * x2) - 50. / 9 - 218. / 9 * x
        + Hm10 * (+32. - 16. / (1. + x) + 16. / x + 16. * x + 16. * x2)
        + H0 * (-100. / 3 + 44. / 3 * x - 176. / 3 * x2)
        + H00 * (8. / (1. + x) + 32. * x - 16. * x2)
        + H10 * (-32. + 16. / x + 16. * x - 16. * x2)
        + H01 * (-32. + 16. / x + 16. * x - 16. * x2)
        + 8. * (gx - g1) * pggsing(x))
    # the last term comes from expanding g(z)[f(z)]_+ = g(1)[f(z)]_+ +
    # (g(z)-g(1))f(z) where (g(z)-g(1))f(z) is regular)

    tmp_CFnf = (-32. + 8. / 3 * 1. / x + 16. * x + 40. / 3 * x2
                      + H0 * (-12. - 20. * x) + H00 * (-8. - 8. * x))

    return tmp_CAnf * CA * nf + tmp_CACA * CA * CA + tmp_CFnf * CF * nf

def Pgg1loc(nf):
    tmp_CAnf = -2. / 3
    tmp_CACA = 8. / 3 + 3. * zeta3
    tmp_CFnf = -1. / 2

    return 4. * (tmp_CAnf * CA * nf + tmp_CACA * CA * CA + tmp_CFnf * CF * nf)

def Pgg1sing(x, nf):
    g1 = 67. / 18 - zeta2

    tmp_CAnf = -10. / 9 * pggsing(x)
    tmp_CACA = 2. * g1 * pggsing(x)

    return 4. * (tmp_CAnf * CA * nf + tmp_CACA * CA * CA)

def Pgg1sing_int(x, nf):
    return - Pgg1sing(0, nf) * np.log(1 - x)

def test_Pgg1():

    sf = ad.SplittingFunction(1, 'g', 'g')

    for x in np.geomspace(1e-5, 1., 10, endpoint=False):
        for nf in range(1, 6 + 1):
            np.testing.assert_allclose(Pgg1reg(x, nf), sf.Regular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(Pgg1sing(x, nf), sf.Singular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(Pgg1sing_int(x, nf), sf.SingularIntegrated(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(Pgg1loc(nf), sf.Local(nf), rtol = 1e-7)


def Pgg0_x_Pgq0(x, nf):
    return ((-4. * CF * nf * (2. + (-2. + x) * x)
            + 2. * CA * CF
                  * (-40. + x * (26. + x * (17. + 8. * x))
                     + 12. * (2. + (-2. + x) * x) * np.log(1. - x)
                     - 24. * (1. + x + x * x) * np.log(x)))
           / 3. / x)

def test_Pgg0_x_Pgq0():

    sf = ad.ConvolutedSplittingFunctions(0, 'g', 'g', 0, 'g', 'q')

    for x in np.geomspace(1e-5, 1., 10, endpoint=False):
        for nf in range(1, 6 + 1):
            np.testing.assert_allclose(Pgg0_x_Pgq0(x, nf), sf.Regular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.Singular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.SingularIntegrated(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.Local(nf), rtol = 1e-7)

def Pqq0_x_Pgq0(x):
    return ((2. * CF * CF
            * (4. * (2. + (-2. + x) * x) * np.log(1. - x)
               - x * (-4. + x + 2. * (-2. + x) * np.log(x))))
           / x)

def test_Pqq0_x_Pgq0():

    sf = ad.ConvolutedSplittingFunctions(0, 'q', 'q', 0, 'g', 'q')

    for x in np.geomspace(1e-5, 1., 10, endpoint=False):
        for nf in range(1, 6 + 1):
            np.testing.assert_allclose(Pqq0_x_Pgq0(x), sf.Regular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.Singular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.SingularIntegrated(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.Local(nf), rtol = 1e-7)

def Pgq0_x_Pqg0(x, nf):
    return (4. * CF * nf
           * (1. + 4. / 3 / x - x - 4. * x * x / 3 + 2. * (1 + x) * np.log(x)))

def test_Pgq0_x_Pqg0():

    sf = ad.ConvolutedSplittingFunctions(0, 'g', 'q', 0, 'q', 'g')

    for x in np.geomspace(1e-5, 1., 10, endpoint=False):
        for nf in range(1, 6 + 1):
            np.testing.assert_allclose(Pgq0_x_Pqg0(x, nf), sf.Regular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.Singular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.SingularIntegrated(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.Local(nf), rtol = 1e-7)
