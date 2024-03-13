import adani as ad
import numpy as np

CF = 4./3
CA = 3

def pgq(x):
    return 2. / x - 2. + x

def pqg(x):
    return 1. - 2. * x + 2. * x * x

def pggreg(x):
    return 1. / x - 2. + x - x * x

def pggsing(x):
    return 1. / (1. - x)

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
