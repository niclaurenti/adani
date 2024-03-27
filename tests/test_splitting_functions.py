import adani as ad
import oldadani as old
import numpy as np

CA = 3
CF = 4./3

def test_Pgq0():

    sf = ad.SplittingFunction(0, 'g', 'q')

    for x in np.geomspace(1e-5, 1., 10, endpoint=True):
        for nf in range(1, 6 + 1):
            np.testing.assert_allclose(old.Pgq0(x), sf.Regular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.Singular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.SingularIntegrated(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.Local(nf), rtol = 1e-7)

def test_Pqg0():

    sf = ad.SplittingFunction(0, 'q', 'g')

    for x in np.geomspace(1e-5, 1., 10, endpoint=True):
        for nf in range(1, 6 + 1):
            np.testing.assert_allclose(old.Pqg0(x, nf), sf.Regular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.Singular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.SingularIntegrated(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.Local(nf), rtol = 1e-7)

def test_Pgg0():

    sf = ad.SplittingFunction(0, 'g', 'g')

    for x in np.geomspace(1e-5, 1., 10, endpoint=False):
        for nf in range(1, 6 + 1):
            np.testing.assert_allclose(old.Pgg0reg(x), sf.Regular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(old.Pgg0sing(x), sf.Singular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(old.Pgg0sing_integrated(x), sf.SingularIntegrated(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(old.Pgg0loc(nf), sf.Local(nf), rtol = 1e-7)

def test_Pqq0():

    sf = ad.SplittingFunction(0, 'q', 'q')

    for x in np.geomspace(1e-5, 1., 10, endpoint=False):
        for nf in range(1, 6 + 1):
            np.testing.assert_allclose(old.Pqq0reg(x), sf.Regular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(old.Pqq0sing(x), sf.Singular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(old.Pqq0sing_integrated(x), sf.SingularIntegrated(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(old.Pqq0loc(), sf.Local(nf), rtol = 1e-7)


def test_Pgq1():

    sf = ad.SplittingFunction(1, 'g', 'q')

    for x in np.geomspace(1e-5, 1., 10, endpoint=False):
        for nf in range(1, 6 + 1):
            np.testing.assert_allclose(old.Pgq1(x, nf), sf.Regular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.Singular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.SingularIntegrated(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(0., sf.Local(nf), rtol = 1e-7)

def test_Pgg1():

    sf = ad.SplittingFunction(1, 'g', 'g')

    for x in np.geomspace(1e-5, 1., 10, endpoint=False):
        for nf in range(1, 6 + 1):
            np.testing.assert_allclose(old.Pgg1reg(x, nf), sf.Regular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(old.Pgg1sing(x, nf), sf.Singular(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(old.Pgg1sing_integrated(x, nf), sf.SingularIntegrated(x, nf), rtol = 1e-7)
            np.testing.assert_allclose(old.Pgg1loc(nf), sf.Local(nf), rtol = 1e-7)


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
