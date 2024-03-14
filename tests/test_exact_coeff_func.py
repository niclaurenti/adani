import adani as ad
import numpy as np

TR = 0.5

def C2_g1(x, m2Q2):
    beta = np.sqrt(1. - 4. * m2Q2 * x / (1 - x))
    x2 = x * x
    m4Q4 = m2Q2 * m2Q2
    L = np.log((1. + beta) / (1. - beta))

    return (4. * TR
           * (L
                  * (-8. * x2 * m4Q4 - 4. * x * m2Q2 * (3. * x - 1) + 2. * x2
                     - 2. * x + 1.)
              + beta * (4. * m2Q2 * x * (x - 1.) - (8. * x2 - 8. * x + 1.)))
    )

def CL_g1(x, m2Q2):
    beta = np.sqrt(1. - 4. * m2Q2 * x / (1 - x))
    x2 = x * x
    L = np.log((1. + beta) / (1. - beta))

    return 16. * TR * (x * (1. - x) * beta - 2. * x2 * m2Q2 * L)

def test_ciao():
    np.testing.assert_allclose(0., 0.)


def test_leading_order():
    nf = 1
    for kind in ['2', 'L']:
        cf = ad.ExactCoefficientFunction(1, kind, 'g', 1e-3, 1e-3, 1000, 1, 25000)
        for xi in np.geomspace(1e-2, 1e4, 10, endpoint=True):
            xmax = 1/(1 + 4/xi)
            for x in np.geomspace(1e-5, xmax, 100, endpoint=False):
                res1 = cf.fx(x, 1./xi, 1./xi, nf)
                res2 = C2_g1(x, 1./xi) if kind == '2' else CL_g1(x, 1./xi)
                np.testing.assert_allclose(res1, res2, rtol = 1e-7)
def eta(x, xi):
    return 0.25 * xi * (1 - x) / x - 1
    
def test_next_leading_order():
    nf = 1
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            if kind == '2' and channel == 'g':
                func = ad.c2nlog_
            elif kind == 'L' and channel == 'g':
                func = ad.clnlog_
            if kind == '2' and channel == 'q':
                func = ad.c2nloq_
            if kind == 'L' and channel == 'q':
                func = ad.clnloq_
            cf = ad.ExactCoefficientFunction(2, kind, channel, 1e-3, 1e-3, 1000, 1, 25000)
            for xi in np.geomspace(1e-2, 1e4, 10, endpoint=True):
                xmax = 1/(1 + 4/xi)
                for x in np.geomspace(1e-5, xmax, 10, endpoint=False):
                    res1 = cf.MuIndependentTerms(x, 1./xi, nf)
                    if eta(x, xi) < 1e6 and eta(x, xi) > 1e-6:
                        res2 = 16 * np.pi * xi * func(eta(x, xi), xi) / x
                    else:
                        res2 = 0
                    np.testing.assert_allclose(res1, res2, rtol = 1e-7)
                    
def test_next_leading_order_mudep():
    nf = 1
    for kind in ['2', 'L']:
        for channel in ['q', 'g']:
            if kind == '2' and channel == 'g':
                func = ad.c2nlobarg_
            elif kind == 'L' and channel == 'g':
                func = ad.clnlobarg_
            if kind == '2' and channel == 'q':
                func = ad.c2nlobarq_
            if kind == 'L' and channel == 'q':
                func = ad.clnlobarq_
            cf = ad.ExactCoefficientFunction(2, kind, channel, 1e-3, 1e-3, 1000, 1, 25000)
            for xi in np.geomspace(1e-2, 1e4, 10, endpoint=True):
                xmax = 1/(1 + 4/xi)
                # we are comparing an exact (but numerical) result with a 2D interpolation:
                # removing very last points in which the interpolation is worse
                for x in np.geomspace(1e-5, xmax*0.9, 10, endpoint=False):
                    res1 = cf.MuDependentTerms(x, 1./xi, 1./xi, nf)
                    if eta(x, xi) < 1e6 and eta(x, xi) > 1e-6:
                        res2 = 16 * np.pi * xi * func(eta(x, xi), xi) / x * np.log(xi)
                    else:
                        continue
                    np.testing.assert_allclose(res1, res2, rtol = 3e-2)

