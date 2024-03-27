import adani as ad
import numpy as np

nf=4

def test_mc_integrals():
    for kind in ['2', 'L']:
        massive_mc = ad.ExactCoefficientFunction(3, kind, 'g', 1e-3, 1e-3, 1000, True, 1000000)
        massive_nomc = ad.ExactCoefficientFunction(3, kind, 'g', 1e-3, 1e-3, 1000, False, 25000)
        for xi in np.geomspace(1e-2, 1e4, 4, endpoint=True):
            m2Q2 = 1/xi
            xmax = 1. / (1. + 4 * m2Q2)
            for x in np.geomspace(1e-5, xmax, 5, endpoint=False):
                    for nf in range(1, 6 + 1):
                        res1 = massive_mc.MuDependentTerms(x, m2Q2, m2Q2, nf)
                        res2 = massive_nomc.MuDependentTerms(x, m2Q2, m2Q2, nf)
                        np.testing.assert_allclose(res1, res2, rtol=5e-3)
