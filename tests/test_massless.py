import adani as ad
import numpy as np

TR = 0.5

def C2_g1(x, nf):
    return (4. * nf * TR
           * (-8. * x * x + 8. * x - 1.
              + np.log((1. - x) / x) * (2. * x * x - 2. * x + 1.)))

def CL_g1(x, nf):
    return 16. * nf * TR * x * (1. - x)

def test_lo():
    for kind in ['2', 'L']:
        ml = ad.MasslessCoefficientFunction(1, kind, 'g')
        for x in np.geomspace(1e-5, 1., 100, endpoint=False):
            for nf in range(6 + 1):
                res1 = ml.MuIndependentTerms(x, nf)
                res2 = C2_g1(x, nf) if kind == '2' else CL_g1(x, nf)
                np.testing.assert_allclose(res1, res2, rtol=1e-7)
