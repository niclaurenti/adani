import adani as ad
import oldadani as old
import numpy as np

def test_as1():
    for kind in ['2', 'L']:
        ml = ad.MasslessCoefficientFunction(1, kind, 'g')
        for x in np.geomspace(1e-5, 1., 100, endpoint=False):
            for nf in range(6 + 1):
                res1 = ml.MuIndependentTerms(x, nf)
                res2 = old.C2_g1_massless(x, nf) if kind == '2' else old.CL_g1_massless(x, nf)
                np.testing.assert_allclose(res1, res2, rtol=1e-7)

def test_as2():
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            ml = ad.MasslessCoefficientFunction(2, kind, channel)
            for x in np.geomspace(1e-5, 1., 100, endpoint=False):
                for nf in range(6 + 1):
                    res1 = ml.MuIndependentTerms(x, nf)
                    if (kind == '2' and channel == 'g'):
                        res2 = old.C2_g2_massless(x, nf)
                    if (kind == '2' and channel == 'q'):
                        res2 = old.C2_ps2_massless(x, nf)
                    if (kind == 'L' and channel == 'g'):
                        res2 = old.CL_g2_massless(x, nf)
                    if (kind == 'L' and channel == 'q'):
                        res2 = old.CL_ps2_massless(x, nf)
                    np.testing.assert_allclose(res1, res2, rtol=1e-7)


def test_as3():
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            ml = ad.MasslessCoefficientFunction(3, kind, channel)
            for x in np.geomspace(1e-5, 1., 100, endpoint=False):
                for nf in range(6 + 1):
                    res1 = ml.MuIndependentTerms(x, nf)
                    if (kind == '2' and channel == 'g'):
                        res2 = old.C2_g3_massless(x, nf)
                    if (kind == '2' and channel == 'q'):
                        res2 = old.C2_ps3_massless(x, nf)
                    if (kind == 'L' and channel == 'g'):
                        res2 = old.CL_g3_massless(x, nf)
                    if (kind == 'L' and channel == 'q'):
                        res2 = old.CL_ps3_massless(x, nf)
                    np.testing.assert_allclose(res1, res2, rtol=1e-7)
