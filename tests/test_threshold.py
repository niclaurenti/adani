import adani as ad
import oldadani as old
import numpy as np

def test_as1():
    for kind in ['2']:
        ml = ad.ThresholdCoefficientFunction(1, kind, 'g')

        for m2Q2 in np.geomspace(1e-2, 1e4, 10):
            xmax = 1. / (1 + 4 * m2Q2)
            for x in np.geomspace(1e-5, xmax, 100, endpoint=False):
                res1 = ml.MuIndependentTerms(x, m2Q2, 1)
                res2 = old.C2_g1_threshold(x, m2Q2)
                np.testing.assert_allclose(res1, res2, rtol=1e-7)

def test_as2():
    for kind in ['2', 'L']:
        thr = ad.ThresholdCoefficientFunction(2, kind, 'g')
        for m2Q2 in np.geomspace(1e-2, 1e4, 10):
            for m2mu2 in np.geomspace(1e-2, 1e4, 10):
                xmax = 1. / (1 + 4 * m2Q2)
                for x in np.geomspace(1e-5, xmax, 100, endpoint=False):
                    res1 = thr.fx(x, m2Q2, m2mu2, 1)
                    if kind == '2':
                        res2 = old.C2_g2_threshold(x, m2Q2, m2mu2)
                    else:
                        res2 = old.CL_g2_threshold(x, m2Q2, m2mu2)
                    np.testing.assert_allclose(res1, res2, rtol=1e-7)

                    res1 = thr.BetaIndependentTerms(x, m2Q2, m2mu2)
                    if kind == '2':
                        res2 = old.C2_g2_threshold_const(x, m2Q2, m2mu2)
                    else:
                        res2 = old.CL_g2_threshold_const(x, m2Q2, m2mu2)
                    np.testing.assert_allclose(res1, res2, rtol=1e-7)

def test_as3():
    for kind in ['2']:
        thr = ad.ThresholdCoefficientFunction(3, kind, 'g')
        for m2Q2 in np.geomspace(1e-2, 1e4, 10):
            for m2mu2 in np.geomspace(1e-2, 1e4, 10):
                xmax = 1. / (1 + 4 * m2Q2)
                for x in np.geomspace(1e-5, xmax, 100, endpoint=False):
                    for nf in range(1, 6 + 1):
                        res1 = thr.fx(x, m2Q2, m2mu2, nf)
                        if kind == '2':
                            res2 = old.C2_g3_threshold(x, m2Q2, m2mu2, nf)
                        else:
                            res2 = old.CL_g3_threshold(x, m2Q2, m2mu2, nf)
                        np.testing.assert_allclose(res1, res2, rtol=1e-7)

                        res1 = thr.BetaIndependentTerms(x, m2Q2, m2mu2)
                        if kind == '2':
                            res2 = old.C2_g3_threshold_const(x, m2Q2, m2mu2)
                        else:
                            res2 = old.CL_g3_threshold_const(x, m2Q2, m2mu2)
                        np.testing.assert_allclose(res1, res2, rtol=1e-7)
