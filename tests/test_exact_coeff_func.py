import adani as ad
import oldadani as old
import numpy as np

def test_as1():
    nf = 1
    for kind in ['2', 'L']:
        cf = ad.ExactCoefficientFunction(1, kind, 'g')
        for xi in np.geomspace(1e-2, 1e4, 10, endpoint=True):
            xmax = 1/(1 + 4/xi)
            for x in np.geomspace(1e-5, xmax, 100, endpoint=False):
                res1 = cf.fx(x, 1./xi, 1./xi, nf)
                res2 = old.C2_g1(x, 1./xi) if kind == '2' else old.CL_g1(x, 1./xi)
                np.testing.assert_allclose(res1, res2, rtol = 1e-7)
def eta(x, xi):
    return 0.25 * xi * (1 - x) / x - 1

def test_as2_muindep():
    nf = 1
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            cf = ad.ExactCoefficientFunction(2, kind, channel, 1e-3, 1e-3, 1000)
            for xi in np.geomspace(1e-2, 1e4, 10, endpoint=True):
                xmax = 1/(1 + 4/xi)
                for x in np.geomspace(1e-5, xmax, 10, endpoint=False):
                    res1 = cf.MuIndependentTerms(x, 1./xi, nf)
                    if eta(x, xi) < 1e6 and eta(x, xi) > 1e-6:
                        if kind == '2' and channel == 'g':
                            res2 = old.C2_g20(x, 1/xi)
                        if kind == 'L' and channel == 'g':
                            res2 = old.CL_g20(x, 1/xi)
                        if kind == '2' and channel == 'q':
                            res2 = old.C2_ps20(x, 1/xi)
                        if kind == 'L' and channel == 'q':
                            res2 = old.CL_ps20(x, 1/xi)
                    else:
                        res2 = 0
                    np.testing.assert_allclose(res1, res2, rtol = 1e-7)

def test_as2_order_mudep():
    nf = 1
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            cf = ad.ExactCoefficientFunction(2, kind, channel, 1e-3, 1e-3, 1000)
            for xi in np.geomspace(1e-2, 1e4, 10, endpoint=True):
                for ximu in np.geomspace(1e-2, 1e4, 10, endpoint=True):
                    xmax = 1/(1 + 4/xi)
                    for x in np.geomspace(1e-5, xmax, 10, endpoint=False):
                        res1 = cf.MuDependentTerms(x, 1./xi, 1./ximu, nf)
                        if kind == '2' and channel == 'g':
                            res2 = old.C2_g21(x, 1/xi)
                        elif kind == 'L' and channel == 'g':
                            res2 = old.CL_g21(x, 1/xi)
                        elif kind == '2' and channel == 'q':
                            res2 = old.C2_ps21(x, 1/xi)
                        else:
                            res2 = old.CL_ps21(x, 1/xi)
                        np.testing.assert_allclose(res1, res2 * np.log(ximu), rtol = 1e-7)

def test_as3_order_mudep():
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            cf = ad.ExactCoefficientFunction(3, kind, channel, 1e-3, 1e-3, 1000, False)
            for xi in np.geomspace(1e-2, 1e4, 5, endpoint=True):
                for ximu in np.geomspace(1e-2, 1e4, 5, endpoint=True):
                    xmax = 1/(1 + 4/xi)
                    for x in np.geomspace(1e-5, xmax, 5, endpoint=False):
                        for nf in range(4, 6 + 1):
                            res1 = cf.MuDependentTerms(x, 1./xi, 1./ximu, nf)
                            if kind == '2' and channel == 'g':
                                res2 = old.C2_g31(x, 1/xi, nf) * np.log(ximu) + old.C2_g32(x, 1/xi, nf) * np.log(ximu)**2
                            elif kind == 'L' and channel == 'g':
                                res2 = old.CL_g31(x, 1/xi, nf) * np.log(ximu) + old.CL_g32(x, 1/xi, nf) * np.log(ximu)**2
                            elif kind == '2' and channel == 'q':
                                res2 = old.C2_ps31(x, 1/xi, nf) * np.log(ximu) + old.C2_ps32(x, 1/xi, nf) * np.log(ximu)**2
                            else:
                                res2 = old.CL_ps31(x, 1/xi, nf) * np.log(ximu) + old.CL_ps32(x, 1/xi, nf) * np.log(ximu)**2
                            np.testing.assert_allclose(res1, res2, rtol = 1e-7)
