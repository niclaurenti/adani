import adani as ad
import oldadani as old
import numpy as np

def test_as2():
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            asy = ad.AsymptoticCoefficientFunction(2, kind, channel)
            for m2Q2 in np.geomspace(1e-2, 1e4, 10):
                for m2mu2 in np.geomspace(1e-2, 1e4, 10):
                    for x in np.geomspace(1e-5, 1., 100, endpoint=False):
                        res1 = asy.fx(x, m2Q2, m2mu2, 1)
                        if kind == '2' and channel == 'g':
                            res2 = old.C2_g2_asymptotic(x, m2Q2, m2mu2)
                        if kind == '2' and channel == 'q':
                            res2 = old.C2_ps2_asymptotic(x, m2Q2, m2mu2)
                        if kind == 'L' and channel == 'g':
                            res2 = old.CL_g2_asymptotic(x, m2Q2, m2mu2)
                        if kind == 'L' and channel == 'q':
                            res2 = old.CL_ps2_asymptotic(x, m2Q2, m2mu2)
                        np.testing.assert_allclose(res1, res2, rtol=1e-7)


def test_as3():
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            hs_version = "exact" if channel == 'q' else "abmp"
            asy = ad.AsymptoticCoefficientFunction(3, kind, channel, True, hs_version)
            for m2Q2 in np.geomspace(1e-2, 1e4, 10):
                for m2mu2 in np.geomspace(1e-2, 1e4, 10):
                    for nf in range(1, 6+1):
                        for x in np.geomspace(1e-5, 1., 100, endpoint=False):
                            res1 = asy.fx(x, m2Q2, m2mu2, nf)
                            if kind == '2' and channel == 'g':
                                res2 = old.C2_g3_asymptotic(x, m2Q2, m2mu2, nf)
                            if kind == '2' and channel == 'q':
                                res2 = old.C2_ps3_asymptotic(x, m2Q2, m2mu2, nf)
                            if kind == 'L' and channel == 'g':
                                res2 = old.CL_g3_asymptotic(x, m2Q2, m2mu2, nf)
                            if kind == 'L' and channel == 'q':
                                res2 = old.CL_ps3_asymptotic(x, m2Q2, m2mu2, nf)
                            np.testing.assert_allclose(res1, res2, rtol=1e-7)

def test_asy():
    for order in range(1, 3 + 1):
        for kind in ['2', 'L']:
            for channel in ['g', 'q']:
                if order == 1 and channel == 'q':
                    continue
                for hs_version in ["exact", "abmp", "klmv"]:
                    if channel == 'q' and hs_version == "abmp":
                        continue
                    if channel == 'g' and hs_version == "exact":
                        continue
                    hs = ad.HighScaleCoefficientFunction(order, kind, channel, hs_version)
                    for nll in [True, False]:
                        pt = ad.PowerTermsCoefficientFunction(order, kind, channel, nll)
                        asy = ad.AsymptoticCoefficientFunction(order, kind, channel, nll, hs_version)
                        for nf in range(1, 6 + 1):
                            for m2mu2 in np.geomspace(1e-4, 1e2, 10):
                                for m2Q2 in np.geomspace(1e-4, 1e2, 10):
                                    for x in np.geomspace(1e-5, 1., 10, endpoint=False):
                                        res1 = pt.fxBand(x, m2Q2, m2mu2, nf) + hs.fxBand(x, m2Q2, m2mu2, nf)
                                        res2 = asy.fxBand(x, m2Q2, m2mu2, nf)

                                        res1_c = pt.fx(x, m2Q2, m2mu2, nf) + hs.fx(x, m2Q2, m2mu2, nf)
                                        res2_c = asy.fx(x, m2Q2, m2mu2, nf)

                                        np.testing.assert_allclose(res1.GetCentral(), res2.GetCentral(), rtol=1e-7)
                                        np.testing.assert_allclose(res1.GetHigher(), res2.GetHigher(), rtol=1e-7)
                                        np.testing.assert_allclose(res1.GetLower(), res2.GetLower(), rtol=1e-7)
                                        np.testing.assert_allclose(res1_c, res1.GetCentral(), rtol=1e-7)
                                        np.testing.assert_allclose(res2_c, res2.GetCentral(), rtol=1e-7)
