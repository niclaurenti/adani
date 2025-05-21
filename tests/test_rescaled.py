import adani as ad
import numpy as np

def test_as1_rescaled():
    for kind in ['2', 'L']:
        hs = ad.HighScaleCoefficientFunction(1, kind, 'g')
        hsr = ad.HighScaleRescaledCoefficientFunction(1, kind, 'g')
        for x in np.geomspace(1e-5, 1., 100, endpoint=False):
            for m2Q2 in np.geomspace(1e-4, 1e-2, 10):
                res1 = hs.fx(x / (1 - 4 * m2Q2 * x), m2Q2, 0., 0)
                res2 = hsr.fx(x, m2Q2, 0., 0)
                np.testing.assert_allclose(res1, res2, rtol=1e-7)

def test_as2_rescaled():
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            hs = ad.HighScaleCoefficientFunction(2, kind, channel)
            hsr = ad.HighScaleRescaledCoefficientFunction(2, kind, channel)
            for x in np.geomspace(1e-5, 1., 100, endpoint=False):
                for m2Q2 in np.geomspace(1e-4, 1e-2, 10):
                    for m2mu2 in np.geomspace(1e-4, 1e-2, 10):
                        res1 = hs.fx(x / (1 - 4 * m2Q2 * x), m2Q2, m2mu2, 0)
                        res2 = hsr.fx(x, m2Q2, m2mu2, 0)
                        np.testing.assert_allclose(res1, res2, rtol=1e-7)

def test_as3_rescaled():
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            for hs_version in ["exact", "abmp", "klmv"]:
                if channel == 'q' and hs_version == "abmp":
                    continue
                hs = ad.HighScaleCoefficientFunction(3, kind, channel, hs_version)
                hsr = ad.HighScaleRescaledCoefficientFunction(3, kind, channel, hs_version)
                for m2Q2 in np.geomspace(1e-2, 1e4, 10):
                    for m2mu2 in np.geomspace(1e-2, 1e4, 10):
                        for nf in range(1, 6+1):
                            for x in np.geomspace(1e-5, 1., 100, endpoint=False):
                                res1 = hs.fxBand(x / (1 - 4 * m2Q2 * x), m2Q2, m2mu2, nf)
                                res2 = hsr.fxBand(x, m2Q2, m2mu2, nf)
                                np.testing.assert_allclose(res1.GetCentral(), res2.GetCentral(), rtol=1e-7)
                                np.testing.assert_allclose(res1.GetHigher(), res2.GetHigher(), rtol=1e-7)
                                np.testing.assert_allclose(res1.GetLower(), res2.GetLower(), rtol=1e-7)
