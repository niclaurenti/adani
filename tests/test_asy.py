import adani as ad
import numpy as np

def test_asy():
    for order in range(1, 3 + 1):
        for kind in ['2', 'L']:
            for channel in ['g', 'q']:
                if order == 1 and channel == 'q':
                    continue
                for exact in [True, False]:
                    if exact and channel == 'g':
                        continue
                    for revised_approx in [True, False]:
                        if exact and revised_approx and channel == 'g':
                            continue
                        hs = ad.HighScaleCoefficientFunction(order, kind, channel, exact, revised_approx)
                        for nll in [True, False]:
                            pt = ad.PowerTermsCoefficientFunction(order, kind, channel, nll)
                            asy = ad.AsymptoticCoefficientFunction(order, kind, channel, nll, exact, revised_approx)
                            for nf in range(1, 6 + 1):
                                for m2mu2 in np.geomspace(1e-4, 1e2, 10):
                                     for m2Q2 in np.geomspace(1e-4, 1e2, 10):
                                        for x in np.geomspace(1e-5, 1., 10, endpoint=False):
                                            res1 = pt.fxBand(x, m2Q2, m2mu2, nf) + hs.fxBand(x, m2Q2, m2mu2, nf)
                                            res2 = asy.fxBand(x, m2Q2, m2mu2, nf)

                                            np.testing.assert_allclose(res1.GetCentral(), res2.GetCentral(), rtol=1e-7)
                                            np.testing.assert_allclose(res1.GetHigher(), res2.GetHigher(), rtol=1e-7)
                                            np.testing.assert_allclose(res1.GetLower(), res2.GetLower(), rtol=1e-7)
