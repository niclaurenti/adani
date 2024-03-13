import adani as ad
import numpy as np

def test_split_logs():
    for kind in ['2', 'L']:
        for channel in ['q', 'g']:
            hs = ad.HighScaleCoefficientFunction(3, kind, channel)
            hs_split = ad.HighScaleSplitLogs(3, kind, channel)
            
            for x in np.geomspace(1e-5, 1, 100, endpoint=False):
                for xi in np.geomspace(1e-2, 1e4, 10, endpoint=True):
                    for nf in range(1, 6 + 1):
                        for v in [0]:
                            res1 = hs.fx(x, 1./xi, 1./xi, nf)
                            res2 = hs_split.fx(x, 1./xi, nf, v)
                            np.testing.assert_allclose(res1, res2, rtol = 1e-7)