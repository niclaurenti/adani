from adani import Exact, GM, ABMP, KLMV
import adani as ad
import numpy as np

def test_splitlogs_vs_highscale():
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            for hs_version in [Exact, ABMP, KLMV]:
                if channel == 'q' and hs_version == ABMP:
                    continue
                if channel == 'g' and hs_version == Exact:
                    continue
                hs = ad.HighScaleCoefficientFunction(3, kind, channel, hs_version)
                hs_split = ad.HighScaleSplitLogs(3, kind, channel, hs_version)
                for x in np.geomspace(1e-5, 1, 100, endpoint=False):
                    for xi in np.geomspace(1e-2, 1e4, 10, endpoint=True):
                        for nf in range(1, 6 + 1):
                            res1 = hs.fxBand(x, 1./xi, 1./xi, nf)
                            res2 = hs_split.fxBand(x, 1./xi, nf)
                            vec1 = [res1.GetCentral(), res1.GetHigher(), res1.GetLower()]
                            vec2 = [res2.GetCentral(), res2.GetHigher(), res2.GetLower()]
                            np.testing.assert_allclose(vec1, vec2, rtol = 1e-4)
                            central1 = hs.fx(x, 1./xi, 1./xi, nf)
                            central2 = hs_split.fx(x, 1./xi, nf)
                            np.testing.assert_allclose(central1, vec1[0], rtol = 1e-4)
                            np.testing.assert_allclose(central2, vec2[0], rtol = 1e-4)

def test_split_logs():
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            for hs_version in [Exact, ABMP, KLMV]:
                if channel == 'q' and hs_version == ABMP:
                    continue
                if channel == 'g' and hs_version == Exact:
                    continue
                hs_split = ad.HighScaleSplitLogs(3, kind, channel, hs_version)
                for x in np.geomspace(1e-5, 1, 100, endpoint=False):
                    for xi in np.geomspace(1e-2, 1e4, 10, endpoint=True):
                        for nf in range(1, 6 + 1):
                            logQ = np.log(1/xi)
                            res1 = (
                                hs_split.N3LL(x, nf)
                                + hs_split.N2LL(x, nf) * logQ
                                + hs_split.NLL(x, nf) * logQ**2
                                + hs_split.LL(x, nf) * logQ**3
                            )
                            res2 = hs_split.fxBand(x, 1./xi, nf)
                            vec1 = [res1.GetCentral(), res1.GetHigher(), res1.GetLower()]
                            vec2 = [res2.GetCentral(), res2.GetHigher(), res2.GetLower()]
                            np.testing.assert_allclose(vec1, vec2, rtol = 1e-7)
                            central1 = (
                                hs_split.LL(x, nf) * logQ**3
                                + hs_split.NLL(x, nf) * logQ**2
                                + hs_split.N2LL(x, nf) * logQ
                                + hs_split.N3LL(x, nf).GetCentral()
                            )
                            central2 = hs_split.fx(x, 1./xi, nf)
                            np.testing.assert_allclose(central1, vec1[0], rtol = 1e-7)
                            np.testing.assert_allclose(central2, vec2[0], rtol = 1e-7)
