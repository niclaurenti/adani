import adani as ad
import oldadani as old
import numpy as np

def test_as1_oldadani():
    for kind in ['2', 'L']:
        hs = ad.HighScaleCoefficientFunction(1, kind, 'g')
        for x in np.geomspace(1e-5, 1., 100, endpoint=False):
            for m2Q2 in np.geomspace(1e-4, 1e-2, 10):
                res1 = hs.fx(x, m2Q2, 0., 0)
                res2 = old.C2_g1_highscale(x, m2Q2) if kind == '2' else old.CL_g1_highscale(x)
                np.testing.assert_allclose(res1, res2, rtol=1e-7)

def test_as2_oldadani():
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            hs = ad.HighScaleCoefficientFunction(2, kind, channel)
            for x in np.geomspace(1e-5, 1., 100, endpoint=False):
                for m2Q2 in np.geomspace(1e-4, 1e-2, 10):
                    for m2mu2 in np.geomspace(1e-4, 1e-2, 10):
                        res1 = hs.fx(x, m2Q2, m2mu2, 0)
                        if (kind == '2' and channel == 'g'):
                            res2 = old.C2_g2_highscale(x, m2Q2, m2mu2)
                        if (kind == '2' and channel == 'q'):
                            res2 = old.C2_ps2_highscale(x, m2Q2, m2mu2)
                        if (kind == 'L' and channel == 'g'):
                            res2 = old.CL_g2_highscale(x, m2Q2, m2mu2)
                        if (kind == 'L' and channel == 'q'):
                            res2 = old.CL_ps2_highscale(x, m2Q2, m2mu2)
                        np.testing.assert_allclose(res1, res2, rtol=1e-7)

def test_as3_oldadani():
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            version = "exact" if channel == 'q' else "abmp"
            hs = ad.HighScaleCoefficientFunction(3, kind, channel, version)
            for x in np.geomspace(1e-5, 1., 100, endpoint=False):
                for m2Q2 in np.geomspace(1e-4, 1e-2, 10):
                    for m2mu2 in np.geomspace(1e-4, 1e-2, 10):
                        for nf in range(6 + 1):
                            res1 = hs.fx(x, m2Q2, m2mu2, nf)
                            if (kind == '2' and channel == 'g'):
                                res2 = old.C2_g3_highscale(x, m2Q2, m2mu2, nf)
                            if (kind == '2' and channel == 'q'):
                                res2 = old.C2_ps3_highscale(x, m2Q2, m2mu2, nf)
                            if (kind == 'L' and channel == 'g'):
                                res2 = old.CL_g3_highscale(x, m2Q2, m2mu2, nf)
                            if (kind == 'L' and channel == 'q'):
                                res2 = old.CL_ps3_highscale(x, m2Q2, m2mu2, nf)
                            np.testing.assert_allclose(res1, res2, rtol=1e-7)

def test_splitlogs_vs_highscale():
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            for hs_version in ["exact", "abmp", "klmv"]:
                if channel == 'q' and hs_version == "abmp":
                    continue
                if channel == 'g' and hs_version == "exact":
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
            for hs_version in ["exact", "abmp", "klmv"]:
                if channel == 'q' and hs_version == "abmp":
                    continue
                if channel == 'g' and hs_version == "exact":
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
