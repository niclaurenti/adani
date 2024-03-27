import adani as ad
import oldadani as old
import numpy as np

def test_highenergy_as2_oldadani():
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            he = ad.HighEnergyCoefficientFunction(2, kind, channel)
            for x in np.geomspace(1e-5, 1., 100, endpoint=False):
                for m2Q2 in np.geomspace(1e-4, 1e-2, 10):
                    for m2mu2 in np.geomspace(1e-4, 1e-2, 10):
                        res1 = he.fx(x, m2Q2, m2mu2, 0)
                        if (kind == '2' and channel == 'g'):
                            res2 = old.C2_g2_highenergy(x, m2Q2, m2mu2)
                        if (kind == '2' and channel == 'q'):
                            res2 = old.C2_ps2_highenergy(x, m2Q2, m2mu2)
                        if (kind == 'L' and channel == 'g'):
                            res2 = old.CL_g2_highenergy(x, m2Q2, m2mu2)
                        if (kind == 'L' and channel == 'q'):
                            res2 = old.CL_ps2_highenergy(x, m2Q2, m2mu2)
                        np.testing.assert_allclose(res1, res2, rtol=1e-7)

def test_highenergy_as3_oldadani():
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            he = ad.HighEnergyCoefficientFunction(3, kind, channel)
            for x in np.geomspace(1e-5, 1., 100, endpoint=False):
                for m2Q2 in np.geomspace(1e-4, 1e-2, 10):
                    for m2mu2 in np.geomspace(1e-4, 1e-2, 10):
                        for nf in range(6 + 1):
                            res1 = he.fx(x, m2Q2, m2mu2, nf)
                            if (kind == '2' and channel == 'g'):
                                res2 = old.C2_g3_highenergy(x, m2Q2, m2mu2, nf, 0)
                            if (kind == '2' and channel == 'q'):
                                res2 = old.C2_ps3_highenergy(x, m2Q2, m2mu2, nf, 0)
                            if (kind == 'L' and channel == 'g'):
                                res2 = old.CL_g3_highenergy(x, m2Q2, m2mu2, nf, 0)
                            if (kind == 'L' and channel == 'q'):
                                res2 = old.CL_ps3_highenergy(x, m2Q2, m2mu2, nf, 0)
                            np.testing.assert_allclose(res1, res2, rtol=1e-7)

def test_highenergyhighscale_as2_oldadani():
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            he = ad.HighEnergyHighScaleCoefficientFunction(2, kind, channel)
            for x in np.geomspace(1e-5, 1., 100, endpoint=False):
                for m2Q2 in np.geomspace(1e-4, 1e-2, 10):
                    for m2mu2 in np.geomspace(1e-4, 1e-2, 10):
                        res1 = he.fx(x, m2Q2, m2mu2, 0)
                        if (kind == '2' and channel == 'g'):
                            res2 = old.C2_g2_highenergy_highscale(x, m2Q2, m2mu2)
                        if (kind == '2' and channel == 'q'):
                            res2 = old.C2_ps2_highenergy_highscale(x, m2Q2, m2mu2)
                        if (kind == 'L' and channel == 'g'):
                            res2 = old.CL_g2_highenergy_highscale(x, m2Q2, m2mu2)
                        if (kind == 'L' and channel == 'q'):
                            res2 = old.CL_ps2_highenergy_highscale(x, m2Q2, m2mu2)
                        np.testing.assert_allclose(res1, res2, rtol=1e-7)

def test_highenergyhighscale_as3_oldadani():
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            he = ad.HighEnergyHighScaleCoefficientFunction(3, kind, channel)
            for x in np.geomspace(1e-5, 1., 100, endpoint=False):
                for m2Q2 in np.geomspace(1e-4, 1e-2, 10):
                    for m2mu2 in np.geomspace(1e-4, 1e-2, 10):
                        for nf in range(6 + 1):
                            res1 = he.fx(x, m2Q2, m2mu2, nf)
                            if (kind == '2' and channel == 'g'):
                                res2 = old.C2_g3_highenergy_highscale(x, m2Q2, m2mu2, nf, 0)
                            if (kind == '2' and channel == 'q'):
                                res2 = old.C2_ps3_highenergy_highscale(x, m2Q2, m2mu2, nf, 0)
                            if (kind == 'L' and channel == 'g'):
                                res2 = old.CL_g3_highenergy_highscale(x, m2Q2, m2mu2, nf, 0)
                            if (kind == 'L' and channel == 'q'):
                                res2 = old.CL_ps3_highenergy_highscale(x, m2Q2, m2mu2, nf, 0)
                            np.testing.assert_allclose(res1, res2, rtol=1e-7)

def test_powerterms_as2_oldadani():
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            he = ad.PowerTermsCoefficientFunction(2, kind, channel)
            for x in np.geomspace(1e-5, 1., 100, endpoint=False):
                for m2Q2 in np.geomspace(1e-4, 1e-2, 10):
                    for m2mu2 in np.geomspace(1e-4, 1e-2, 10):
                        res1 = he.fx(x, m2Q2, m2mu2, 0)
                        if (kind == '2' and channel == 'g'):
                            res2 = old.C2_g2_power_terms(x, m2Q2, m2mu2)
                        if (kind == '2' and channel == 'q'):
                            res2 = old.C2_ps2_power_terms(x, m2Q2, m2mu2)
                        if (kind == 'L' and channel == 'g'):
                            res2 = old.CL_g2_power_terms(x, m2Q2, m2mu2)
                        if (kind == 'L' and channel == 'q'):
                            res2 = old.CL_ps2_power_terms(x, m2Q2, m2mu2)
                        np.testing.assert_allclose(res1, res2, rtol=1e-7)

def test_powerterms_as3_oldadani():
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            he = ad.PowerTermsCoefficientFunction(3, kind, channel)
            for x in np.geomspace(1e-5, 1., 100, endpoint=False):
                for m2Q2 in np.geomspace(1e-4, 1e-2, 10):
                    for m2mu2 in np.geomspace(1e-4, 1e-2, 10):
                        for nf in range(6 + 1):
                            res1 = he.fx(x, m2Q2, m2mu2, nf)
                            if (kind == '2' and channel == 'g'):
                                res2 = old.C2_g3_power_terms(x, m2Q2, m2mu2, nf, 0)
                            if (kind == '2' and channel == 'q'):
                                res2 = old.C2_ps3_power_terms(x, m2Q2, m2mu2, nf, 0)
                            if (kind == 'L' and channel == 'g'):
                                res2 = old.CL_g3_power_terms(x, m2Q2, m2mu2, nf, 0)
                            if (kind == 'L' and channel == 'q'):
                                res2 = old.CL_ps3_power_terms(x, m2Q2, m2mu2, nf, 0)
                            np.testing.assert_allclose(res1, res2, rtol=1e-7)

def test_powerterms():
    for order in range(1, 3 + 1):
        for kind in ['2', 'L']:
            for channel in ['g', 'q']:
                for NLL in [True, False]:
                    he = ad.HighEnergyCoefficientFunction(order, kind, channel, NLL)
                    hehs = ad.HighEnergyHighScaleCoefficientFunction(order, kind, channel, NLL)
                    pt = ad.PowerTermsCoefficientFunction(order, kind, channel, NLL)
                    for nf in range(1, 6 + 1):
                        for m2mu2 in np.geomspace(1e-4, 1e2, 10):
                            for m2Q2 in np.geomspace(1e-4, 1e2, 10):
                                for x in np.geomspace(1e-5, 1., 10, endpoint=False):
                                    res1 = pt.fxBand(x, m2Q2, m2mu2, nf)

                                    res2_c = he.fxBand(x, m2Q2, m2mu2, nf).GetCentral() - hehs.fxBand(x, m2Q2, m2mu2, nf).GetCentral()
                                    tmp1 = he.fxBand(x, m2Q2, m2mu2, nf).GetHigher() - hehs.fxBand(x, m2Q2, m2mu2, nf).GetHigher()
                                    tmp2 = he.fxBand(x, m2Q2, m2mu2, nf).GetLower() - hehs.fxBand(x, m2Q2, m2mu2, nf).GetLower()
                                    if tmp1 > tmp2:
                                        res2_h = tmp1
                                        res2_l = tmp2
                                    else:
                                        res2_h = tmp2
                                        res2_l = tmp1

                                    np.testing.assert_allclose(res1.GetCentral(), res2_c, rtol=1e-7)
                                    np.testing.assert_allclose(res1.GetHigher(), res2_h, rtol=1e-7)
                                    np.testing.assert_allclose(res1.GetLower(), res2_l, rtol=1e-7)
                    del he
                    del pt
                    del hehs
