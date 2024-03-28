import adani as ad
import oldadani as oldad
import numpy as np

def test_mudependent_terms():
    for order in [2, 3]:
        for channel in ['g', 'q']:
            for kind in ['2', 'L']:
                for dim in [1000]:
                    for MCcalls in [20000]:
                        for mf in [False, True]:
                            for abserr in [1e-3]:
                                relerr = abserr
                                massive = ad.ExactCoefficientFunction(order, kind, channel, abserr, relerr, dim, mf, MCcalls)
                                app = ad.ApproximateCoefficientFunction(order, kind, channel, True, "klmv", abserr, relerr, dim, mf, MCcalls)
                                x = np.geomspace(1e-5, 1., 5, endpoint=True)
                                for xi in np.geomspace(1e-2, 1e4, 4, endpoint=True):
                                    for nf in [4, 5]:
                                        res1 = [massive.MuDependentTerms(x_, 1/xi, 1/xi, nf) for x_ in x]
                                        res2 = [app.MuDependentTerms(x_, 1/xi, 1/xi, nf) for x_ in x]
                                        np.testing.assert_allclose(res1, res2, rtol=1e-7)
                                del massive
                                del app

def test_as2_muindep_oldversion():
    for channel in ['g', 'q']:
        for kind in ['2', 'L']:
            app = ad.ApproximateCoefficientFunction(2, kind, channel)
            for xi in np.geomspace(1e-2, 1e2, 10):
                m2Q2 = 1/xi
                for x in np.geomspace(1e-5, 1, 10):
                    for v in [0, 1, -1]:
                        if kind == '2':
                            if channel == 'g':
                                res_old = oldad.C2_g2_approximation(x, m2Q2, 1., v)
                            if channel == 'q':
                                res_old = oldad.C2_ps2_approximation(x, m2Q2, 1., v)
                        if kind == 'L':
                            if channel == 'g':
                                res_old = oldad.CL_g2_approximation(x, m2Q2, 1., v)
                            if channel == 'q':
                                res_old = oldad.CL_ps2_approximation(x, m2Q2, 1., v)
                        tmp = app.MuIndependentTermsBand(x, m2Q2, 1)
                        if v == 0:
                            res_new = tmp.GetCentral()
                        if v == -1:
                            res_new = tmp.GetLower()
                        if v == 1:
                            res_new = tmp.GetHigher()
                        np.testing.assert_allclose(res_old, res_new, rtol=1e-7)

def test_as3_muindep_oldversion():
    for channel in ['g', 'q']:
        for kind in ['2', 'L']:
            highscale_version = "exact" if channel == 'q' else "abmp"
            app = ad.ApproximateCoefficientFunction(3, kind, channel,True, highscale_version)
            for xi in np.geomspace(1e-2, 1e2, 10):
                m2Q2 = 1/xi
                for x in np.geomspace(1e-5, 1, 10):
                    for nf in range(1, 6 + 1):
                        for v in [0, 1, -1]:
                            if kind == '2':
                                if channel == 'g':
                                    res_old = oldad.C2_g3_approximation(x, m2Q2, 1., nf, v)
                                if channel == 'q':
                                    res_old = oldad.C2_ps3_approximation(x, m2Q2, 1., nf, v)
                            if kind == 'L':
                                if channel == 'g':
                                    res_old = oldad.CL_g3_approximation(x, m2Q2, 1., nf, v)
                                if channel == 'q':
                                    res_old = oldad.CL_ps3_approximation(x, m2Q2, 1., nf, v)
                            tmp = app.MuIndependentTermsBand(x, m2Q2, nf)
                            if v == 0:
                                res_new = tmp.GetCentral()
                            if v == -1:
                                res_new = tmp.GetLower()
                            if v == 1:
                                res_new = tmp.GetHigher()
                            np.testing.assert_allclose(res_old, res_new, rtol=1e-7)

def test_mudep_as2_oldversion():
    for channel in ['g', 'q']:
        for kind in ['2', 'L']:
            app = ad.ApproximateCoefficientFunction(2, kind, channel)
            for xi in np.geomspace(1e-2, 1e2, 10):
                m2Q2 = 1/xi
                xmax = 1/(1 + 4*m2Q2)
                Lmu = -np.log(m2Q2)
                for x in np.geomspace(1e-5, 1, 10):
                    for nf in [3, 4]:
                        if kind == '2':
                            if channel == 'g':
                                res_old = oldad.C2_g21(x, m2Q2) * Lmu
                            if channel == 'q':
                                res_old = oldad.C2_ps21(x, m2Q2) * Lmu
                        if kind == 'L':
                            if channel == 'g':
                                res_old = oldad.CL_g21(x, m2Q2) * Lmu
                            if channel == 'q':
                                res_old = oldad.CL_ps21(x, m2Q2) * Lmu
                        if x > xmax:
                            res_old = 0.
                        res_new = app.MuDependentTerms(x, m2Q2, m2Q2, nf)
                        np.testing.assert_allclose(res_old, res_new, rtol=1e-7)
            del app

def test_mudep_as3_oldversion():
    for channel in ['g', 'q']:
        for kind in ['2', 'L']:
            highscale_version = "exact" if channel == 'q' else "abmp"
            for mf in [False]:
                app = ad.ApproximateCoefficientFunction(3, kind, channel, True, highscale_version, 1e-3, 1e-3, 1000, mf, 25000)
                for xi in np.geomspace(1e-2, 1e2, 10):
                    m2Q2 = 1/xi
                    xmax = 1/(1 + 4*m2Q2)
                    Lmu = -np.log(m2Q2)
                    for x in np.geomspace(1e-5, 1, 10):
                        for nf in [3, 4]:
                            if kind == '2':
                                if channel == 'g':
                                    res_old = oldad.C2_g31(x, m2Q2, nf) * Lmu + oldad.C2_g32(x, m2Q2, nf, mf) * Lmu**2
                                if channel == 'q':
                                    res_old = oldad.C2_ps31(x, m2Q2, nf) * Lmu + oldad.C2_ps32(x, m2Q2, nf) * Lmu**2
                            if kind == 'L':
                                if channel == 'g':
                                    res_old = oldad.CL_g31(x, m2Q2, nf) * Lmu + oldad.CL_g32(x, m2Q2, nf, mf) * Lmu**2
                                if channel == 'q':
                                    res_old = oldad.CL_ps31(x, m2Q2, nf) * Lmu + oldad.CL_ps32(x, m2Q2, nf) * Lmu**2
                            if x > xmax:
                                res_old = 0.
                            res_new = app.MuDependentTerms(x, m2Q2, m2Q2, nf)
                            np.testing.assert_allclose(res_old, res_new, rtol=1e-7)
                del app

def test_klmv_as2():
    for channel in ['g', 'q']:
        app = ad.ApproximateCoefficientFunctionKLMV(2, '2', channel, "abmp")
        for xi in np.geomspace(1e-2, 1e2, 10):
            m2Q2 = 1/xi
            for x in np.geomspace(1e-5, 1, 10):
                if channel == 'g':
                    tmpA = oldad.C2_g2_approximationA_klmv(x, m2Q2, 1.)
                    tmpB = oldad.C2_g2_approximationB_klmv(x, m2Q2, 1.)
                if channel == 'q':
                    tmpA = oldad.C2_ps2_approximationA_klmv(x, m2Q2, 1.)
                    tmpB = oldad.C2_ps2_approximationB_klmv(x, m2Q2, 1.)
                if tmpA > tmpB:
                    higher = tmpA
                    lower = tmpB
                else:
                    higher = tmpB
                    lower = tmpA
                tmp = app.MuIndependentTermsBand(x, m2Q2, 1)
                np.testing.assert_allclose(higher, tmp.GetHigher(), rtol=1e-7)
                np.testing.assert_allclose(lower, tmp.GetLower(), rtol=1e-7)

def test_klmv_as3():
    for channel in ['g', 'q']:
        hs_version = "abmp" if channel == 'g' else 'exact'
        app = ad.ApproximateCoefficientFunctionKLMV(3, '2', channel, hs_version)
        for xi in np.geomspace(1e-2, 1e2, 10):
            m2Q2 = 1/xi
            for x in np.geomspace(1e-5, 1, 10):
                for nf in [4, 5, 6]:
                    if channel == 'g':
                        tmpA = oldad.C2_g3_approximationA_klmv(x, m2Q2, 1., nf, 0)
                        tmpB = oldad.C2_g3_approximationB_klmv(x, m2Q2, 1., nf, 0)
                    if channel == 'q':
                        tmpA = oldad.C2_ps3_approximationA_klmv(x, m2Q2, 1., nf)
                        tmpB = oldad.C2_ps3_approximationB_klmv(x, m2Q2, 1., nf)
                    if tmpA > tmpB:
                        higher = tmpA
                        lower = tmpB
                    else:
                        higher = tmpB
                        lower = tmpA
                    tmp = app.MuIndependentTermsBand(x, m2Q2, nf)
                    np.testing.assert_allclose(higher, tmp.GetHigher(), rtol=1e-7)
                    np.testing.assert_allclose(lower, tmp.GetLower(), rtol=1e-7)

def test_klmv_paper_as3():
    for channel in ['g', 'q']:
        app = ad.ApproximateCoefficientFunctionKLMV(3, '2', channel, "klmv")
        for xi in np.geomspace(1e-2, 1e2, 10):
            m2Q2 = 1/xi
            for x in np.geomspace(1e-5, 1, 10):
                for nf in [4, 5, 6]:
                    if channel == 'g':
                        tmpA = oldad.C2_g3_approximationA_klmv(x, m2Q2, 1., nf, 0)
                        tmpB = oldad.C2_g3_approximationB_klmv_paper(x, m2Q2, 1., nf, 0)
                    if channel == 'q':
                        tmpA = oldad.C2_ps3_approximationA_klmv_paper(x, m2Q2, 1., nf)
                        tmpB = oldad.C2_ps3_approximationB_klmv_paper(x, m2Q2, 1., nf)
                    if tmpA > tmpB:
                        higher = tmpA
                        lower = tmpB
                    else:
                        higher = tmpB
                        lower = tmpA
                    tmp = app.MuIndependentTermsBand(x, m2Q2, nf)
                    np.testing.assert_allclose(higher, tmp.GetHigher(), rtol=1e-7)
                    np.testing.assert_allclose(lower, tmp.GetLower(), rtol=1e-7)

def test_klmv_lowxi_as3():
    app = ad.ApproximateCoefficientFunctionKLMV(3, '2', 'g', "klmv", True)
    for xi in np.geomspace(1e-2, 1e2, 10):
        m2Q2 = 1/xi
        for x in np.geomspace(1e-5, 1, 10):
            for nf in [4, 5, 6]:
                tmpA = oldad.C2_g3_approximationA_klmv(x, m2Q2, 1., nf, 0)
                tmpB = oldad.C2_g3_approximationBlowxi_klmv(x, m2Q2, 1., nf, 0)
                if tmpA > tmpB:
                    higher = tmpA
                    lower = tmpB
                else:
                    higher = tmpB
                    lower = tmpA
                tmp = app.MuIndependentTermsBand(x, m2Q2, nf)
                np.testing.assert_allclose(higher, tmp.GetHigher(), rtol=1e-7)
                np.testing.assert_allclose(lower, tmp.GetLower(), rtol=1e-7)
