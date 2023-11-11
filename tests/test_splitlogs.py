import adani as ad
import numpy as np

nf = 4
m = 5

xmin = 1e-6
qmin = 1.
qmax = 1000
npoints = 100

def test_C2_g3():
    for x in np.geomspace(xmin, 1., npoints):
        for q in np.geomspace(qmin, qmax, npoints):
            m2Q2 = m**2/q**2
            log = np.log(m2Q2)
            total = ad.C2_g3_highscale(x, m2Q2, m2Q2, nf)
            splitted_logs = (
                ad.C2_g3_highscale_LL(x, nf) * log**3
                + ad.C2_g3_highscale_NLL(x, nf) * log**2
                + ad.C2_g3_highscale_N2LL(x, nf) * log**1
                + ad.C2_g3_highscale_N3LL(x, nf, 0) * log**0
            )
            np.testing.assert_allclose(total, splitted_logs, rtol=1e-8)

def test_C2_ps3():
    for x in np.geomspace(xmin, 1., npoints):
        for q in np.geomspace(qmin, qmax, npoints):
            m2Q2 = m**2/q**2
            log = np.log(m2Q2)
            total = ad.C2_ps3_highscale(x, m2Q2, m2Q2, nf)
            splitted_logs = (
                ad.C2_ps3_highscale_LL(x, nf) * log**3
                + ad.C2_ps3_highscale_NLL(x, nf) * log**2
                + ad.C2_ps3_highscale_N2LL(x, nf) * log**1
                + ad.C2_ps3_highscale_N3LL(x, nf) * log**0
            )
            np.testing.assert_allclose(total, splitted_logs, rtol=1e-8)

def test_CL_g3():
    for x in np.geomspace(xmin, 1., npoints):
        for q in np.geomspace(qmin, qmax, npoints):
            m2Q2 = m**2/q**2
            log = np.log(m2Q2)
            total = ad.CL_g3_highscale(x, m2Q2, m2Q2, nf)
            splitted_logs = (
                + ad.CL_g3_highscale_NLL(x) * log**2
                + ad.CL_g3_highscale_N2LL(x, nf) * log**1
                + ad.CL_g3_highscale_N3LL(x, nf) * log**0
            )
            np.testing.assert_allclose(total, splitted_logs, rtol=1e-8)

def test_CL_ps3():
    for x in np.geomspace(xmin, 1., npoints):
        for q in np.geomspace(qmin, qmax, npoints):
            m2Q2 = m**2/q**2
            log = np.log(m2Q2)
            total = ad.CL_ps3_highscale(x, m2Q2, m2Q2, nf)
            splitted_logs = (
                + ad.CL_ps3_highscale_NLL(x) * log**2
                + ad.CL_ps3_highscale_N2LL(x) * log**1
                + ad.CL_ps3_highscale_N3LL(x, nf) * log**0
            )
            np.testing.assert_allclose(total, splitted_logs, rtol=1e-8)
