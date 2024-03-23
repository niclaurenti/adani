import adani as ad
import oldadani as old
import numpy as np
from scipy.special import spence
from mpmath import polylog, nsum, sin

TR = 0.5
CF = 4./3
CA = 3

ln2 = 0.6931471805599453
Li2_1_2 = 0.5822405265
Li3_1_2 = 0.5372131936

zeta2 = 1.6449340668482264
zeta3 = 1.2020569031595942


def Li2(x):
    return spence(1-x)

def Li3(x):
    return float(polylog(3, x).real)

def S12(x):
    if (x > 0 and x < 1):
        return -Li3(1. - x) + zeta3 + np.log(1. - x) * Li2(1. - x) + 0.5 * np.log(x) * np.log(1. - x) * np.log(1. - x)
    elif x < 0:
        c = 1. / (1. - x)
        return (
            -Li3(c) + zeta3 + np.log(c) * Li2(c)
            + 0.5 * np.log(c) * np.log(c) * np.log(1. - c)
            - 1. / 6. * np.log(c) * np.log(c) * np.log(c)
        )
    elif x == 0:
        return 0.
    elif x ==1:
        return zeta3
    else:
        raise ValueError("x > 1 in S12 not implemented!")
def H_1(x, i):

    if (i == 1):
        return -np.log(1. - x)
    if (i == 0):
        return np.log(x)
    if (i == -1):
        return np.log(1. + x)

def H_2(x, k):

    (i,j) = k

    if (i == 0 and j == 0):
        L = np.log(x)
        return 0.5 * L * L
    if (i == 0 and j == 1) :
        return Li2(x)
    if (i == 0 and j == -1):
        return -Li2(-x)
    if (i == 1 and j == 0):
        return -np.log(x) * np.log(1. - x) - Li2(x)
    # there's a typo in the paper: +Li2(x)
    if (i == 1 and j == 1):
        Lm = np.log(1. - x)
        return 0.5 * Lm * Lm
    if (i == 1 and j == -1):
        xm = 1. - x
        return Li2(xm / 2.) - ln2 * np.log(xm) - Li2_1_2
    if (i == -1 and j == 0):
        return np.log(x) * np.log(1. + x) + Li2(-x)
    if (i == -1 and j == 1) :
        xp = 1. + x
        return Li2(xp / 2) - ln2 * np.log(xp) - Li2_1_2
    if (i == -1 and j == -1):
        Lp = np.log(1. + x)
        return 0.5 * Lp * Lp

def H_3(x, l):
    (i, j, k) = l
    if (i == 1 and j == 1 and k == 1) :
        Lm = np.log(1. - x)
        return -1. / 6. * Lm * Lm * Lm

    if (i == 1 and j == 1 and k == 0) :
        xm = 1. - x
        return zeta2 * np.log(xm) - Li3(xm) + zeta3

    if (i == 1 and j == 1 and k == -1):
        return 0.5 * H_1(x, 1) * H_2(x, (1, -1)) - 0.5 * H_3(x, (1, -1, 1))
    if (i == 1 and j == 0 and k == 1):
        return -2. * S12(x) - np.log(1. - x) * Li2(x)
    if (i == 1 and j == 0 and k == 0):
        L = np.log(x)
        return -1. / 2 * L * L * np.log(1. - x) - L * Li2(x) + Li3(x)

    if (i == 1 and j == 0 and k == -1):
        return H_1(x, 1) * H_2(x, (0, -1)) - H_3(x, (0, 1, -1)) - H_3(x, (0, -1, 1))
    if (i == 1 and j == -1 and k == 1):
        return H_1(x, 1) * H_2(x, (-1, 1)) - 2. * H_3(x, (-1, 1, 1))
    if (i == 1 and j == -1 and k == 0):
        return H_1(x, 0) * H_2(x, (1, -1)) - H_3(x, (0, 1, -1)) - H_3(x, (1, 0, -1))
    if (i == 1 and j == -1 and k == -1):
        xp = 1. + x
        Lp = np.log(xp)
        return (-0.5 * np.log((1. - x) / 2.) * Lp * Lp - Li3_1_2 - Lp * Li2(xp / 2.)
               + Li3(xp / 2.))


    if (i == 0 and j == 1 and k == 1):
        return S12(x)
    if (i == 0 and j == 1 and k == 0):
        return np.log(x) * Li2(x) - 2. * Li3(x)
    if (i == 0 and j == 1 and k == -1):
        xp = 1. + x
        Lp = np.log(xp)
        return (Li3(2. * x / xp) - Li3(x / xp) - Li3(xp / 2.) - Li3(x)
               + Lp * Li2(0.5) + Lp * Li2(x) + 0.5 * ln2 * Lp * Lp + Li3_1_2)

    if (i == 0 and j == 0 and k == 1):
        return Li3(x)
    if (i == 0 and j == 0 and k == 0):
        L = np.log(x)
        return 1. / 6. * L * L * L

    if (i == 0 and j == 0 and k == -1):
        return -Li3(-x)
    if (i == 0 and j == -1 and k == 1):
        xm = 1. - x
        Lm = np.log(xm)
        Lm2 = Lm * Lm
        Lm3 = Lm2 * Lm
        return (-S12(x) + Li3(-2. * x / xm) - Li3(xm / 2) - Li3(-x) + Li3_1_2
               + Li3(x) + Lm * Li2(-x) + Lm * Li2_1_2 - Lm * Li2(x)
               + 0.5 * ln2 * Lm2 - 1. / 6. * Lm3)

    if (i == 0 and j == -1 and k == 0):
        return -np.log(x) * Li2(-x) + 2. * Li3(-x)
    if (i == 0 and j == -1 and k == -1):
        return S12(-x)

    if (i == -1 and j == 1 and k == 1):
        xm = 1. - x
        Lm = np.log(xm)
        return (0.5 * np.log((1. + x) / 2.) * Lm * Lm + Li3(0.5) + Lm * Li2(xm / 2.)
               - Li3(xm / 2))

    if (i == -1 and j == 1 and k == 0):
        return H_1(x, 0) * H_2(x, (-1, 1)) - H_3(x, (0, -1, 1)) - H_3(x, (-1, 0, 1))
    if (i == -1 and j == 1 and k == -1):
        return H_1(x, -1) * H_2(x, (1, -1)) - 2. * H_3(x, (1, -1, -1))
    if (i == -1 and j == 0 and k == 1):
        return H_1(x, -1) * H_2(x, (0, 1)) - H_3(x, (0, -1, 1)) - H_3(x, (0, 1, -1))
    if (i == -1 and j == 0 and k == 0):
        L = np.log(x)
        return 0.5 * L * L * np.log(1. + x) + L * Li2(-x) - Li3(-x)

    if (i == -1 and j == 0 and k == -1):
        return H_2(x, (0, -1)) * H_1(x, -1) - 2. * H_3(x, (0, -1, -1))
    if (i == -1 and j == -1 and k == 1):
        return 0.5 * H_1(x, -1) * H_2(x, (-1, 1)) - 0.5 * H_3(x, (-1, 1, -1))
    if (i == -1 and j == -1 and k == 0):
        return H_1(x, 0) * H_2(x, (-1, -1)) - H_3(x, (0, -1, -1)) - H_3(x, (-1, 0, -1))
    if (i == -1 and j == -1 and k == -1):
        Lp = np.log(1. + x)
        return 1. / 6. * Lp * Lp * Lp


def C2_g1(x, nf):
    return (4. * nf * TR
           * (-8. * x * x + 8. * x - 1.
              + np.log((1. - x) / x) * (2. * x * x - 2. * x + 1.)))

def CL_g1(x, nf):
    return 16. * nf * TR * x * (1. - x)

def test_lo():
    for kind in ['2', 'L']:
        ml = ad.MasslessCoefficientFunction(1, kind, 'g')
        for x in np.geomspace(1e-5, 1., 100, endpoint=False):
            for nf in range(6 + 1):
                res1 = ml.MuIndependentTerms(x, nf)
                res2 = C2_g1(x, nf) if kind == '2' else CL_g1(x, nf)
                np.testing.assert_allclose(res1, res2, rtol=1e-7)

def C2_g2(x, nf):
    x2 = x * x
    x3 = x2 * x

    Hm1 = H_1(x, -1)
    H0 = H_1(x, 0)
    H1 = H_1(x, 1)

    H00 = H_2(x, (0, 0))
    H10 = H_2(x, (1, 0))
    Hm10 = H_2(x, (-1, 0))
    H01 = H_2(x, (0, 1))
    H11 = H_2(x, (1, 1))

    Hm1m10 = H_3(x, (-1, -1, 0))
    H0m10 = H_3(x, (0, -1, 0))
    Hm100 = H_3(x, (-1, 0, 0))
    H000 = H_3(x, (0, 0, 0))
    H100 = H_3(x, (1, 0, 0))
    H010 = H_3(x, (0, 1, 0))
    H110 = H_3(x, (1, 1, 0))
    Hm101 = H_3(x, (-1, 0, 1))
    H001 = H_3(x, (0, 0, 1))
    H101 = H_3(x, (1, 0, 1))
    H011 = H_3(x, (0, 1, 1))
    H111 = H_3(x, (1, 1, 1))

    return (nf * CF
               * (-647. / 15 - 104. / 3 * zeta2 * x + 72 * zeta2 * x2
                  + 96. / 5 * zeta2 * x3 + 16 * zeta2 + 72 * zeta3 * x2
                  + 32 * zeta3 + 8. / 15 / x + 239. / 5 * x - 36. / 5 * x2
                  + 32 * H0m10 + 32 * H0m10 * x2 - 32 * Hm1 * zeta2 * x
                  - 16 * Hm1 * zeta2 * x2 - 16 * Hm1 * zeta2 - 32 * Hm1m10
                  - 64 * Hm1m10 * x - 32 * Hm1m10 * x2 + 48 * Hm10
                  + 8. / 15 * Hm10 / x2 + 64. / 3 * Hm10 * x
                  + 96. / 5 * Hm10 * x3 + 16 * Hm100 + 32 * Hm100 * x
                  + 16 * Hm100 * x2 - 236. / 15 * H0 - 32 * H0 * zeta2 * x
                  + 48 * H0 * zeta2 * x2 + 16 * H0 * zeta2 - 8. / 15 * H0 / x
                  + 113. / 5 * H0 * x - 216. / 5 * H0 * x2 - 3 * H00
                  + 44. / 3 * H00 * x - 72 * H00 * x2 - 96. / 5 * H00 * x3
                  - 10 * H000 + 20 * H000 * x - 40 * H000 * x2 - 14 * H1
                  - 16 * H1 * zeta2 * x + 32 * H1 * zeta2 * x2 + 8 * H1 * zeta2
                  + 40 * H1 * x - 24 * H1 * x2 - 26 * H10 + 80 * H10 * x
                  - 72 * H10 * x2 - 4 * H100 + 8 * H100 * x - 24 * H100 * x2
                  - 26 * H11 + 80 * H11 * x - 72 * H11 * x2 - 16 * H110
                  + 32 * H110 * x - 32 * H110 * x2 - 20 * H111 + 40 * H111 * x
                  - 40 * H111 * x2 - 24 * H101 + 48 * H101 * x - 48 * H101 * x2
                  - 16 * H01 + 56 * H01 * x - 72 * H01 * x2 - 12 * H010
                  + 24 * H010 * x - 32 * H010 * x2 - 16 * H011 + 32 * H011 * x
                  - 40 * H011 * x2 - 16 * H001 + 32 * H001 * x - 48 * H001 * x2)
           + nf * CA
                 * (+239. / 9 - 16. / 3 * zeta2 / x - 144 * zeta2 * x
                    + 148 * zeta2 * x2 + 8 * zeta2 - 48 * zeta3 * x
                    + 24 * zeta3 * x2 + 4 * zeta3 + 344. / 27 / x
                    + 1072. / 9 * x - 4493. / 27 * x2 + 16 * H0m10 * x2
                    - 8 * Hm1 * zeta2 * x - 16 * Hm1 * zeta2 * x2
                    - 4 * Hm1 * zeta2 + 8 * Hm1m10 + 16 * Hm1m10 * x - 24 * Hm10
                    - 16. / 3 * Hm10 / x + 80. / 3 * Hm10 * x2 + 8 * Hm100
                    + 16 * Hm100 * x + 24 * Hm100 * x2 + 8 * Hm101
                    + 16 * Hm101 * x + 16 * Hm101 * x2 + 58 * H0
                    - 64 * H0 * zeta2 * x + 16 * H0 * zeta2 * x2
                    - 8 * H0 * zeta2 + 584. / 3 * H0 * x - 2090. / 9 * H0 * x2
                    - 2 * H00 + 176 * H00 * x - 388. / 3 * H00 * x2 + 20 * H000
                    + 56 * H000 * x + 62. / 3 * H1 - 16 * H1 * zeta2 * x
                    + 8 * H1 * zeta2 * x2 + 8 * H1 * zeta2 - 104. / 9 * H1 / x
                    + 454. / 3 * H1 * x - 1570. / 9 * H1 * x2 - 4 * H10
                    + 16. / 3 * H10 / x + 80 * H10 * x - 268. / 3 * H10 * x2
                    - 12 * H100 + 24 * H100 * x - 16 * H100 * x2 - 4 * H11
                    + 16. / 3 * H11 / x + 72 * H11 * x - 244. / 3 * H11 * x2
                    - 12 * H110 + 24 * H110 * x - 24 * H110 * x2 - 4 * H111
                    + 8 * H111 * x - 8 * H111 * x2 - 4 * H101 + 8 * H101 * x
                    - 8 * H101 * x2 - 8 * H01 + 144 * H01 * x - 148 * H01 * x2
                    + 48 * H010 * x - 16 * H010 * x2 + 48 * H011 * x
                    - 16 * H011 * x2 + 8 * H001 + 64 * H001 * x
                    - 16 * H001 * x2))

def C2_ps2(x, nf):
    x2 = x * x

    H0 = H_1(x, 0)
    H1 = H_1(x, 1)

    H00 = H_2(x, (0, 0))
    H10 = H_2(x, (1, 0))
    Hm10 = H_2(x, (-1, 0))
    H01 = H_2(x, (0, 1))
    H11 = H_2(x, (1, 1))

    H000 = H_3(x, (0, 0, 0))
    H010 = H_3(x, (0, 1, 0))
    H001 = H_3(x, (0, 0, 1))
    H011 = H_3(x, (0, 1, 1))

    return (nf * CF
           * (+158. / 9 - 16. / 3 * zeta2 / x - 16 * zeta2 * x + 16 * zeta2 * x2
              - 8 * zeta3 * x - 8 * zeta3 + 344. / 27 / x - 422. / 9 * x
              + 448. / 27 * x2 - 16 * Hm10 - 16. / 3 * Hm10 / x - 16 * Hm10 * x
              - 16. / 3 * Hm10 * x2 + 56 * H0 - 16 * H0 * zeta2 * x
              - 16 * H0 * zeta2 - 88. / 3 * H0 * x - 128. / 9 * H0 * x2
              - 2 * H00 + 30 * H00 * x - 64. / 3 * H00 * x2 + 20 * H000
              + 20 * H000 * x + 104. / 3 * H1 - 104. / 9 * H1 / x
              - 80. / 3 * H1 * x + 32. / 9 * H1 * x2 + 4 * H10
              + 16. / 3 * H10 / x - 4 * H10 * x - 16. / 3 * H10 * x2 + 4 * H11
              + 16. / 3 * H11 / x - 4 * H11 * x - 16. / 3 * H11 * x2
              - 16 * H01 * x2 + 8 * H010 + 8 * H010 * x + 8 * H011
              + 8 * H011 * x + 16 * H001 + 16 * H001 * x))


def CL_g2(x, nf):
    x2 = x * x
    x3 = x2 * x

    H0 = H_1(x, 0)
    H1 = H_1(x, 1)

    H00 = H_2(x, (0, 0))
    H10 = H_2(x, (1, 0))
    Hm10 = H_2(x, (-1, 0))
    H01 = H_2(x, (0, 1))
    H11 = H_2(x, (1, 1))

    return (nf * CF
               * (-128. / 15 + 16. / 3 * zeta2 * x + 64. / 5 * zeta2 * x3
                  + 32. / 15 / x - 304. / 5 * x + 336. / 5 * x2
                  + 32. / 15 * Hm10 / x2 - 32. / 3 * Hm10 * x
                  + 64. / 5 * Hm10 * x3 - 104. / 15 * H0 - 32. / 15 * H0 / x
                  - 208. / 5 * H0 * x + 96. / 5 * H0 * x2 - 64. / 3 * H00 * x
                  - 64. / 5 * H00 * x3 - 8 * H1 - 24 * H1 * x + 32 * H1 * x2
                  - 16 * H01 * x)
           + nf * CA
                 * (+16. / 3 - 64 * zeta2 * x + 32 * zeta2 * x2 - 16. / 9 / x
                    + 272. / 3 * x - 848. / 9 * x2 + 32 * Hm10 * x
                    + 32 * Hm10 * x2 + 16 * H0 + 128 * H0 * x - 208 * H0 * x2
                    + 96 * H00 * x + 16 * H1 - 16. / 3 * H1 / x + 144 * H1 * x
                    - 464. / 3 * H1 * x2 + 32 * H10 * x - 32 * H10 * x2
                    + 32 * H11 * x - 32 * H11 * x2 + 96 * H01 * x
                    - 32 * H01 * x2))

def CL_ps2(x, nf):
    x2 = x * x

    H0 = H_1(x, 0)
    H1 = H_1(x, 1)

    H00 = H_2(x, (0, 0))
    H01 = H_2(x, (0, 1))

    return (nf * CF
           * (+16. / 3 - 16 * zeta2 * x - 16. / 9 / x - 64. / 3 * x
              + 160. / 9 * x2 + 16 * H0 - 16 * H0 * x - 32 * H0 * x2
              + 32 * H00 * x + 16 * H1 - 16. / 3 * H1 / x - 32. / 3 * H1 * x2
              + 16 * H01 * x))

def test_nnlo():
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            ml = ad.MasslessCoefficientFunction(2, kind, channel)
            for x in np.geomspace(1e-5, 1., 100, endpoint=False):
                for nf in range(6 + 1):
                    res1 = ml.MuIndependentTerms(x, nf)
                    if (kind == '2' and channel == 'g'):
                        res2 = C2_g2(x, nf)
                    if (kind == '2' and channel == 'q'):
                        res2 = C2_ps2(x, nf)
                    if (kind == 'L' and channel == 'g'):
                        res2 = CL_g2(x, nf)
                    if (kind == 'L' and channel == 'q'):
                        res2 = CL_ps2(x, nf)
                    np.testing.assert_allclose(res1, res2, rtol=1e-7)


def test_n3lo():
    for kind in ['2', 'L']:
        for channel in ['g', 'q']:
            ml = ad.MasslessCoefficientFunction(3, kind, channel)
            for x in np.geomspace(1e-5, 1., 100, endpoint=False):
                for nf in range(6 + 1):
                    res1 = ml.MuIndependentTerms(x, nf)
                    if (kind == '2' and channel == 'g'):
                        res2 = old.C2_g3_massless(x, nf)
                    if (kind == '2' and channel == 'q'):
                        res2 = old.C2_ps3_massless(x, nf)
                    if (kind == 'L' and channel == 'g'):
                        res2 = old.CL_g3_massless(x, nf)
                    if (kind == 'L' and channel == 'q'):
                        res2 = old.CL_ps3_massless(x, nf)
                    np.testing.assert_allclose(res1, res2, rtol=1e-7)
