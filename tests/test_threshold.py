import adani as ad
import numpy as np
from test_exact_coeff_func import C2_g1, CL_g1
from scipy.special import spence

TR = 0.5
CA = 3.
CF = 4. / 3
ln2 = np.log(2)
zeta2 = 1.6449340668482264
zeta3 = 1.2020569031595942

def C2_g1_threshold(x, m2Q2):
    beta = np.sqrt(1. - 4. * m2Q2 * x / (1. - x))
    xi = 1/m2Q2
    return xi * TR * beta / (1. + xi / 4.) / x

def Li2(x):
    return spence(1-x)

def c0(xi):
    y = np.sqrt(1. + 4. / xi)

    L1 = np.log(1. + xi / 2)
    L2 = np.log(2. + xi / 2)
    L3 = np.log(np.sqrt(xi) * (y - 1.) / 2)

    xip2 = 2. + xi
    xip4 = 4. + xi

    Li_2 = Li2(-2. / xip2)
    pi2 = np.pi * np.pi

    c_CA = (50. - pi2 + 12. * L3 / y + 4. * L3 * L3 + L1 * L1 + 6. * L2
                  - 4. * L2 * L2 + 2 * Li_2 + 48. / xip2 - 4. * L2 / xip2
                  + 64. * L2 / xip2 / xip2 - 128. * L2 / (xip2 * xip2 * xip4)
                  - 160. / xip2 / xip4 - 64. * L2 / xip2 / xip4
                  + 128. / (xip2 * xip4 * xip4) - 12. * (4. + zeta2) / xip4
                  - 8. * L3 * L3 / xip4 + 64. / xip4 / xip4)

    c_CF = (-18. - 2. / 3. * pi2 - 24. * L3 / y - 8. * L3 * L3
                  + 2. * L1 * L1 - 6. * L2 + 4. * Li_2 - 48. / xip2
                  + 8. * L2 / xip2 + 360. / xip2 / xip4
                  + 128. * L2 / xip2 / xip4 - 544. / (xip2 * xip4 * xip4)
                  + 48. * L3 * L3 / xip4 - 8. * L1 * L1 / xip4
                  + (44. + 40. * zeta2) / xip4 - 120. * L2 / xip2 / xip2
                  + 256. * L2 / (xip2 * xip2 * xip4) - 16. * Li_2 / xip4
                  - 272. / xip4 / xip4)

    return CA * c_CA + CF * c_CF

def c0_bar(xi):
    return 4. * CA * (2. + np.log(1. + xi / 4.)) - 4. / 3. * TR

def C_g2_threshold(x, m2Q2, m2mu2, kind):
    beta = np.sqrt(1. - 4. * m2Q2 * x / (1. - x))
    logb = np.log(beta)
    log2b = logb**2
    xi = 1/m2Q2
    exp = (16. * CA * log2b + (48. * CA * ln2 - 40. * CA) * logb
           + (2 * CF - CA) * np.pi**2 / beta + 8. * CA * np.log(m2mu2) * logb)
    exp_const = (c0(xi) + 36. * CA * ln2 * ln2 - 60 * CA * ln2
           + np.log(m2mu2) * (8. * CA * ln2 - c0_bar(xi)))
    lo = C2_g1(x, m2Q2) if kind == '2' else CL_g1(x, m2Q2)
    return lo * (exp + exp_const)

def C_g2_threshold_const(x, m2Q2, m2mu2, kind):
    xi = 1/m2Q2
    exp_const = (c0(xi) + 36. * CA * ln2 * ln2 - 60 * CA * ln2
           + np.log(m2mu2) * (8. * CA * ln2 - c0_bar(xi)))
    lo = C2_g1(x, m2Q2) if kind == '2' else CL_g1(x, m2Q2)
    return lo * (exp_const)

def C_g3_threshold(x, m2Q2, m2mu2, nf, kind):
    
    xi = 1. / m2Q2
    beta = np.sqrt(1. - 4. * m2Q2 * x / (1. - x))

    Lm = np.log(m2mu2)
    Lm2 = Lm * Lm
    l = np.log(beta)
    l2 = l * l
    l3 = l2 * l
    l4 = l3 * l

    ln2_2 = ln2 * ln2
    ln2_3 = ln2_2 * ln2

    pi2 = np.pi * np.pi
    pi4 = pi2 * pi2

    c_log4 = 128. * CA * CA

    c_log3 = ((768. * ln2 - 6464. / 9.) * CA * CA + 128. / 9. * CA * nf
                    + 128. * CA * CA * Lm)

    c_log2 = ((1728. * ln2_2 - 3232. * ln2 - 208. / 3. * pi2 + 15520. / 9.) * CA * CA
        + (64. * ln2 - 640. / 9.) * CA * nf + 16. * CA * c0(xi)
        + 32. * CA * (CF - CA / 2) * pi2 / beta
        - ((-512 * ln2 + 1136. / 3.) * CA * CA - 32. / 3. * CA * nf
           + 16 * CA * c0_bar(xi))
              * Lm
        + 32 * CA * CA * Lm2)

    c_log_const = ((1728. * ln2_3 - 4848 * ln2_2 + 15520. / 3. * ln2 - 208 * pi2 * ln2
         + 936 * zeta3 + 608. / 3. * pi2 - 88856. / 27.)
            * CA * CA
        + (96. * ln2_2 - 640. / 3. * ln2 - 16. / 3. * pi2 + 4592. / 27.) * CA
              * nf
        - 32. * CF * (CF - CA / 2) * pi2 + (48. * ln2 - 40.) * CA * c0(xi))

    c_log_fracbeta = ((-92. / 3. + 32. * ln2) * CA + 8. / 3. * nf) * (CF - CA / 2) * pi2

    c_log_Lm = -((-672. * ln2_2 + 976 * ln2 + 104. / 3. * pi2 - 4160. / 9.) * CA * CA
          + (-32. * ln2 + 320. / 9.) * CA * nf
          + (48. * ln2 - 40.) * CA * c0_bar(xi) - 8. * CA * c0(xi)
          - 16. * CA * (CF - CA / 2) * pi2 / beta)

    c_log_Lm2 = ((64. * ln2 - 44. / 3.) * CA * CA + 8. / 3. * CA * nf
                       - 8. * CA * c0_bar(xi))

    c_log = c_log_const + c_log_fracbeta / beta + c_log_Lm * Lm + c_log_Lm2 * Lm2

    c_fracbeta = (((8. * ln2_2 - 68. / 3. * ln2 + 8. / 3. * pi2 - 658. / 9.) * CA
         + (8. / 3. * ln2 - 20. / 9.) * nf + 2 * c0(xi)
         + (26. / 3. * CA + 4. / 3. * nf - 2 * c0_bar(xi)) * Lm)
        * (CF - CA / 2) * pi2)

    c_fracbeta2 = 4. / 3. * (CF - CA / 2) * (CF - CA / 2) * pi4

    exp = (c_log4 * l4 + c_log3 * l3 + c_log2 * l2 + c_log * l
           + c_fracbeta / beta + c_fracbeta2 / beta / beta)
    
    exp_const = (c0(xi) + 36. * CA * ln2 * ln2 - 60. * CA * ln2
                          + Lm * (8. * CA * ln2 - c0_bar(xi)))**2
    
    lo = C2_g1(x, m2Q2) if kind == '2' else CL_g1(x, m2Q2)
    return lo * (exp + exp_const)

def C_g3_threshold_const(x, m2Q2, m2mu2, kind):
    xi = 1/m2Q2
    exp_const = (c0(xi) + 36. * CA * ln2 * ln2 - 60 * CA * ln2
           + np.log(m2mu2) * (8. * CA * ln2 - c0_bar(xi)))**2
    lo = C2_g1(x, m2Q2) if kind == '2' else CL_g1(x, m2Q2)
    return lo * (exp_const)

def test_lo():
    for kind in ['2']:
        ml = ad.ThresholdCoefficientFunction(1, kind, 'g')
        
        for m2Q2 in np.geomspace(1e-2, 1e4, 10):
            xmax = 1. / (1 + 4 * m2Q2)
            for x in np.geomspace(1e-5, xmax, 100, endpoint=False):
                res1 = ml.MuIndependentTerms(x, m2Q2, 1)
                res2 = C2_g1_threshold(x, m2Q2)
                np.testing.assert_allclose(res1, res2, rtol=1e-7)

def test_nlo():
    for kind in ['2', 'L']:
        thr = ad.ThresholdCoefficientFunction(2, kind, 'g')
        for m2Q2 in np.geomspace(1e-2, 1e4, 10):
            for m2mu2 in np.geomspace(1e-2, 1e4, 10):
                xmax = 1. / (1 + 4 * m2Q2)
                for x in np.geomspace(1e-5, xmax, 100, endpoint=False):
                    res1 = thr.fx(x, m2Q2, m2mu2, 1)
                    res2 = C_g2_threshold(x, m2Q2, m2mu2, kind)
                    np.testing.assert_allclose(res1, res2, rtol=1e-7)
                    
                    res1 = thr.BetaIndependentTerms(x, m2Q2, m2mu2)
                    res2 = C_g2_threshold_const(x, m2Q2, m2mu2, kind)
                    np.testing.assert_allclose(res1, res2, rtol=1e-7)
                    
def test_nnlo():
    for kind in ['2']:
        thr = ad.ThresholdCoefficientFunction(3, kind, 'g')
        for m2Q2 in np.geomspace(1e-2, 1e4, 10):
            for m2mu2 in np.geomspace(1e-2, 1e4, 10):
                xmax = 1. / (1 + 4 * m2Q2)
                for x in np.geomspace(1e-5, xmax, 100, endpoint=False):
                    for nf in range(1, 6 + 1):
                        res1 = thr.fx(x, m2Q2, m2mu2, nf)
                        res2 = C_g3_threshold(x, m2Q2, m2mu2, nf, kind)
                        np.testing.assert_allclose(res1, res2, rtol=1e-7)
                        
                        res1 = thr.BetaIndependentTerms(x, m2Q2, m2mu2)
                        res2 = C_g3_threshold_const(x, m2Q2, m2mu2, kind)
                        np.testing.assert_allclose(res1, res2, rtol=1e-7)

