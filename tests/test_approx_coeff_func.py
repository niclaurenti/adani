from adani import Analytical, DoubleNumerical, MonteCarlo
from adani import Exact, GM, ABMP, KLMV
import adani as ad
import numpy as np

def test_coeff_func():
    for order in range(1, 3 + 1):
        for kind in ["2", "L"]:
            for channel in ["q", "g"]:
                if order == 1 and channel == "q":
                    continue
                app = ad.ApproximateCoefficientFunction(order, kind, channel)
                assert app.GetOrder() == order
                assert app.GetKind() == kind
                assert app.GetChannel() == channel

def test_mudependent_terms():
    for order in [2, 3]:
        for channel in ['g', 'q']:
            for kind in ['2', 'L']:
                for dim in [1000]:
                    for MCcalls in [20000]:
                        for mf in [Analytical, DoubleNumerical, MonteCarlo]:
                            for abserr in [1e-3]:
                                hscale = KLMV if order == 3 else Exact
                                relerr = abserr
                                massive = ad.ExactCoefficientFunction(order, kind, channel, abserr, relerr, dim)
                                app = ad.ApproximateCoefficientFunction(order, kind, channel, True, hscale, abserr, relerr, dim)
                                if mf in [DoubleNumerical, MonteCarlo]:
                                    massive.SetDoubleIntegralMethod(mf, abserr, relerr, dim, MCcalls)
                                    app.SetDoubleIntegralMethod(mf, abserr, relerr, dim, MCcalls)
                                x = np.geomspace(1e-5, 1., 5, endpoint=True)
                                for xi in np.geomspace(1e-2, 1e4, 4, endpoint=True):
                                    for nf in [4, 5]:
                                        res1 = [massive.MuDependentTerms(x_, 1/xi, 1/xi, nf) for x_ in x]
                                        res2 = [app.MuDependentTerms(x_, 1/xi, 1/xi, nf) for x_ in x]
                                        np.testing.assert_allclose(res1, res2, rtol=1e-7)
