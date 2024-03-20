import adani as ad
import numpy as np

def test_C2_g3_old_version():
    x = np.geomspace(1e-5, 1., 10, endpoint=False)
    results = {
        (1., 4): [145471440.12980354, 33179545.935143087, 6929462.250731408, 1318599.6838400103, 255384.0801166295, 69740.70469675733, 24844.823116321517, 4942.881318050308, -362.8192245755473, 0.0],
        (5., 4): [743188051.6500003, 179481168.73378846, 40035815.07912805, 7854244.375811367, 1219413.3148602317, 108207.1999990224, -3050.0624306276713, 3100.9036119199877, 4377.376650681199, 971.523210850801],
        (10., 5): [1334772280.7332137, 328299400.3694369, 75154548.09314843, 15292451.546245273, 2484227.8197733606, 209508.06920325314, -39338.199972231814, -16401.85959599065, 1393.592019854113, 2372.2634137450464],
        (100., 5): [5454248887.211178, 1375742902.3787572, 326416088.93319005, 69974561.47478054, 12282692.537894182, 1152922.8952283196, -295711.53190138395, -204065.27770961143, -59194.17871140024, -4571.839774559608],
        (10000, 5): [28176049877.850338, 7210310885.839126, 1746703497.5839772, 386475233.25839806, 71595674.99450848, 7793459.921879084, -1534915.8803332832, -1358662.4950624409, -504689.56072508526, -104087.00690700817],
    }
    approx = ad.ApproximateCoefficientFunction(3, '2', 'g', True, False, True, 1e-3, 1e-3, 1000, 1, 25000)

    for (key, value) in results.items():
        (xi, nf) = key
        m2Q2 = 1/xi
        newres = [approx.MuIndependentTerms(x_, m2Q2, nf) for x_ in x]
        np.testing.assert_allclose(newres, value, 1e-5)

def test_mudependent_terms():
    for order in [3]:
        for channel in ['g', 'q']:
            for kind in ['2', 'L']:
                for dim in [100, 1000]:
                    for MCcalls in [10000, 20000]:
                        for mf in [0, 1]:
                            for abserr in [1e-2, 1e-3]:
                                relerr = abserr
                                massive = ad.ExactCoefficientFunction(order, kind, channel, abserr, relerr, dim, mf, MCcalls)
                                app = ad.ApproximateCoefficientFunction(order, kind, channel,True, False, True, abserr, relerr, dim, mf, MCcalls)
                                x = np.geomspace(1e-5, 1., 10, endpoint=True)
                                for xi in np.geomspace(1e-2, 1e4, 4, endpoint=True):
                                    for nf in range(1, 6 + 1):
                                        res1 = [massive.MuDependentTerms(x_, 1/xi, 1/xi, nf) for x_ in x]
                                        res2 = [app.MuDependentTerms(x_, 1/xi, 1/xi, nf) for x_ in x]
                                        np.testing.assert_allclose(res1, res2, rtol=1e-7)
                                del massive
                                del app
