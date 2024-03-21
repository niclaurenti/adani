import adani as ad
import oldadani as oldad
import numpy as np

# def test_C2_g3_old_version():
#     x = np.geomspace(1e-5, 1., 10, endpoint=False)
#     # results from adani 0.14.3
#     results = {
#         (1., 4): [145471440.12980354, 33179545.935143087, 6929462.250731408, 1318599.6838400103, 255384.0801166295, 69740.70469675733, 24844.823116321517, 4942.881318050308, -362.8192245755473, 0.0],
#         (5., 4): [743188051.6500003, 179481168.73378846, 40035815.07912805, 7854244.375811367, 1219413.3148602317, 108207.1999990224, -3050.0624306276713, 3100.9036119199877, 4377.376650681199, 971.523210850801],
#         (10., 5): [1334772280.7332137, 328299400.3694369, 75154548.09314843, 15292451.546245273, 2484227.8197733606, 209508.06920325314, -39338.199972231814, -16401.85959599065, 1393.592019854113, 2372.2634137450464],
#         (100., 5): [5454248887.211178, 1375742902.3787572, 326416088.93319005, 69974561.47478054, 12282692.537894182, 1152922.8952283196, -295711.53190138395, -204065.27770961143, -59194.17871140024, -4571.839774559608],
#         (10000, 5): [28176049877.850338, 7210310885.839126, 1746703497.5839772, 386475233.25839806, 71595674.99450848, 7793459.921879084, -1534915.8803332832, -1358662.4950624409, -504689.56072508526, -104087.00690700817],
#     }
#     approx = ad.ApproximateCoefficientFunction(3, '2', 'g', True, False, True, 1e-3, 1e-3, 1000, 1, 25000)

#     for (key, value) in results.items():
#         (xi, nf) = key
#         m2Q2 = 1/xi
#         newres = [approx.MuIndependentTerms(x_, m2Q2, nf) for x_ in x]
#         np.testing.assert_allclose(newres, value, 1e-7)

# def test_C2_ps3_old_version():
#     x = np.geomspace(1e-5, 1., 10, endpoint=False)
#     # results from adani 0.14.3
#     results = {
#         (1., 4): [81839166.40060224, 19983631.262041003, 4621934.937475589, 1012838.894594633, 219236.1322225202, 52411.4218696204, 13933.827317278641, 2179.65043103194, -404.9223847223813, 0.],
#         (5., 4): [370433204.98138374, 92312106.96733469, 21677300.29729032, 4675689.2469764855, 896142.4651095183, 151217.03489076335, 27461.923963425932, 8469.665117746337, 2759.8739119587945, 225.60106649013127],
#         (10., 5): [650738512.3677523, 164007577.59588692, 39076424.97397647, 8566796.358678792, 1652696.3401824397, 261530.10191432017, 32908.712435241105, 6335.607095837572, 2964.9567632109183, 722.7467948775832],
#         (100., 5): [2575301519.877781, 659311658.3163434, 160259607.07747307, 35936597.2288462, 7010522.604863436, 1014396.8256257733, 30625.973310558875, -41151.648587472344, -14184.364524809693, -1530.8517412774409],
#         (10000, 5): [13106578380.728123, 3389774971.6744285, 835258542.4674901, 190660545.97906336, 37951895.90912031, 5490598.134831231, -291.29541339835316, -377431.45330583374, -158766.35333458823, -34205.82581016607],
#     }
#     approx = ad.ApproximateCoefficientFunction(3, '2', 'q', True, True)

#     for (key, value) in results.items():
#         (xi, nf) = key
#         m2Q2 = 1/xi
#         newres = [approx.MuIndependentTerms(x_, m2Q2, nf) for x_ in x]
#         np.testing.assert_allclose(newres, value, 1e-7)

# def test_CL_g3_old_version():
#     x = np.geomspace(1e-5, 1., 10, endpoint=False)
#     # results from adani 0.14.3
#     results = {
#         (1., 4): [34898300.0009964, 10602834.023263108, 3262350.460220233, 1034087.369598095, 346930.64915267046, 125205.18438577586, 42406.44435516963, 5850.420132240786, 143.01488937933098, 0.0],
#         (5., 4): [102671099.55566017, 26247870.32068742, 6400147.298576518, 1470243.4047489867, 321810.9937162599, 76803.03951467396, 27946.654766119653, 14624.147780325871, 4458.363807712411, 123.41372384391079],
#         (10., 5): [223551430.11224937, 56850702.63623221, 13672875.726310845, 3015035.831987956, 576309.5624274018, 87272.80988683799, 12342.26016424239, 6778.871070476519, 5284.620133376426, 653.6133370305272],
#         (100., 5): [1212706459.2265527, 313445500.42034954, 77076002.982617, 17482900.176525474, 3408500.8550913734, 449111.71224176907, -29035.195438674924, -46746.07742777833, -18226.57667534148, -531.6644084560719],
#         (10000, 5): [5229207028.838566, 1364581011.5641327, 340307696.5436845, 78922913.9855198, 16002043.18906737, 2317376.5350905857, -69847.11917414087, -233526.22458851323, -118252.60468366444, -30733.81096799032],
#     }
#     approx = ad.ApproximateCoefficientFunction(3, 'L', 'g', True, False, True, 1e-3, 1e-3, 1000, 1, 25000)

#     for (key, value) in results.items():
#         (xi, nf) = key
#         m2Q2 = 1/xi
#         newres = [approx.MuIndependentTerms(x_, m2Q2, nf) for x_ in x]
#         np.testing.assert_allclose(newres, value, 1e-7)
oldad
# def test_CL_ps3_old_version():
#     x = np.geomspace(1e-5, 1., 10, endpoint=False)
#     # results from adani 0.14.3
#     results = {
#         (1., 4): [18264900.065671325, 5591545.049142779, 1733424.9872403524, 552680.6241826949, 185541.3541547223, 66311.29117515919, 21919.537089617108, 2880.1691002468774, 51.25789565581389, 0.],
#         (5., 4): [53316533.64214351, 14107092.196654607, 3624404.0989812785, 905296.4869830979, 225878.65406261428, 62063.4502627242, 21955.383406969613, 9442.354828807334, 2309.5435010836586, 48.52355451171832],
#         (10., 5): [110769910.10295145, 28892841.701933067, 7235023.398763551, 1713975.0243687262, 379150.3434254426, 80320.6560868768, 19844.226358982785, 7878.731021058479, 3444.4867439240484, 287.1410019321237],
#         (100., 5): [568176213.7950808, 148566058.93647945, 37201210.981392816, 8714077.753391517, 1821674.7528531721, 301718.89800044947, 22312.082572018026, -7945.338782486486, -3230.1188120387897, 108.24701694208481],
#         (10000, 5): [2416508349.1613946, 635750606.6713928, 160538969.44249067, 38041116.7236875, 8068124.978092233, 1345060.6634551995, 77162.87091539113, -63049.750050451476, -33836.14451590534, -8501.133810409003],
#     }
#     approx = ad.ApproximateCoefficientFunction(3, 'L', 'q', True, True)

#     for (key, value) in results.items():
#         (xi, nf) = key
#         m2Q2 = 1/xi
#         newres = [approx.MuIndependentTerms(x_, m2Q2, nf) for x_ in x]
#         np.testing.assert_allclose(newres, value, 1e-7)

# def test_mudependent_terms():
#     for order in [2, 3]:
#         for channel in ['g', 'q']:
#             for kind in ['2', 'L']:
#                 for dim in [100, 1000]:
#                     for MCcalls in [10000, 20000]:
#                         for mf in [0, 1]:
#                             for abserr in [1e-2, 1e-3]:
#                                 relerr = abserr
#                                 massive = ad.ExactCoefficientFunction(order, kind, channel, abserr, relerr, dim, mf, MCcalls)
#                                 app = ad.ApproximateCoefficientFunction(order, kind, channel,True, False, True, abserr, relerr, dim, mf, MCcalls)
#                                 x = np.geomspace(1e-5, 1., 10, endpoint=True)
#                                 for xi in np.geomspace(1e-2, 1e4, 4, endpoint=True):
#                                     for nf in range(1, 6 + 1):
#                                         res1 = [massive.MuDependentTerms(x_, 1/xi, 1/xi, nf) for x_ in x]
#                                         res2 = [app.MuDependentTerms(x_, 1/xi, 1/xi, nf) for x_ in x]
#                                         np.testing.assert_allclose(res1, res2, rtol=1e-7)
#                                 del massive
#                                 del app

# def test_C30_oldversion():
#     for channel in ['g', 'q']:
#         for kind in ['2', 'L']:
#             exact = True if channel == 'q' else False
#             hs_appr = True if channel == 'g' else False
#             app = ad.ApproximateCoefficientFunction(3, kind, channel,True, exact, hs_appr)
#             for xi in np.geomspace(1e-2, 1e2, 10):
#                 m2Q2 = 1/xi
#                 for x in np.geomspace(1e-5, 1, 10):
#                     for nf in range(1, 6 + 1):
#                         for v in [0, 1, -1]:
#                             if kind == '2':
#                                 if channel == 'g':
#                                     res_old = oldad.C2_g3_approximation(x, m2Q2, 1., nf, v)
#                                 if channel == 'q':
#                                     res_old = oldad.C2_ps3_approximation(x, m2Q2, 1., nf, v)
#                             if kind == 'L':
#                                 if channel == 'g':
#                                     res_old = oldad.CL_g3_approximation(x, m2Q2, 1., nf, v)
#                                 if channel == 'q':
#                                     res_old = oldad.CL_ps3_approximation(x, m2Q2, 1., nf, v)
#                             tmp = app.MuIndependentTermsBand(x, m2Q2, nf)
#                             if v == 0:
#                                 res_new = tmp.GetCentral()
#                             if v == -1:
#                                 res_new = tmp.GetLower()
#                             if v == 1:
#                                 res_new = tmp.GetHigher()
#                             np.testing.assert_allclose(res_old, res_new, rtol=1e-7)

def test_mudep_oldversion():
    for channel in ['g', 'q']:
        for kind in ['2', 'L']:
            exact = True if channel == 'q' else False
            hs_appr = True if channel == 'g' else False
            for mf in [0]:
                app = ad.ApproximateCoefficientFunction(3, kind, channel,True, exact, hs_appr, 1e-3, 1e-3, 1000, mf, 25000)
                for xi in np.geomspace(1e-2, 1e2, 10):
                    m2Q2 = 1/xi
                    xmax = 1/(1 + 4*m2Q2)
                    Lmu = -np.log(m2Q2)
                    for x in np.geomspace(1e-5, 1, 10):
                        for nf in range(1, 6 + 1):
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
