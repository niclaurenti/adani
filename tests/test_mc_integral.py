import adani
import numpy as np

mb=4.92
nf=4

def test_mc_integrals():
    for q in np.geomspace(2 * mb, 100., 10):
        m2Q2 = (mb / q)**2
        xmax = 1. / (1. + 4 * m2Q2)
        for x in np.geomspace(1e-4, xmax, 10, endpoint=False):
            mc = adani.C2_g1_x_Pgg0_x_Pgg0_MC(x, m2Q2, nf)
            no_mc = adani.C2_g1_x_Pgg0_x_Pgg0(x, m2Q2, nf)
            # np.testing.assert_allclose(mc, no_mc, rtol=1e-2)
