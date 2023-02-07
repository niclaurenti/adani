import adani
import numpy as np

nf = 4
m = 4.92

for x in np.geomspace(1e-4, 1, 5):
    for Q in np.geomspace(10, 100, 5):
        mQ = m**2 / Q**2
        mMu = mQ
        print("C2_g^(3)(x=", x, ", Q=", Q, ") = ", adani.C2_g3_approximation(x, mQ, mMu, nf))
        print("C2_ps^(3)(x=", x, ", Q=", Q, ") = ", adani.C2_ps3_approximation(x, mQ, mMu, nf))
        print("CL_g^(3)(x=", x, ", Q=", Q, ") = ", adani.CL_g3_approximation(x, mQ, mMu, nf))
        print("CL_ps^(3)(x=", x, ", Q=", Q, ") = ", adani.CL_ps3_approximation(x, mQ, mMu, nf))
