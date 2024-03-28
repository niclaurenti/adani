import adani as ad
import numpy as np

nf = 4
m = 4.92

F2g = ad.ApproximateCoefficientFunction(3, '2', 'g', True, "abmp", 1e-3, 1e-3, 1000, False, 25000)
FLg = ad.ApproximateCoefficientFunction(3, 'L', 'g', True, "abmp", 1e-3, 1e-3, 1000, False, 25000)
F2q = ad.ApproximateCoefficientFunction(3, '2', 'q', True, "exact", 1e-3, 1e-3, 1000, False, 25000)
FLq = ad.ApproximateCoefficientFunction(3, 'L', 'q', True, "exact", 1e-3, 1e-3, 1000, False, 25000)

for x in np.geomspace(1e-4, 1, 5):
    for Q in np.geomspace(10, 100, 5):
        m2Q2 = m**2 / Q**2
        m2mu2 = m2Q2
        print("C2_g^(3)(x=", x, ", Q=", Q, ") = ", F2g.fx(x, m2Q2, m2mu2, nf))
        print("C2_ps^(3)(x=", x, ", Q=", Q, ") = ", F2g.fx(x, m2Q2, m2mu2, nf))
        print("CL_g^(3)(x=", x, ", Q=", Q, ") = ", F2g.fx(x, m2Q2, m2mu2, nf))
        print("CL_ps^(3)(x=", x, ", Q=", Q, ") = ", F2g.fx(x, m2Q2, m2mu2, nf))
