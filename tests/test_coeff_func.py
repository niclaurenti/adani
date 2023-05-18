import adani

def test_C2_g3_approx():
    m2Q2 = 1/10**2
    nf=4
    for x in [1e-4, 1e-3, 1e-2, 1e-1, 0.5]:
        adani.C2_g3_approximation(x, m2Q2, m2Q2, nf, v=0, method_flag=1)
