import adani as ad
import numpy as np

def test_nll():
    for nll in [True, False]:
        he = ad.HighEnergyCoefficientFunction(3, "2", "g", NLL=nll)
        hehs = ad.HighEnergyHighScaleCoefficientFunction(3, "2", "g", NLL=nll)
        assert he.GetNLL() == nll
        assert hehs.GetNLL() == nll
