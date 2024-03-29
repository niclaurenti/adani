import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import yaml

import adani

runcard = sys.argv[1]

with open(runcard, "r") as stream:
    try:
        runcard = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

m = runcard["m"]
nf = runcard["nf"]
mufrac = runcard["mufrac"]

hs_version = "exact" if runcard["channel"][1] == "q" else "abmp"
appr = adani.ApproximateCoefficientFunction(
    3,
    runcard["channel"][0],
    runcard["channel"][1],
    True,
    hs_version,
    1e-3,
    1e-3,
    1000,
    False,
    25000,
)


def function_to_exe_in_parallel(pair):
    x, q = pair
    mu = mufrac * q
    m2Q2 = m**2 / q**2
    m2mu2 = m**2 / mu**2

    res = appr.fxBand(x, m2Q2, m2mu2, nf)

    return [res.GetCentral(), res.GetHigher(), res.GetLower()]


def run(n_threads, x_grid, Q_grid):
    grid = []
    for q in Q_grid:
        for x in x_grid:
            grid.append((x, q))
    args = (function_to_exe_in_parallel, grid)
    with Pool(n_threads) as pool:
        result = pool.map(*args)
    return result


if __name__ == "__main__":

    here = Path(os.path.dirname(os.path.realpath(__file__)))

    xfname = "x.txt"
    xlines = [l.strip() for l in open(here / xfname)]
    xdata = np.array([l.split() for l in xlines])
    x_grid = np.array([float(i) for i in (xdata[0, :])])

    Qfname = "Q.txt"
    Qlines = [l.strip() for l in open(here / Qfname)]
    Qdata = np.array([l.split() for l in Qlines])
    Q_grid = np.array([float(i) for i in (Qdata[0, :])])

    if runcard["debug"]:
        x_grid = np.geomspace(1e-4, 1, 10)
        Q_grid = np.geomspace(1.0, 100, 10)

    if runcard["verbose"]:
        print(
            f"Computation of the grid for the coefficient function C{runcard['channel']} for m = {m} GeV, nf = {nf}, and µ/Q = {mufrac}"
        )
        print(f"Size of the grid (x,Q) = ({len(x_grid)},{len(Q_grid)})")
        print(
            "This may take a while (depending on the number of threads you choose). In order to spend this time, I would suggest you this interesting view:"
        )
        print("https://www.youtube.com/watch?v=53pG68KCUMI")

    start = time.perf_counter()
    res_vec = np.array(run(runcard["n_threads"], x_grid, Q_grid))
    if runcard["verbose"]:
        print("total running time: ", time.perf_counter() - start)

    res_mat = res_vec.reshape(len(Q_grid), len(x_grid), 3)

    output_name = f"C_{runcard['channel']}_nf{nf}"

    names = ["central", "higher", "lower"]

    os.system(f"mkdir -p {here / 'results'}")

    for i in range(3):
        if runcard["verbose"]:
            print(
                f"Saving {names[i]} grid in: {here / 'results'}/{output_name}_{names[i]}.dat"
            )
        with open(f"{here / 'results'}/" + output_name + names[i] + ".dat", "w") as f:
            np.savetxt(f, res_mat[:, :, i])
