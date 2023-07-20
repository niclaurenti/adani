import adani
from multiprocessing import Pool
import yaml
import numpy as np
import time
import sys

runcard = sys.argv[1]

with open(runcard, 'r') as stream:
    try:
        runcard=yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)

m = runcard["m"]
nf = runcard["nf"]
mufrac = runcard["mufrac"]
v = runcard.get("v", 0)

def function_to_exe_in_parallel(pair):
    x, q = pair
    mu = mufrac * q
    m2Q2 = m**2 / q**2
    m2mu2 = m**2 / mu**2
    if runcard["channel"] == "2g":
        return adani.C2_g3_approximation(x, m2Q2, m2mu2, nf, v, method_flag=1)
    elif runcard["channel"] == "2q":
        return adani.C2_ps3_approximation(x, m2Q2, m2mu2, nf, v)
    elif runcard["channel"] == "Lg":
        return adani.CL_g3_approximation(x, m2Q2, m2mu2, nf, v, method_flag=1)
    elif runcard["channel"] == "Lq":
        return adani.CL_ps3_approximation(x, m2Q2, m2mu2, nf, v)
    else:
        raise ValueError("Set channel to one of these: 2g 2q Lg Lq")


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
    xfname = "x.txt"
    xlines = [l.strip() for l in open(xfname)]
    xdata = np.array([ l.split() for l in xlines ])
    x_grid = np.array([float(i) for i in (xdata[0,:])])
    Qfname = "Q.txt"
    Qlines = [l.strip() for l in open(Qfname)]
    Qdata = np.array([ l.split() for l in Qlines ])
    Q_grid = np.array([float(i) for i in (Qdata[0,:])])
    if runcard["debug"]:
        x_grid = np.geomspace(1e-4, 1, 10)
        Q_grid = np.geomspace(1., 100, 10)

    if runcard["verbose"]:
        print(f"Computation of the grid for the coefficient function C{runcard['channel']} for m = {m} GeV, nf = {nf}, and µ/Q = {mufrac}")
        print(f"Size of the grid (x,Q) = ({len(x_grid)},{len(Q_grid)})")
        print("This may take a while (depending on the number of threads you choose). In order to spend this time, I would suggest you this interesting view:")
        print("https://www.youtube.com/watch?v=53pG68KCUMI")

    start = time.perf_counter()
    res_vec=np.array(run(runcard["n_threads"], x_grid, Q_grid))
    if runcard["verbose"]:
        print("total running time: ", time.perf_counter() - start)

    res_mat = res_vec.reshape(len(Q_grid), len(x_grid))
    if runcard["verbose"]:
        print("Saving grid in ", runcard["output_file"])
    with open(runcard["output_file"], 'w') as f:
        np.savetxt(f, res_mat)
