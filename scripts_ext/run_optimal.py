"""Optimal recovery (SDP, optimal_recovery.solve_sdp) for the cat and Fock codes on the E1 and E2 grids.

    CVXPY_PATH=<dir with cvxpy + clarabel> conda run -n qcsim python scripts_ext/run_optimal.py
GKP is not included: its output support (d = 70-140 Fock directions) makes the SDP too slow for Clarabel here.
For GKP we use the Barnum-Knill bracket instead (make_ext_figures.py): F_e(Petz) >= F_e(opt)^2.
Writes results_ext/optimal.csv/.json; each row stores code, channel and the certified dual upper bound.
"""
import json
import time
import numpy as np
from cvqec.circuits.teleportation import TeleportationParams, effective_channel
from cvqec.experiments.runner import save_rows
from cvqec.noise import Dephasing
from optimal_recovery import choi_X, solve_sdp
from run_ext import N, OUT, code_of, lineup

prev = OUT / "optimal.json"
rows = json.loads(prev.read_text()) if prev.exists() else []      # resume: keep finished points
done = {(r["experiment"], r["label"], r.get("sigma_phi", r.get("db"))) for r in rows}
codes = [c for c in lineup() if c[1] in ("cat4", "fock")]
jobs = [("E1", "sigma_phi", s, Dephasing(s)) for s in np.round(np.arange(0.02, 0.401, 0.02), 4).tolist()]
# E2 on a coarser grid: the support (and the SDP cost) grows fast at low squeezing (cat: d = 59 at 2 dB).
jobs += [("E2", "db", db, effective_channel(TeleportationParams.from_db(db))) for db in range(2, 17, 2)]
E2_MIN_DB = {"cat4": 8, "fock": 2}
for label, key, params in codes:
    code = code_of(key, params)
    for exp, xk, x, ch in jobs:
        if (exp, label, x) in done or (exp == "E2" and x < E2_MIN_DB[key]):
            continue
        t = time.time()
        tol = 1e-13 if exp == "E1" else 1e-11      # the neglected weight is added to the upper bound
        r = solve_sdp(choi_X(code, ch, N)[0], N, support_tol=tol)
        r["support_tol"] = tol
        row = {"experiment": exp, "label": label, "code": key, **{f"code.{k}": v for k, v in params.items()},
               "N": N, xk: x, "channel.repr": ch.name, **r,
               "F_avg_opt_upper": (2 * r["F_e_opt_upper"] + 1) / 3, "seconds": round(time.time() - t, 2)}
        rows.append(row)
        print(f"{exp} {label} {xk}={x}: 1-F_avg opt {1 - r['F_avg_opt']:.4e}  gap {r['F_e_opt_upper'] - r['F_e_opt']:.1e}"
              f"  d={r['support_dim']}  {row['seconds']} s", flush=True)
        save_rows(rows, OUT / "optimal.csv")
print("== saved optimal", len(rows))
