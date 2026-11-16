"""Compare the N = 140 and N = 180 rows of results_ext/conv.json; writes results_ext/convergence_summary.json."""
import json
from pathlib import Path
RES = Path(__file__).resolve().parent.parent / "results_ext"
rows = json.loads((RES / "conv.json").read_text())
key = lambda r: (r["experiment"], r["label"], r["decoder"], r.get("channel.sigma_phi"), r.get("channel.db"), r.get("channel.gain"))
by = {}
for r in rows:
    by.setdefault(key(r), {})[r["N"]] = 1 - r["average_fidelity"]
out = []
for k, v in by.items():
    a, b = v[140], v[180]
    out.append({"experiment": k[0], "label": k[1], "decoder": k[2], "sigma_phi": k[3], "db": k[4], "gain": k[5],
                "infid_N140": a, "infid_N180": b, "abs_diff": abs(a - b), "rel_diff": abs(a - b) / b})
out.sort(key=lambda o: -o["rel_diff"])
(RES / "convergence_summary.json").write_text(json.dumps(out, indent=1))
print("max abs diff:", max(o["abs_diff"] for o in out), " max rel diff:", out[0]["rel_diff"])
for o in out[:8]:
    print(o)
