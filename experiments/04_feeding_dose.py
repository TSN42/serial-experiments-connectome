"""摂食回路（T1）の用量反応。局在域の重み倍率で、糖の刺激頻度を振り、苦味の同時刺激で抑制されるかを見る。

    python experiments/04_feeding_dose.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from secon.data import load  # noqa: E402
from secon.lif import LIFParams, Network, simulate  # noqa: E402

OUT = ROOT / "results" / "lif"


def main():
    c = load()
    sugar = c.select(type=["LB3b", "LB3c"])
    bitter = c.select(type_regex=r"^LB1[a-d]$")
    water = c.select(type="LB3a")
    mn9 = c.select(type="MN9", side="L")
    rows = []
    for w in (0.40, 0.44, 0.46):
        p = LIFParams(w_scale=w)
        net = Network.from_connectome(c, p)
        conds = [(f"sugar {f}Hz", [(sugar, f)]) for f in (100, 125, 150, 175, 200, 250)]
        conds += [(f"sugar 200Hz + bitter {f}Hz", [(sugar, 200), (bitter, f)]) for f in (25, 50, 100)]
        conds += [("bitter 200Hz", [(bitter, 200)]), ("water 200Hz", [(water, 200)])]
        for name, stim in conds:
            r = simulate(net, stim, p, n_trials=3, seed=11).rates
            rows.append({"w_scale": w, "condition": name, "MN9_L_hz": float(r[mn9].mean()), "active_neurons": int((r > 1).sum())})
            print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(OUT / "feeding_dose.csv", index=False)


if __name__ == "__main__":
    main()
