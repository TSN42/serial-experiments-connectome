"""MaleCNS v1.0 の構造解析。結果は results/structure/ に、図は docs/figures/ に出す。

    python experiments/01_structure.py
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from secon import structure as S  # noqa: E402
from secon.data import load  # noqa: E402
from secon.groups import effector_groups, sensory_groups  # noqa: E402

OUT = ROOT / "results" / "structure"
OUT.mkdir(parents=True, exist_ok=True)


def main():
    t0 = time.time()
    c = load()
    nr = c.neurons
    summary = {"neurons": int(c.n), "edges": int(len(c.pre)), "synapses": int(c.count.sum())}

    # --- 1. 組成と伝達物質 -------------------------------------------------
    comp = pd.crosstab(nr.superclass, nr.nt)
    comp.to_csv(OUT / "superclass_by_nt.csv")
    summary["nt_counts"] = nr.nt.value_counts().to_dict()
    summary["nt_source"] = nr.nt_source.value_counts().to_dict()
    w = c.signed_weights("fast_only")
    tot_in = np.bincount(c.post, weights=c.count, minlength=c.n)
    inh_in = np.bincount(c.post, weights=np.where(w < 0, c.count, 0), minlength=c.n)
    zero_in = np.bincount(c.post, weights=np.where(w == 0, c.count, 0), minlength=c.n)
    ei = pd.DataFrame({"superclass": nr.superclass, "inh_frac": inh_in / np.maximum(tot_in, 1),
                       "unsigned_frac": zero_in / np.maximum(tot_in, 1), "in_syn": tot_in})
    ei.groupby("superclass")[["inh_frac", "unsigned_frac", "in_syn"]].median().to_csv(OUT / "input_ei_by_superclass.csv")
    summary["median_inh_input_fraction"] = float(np.median(ei.inh_frac[tot_in > 0]))
    summary["synapses_per_edge"] = {q: float(np.percentile(c.count, q)) for q in (50, 90, 99)}
    summary["edges_ge5_fraction"] = float((c.count >= 5).mean())
    summary["synapses_on_edges_ge5_fraction"] = float(c.count[c.count >= 5].sum() / c.count.sum())

    # --- 2. 左右差（再構成の欠けの検出） ---------------------------------------
    lr = S.lr_asymmetry(c)
    lr.to_csv(OUT / "lr_asymmetry.csv")
    flag = lr.log2_inp_ratio.abs() > 1
    summary["lr_types"] = int(len(lr))
    summary["lr_flagged_2x"] = int(flag.sum())
    summary["lr_flagged_by_superclass"] = (lr[flag].groupby("superclass").size()
                                           / lr.groupby("superclass").size()).fillna(0).round(3).to_dict()
    print(f"[{time.time()-t0:.0f}s] composition + LR done", flush=True)

    # --- 3. 感覚 → 効果器: 確率的走査の層 ------------------------------------
    sens, eff = sensory_groups(c), effector_groups(c)
    eff_all = np.concatenate(list(eff.values()))
    summary["sensory_group_sizes"] = {k: int(len(v)) for k, v in sens.items()}
    summary["effector_group_sizes"] = {k: int(len(v)) for k, v in eff.items()}
    rows = []
    layers = {}
    for name, seeds in sens.items():
        lay = S.traversal_layers(c, seeds, n_runs=64)
        layers[name] = lay
        for ename, idx in eff.items():
            l = lay[idx]
            rows.append({"sensory": name, "effector": ename, "median_layer": float(np.nanmedian(l)) if np.isfinite(l).any() else np.nan,
                         "reached_frac": float(np.isfinite(l).mean())})
        print(f"[{time.time()-t0:.0f}s] traversal {name}: reached {np.isfinite(lay).mean():.3f} of CNS", flush=True)
    pd.DataFrame(rows).to_csv(OUT / "traversal_sensory_to_effector.csv", index=False)
    pd.DataFrame(layers, index=nr.bodyId).to_parquet(OUT / "traversal_layers.parquet")

    # --- 4. 線形影響度（BANC 方式、非符号・符号付き、λ = 0.99 と 0.5） ---------------
    # λ=0.99 だと定常解は先頭固有モード（利得 1/(1-λ) = 100 倍）に支配される。MaleCNS ではその先頭モードが
    # 触角葉の局所ニューロン塊に局在しているので、λ=0.5 と並べて比べる。値は「種 1 細胞あたり」に揃える。
    rows = []
    for signed in (False, True):
        key = "signed" if signed else "unsigned"
        W = S.signed_matrix(c, "fast_only", min_syn=5)
        if not signed:
            W = abs(W)
        lead, vec = S.leading_mode(W)
        mass = vec / vec.sum()
        top = np.argsort(-mass)[:20]
        summary[f"leading_mode_{key}_min5"] = {
            "eigenvalue": [lead.real, lead.imag],
            "participation_ratio": float(1 / (mass ** 2).sum()),
            "top_types": nr.type.iloc[top].value_counts().to_dict(),
            "mass_by_superclass": pd.Series(mass).groupby(nr.superclass.to_numpy()).sum().round(4).sort_values(ascending=False).head(5).to_dict(),
        }
        for lam in (0.99, 0.5):
            profiles = {}
            for name, seeds in sens.items():
                r = S.linear_influence(W, seeds, scale=abs(lead), lam=lam) / len(seeds)
                fname = name.replace(" ", "_").replace("/", "-").replace(":", "")
                np.save(OUT / f"influence_{key}_lam{lam}_{fname}.npy", r.astype(np.float32))
                for ename, idx in eff.items():
                    rows.append({"signed": signed, "lam": lam, "sensory": name, "effector": ename,
                                 "mean": float(r[idx].mean()), "mean_abs": float(np.abs(r[idx]).mean())})
                profiles[name] = np.log10(np.abs(r[eff_all]) + 1e-30)
            # 感覚群ごとの効果器プロファイルがどれだけ似てしまうか（先頭モード支配の指標）
            cc = np.corrcoef(np.array(list(profiles.values())))
            summary[f"influence_profile_mean_corr_{key}_lam{lam}"] = float(cc[np.triu_indices_from(cc, 1)].mean())
            print(f"[{time.time()-t0:.0f}s] influence {key} lam={lam} done", flush=True)
    pd.DataFrame(rows).to_csv(OUT / "influence_sensory_to_effector.csv", index=False)

    # --- 5. 効果器はどこから入力を受けているか -------------------------------
    src = {name: S.input_sources(c, idx) for name, idx in eff.items()}
    pd.DataFrame(src).fillna(0).to_csv(OUT / "effector_input_sources.csv")

    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False, default=float))
    print(f"done in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
