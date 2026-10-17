"""MaleCNS 全中枢 LIF の検証スイート。重みスケール w_scale を振り、各回路テストの指標を並べる。

    python experiments/02_lif.py                 # 重み = シナプス数 × 0.275 mV × w_scale（Shiu 型）を振る
    python experiments/02_lif.py --alpha 0.5     # さらに受け手の大きさで補正: × (346 / 総入力)^alpha
    python experiments/02_lif.py --adapt-b 3     # スパイク頻度順応（発火ごとに 3 mV、τ = 100 ms）
    python experiments/02_lif.py --quick         # 動作確認用

テストは「配線のどこを刺激し、どこが応答すべきか」を文献から決めたもの。
  T1 摂食   : 糖 GRN (LB3b/c) → MN9。苦味 (LB1a-d) は駆動せず、糖+苦味で抑制（Shiu 2024; Tastekin 2026）
  T2 グルーミング: JO-C/E → aDN1 (DNg62), aDN2 (DNge078)（Hampel 2015; Shiu 2024）
  T3 逃避   : LPLC2 → 巨大線維 DNp01 → TTMn。GF→TTMn は電気シナプスなので、無い/有りで比較（Allen 2006 ほか）
  T4 歩行   : DNg100 (BDN2) → 脚 MN のリズム（Bidaye 2020; Pugliese 2025）
  T5 求愛   : P1 相当（オス特異的 pC1, fru+dsx 共発現）→ pIP10（Ding 2019 ほか）
  R  対照   : 無作為に選んだ中枢 ACh ニューロン 30 個 → 活動がどこまで広がるか（暴走の指標）
"""
import argparse
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from secon.data import load  # noqa: E402
from secon.lif import LIFParams, Network, simulate  # noqa: E402

OUT = ROOT / "results" / "lif"
SCALES = [0.30, 0.35, 0.40, 0.42, 0.44, 0.46, 0.50, 0.55, 0.60, 0.70, 1.00]
GAP_KICK_MV = 9.0   # 1 回の GF スパイクで静止中の TTMn/PSI を確実に閾値超えさせる（実物は 1:1 伝達）


def neuron_sets(c):
    nr = c.neurons
    s = {
        "sugar": c.select(type=["LB3b", "LB3c"]),
        "bitter": c.select(type_regex=r"^LB1[a-d]$"),
        "mn9_L": c.select(type="MN9", side="L"),
        "jo_ce": c.select(type_regex=r"^JO-[CE]"),
        "jo_ab": c.select(type_regex=r"^JO-[AB]"),
        "groom_dn": c.select(type=["DNg62", "DNge078"]),
        "lplc2": c.select(type="LPLC2"),
        "gf": c.select(type="DNp01"),
        "ttmn": c.select(type="TTMn"),
        "psi": c.select(type="PSI"),
        "dng100": c.select(type="DNg100"),
        "leg_mn": np.flatnonzero(nr.superclass.eq("vnc_motor") & nr.subclass.isin(["fl", "ml", "hl"])),
        "leg_flexor": np.flatnonzero(nr.superclass.eq("vnc_motor") & nr.type.fillna("").str.contains(r"(?:Ti|Tr) flexor MN$")),
        "leg_extensor": np.flatnonzero(nr.superclass.eq("vnc_motor") & nr.type.fillna("").str.contains(r"(?:Ti|Tr) extensor MN$")),
        "dn": c.select(superclass="descending_neuron"),
        "pip10": c.select(type="pIP10"),
        "p1_like": np.flatnonzero(nr.type.fillna("").str.startswith("pC1") & nr.dimorphism.eq("male-specific")
                                  & nr.fruDsx.fillna("").str.startswith("coexpress")),
        "wing_mn": np.flatnonzero(nr.superclass.eq("vnc_motor") & nr.subclass.eq("wm")),
    }
    rng = np.random.default_rng(42)
    pool = np.flatnonzero(nr.superclass.eq("cb_intrinsic") & nr.nt.eq("acetylcholine"))
    s["random"] = [rng.choice(pool, 30, replace=False) for _ in range(3)]
    return s


def gap_pairs(c, sets):
    """GF → TTMn / PSI の電気シナプス（左右 2 対 2）。組み合わせは化学シナプス数の合計が最大になる対応で決める。"""
    pairs = []
    g0, g1 = sets["gf"]
    for tgt in ("ttmn", "psi"):
        t0, t1 = sets[tgt]
        syn = lambda a, b: c.count[(c.pre == a) & (c.post == b)].sum()
        if syn(g0, t0) + syn(g1, t1) >= syn(g0, t1) + syn(g1, t0):
            pairs += [(int(g0), int(t0), GAP_KICK_MV), (int(g1), int(t1), GAP_KICK_MV)]
        else:
            pairs += [(int(g0), int(t1), GAP_KICK_MV), (int(g1), int(t0), GAP_KICK_MV)]
    return pairs


def follow_metrics(spikes, pairs):
    """GF スパイクの直後 5 ms 以内に相手 (TTMn) が発火した割合と、その遅れの中央値。"""
    if len(spikes) == 0:
        return np.nan, np.nan
    t, i = spikes[:, 0], spikes[:, 1].astype(int)
    lags, followed, total = [], 0, 0
    for a, b, _ in pairs:
        ta, tb = t[i == a], t[i == b]
        for x in ta:
            total += 1
            later = tb[(tb >= x) & (tb < x + 5.0)]
            if len(later):
                followed += 1
                lags.append(later[0] - x)
    return (followed / total if total else np.nan), (float(np.median(lags)) if lags else np.nan)


def binned(spikes, idx, t_run, dt_bin, t_skip):
    t, i = spikes[:, 0], spikes[:, 1].astype(int)
    keep = (t >= t_skip) & np.isin(i, idx)
    return np.histogram(t[keep], np.arange(t_skip, t_run + dt_bin / 2, dt_bin))[0].astype(float)


def rhythm_metrics(spikes, idx, flex, ext, t_run, dt_bin=5.0, t_skip=500.0):
    """脚 MN 集団発火の周期性と、屈筋/伸筋 MN の交代。

    peak_ratio: 2–20 Hz 帯の PSD ピーク ÷ 同じ帯の中央値（平坦な雑音なら 1 前後、明確なリズムなら ≫ 1）。
    flex_ext_corr: 屈筋 MN 群と伸筋 MN 群の発火（10 ms ビン）の相関。歩行なら負（交代）になるはず。
    """
    out = {"leg_mn_active": 0, "rhythm_peak_hz": np.nan, "rhythm_peak_ratio": np.nan, "flex_ext_corr": np.nan}
    if len(spikes) == 0:
        return out
    t, i = spikes[:, 0], spikes[:, 1].astype(int)
    out["leg_mn_active"] = len(np.unique(i[(t >= t_skip) & np.isin(i, idx)]))
    pop = binned(spikes, idx, t_run, dt_bin, t_skip)
    if pop.sum() >= 20:
        pop -= pop.mean()
        f = np.fft.rfftfreq(len(pop), dt_bin / 1000.0)
        psd = np.abs(np.fft.rfft(pop)) ** 2
        band = (f >= 2) & (f <= 20)
        k = np.argmax(psd[band])
        out["rhythm_peak_hz"] = float(f[band][k])
        out["rhythm_peak_ratio"] = float(psd[band][k] / (np.median(psd[band]) + 1e-12))
    a, b = binned(spikes, flex, t_run, 10.0, t_skip), binned(spikes, ext, t_run, 10.0, t_skip)
    if a.std() > 0 and b.std() > 0:
        out["flex_ext_corr"] = float(np.corrcoef(a, b)[0, 1])
    return out


def params(alpha, adapt_b, w_scale, t_run):
    return LIFParams(w_scale=w_scale, size_alpha=alpha, adapt_b=adapt_b, t_run=t_run)


def variant_tag(alpha, adapt_b):
    return ("" if alpha == 0 else f"_alpha{alpha:g}") + ("" if adapt_b == 0 else f"_adapt{adapt_b:g}")


def run_scale(args):
    alpha, adapt_b, w_scale, n_trials, t_run = args
    tag = variant_tag(alpha, adapt_b).lstrip("_") + ("_" if (alpha or adapt_b) else "")
    c = load()
    sets = neuron_sets(c)
    p = params(alpha, adapt_b, w_scale, t_run)
    net = Network.from_connectome(c, p)
    net_gap = net.with_gaps(gap_pairs(c, sets))
    rows = []

    def add(test, cond, metric, value):
        rows.append({"w_scale": w_scale, "test": test, "condition": cond, "metric": metric, "value": float(value)})

    def run(stim, network=net, **kw):
        return simulate(network, stim, p, n_trials=n_trials, seed=int(w_scale * 1000), **kw)

    # R: 無作為刺激の広がり
    for k, rnd in enumerate(sets["random"]):
        r = run([(rnd, 100.0)]).rates
        add("R", f"random{k}", "active_neurons", (r > 1).sum())
        add("R", f"random{k}", "saturated_neurons", (r > 300).sum())
    # T1 摂食
    for cond, stim in [("sugar", [(sets["sugar"], 100.0)]), ("sugar200", [(sets["sugar"], 200.0)]),
                       ("bitter", [(sets["bitter"], 100.0)]),
                       ("sugar+bitter", [(sets["sugar"], 100.0), (sets["bitter"], 100.0)])]:
        r = run(stim).rates
        add("T1", cond, "MN9_L_hz", r[sets["mn9_L"]].mean())
        add("T1", cond, "active_neurons", (r > 1).sum())
        if cond == "sugar":
            np.save(OUT / f"{tag}rates_sugar_w{w_scale:.2f}.npy", r.astype(np.float32))
    # T2 グルーミング
    for cond in ("jo_ce", "jo_ab"):
        r = run([(sets[cond], 100.0)]).rates
        add("T2", cond, "groom_DN_hz", r[sets["groom_dn"]].mean())
        dn_rates = r[sets["dn"]]
        best = max(r[g] for g in sets["groom_dn"])
        add("T2", cond, "groom_DN_rank_among_DNs", (dn_rates > best).sum() + 1 if best > 0 else np.nan)
        add("T2", cond, "DNs_active", (dn_rates > 1).sum())
        add("T2", cond, "active_neurons", (r > 1).sum())
    # T3 逃避（電気シナプス無/有）
    gp = [x for x in gap_pairs(c, sets) if x[1] in set(sets["ttmn"])]
    rec = np.concatenate([sets["gf"], sets["ttmn"]])
    for cond, network in [("chemical_only", net), ("with_gap", net_gap)]:
        for rate in (100.0, 30.0):
            res = run([(sets["lplc2"], rate)], network=network, record=rec)
            r = res.rates
            ctag = f"{cond}@{int(rate)}Hz"
            add("T3", ctag, "GF_hz", r[sets["gf"]].mean())
            add("T3", ctag, "TTMn_hz", r[sets["ttmn"]].mean())
            frac, lag = follow_metrics(res.spikes[0], gp)
            add("T3", ctag, "TTMn_follow_frac", frac)
            add("T3", ctag, "GF_TTMn_lag_ms", lag)
            add("T3", ctag, "active_neurons", (r > 1).sum())
    # T4 歩行リズム
    p_long = params(alpha, adapt_b, w_scale, 3 * t_run)
    res = simulate(net, [(sets["dng100"], 100.0)], p_long, n_trials=1, seed=7, record=sets["leg_mn"])
    for key, val in rhythm_metrics(res.spikes[0], sets["leg_mn"], sets["leg_flexor"], sets["leg_extensor"], p_long.t_run).items():
        add("T4", "dng100", key, val)
    add("T4", "dng100", "leg_MN_mean_hz", res.rates[sets["leg_mn"]].mean())
    add("T4", "dng100", "active_neurons", (res.rates > 1).sum())
    np.save(OUT / f"{tag}legmn_spikes_dng100_w{w_scale:.2f}.npy", res.spikes[0])
    # T5 求愛
    r = run([(sets["p1_like"], 100.0)]).rates
    add("T5", "p1_like", "pIP10_hz", r[sets["pip10"]].mean())
    add("T5", "p1_like", "wing_MN_active", (r[sets["wing_mn"]] > 1).sum())
    add("T5", "p1_like", "active_neurons", (r > 1).sum())
    print(f"alpha={alpha} adapt_b={adapt_b} w_scale={w_scale} done", flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--alpha", type=float, default=0.0)
    ap.add_argument("--adapt-b", type=float, default=0.0)
    ap.add_argument("--scales", type=float, nargs="*")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    scales = [0.44] if args.quick else (args.scales or SCALES)
    n_trials, t_run = (1, 500.0) if args.quick else (3, 1000.0)
    t0 = time.time()
    c = load()
    sets = neuron_sets(c)
    meta = {k: (len(v) if not isinstance(v, list) else [len(x) for x in v]) for k, v in sets.items()}
    meta["gap_pairs"] = gap_pairs(c, sets)
    meta["gap_pairs_types"] = [[c.neurons.instance[a], c.neurons.instance[b], k] for a, b, k in meta["gap_pairs"]]
    meta.update({"scales": scales, "n_trials": n_trials, "t_run_ms": t_run})
    (OUT / f"meta{variant_tag(args.alpha, args.adapt_b) or '_base'}.json").write_text(json.dumps(meta, indent=2, default=int))
    with ProcessPoolExecutor(args.workers) as ex:
        rows = [r for chunk in ex.map(run_scale, [(args.alpha, args.adapt_b, s, n_trials, t_run) for s in scales])
                for r in chunk]
    name = "sweep" + variant_tag(args.alpha, args.adapt_b) + ("_quick" if args.quick else "")
    pd.DataFrame(rows).to_csv(OUT / f"{name}.csv", index=False)
    print(f"done in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
