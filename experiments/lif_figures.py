"""LIF スイープ（results/lif/sweep*.csv）の図と合否表。03_figures.py から呼ばれる。"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from secon import plotstyle as ps  # noqa: E402
from secon.plotstyle import plt  # noqa: E402

RL, RS, FIG = ROOT / "results" / "lif", ROOT / "results" / "structure", ROOT / "docs" / "figures"
AVALANCHE = 1000   # 無作為 30 細胞の刺激で活動細胞がこれを超えたら「暴走」


def val(d, test, cond, metric):
    s = d[(d.test == test) & (d.condition == cond) & (d.metric == metric)].set_index("w_scale").value
    return s.sort_index()


def verdicts(d):
    """スケールごとの合否表。基準は docs/02 に記載。"""
    rows = []
    for w in sorted(d.w_scale.unique()):
        x = d[d.w_scale == w]
        g = lambda t, c, m: float(x[(x.test == t) & (x.condition == c) & (x.metric == m)].value.iloc[0])
        rnd = np.median([g("R", f"random{k}", "active_neurons") for k in range(3)])
        su, bi, sb = g("T1", "sugar", "MN9_L_hz"), g("T1", "bitter", "MN9_L_hz"), g("T1", "sugar+bitter", "MN9_L_hz")
        ce, ab, rank = g("T2", "jo_ce", "groom_DN_hz"), g("T2", "jo_ab", "groom_DN_hz"), g("T2", "jo_ce", "groom_DN_rank_among_DNs")
        gap_f, gap_l = g("T3", "with_gap@30Hz", "TTMn_follow_frac"), g("T3", "with_gap@30Hz", "GF_TTMn_lag_ms")
        ch_f, ch_l = g("T3", "chemical_only@30Hz", "TTMn_follow_frac"), g("T3", "chemical_only@30Hz", "GF_TTMn_lag_ms")
        gf30 = g("T3", "chemical_only@30Hz", "GF_hz")
        ratio, peak, fe = g("T4", "dng100", "rhythm_peak_ratio"), g("T4", "dng100", "rhythm_peak_hz"), g("T4", "dng100", "flex_ext_corr")
        pip = g("T5", "p1_like", "pIP10_hz")
        rows.append({
            "w_scale": w,
            "R 局在（無作為刺激の活動 < 1000）": rnd < AVALANCHE,
            "T1 糖→MN9 / 苦味で抑制": (su > 10) and (bi < 1) and (sb < 0.5 * su),
            "T2 JO-C/E→グルーミング DN（選択的）": (ce > 10) and (ab < 0.5 * ce) and (rank <= 10),
            "T3a GF→TTMn 1:1・<1.5 ms（ギャップ有）": (gap_f >= 0.9) and (gap_l < 1.5),
            "T3b 化学シナプスだけでは速い 1:1 にならない": not ((ch_f >= 0.9) and (ch_l < 1.5)),
            "T3c GF が飽和しない（LPLC2 30 Hz で < 100 Hz）": gf30 < 100,
            "T4 DNg100→脚 MN リズム（4–20 Hz・交代あり）": (ratio > 5) and (peak >= 4) and (fe < -0.2),
            "T5 P1 相当→pIP10": pip > 10,
            "_random_active": rnd, "_MN9_sugar": su, "_GF30": gf30, "_flex_ext": fe, "_peak_ratio": ratio, "_peak_hz": peak,
        })
    return pd.DataFrame(rows)


def fig_phase(d, name, suptitle=None):
    fig, axes = plt.subplots(2, 3, figsize=(13, 7))
    if suptitle:
        fig.suptitle(suptitle, x=0.01, ha="left", fontsize=12, fontweight="bold")
    rnd = pd.concat([val(d, "R", f"random{k}", "active_neurons") for k in range(3)], axis=1)
    aval = rnd.median(axis=1) > AVALANCHE
    x = rnd.index.to_numpy()

    def shade(ax):
        if aval.any():
            ax.axvspan(x[aval.to_numpy()].min() - 0.01, x.max() + 0.02, color=ps.GRID, alpha=0.6, lw=0)

    def lines(ax, series, ylabel, title, logy=False):
        shade(ax)
        for (label, s), col in zip(series, ps.SERIES):
            ax.plot(s.index, s.values, color=col, marker="o", ms=4, label=label)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.set_xlabel("w_scale（1 シナプスの重み倍率）")
        if logy:
            ax.set_yscale("log")
        if len(series) > 1:
            ax.legend()

    lines(axes[0, 0], [("無作為 30 細胞（中央値）", rnd.median(axis=1)),
                       ("糖 GRN", val(d, "T1", "sugar", "active_neurons"))],
          "活動ニューロン数 (>1 Hz)", "R 活動の広がり（灰色 = 暴走域）", logy=True)
    lines(axes[0, 1], [("糖", val(d, "T1", "sugar", "MN9_L_hz")), ("糖+苦味", val(d, "T1", "sugar+bitter", "MN9_L_hz")),
                       ("苦味", val(d, "T1", "bitter", "MN9_L_hz"))], "左 MN9 [Hz]", "T1 摂食：MN9")
    lines(axes[0, 2], [("JO-C/E（風・重力）", val(d, "T2", "jo_ce", "groom_DN_hz")),
                       ("JO-A/B（聴覚）", val(d, "T2", "jo_ab", "groom_DN_hz"))], "aDN1/aDN2 [Hz]", "T2 触角グルーミング DN")
    lines(axes[1, 0], [("LPLC2 100 Hz", val(d, "T3", "chemical_only@100Hz", "GF_hz")),
                       ("LPLC2 30 Hz", val(d, "T3", "chemical_only@30Hz", "GF_hz"))], "巨大線維 GF [Hz]", "T3 逃避：GF の発火")
    lines(axes[1, 1], [("屈筋−伸筋の相関", val(d, "T4", "dng100", "flex_ext_corr"))], "相関係数",
          "T4 DNg100→脚 MN：交代（負なら交代）")
    axes[1, 1].axhline(0, color=ps.TEXT2, lw=1)
    lines(axes[1, 2], [("pIP10", val(d, "T5", "p1_like", "pIP10_hz"))], "pIP10 [Hz]", "T5 P1 相当→pIP10")
    fig.tight_layout()
    ps.save(fig, FIG / name)


def fig_rhythm(w):
    path = RL / f"legmn_spikes_dng100_w{w:.2f}.npy"
    if not path.exists():
        return
    sys.path.insert(0, str(ROOT))
    from secon.data import load
    c = load()
    sp = np.load(path)
    sub = c.neurons.subclass.to_numpy()
    fig, ax = plt.subplots(figsize=(10, 3.6))
    t0 = 500.0
    m = (sp[:, 0] >= t0) & (sp[:, 0] < t0 + 1500)
    idx = sp[m, 1].astype(int)
    order = {v: k for k, v in enumerate(np.unique(idx))}
    for code, label, col in (("fl", "前脚 MN", ps.SERIES[0]), ("ml", "中脚 MN", ps.SERIES[1]), ("hl", "後脚 MN", ps.SERIES[2])):
        k = sub[idx] == code
        ax.scatter(sp[m, 0][k], [order[i] for i in idx[k]], s=2, color=col, label=label, linewidths=0)
    ax.set_xlabel("時間 [ms]")
    ax.set_ylabel("発火した脚 MN")
    ax.set_title(f"DNg100 (BDN2) 刺激中の脚運動ニューロン（w_scale={w}）")
    ax.legend(markerscale=4, loc="upper right")
    ps.save(fig, FIG / "fig5_dng100_raster.png")


def auc_active(score, active):
    """score が LIF で活動する細胞を見分けられるか（Mann–Whitney の AUC）。"""
    ranks = pd.Series(score).rank().to_numpy()
    n1, n0 = active.sum(), (~active).sum()
    return float((ranks[active].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def fig_influence_vs_lif(w, lam_plot=0.5):
    lif_path = RL / f"rates_sugar_w{w:.2f}.npy"
    if not lif_path.exists():
        return None
    r = np.load(lif_path)
    from secon.data import load
    from secon.groups import sensory_groups
    seeds = np.zeros(len(r), bool)
    seeds[sensory_groups(load())["taste: sugar"]] = True     # 種（糖 GRN）は除いて評価
    aucs = {}
    for lam in (0.99, 0.5):
        path = RS / f"influence_signed_lam{lam}_taste_sugar.npy"
        if not path.exists():
            continue
        inf = np.load(path)
        lab = (r > 1)[~seeds]
        aucs[f"lam{lam}"] = auc_active(inf[~seeds], lab)
        if lam == lam_plot:
            score = inf[~seeds]
            xs = np.sign(score) * np.log10(1 + np.abs(score) * 1e6)
            fig, ax = plt.subplots(figsize=(6, 4.2))
            jitter = np.random.default_rng(0).uniform(-0.3, 0.3, (~lab).sum())
            ax.scatter(xs[~lab], jitter, s=2, color=ps.NEUTRAL, linewidths=0, label="LIF で不活動")
            ax.scatter(xs[lab], np.log10(1 + r[~seeds][lab]), s=6, color=ps.SERIES[0], linewidths=0, label="LIF で活動 (>1 Hz)")
            ax.set_xlabel("符号付き線形影響度（sign·log10(1+|x|·1e6)）")
            ax.set_ylabel("LIF 発火率 log10(1+Hz)")
            ax.set_title(f"糖 GRN 刺激：線形影響度(λ={lam}) vs LIF w={w}（AUC={aucs[f'lam{lam}']:.2f}）")
            ax.legend(markerscale=3)
            ps.save(fig, FIG / "fig6_influence_vs_lif.png")
    return aucs


def main(work_w=0.44):
    out = {}
    for mode_file, fig_name in (("sweep.csv", "fig4_lif_phase.png"), ("sweep_alpha0.5.csv", "fig7_lif_phase_alpha0.5.png"),
                                ("sweep_adapt1.csv", "fig8_lif_phase_adapt1.png")):
        if (RL / mode_file).exists():
            d = pd.read_csv(RL / mode_file)
            fig_phase(d, fig_name, {"sweep.csv": None, "sweep_alpha0.5.csv": "細胞の大きさ補正 α = 0.5",
                                    "sweep_adapt1.csv": "スパイク頻度順応 b = 1 mV, τ = 100 ms"}[mode_file])
            v = verdicts(d)
            v.to_csv(RL / mode_file.replace("sweep", "verdicts"), index=False)
            out[mode_file] = v
    fig_rhythm(work_w)
    aucs = fig_influence_vs_lif(work_w)
    if aucs:
        (RL / "influence_vs_lif_auc.json").write_text(__import__("json").dumps(aucs, indent=2))
    return out


if __name__ == "__main__":
    for k, v in main().items():
        print(k)
        print(v.to_string())
