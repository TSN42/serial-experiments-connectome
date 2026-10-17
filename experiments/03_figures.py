"""results/ の数値から docs/figures/ の図を作る。

    python experiments/03_figures.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from secon import plotstyle as ps  # noqa: E402
from secon.plotstyle import plt  # noqa: E402

RS, RL, FIG = ROOT / "results" / "structure", ROOT / "results" / "lif", ROOT / "docs" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

NT_GROUPS = [("acetylcholine", "ACh（興奮）"), ("glutamate", "Glu（抑制扱い）"), ("gaba", "GABA（抑制）"),
             ("histamine", "ヒスタミン（抑制）"), ("monoamine", "モノアミン（符号なし）"), ("unknown", "不明（符号なし）")]


def fig_composition():
    df = pd.read_csv(RS / "superclass_by_nt.csv", index_col=0)
    df["monoamine"] = df[[c for c in ("dopamine", "octopamine", "serotonin") if c in df]].sum(axis=1)
    df = df[[k for k, _ in NT_GROUPS]]
    total = df.sum(axis=1)
    df = df.loc[total[total >= 40].sort_values().index]
    frac = df.div(df.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(8, 5.2))
    left = np.zeros(len(frac))
    colors = ps.SERIES[:5] + [ps.NEUTRAL]
    for (k, label), col in zip(NT_GROUPS, colors):
        ax.barh(frac.index, frac[k], left=left, height=0.62, color=col, label=label, edgecolor=ps.SURFACE, linewidth=1)
        left += frac[k].to_numpy()
    for y, n in enumerate(df.sum(axis=1)):
        ax.text(1.01, y, f"{int(n):,}", va="center", fontsize=8, color=ps.TEXT2)
    ax.set_xlim(0, 1)
    ax.set_xlabel("ニューロンの割合（右の数字は細胞数）")
    ax.grid(axis="y", visible=False)
    ax.set_title("MaleCNS v1.0：スーパークラスごとの伝達物質")
    ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.45, -0.1))
    ps.save(fig, FIG / "fig1_composition.png")


def fig_lr():
    lr = pd.read_csv(RS / "lr_asymmetry.csv", index_col=0)
    one = lr[(lr.cells_L == 1) & (lr.cells_R == 1)]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), gridspec_kw={"width_ratios": [1.3, 1]})
    ax = axes[0]
    ax.hist(one.log2_inp_ratio.clip(-5, 5), bins=np.linspace(-5, 5, 81), color=ps.SERIES[0], edgecolor=ps.SURFACE, linewidth=0.5)
    ax.set_yscale("log")
    for x in (-1, 1):
        ax.axvline(x, color=ps.TEXT2, linewidth=1)
    for name, ytext, dx in (("MN9", 30, 0.6), ("DNb06", 6, -0.6), ("MNhl02", 8, 0.0), ("GNG494", 30, -0.4)):
        if name in one.index:
            x = float(np.clip(one.loc[name, "log2_inp_ratio"], -5, 5))
            ax.annotate(name, (x, 1.2), xytext=(x + dx, ytext), ha="center", fontsize=8, color=ps.TEXT,
                        arrowprops=dict(arrowstyle="-", color=ps.TEXT2, lw=0.8))
    ax.set_xlabel("log2(左の入力シナプス数 / 右)　※左右1細胞ずつの型")
    ax.set_ylabel("型の数")
    ax.set_title(f"左右差の分布（{len(one):,} 型）")
    ax = axes[1]
    flag = lr.log2_inp_ratio.abs() > 1
    g = (lr[flag].groupby("superclass").size() / lr.groupby("superclass").size()).fillna(0)
    n = lr.groupby("superclass").size()
    g = g[n >= 5].sort_values()
    ax.barh(g.index, g.values, height=0.6, color=ps.SERIES[1])
    for y, v in enumerate(g.values):
        ax.text(v + 0.005, y, f"{v:.0%}", va="center", fontsize=8, color=ps.TEXT2)
    ax.set_xlabel("左右で入力が2倍以上違う型の割合")
    ax.grid(axis="y", visible=False)
    ax.set_title("スーパークラス別")
    fig.tight_layout()
    ps.save(fig, FIG / "fig2_lr_asymmetry.png")


def heatmap(ax, mat, fmt, cmap, label, reverse=False):
    im = ax.imshow(mat.to_numpy(), cmap=cmap.reversed() if reverse else cmap, aspect="auto")
    ax.set_xticks(range(mat.shape[1]), mat.columns, rotation=40, ha="right")
    ax.set_yticks(range(mat.shape[0]), mat.index)
    ax.grid(False)
    vals = mat.to_numpy()
    lo, hi = np.nanmin(vals), np.nanmax(vals)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = vals[i, j]
            if np.isfinite(v):
                x = (v - lo) / (hi - lo + 1e-12)
                dark = (1 - x if reverse else x) > 0.55
                ax.text(j, i, fmt(v), ha="center", va="center", fontsize=7, color="white" if dark else ps.TEXT)
    cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label(label)
    cb.outline.set_visible(False)


def fig_sensory_effector():
    tr = pd.read_csv(RS / "traversal_sensory_to_effector.csv")
    inf = pd.read_csv(RS / "influence_sensory_to_effector.csv")
    order_s = list(dict.fromkeys(tr.sensory))
    order_e = list(dict.fromkeys(tr.effector))
    m1 = tr.pivot(index="sensory", columns="effector", values="median_layer").loc[order_s, order_e]
    u = inf[~inf.signed & (inf.lam == 0.5)].pivot(index="sensory", columns="effector", values="mean").loc[order_s, order_e]
    m2 = np.log10(u.clip(lower=1e-30))
    fig, axes = plt.subplots(1, 2, figsize=(15, 5.4))
    heatmap(axes[0], m1, lambda v: f"{v:.1f}", ps.SEQ, "到達ステップ（中央値）", reverse=True)
    axes[0].set_title("確率的走査：感覚から何ステップで届くか（濃いほど近い）")
    heatmap(axes[1], m2, lambda v: f"{v:.1f}", ps.SEQ, "log10 影響度（種 1 細胞あたり、群平均）")
    axes[1].set_title("線形影響度（BANC 方式、符号なし、λ=0.5、濃いほど強い）")
    axes[1].set_yticklabels([])
    fig.tight_layout()
    ps.save(fig, FIG / "fig3_sensory_to_effector.png")


def main():
    ps.setup()
    if (RS / "superclass_by_nt.csv").exists():
        fig_composition()
        fig_lr()
    if (RS / "influence_sensory_to_effector.csv").exists():
        fig_sensory_effector()
    if (RL / "sweep.csv").exists():
        import lif_figures  # noqa: F401  （LIF 側の図は別ファイル）
        lif_figures.main()


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).parent))
    main()
