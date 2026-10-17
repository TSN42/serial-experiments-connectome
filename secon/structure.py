"""MaleCNS の構造解析。動的モデルを回す前に「配線そのもの」が何を言っているかを見る道具。

- lr_asymmetry      : 左右一対の型で入力シナプス数を比べ、再構成の欠け（片側だけ極端に少ない）を検出する
- traversal_layers  : Schlegel et al. (2021, 2024) 型の確率的走査。感覚ニューロンから何層目に到達するか
- linear_influence  : Bates et al. (2026, BANC) の影響度スコア。tau dr/dt = -r + W r + s の定常解
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .data import Connectome


def input_totals(conn: Connectome) -> np.ndarray:
    return np.bincount(conn.post, weights=conn.count, minlength=conn.n)


def output_totals(conn: Connectome) -> np.ndarray:
    return np.bincount(conn.pre, weights=conn.count, minlength=conn.n)


def lr_asymmetry(conn: Connectome, min_cells_per_side: int = 1) -> pd.DataFrame:
    """型ごとに左右の1細胞あたり入力シナプス数を比べる。

    同じ型の左右は形態も入力数もほぼ等しいはずなので、比が大きく外れる型は
    片側の再構成が欠けている可能性が高い（MN9 の右側など）。
    """
    nr = conn.neurons
    side = nr.somaSide.fillna(nr.rootSide)
    df = pd.DataFrame({"type": nr.type, "side": side, "superclass": nr.superclass,
                       "inp": input_totals(conn), "out": output_totals(conn)})
    df = df[df.type.notna() & df.side.isin(["L", "R"])]
    g = df.groupby(["type", "side"]).agg(cells=("inp", "size"), inp=("inp", "mean"), out=("out", "mean"),
                                        superclass=("superclass", "first")).unstack("side")
    g = g[(g[("cells", "L")] >= min_cells_per_side) & (g[("cells", "R")] >= min_cells_per_side)]
    out = pd.DataFrame({
        "superclass": g[("superclass", "L")],
        "cells_L": g[("cells", "L")].astype(int), "cells_R": g[("cells", "R")].astype(int),
        "inp_L": g[("inp", "L")], "inp_R": g[("inp", "R")],
        "out_L": g[("out", "L")], "out_R": g[("out", "R")],
    })
    out["log2_inp_ratio"] = np.log2((out.inp_L + 1) / (out.inp_R + 1))
    out["log2_out_ratio"] = np.log2((out.out_L + 1) / (out.out_R + 1))
    return out.sort_values("log2_inp_ratio", key=np.abs, ascending=False)


def input_fraction_matrix(conn: Connectome) -> sp.csr_matrix:
    """A[post, pre] = 入力シナプスに占める割合（列和ではなく行和が 1）。"""
    tot = input_totals(conn)
    vals = conn.count / np.maximum(tot[conn.post], 1)
    return sp.csr_matrix((vals.astype(np.float32), (conn.post, conn.pre)), shape=(conn.n, conn.n))


def traversal_layers(conn: Connectome, seeds, n_runs: int = 64, threshold: float = 0.3,
                     max_steps: int = 20, seed: int = 0) -> np.ndarray:
    """確率的走査で各ニューロンに到達する平均ステップ（seed = 0）。到達しなければ NaN。

    各ステップで、未到達ニューロンは「既到達ニューロンからの入力割合 / threshold」の確率（上限 1）で到達する。
    Schlegel et al. 2021 (hemibrain) / 2024 (FlyWire) の information flow と同じ考え方。
    """
    rng = np.random.default_rng(seed)
    A = input_fraction_matrix(conn)
    reached = np.zeros((conn.n, n_runs), bool)
    layer = np.full((conn.n, n_runs), np.nan, np.float32)
    reached[seeds] = True
    layer[seeds] = 0
    for step in range(1, max_steps + 1):
        frac = A @ reached.astype(np.float32)
        p = np.minimum(frac / threshold, 1.0)
        new = (~reached) & (rng.random(p.shape) < p)
        if not new.any():
            break
        reached |= new
        layer[new] = step
    with np.errstate(all="ignore"):
        return np.nanmean(layer, axis=1)


def signed_matrix(conn: Connectome, preset: str = "fast_only", min_syn: int = 5) -> sp.csr_matrix:
    c = conn.filtered(min_syn) if min_syn > 1 else conn
    w = c.signed_weights(preset)
    nz = w != 0
    return sp.csr_matrix((w[nz], (c.post[nz], c.pre[nz])), shape=(c.n, c.n))


def leading_mode(W: sp.csr_matrix):
    """絶対値最大の固有値と、その固有ベクトルの |v|^2（ARPACK）。"""
    vals, vecs = spla.eigs(W.astype(np.float64), k=1, which="LM", tol=1e-6, maxiter=5000)
    return complex(vals[0]), np.abs(vecs[:, 0]) ** 2


def linear_influence(W: sp.csr_matrix, seeds, scale: float, lam: float = 0.99,
                     tol: float = 1e-10, max_iter: int = 2000) -> np.ndarray:
    """定常状態 r = (I - W~)^-1 s、W~ = lam / scale × W（scale は W のスペクトル半径）。

    ConnectomeInfluenceCalculator と同じ定式化。疎な LU は埋め込みが大きいので Neumann 級数で解く。
    """
    Wt = W * (lam / scale)
    s = np.zeros(W.shape[0])
    s[seeds] = 1.0
    r, term = s.copy(), s.copy()
    for _ in range(max_iter):
        term = Wt @ term
        r += term
        if np.abs(term).max() < tol:
            break
    return r


def input_sources(conn: Connectome, targets, by: str = "superclass") -> pd.Series:
    """targets 群が受ける入力シナプスの、送り手カテゴリ別の割合。"""
    m = np.isin(conn.post, targets)
    lab = conn.neurons[by].fillna("(none)").to_numpy()[conn.pre[m]]
    s = pd.Series(conn.count[m]).groupby(lab).sum()
    return (s / s.sum()).sort_values(ascending=False)
