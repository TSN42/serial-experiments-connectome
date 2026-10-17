import numpy as np
import pandas as pd

from secon import structure as S
from secon.data import Connectome


def toy():
    """0 → 1 → 2 の鎖と、左右ペア（型 A の L/R）。A_R だけ入力が少ない。"""
    nr = pd.DataFrame({
        "idx": np.arange(5), "bodyId": np.arange(5) + 100,
        "type": ["S", "A", "A", "M", "M"], "somaSide": ["L", "L", "R", "L", "R"], "rootSide": [None] * 5,
        "superclass": ["cb_sensory", "cb_intrinsic", "cb_intrinsic", "cb_motor", "cb_motor"],
        "class": [None] * 5, "subclass": [None] * 5,
        "nt": ["acetylcholine", "acetylcholine", "gaba", "glutamate", "glutamate"],
    })
    pre = np.array([0, 0, 1, 2], np.int32)
    post = np.array([1, 2, 3, 4], np.int32)
    count = np.array([40, 5, 10, 10], np.int32)
    return Connectome(nr, pre, post, count)


def test_lr_asymmetry_flags_incomplete_side():
    a = S.lr_asymmetry(toy())
    assert a.loc["A", "log2_inp_ratio"] > 2      # A_L は 40、A_R は 5
    assert abs(a.loc["M", "log2_inp_ratio"]) < 1e-9


def test_traversal_reaches_chain_in_order():
    lay = S.traversal_layers(toy(), [0], n_runs=8)
    assert lay[0] == 0 and lay[1] == 1 and lay[3] == 2


def test_linear_influence_matches_closed_form():
    c = toy()
    W = S.signed_matrix(c, min_syn=1).toarray()
    lead = 1.0  # 冪零なので固有値 0。scale は任意
    r = S.linear_influence(S.signed_matrix(c, min_syn=1), [0], scale=lead, lam=0.01)
    Wt = W * 0.01
    expect = np.linalg.solve(np.eye(5) - Wt, np.eye(5)[0])
    assert np.allclose(r, expect)


def test_input_sources_sum_to_one():
    s = S.input_sources(toy(), [3, 4])
    assert np.isclose(s.sum(), 1.0)
