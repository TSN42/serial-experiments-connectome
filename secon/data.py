"""MaleCNS v1.0 の取得・読み込み。

ニューロン表（型・部位・伝達物質・符号）と、ニューロン間の結合（pre, post, シナプス数）を返す。
大きい結合表 (1.1 GB) は初回にニューロン間だけへ絞って npz にキャッシュする。
"""
from __future__ import annotations

import urllib.request
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.compute as pc
import pyarrow.feather as feather

BASE = "https://storage.googleapis.com/flyem-male-cns/v1.0/connectome-data/flat-connectome/"
FILES = {
    "ann": "body-annotations-male-cns-v1.0-minconf-0.5.feather",
    "nt": "body-neurotransmitters-male-cns-v1.0.feather",
    "w": "connectome-weights-male-cns-v1.0-minconf-0.5.feather",
}
DEFAULT_DIR = Path(__file__).resolve().parent.parent / "data" / "malecns-v1.0"

# 伝達物質 → シナプス符号。
#   fast_only : 速い伝達物質だけ符号を付け、モノアミンと不明は 0（駆動しない）。既定。
#   shiu2024  : Shiu et al. 2024 に近い扱い（ACh 以外の速い伝達物質を抑制、他は興奮として扱う近似）。
# ヒスタミンは視細胞などの伝達物質で、受け手の ort/HisCl1（Cl- チャネル）を開くので抑制とする。
SIGN_PRESETS = {
    "fast_only": {"acetylcholine": 1, "gaba": -1, "glutamate": -1, "histamine": -1},
    "shiu2024": {"acetylcholine": 1, "gaba": -1, "glutamate": -1, "histamine": -1,
                 "dopamine": 1, "octopamine": 1, "serotonin": 1, "tyramine": 1, "unknown": 1},
}


def fetch(key: str, data_dir: Path = DEFAULT_DIR) -> Path:
    path = Path(data_dir) / FILES[key]
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        print(f"downloading {FILES[key]} ...")
        urllib.request.urlretrieve(BASE + FILES[key], path)
    return path


def resolve_nt(nt: pd.DataFrame, min_conf: float = 0.5) -> pd.DataFrame:
    """ニューロンごとの伝達物質を1つに決める。

    consensus（文献の正解 + 型単位の予測）を優先し、unclear なら型単位の予測、
    それも低信頼なら個体単位の予測、どれも閾値未満なら unknown。
    """
    out = nt[["body"]].copy()
    out["nt"] = "unknown"
    out["nt_source"] = "none"
    for col, conf, label in [
        ("predicted_nt", "predicted_nt_confidence", "neuron"),
        ("celltype_predicted_nt", "celltype_predicted_nt_confidence", "celltype"),
    ]:
        ok = nt[col].notna() & (nt[conf] >= min_conf)
        out.loc[ok, "nt"] = nt.loc[ok, col]
        out.loc[ok, "nt_source"] = label
    ok = nt.consensus_nt.notna() & (nt.consensus_nt != "unclear")
    out.loc[ok, "nt"] = nt.loc[ok, "consensus_nt"]
    out.loc[ok, "nt_source"] = np.where(nt.loc[ok, "ground_truth"].notna(), "ground_truth", "consensus")
    unclear = out.nt == "unclear"
    out.loc[unclear, ["nt", "nt_source"]] = ["unknown", "none"]
    return out


@dataclass
class Connectome:
    neurons: pd.DataFrame   # idx 順。bodyId, type, superclass, class, subclass, nt, sign, ...
    pre: np.ndarray         # int32, ニューロン idx
    post: np.ndarray        # int32
    count: np.ndarray       # int32, シナプス数

    @property
    def n(self) -> int:
        return len(self.neurons)

    def select(self, type=None, type_regex=None, superclass=None, cls=None, subclass=None, side=None):
        """条件に合うニューロンの idx を返す（すべて AND）。"""
        m = np.ones(self.n, bool)
        nr = self.neurons
        if type is not None:
            m &= nr.type.isin([type] if isinstance(type, str) else type).to_numpy()
        if type_regex is not None:
            m &= nr.type.fillna("").str.match(type_regex).to_numpy()
        if superclass is not None:
            m &= nr.superclass.isin([superclass] if isinstance(superclass, str) else superclass).to_numpy()
        if cls is not None:
            m &= (nr["class"] == cls).to_numpy()
        if subclass is not None:
            m &= nr.subclass.isin([subclass] if isinstance(subclass, str) else subclass).to_numpy()
        if side is not None:
            s = nr.somaSide.fillna(nr.rootSide)
            m &= (s == side).to_numpy()
        return np.flatnonzero(m)

    def signed_weights(self, preset: str = "fast_only") -> np.ndarray:
        sign = self.neurons.nt.map(SIGN_PRESETS[preset]).fillna(0).to_numpy(np.int8)
        return sign[self.pre].astype(np.float32) * self.count

    def filtered(self, min_syn: int) -> "Connectome":
        keep = self.count >= min_syn
        return Connectome(self.neurons, self.pre[keep], self.post[keep], self.count[keep])


def load(data_dir: Path = DEFAULT_DIR, min_syn: int = 1) -> Connectome:
    data_dir = Path(data_dir)
    ann = pd.read_feather(fetch("ann", data_dir))
    neurons = ann[ann.superclass.notna()].sort_values("bodyId").reset_index(drop=True)
    nt = resolve_nt(pd.read_feather(fetch("nt", data_dir)))
    neurons = neurons.merge(nt, left_on="bodyId", right_on="body", how="left").drop(columns="body")
    neurons["nt"] = neurons.nt.fillna("unknown")
    neurons["nt_source"] = neurons.nt_source.fillna("none")
    neurons.insert(0, "idx", np.arange(len(neurons), dtype=np.int32))

    cache = data_dir / "neuron-edges-v1.0.npz"
    if not cache.exists():
        print("building neuron-neuron edge cache (once) ...")
        table = feather.read_table(fetch("w", data_dir))
        ids = neurons.bodyId.to_numpy()
        keep = pc.and_(pc.is_in(table["body_pre"], value_set=pa_array(ids)),
                       pc.is_in(table["body_post"], value_set=pa_array(ids)))
        table = table.filter(keep)
        pre = np.searchsorted(ids, table["body_pre"].to_numpy()).astype(np.int32)
        post = np.searchsorted(ids, table["body_post"].to_numpy()).astype(np.int32)
        count = table["weight"].to_numpy().astype(np.int32)
        np.savez(cache, pre=pre, post=post, count=count, ids=ids)
    z = np.load(cache)
    assert np.array_equal(z["ids"], neurons.bodyId.to_numpy()), "cache と注釈の bodyId が一致しない"
    conn = Connectome(neurons, z["pre"], z["post"], z["count"])
    return conn.filtered(min_syn) if min_syn > 1 else conn


def pa_array(x):
    import pyarrow as pa
    return pa.array(x)
