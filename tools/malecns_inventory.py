"""MaleCNS v1.0 の中身を実ファイルから棚卸しする。

使い方:
    python tools/malecns_inventory.py               # アノテーション + 伝達物質予測 (~57 MB)
    python tools/malecns_inventory.py --weights     # 全結合表 (~1.1 GB) も読む

データは CC-BY 4.0 (https://male-cns.janelia.org/download/)。
"""
import argparse
import re
import urllib.request
from pathlib import Path

import pandas as pd

BASE = "https://storage.googleapis.com/flyem-male-cns/v1.0/connectome-data/flat-connectome/"
FILES = {
    "ann": "body-annotations-male-cns-v1.0-minconf-0.5.feather",
    "nt": "body-neurotransmitters-male-cns-v1.0.feather",
    "w": "connectome-weights-male-cns-v1.0-minconf-0.5.feather",
}
# 筋名で命名された運動ニューロン (例: "Ti flexor MN", "DLMn c-f", "b1 MN")。"MNhl64" のような系統名は除外
MUSCLE_NAMED = re.compile(r"(?:MN$|MN |Mn|MNs)")
SYSTEMATIC = re.compile(r"^MN[a-z]{2}\d+")


def fetch(key, data_dir):
    path = data_dir / FILES[key]
    if not path.exists():
        data_dir.mkdir(parents=True, exist_ok=True)
        print(f"downloading {FILES[key]} ...")
        urllib.request.urlretrieve(BASE + FILES[key], path)
    return path


def section(title):
    print(f"\n## {title}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", type=Path, default=Path("data/malecns-v1.0"))
    ap.add_argument("--weights", action="store_true", help="1.1 GB の結合表も解析する")
    args = ap.parse_args()

    ann = pd.read_feather(fetch("ann", args.data_dir))
    neurons = ann[ann.superclass.notna()]

    section("neurons")
    print(f"neurons (superclass 付き): {len(neurons):,}  typed: {neurons.type.notna().sum():,}  "
          f"types: {neurons.type.nunique():,}")
    print(neurons.superclass.value_counts().to_string())

    section("sensory neurons by modality (class)")
    sensory = neurons[neurons.superclass.str.contains("sensory")]
    print(sensory["class"].value_counts(dropna=False).to_string())
    print(sensory.subclass.value_counts().head(20).to_string())

    section("photoreceptors / lamina coverage")
    for pat in [r"^R1-R6$", r"^R[78]", r"^L[1-5]$"]:
        sel = neurons[neurons.type.fillna("").str.match(pat)]
        print(f"{pat:10s} n={len(sel):5d}  rootSide={sel.rootSide.value_counts().to_dict()}  "
              f"somaSide={sel.somaSide.value_counts().to_dict()}")
    cols = neurons[neurons.assignedOlHex1.notna()]
    per_side = cols.groupby("somaSide")[["assignedOlHex1", "assignedOlHex2"]].apply(
        lambda d: len(d.drop_duplicates()))
    print("hex columns per side:", per_side.to_dict())

    section("motor / efferent")
    motor = neurons[neurons.superclass.isin(["vnc_motor", "cb_motor"])]
    print(f"motor neurons: {len(motor)}")
    print("subclass:", motor.subclass.value_counts(dropna=False).to_dict())
    t = motor.type.fillna("")
    named = t.str.contains(MUSCLE_NAMED) & ~t.str.match(SYSTEMATIC)
    print(f"muscle-named: {named.sum()} / {len(motor)}")
    print(motor[named].type.value_counts().head(25).to_string())
    for sc in ["cb_endocrine", "vnc_endocrine", "vnc_efferent", "ENS"]:
        sel = neurons[neurons.superclass == sc]
        print(f"{sc}: {len(sel)} {sel.type.value_counts().head(6).to_dict()}")

    section("descending / ascending")
    for sc in ["descending_neuron", "ascending_neuron"]:
        sel = neurons[neurons.superclass == sc]
        print(f"{sc}: {len(sel)} neurons, {sel.type.nunique()} types, MANC-matched {sel.mancType.notna().sum()}")

    section("sexual dimorphism / fru-dsx")
    print(neurons.dimorphism.value_counts().to_dict())
    print(neurons.fruDsx.value_counts().to_dict())

    section("neurotransmitter (consensus)")
    nt = pd.read_feather(fetch("nt", args.data_dir))
    nt = nt[nt.body.isin(neurons.bodyId)]
    print(nt.consensus_nt.value_counts(dropna=False).to_string())

    if args.weights:
        section("connectivity")
        w = pd.read_feather(fetch("w", args.data_dir))
        print(f"segment-pair rows: {len(w):,}  total synapses (PSD): {w.weight.sum():,}")
        ids = set(neurons.bodyId)
        nn = w[w.body_pre.isin(ids) & w.body_post.isin(ids)]
        print(f"neuron-neuron edges: {len(nn):,}  synapses: {nn.weight.sum():,}")
        for th in (3, 5):
            print(f"  edges with weight >= {th}: {(nn.weight >= th).sum():,}")
        print(f"CSR (fp32 weight + int32 index): {len(nn) * 8 / 1e9:.2f} GB")


if __name__ == "__main__":
    main()
