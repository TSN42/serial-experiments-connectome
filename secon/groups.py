"""感覚モダリティと効果器（運動ニューロンなど）のグループ定義。MaleCNS v1.0 の注釈名に合わせてある。

味覚の対応は Tastekin et al. (Cell 2026) の受容体ドライバー照合に従う：
LB1a–d = 苦味 (Gr33a)、LB3a = 水 (ppk28)、LB3b–c = 糖 (Gr64f)。
（Shiu 2024 が使った旧「sugar/water」クラスには水・高塩の細胞が混ざっていた。）
"""
from __future__ import annotations

import numpy as np

from .data import Connectome

LEG_NERVES = ["ProLN", "MesoLN", "MetaLN"]


def sensory_groups(c: Connectome) -> dict[str, np.ndarray]:
    nr = c.neurons
    leg = nr.entryNerve.isin(LEG_NERVES).to_numpy()
    cls = nr["class"].to_numpy()
    sub = nr.subclass.to_numpy()
    g = {
        "photoreceptors": c.select(type_regex=r"^(R1-R6|R7|R8)"),
        "olfactory ORN": c.select(type_regex=r"^ORN_"),
        "taste: sugar": c.select(type=["LB3b", "LB3c"]),
        "taste: bitter": c.select(type_regex=r"^LB1[a-d]$"),
        "taste: water": c.select(type="LB3a"),
        "JO auditory (A/B)": c.select(type_regex=r"^JO-[AB]"),
        "JO wind/gravity (C/E)": c.select(type_regex=r"^JO-[CE]"),
        "leg chordotonal": np.flatnonzero((sub == "chordotonal organ") & leg),
        "campaniform (all)": np.flatnonzero(sub == "campaniform sensilla"),
        "leg touch bristles": np.flatnonzero((cls == "mechanosensory_tactile") & leg),
        "haltere": np.flatnonzero(sub == "haltere"),
        "thermo": np.flatnonzero(cls == "thermosensory"),
        "hygro": np.flatnonzero(cls == "hygrosensory"),
    }
    return {k: v for k, v in g.items() if len(v)}


MN_SUBCLASS = {"fl": "MN front leg", "ml": "MN mid leg", "hl": "MN hind leg", "wm": "MN wing",
               "hm": "MN haltere", "nm": "MN neck", "pm": "MN proboscis/pharynx", "am": "MN antenna",
               "ad": "MN abdomen"}


def effector_groups(c: Connectome) -> dict[str, np.ndarray]:
    mn = c.neurons.superclass.isin(["vnc_motor", "cb_motor"]).to_numpy()
    sub = c.neurons.subclass.to_numpy()
    g = {name: np.flatnonzero(mn & (sub == code)) for code, name in MN_SUBCLASS.items()}
    g["descending neurons"] = c.select(superclass="descending_neuron")
    g["endocrine"] = c.select(superclass=["cb_endocrine", "vnc_endocrine"])
    return g
