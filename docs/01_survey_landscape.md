# 調査報告 #01：コネクトームから「放てるハエ」までの全体像

調査日: 2026-10-13 ／ 対象: MaleCNS を中心としたショウジョウバエの全脳・全身シミュレーションのデータ・モデル・ライブラリ

> この報告の数値には出典が3種類ある。
> - **[実測]** MaleCNS v1.0 の公開ファイルをこのリポジトリの `tools/malecns_inventory.py` で集計した値。
> - **[論文]** 一次文献（論文本体・公式ページ）から確認した値。
> - **[二次]** ニュース記事や README しか確認できていない値。裏取りが必要。

---

## 0. 結論（先に）

1. **土台は MaleCNS v1.0 で決まり。** 1匹のオスの脳・両視葉・腹側神経索 (VNC) を、首の結合 (neck connective) が切れていない状態でまるごと再構成した、現時点で唯一の「感覚ニューロン → 中枢 → 運動ニューロン」が同じ個体内でつながったデータ。ニューロン 166,700 個、ニューロン間シナプス 1.24 億、CC-BY 4.0 [実測/論文]。
2. **配線の正しさは十分に高い。足りないのは配線以外。** 弱点は順に、(a) ニューロン単体の生理（全部同じ LIF で近似されている）、(b) 神経修飾・ペプチド・ギャップ結合が見えない、(c) **運動ニューロン → 筋 → 関節トルクの変換**、(d) 感覚器の物理と変換、の4つ。
3. **全中枢の動的モデルで検証済みのものは存在しない。** Shiu ら (2024) の LIF 全脳モデルは FlyWire の脳だけで VNC が無い。視覚 (flyvis)、歩行 CPG (VNC)、摂食・グルーミングといった個別回路はそれぞれ検証されている。MaleCNS 全体のシミュレーションは自分たちで組むことになるが、計算量は小さい（結合行列は 1匹あたり約 0.2 GB）。
4. **身体は NeuroMechFly/FlyGym 2.x（歩行・感覚）と flybody（飛行）の2択。** どちらもメスの個体から作られ、関節を位置制御で動かす（筋は無い）。筋骨格モデルは前脚と翅ヒンジに部分的にあるだけ。
5. **「コネクトームで体を動かす」先行例はすでにある**（Eon Systems の embodied fly、Doom を遊ぶハエなど）。ただし全例が「少数の下行ニューロン (DN) を読んで手設計の歩行コントローラに渡す」近道を使っており、VNC も運動ニューロンも通していない。MaleCNS を使えば **運動ニューロンの活動で体を直接動かす** ことが原理的に初めて可能になる。ここが本プロジェクトで最初に狙う価値のある挑戦。
6. プシュケー・チップ（読み出し → 送信 → 受信 → 書き込み）は、計算上はほぼタダで載る。意味のある遊びになるかは「どこを読み、どこに書くか」の部位が **生理的に検証された回路** かどうかで決まる（§9）。

---

## 1. 全体像：何をつなぐ必要があるか

```
 環境(物理・匂い・光・音・他個体)
   │  ①感覚器の物理（個眼の光学、匂い分子の拡散、接触、関節角）
   ▼
 感覚ニューロン(求心性)  ── ②変換モデル（光→電位、匂い→受容体応答…）
   │
   ▼
 中枢(脳 + 視葉 + VNC)    ── ③ダイナミクス（ニューロンモデル、シナプス符号・強度、修飾）
   │                         ← ここにプシュケー・チップを挿す
   ▼
 運動ニューロン(遠心性)   ── ④MN → 筋活性 → 力 → 関節トルク
   │
   ▼
 身体(剛体・関節・接着・翅) ── ⑤物理エンジン(MuJoCo)
   │
   └──→ 環境へ（→①に戻る）
```

| 区間 | 使えるもの | 解剖学的な正しさ（筆者評価） |
|---|---|---|
| 中枢の配線 | MaleCNS / BANC / FlyWire / MANC | ★★★★★ |
| シナプス符号（興奮/抑制） | 伝達物質予測（EM画像からCNNで推定） | ★★★★ |
| シナプス強度 | シナプス数で代用 | ★★ |
| ニューロン単体の生理 | 一律 LIF、flyvis のような学習済みパラメータ | ★ 〜 ★★ |
| 神経修飾・ホルモン・ギャップ結合 | ほぼ無い（内分泌細胞の位置だけ分かる） | ☆ |
| 感覚器 → 感覚ニューロン | 視覚は個眼↔カラム対応あり、嗅覚は受容体IDあり、他は粗い | ★★ 〜 ★★★ |
| MN → 筋 → 関節 | MN の筋同定は約半数、筋骨格モデルは前脚と翅ヒンジのみ | ★ 〜 ★★ |
| 身体形状・物理 | micro-CT・共焦点ベース、MuJoCo | ★★★★（ただしメス） |

---

## 2. コネクトーム・データ

### 2.1 主要データセット比較

| データ | 個体 | 範囲 | ニューロン数 | 手法 | ライセンス/アクセス |
|---|---|---|---|---|---|
| **MaleCNS v1.0** (Berg et al., *Cell* 2026) | オス1匹 | 脳 + 両視葉 + VNC（首つながったまま） | 166,700 [実測] | FIB-SEM 8 nm 等方 + Google FFN 分割 | **CC-BY 4.0**、neuPrint / Codex / GCS 一括DL |
| BANC (Bates et al., *Nature* 2026 [二次: 掲載日]) | メス1匹 | 脳 + VNC（視細胞は含まず） | 約16万（図中の proofread 済みは 114,518）[論文] | GridTape TEM 4×4×45 nm | CAVE / Codex（ライセンス要確認） |
| FlyWire / FAFB v783 (Dorkenwald 2024; Schlegel 2024) | メス1匹 | 脳 + 両視葉（VNC 無し） | 139,255 | ssTEM + 市民科学で校正 | CC BY-NC 4.0 [二次]、Codex / CAVE |
| optic-lobe v1.0 (Nern et al. 2025) | オス | 右視葉 | 約5万 | FIB-SEM | neuPrint |
| MANC v1.2 (Takemura / Marin / Cheong et al. 2024) | オス | VNC のみ | 約2.3万 | FIB-SEM | neuPrint |
| FANC (Azevedo et al. 2024) | メス | VNC のみ | 要確認 | ssTEM（MN と筋の対応を X 線で補完） | CAVE |
| hemibrain v1.2.1 (Scheffer 2020) | メス | 中枢脳の約半分 | 約2.5万 | FIB-SEM | neuPrint |

**MaleCNS が土台に向く理由**
- 脳と VNC が**同一個体**。FlyWire（脳）と MANC（VNC）を DN の型マッチングでつなぐ方式だと、接合部が推測になる。MaleCNS では DN 1,314 個・AN 1,846 個がそのまま両側につながっている [実測]。
- **CC-BY**。FlyWire の NC（非営利）条件が無いので、遊びを公開・配布しやすい。
- 感覚ニューロンに受容体型・入力神経、運動ニューロンに出口神経・支配筋の注釈が付いている。
- 型名は FlyWire / hemibrain / MANC と 97.5% 対応づけ済み [論文]。他データで得られた生理の知見を流用できる。

**BANC を併用すべき場面**
- 末梢側の注釈が厚い（感覚 14,410、MN 833、内分泌 193、内臓遠心性など）。全身の「局所フィードバックループ」と、DN/AN による行動モジュール（15 スーパークラスター）の解析がある [論文]。
- **メスの全中枢**として、オス・メス混在の個体群を作るときに必要になる。

### 2.2 MaleCNS v1.0 の中身（実ファイルで確認）

`python tools/malecns_inventory.py --weights` の出力から抜粋 [実測]：

**規模**
- アノテーション付きニューロン 166,700。型付き 164,506、型ラベル 11,751 種（論文値は約 11,700）。
- 結合表：セグメント対 1.52 億行、総シナプス後部 (PSD) 3.12 億。
  - **ニューロン ↔ ニューロン**：エッジ 25,582,938 本、シナプス 124,177,617 個（報道の「1.25 億シナプス」はこの数）。
  - 重み ≥3 のエッジ 1,052 万本、≥5 のエッジ 624 万本。
- スーパークラス：視葉内在 89,403、中枢脳内在 32,164、VNC 内在 13,161、視覚投射 9,201、感覚（VNC 6,370／視葉 6,098／脳 4,868）、AN 1,846、DN 1,314、運動 815、内分泌 94。

**感覚ニューロン（modality 別）**

| class | 数 | 主な内訳 |
|---|---|---|
| visual | 6,091 | R1–R6 3,377、R7/R8 2,714（p/y/d サブタイプ付き） |
| olfactory | 2,639 | ORN 2,635、53 糸球体型 |
| mechanosensory_tactile | 2,558 | 剛毛 2,206（脚剛毛 768、翅剛毛 385、唇弁 163…） |
| mechanosensory | 1,733 | ジョンストン器官 (JO)：聴覚 114 + 風・重力 475（JO-A/B/C/D/E の型付き） |
| mechanosensory_proprioceptive | 1,454 | 弦音器官 425、鐘状感覚子 426、毛板 113、平均棍 205 |
| gustatory | 1,428 | 口器・脚・翅。受容体は putative（ppk23 269、ppk25 257、IR52b 226） |
| unknown_sensory | 1,712 | 腹部 1,143 など |
| hygro / thermo | 66 / 25 | VP 糸球体系 |

**視覚入力に関する重要な欠損**
- 視葉のカラムは左 879・右 892 に六角座標が振られ、L1–L5 はほぼ全カラム揃っている（各型 約 1,775）。
- 一方、**R1–R6 は 3,377 しかない**。期待値の約 10,600（≒885 カラム × 6 × 両眼）の約 1/3 で、特に左が少ない（左 1,112、右 2,265）。論文も、試料端のアーティファクトで一部の R1–6 が再構成されていないと明記している [論文]。R7 にも欠損がある。
- → 視覚入力は**ラミナ単極細胞 L1–L3 にカラム単位で注入する**か、欠けた R1–6 → L 結合を神経重畳の規則で合成する必要がある。

**運動・遠心性**

| subclass | 数 | 意味（出口神経から推定） | 筋名付きの例 |
|---|---|---|---|
| fl / ml / hl | 135 / 116 / 130 | 前・中・後脚 | Ti flexor 37、Acc. ti flexor 47、Tr flexor 35、Fe reductor 20、Ta levator/depressor、sternal rotator… |
| wm | 67 | 翅 | DLMn・DVMn（パワー筋）、b1–b3・i1–i2・iii1/3・hg1–4・tp1–2・ps1–2（舵取り筋）、TTMn |
| hm | 16 | 平均棍 | hi1/hi2/hiii2、hDVM |
| nm | 44 | 首（CvN, DProN） | |
| pm | 67 | 口吻・咽頭（PhN, MxLbN） | MN1–MN13（MN9 を含む） |
| am / rm / xm | 13 / 7 / 6 | 触角 (AN) / 眼 (ON) / 不明 | |
| ad | 214 | 腹部 | |

- **815 個中 400 個は支配筋の名前で型付けされている** [実測]。残りは MNhl64 のような系統名。
- 内分泌：IPC（インスリン産生細胞）16、DH44、ITP、Hugin、LK、CAPA など 94 個。ホルモン状態をモデル化する足場はある。

**伝達物質（コンセンサス予測）**
- ACh 103,720、Glu 29,302、GABA 22,069、ヒスタミン 7,891、DA 392、OA 101、5-HT 48、unclear 2,999 [実測]。
- 手法は T-bar 周辺画像を ResNet50 で分類。オクトパミンとセロトニンは検証データが乏しく、慎重に扱われている [論文]。
- 参考：同系統の手法で FlyWire 全脳を評価した Eckstein et al. (2024) の正答率は、シナプス単位 87%・ニューロン単位 94% [論文]。

**再構成の品質**
- シナプス検出の precision/recall は 0.82/0.81。ニューロピル内の完了率は前シナプス 94%、後シナプス 42%。**両側とも校正済みの結合は 40.1%** [論文]。
  - → 重みは「本物の結合のうち約 4 割を数えた値」。ニューロン間の相対比は保たれていると仮定することになる。
- 核 141,780 個の 98.9% が校正済みニューロンに属する [論文]。
- 個体: 5日齢オス（Canton S G1 × w1118）[論文]。
- ORN の完全再構成は約 70%。JO と一部 MN には分割上の問題がある [論文]。

**性差**
- male-specific 1,258、sexually dimorphic 771（+ potentially 339）。fru のみ 4,600、dsx のみ 154、共発現 258（いずれも high+low）[実測]。
- 求愛歌 (pulse/sine song) の回路がそのまま入っている。→ §9 の「本物の通信チャネルとの比較」に使える。

### 2.3 コネクトームに**入っていない**もの

| 欠けているもの | 影響 | 回避策 |
|---|---|---|
| **ギャップ結合（電気シナプス）** | EM の解像度では見えない。例：**巨大線維 (GF, DNp01) → TTMn の逃避経路は shakB ギャップ結合が主役**で、化学シナプスだけのモデルでは逃避が遅れたり出なかったりする | 文献既知のギャップ結合を手で追加する |
| 神経ペプチド・モノアミンの体積伝達 | 空腹・覚醒・求愛状態などの「内部状態」が表せない | 内分泌・修飾ニューロンを「グローバル状態変数」として扱う |
| 受容体型・時定数・Glu の符号 | Glu は GluCl を介すと抑制だが、場所によっては興奮 | 型ごとに文献・トランスクリプトーム（Fly Cell Atlas）で上書き |
| シナプス強度 | 数 ≠ 生理的強度 | flyvis 式に「型ごとのパラメータを学習」 |
| 非スパイクニューロン | 視葉の多くは段階電位。LIF は不適 | 視葉は rate / graded モデルにする |
| 可塑性・学習状態 | キノコ体学習が初期状態で凍結 | DAN → KC-MBON の可塑則を後付け |
| 個体差 | オス1匹分の n=1 | BANC / FlyWire / MANC との型レベル比較で頑健性を見る |

---

## 3. 脳ダイナミクスのモデル

### 3.1 既存モデル

| モデル | 範囲 | 方式 | 検証されたこと | 実装 |
|---|---|---|---|---|
| **Shiu et al. 2024** (*Nature*) | FlyWire 脳（約 12.7–14 万） | LIF。重みはシナプス数 × 符号（NT 予測）。ギャップ結合・修飾なし | 糖/水の GRN 刺激 → 摂食開始ニューロン・MN9 発火を予測し、光遺伝学で確認。機械感覚 → 触角グルーミング回路 | Brian2（MIT）: [philshiu/Drosophila_brain_model](https://github.com/philshiu/Drosophila_brain_model) |
| Eon fork | 同上 v783 | 同上。Brian2 / Brian2CUDA / PyTorch / NEST GPU / GeNN の多バックエンド比較 | 実装間の一致検証 | GPL-2.0+: [eonsystemspbc/fly-brain](https://github.com/eonsystemspbc/fly-brain) |
| **flyvis** (Lappalainen et al. 2024, *Nature*) | 視葉の運動視経路 64 型、六角格子カラム | 配線は固定、単一ニューロン・シナプスのパラメータをオプティカルフロー課題で学習（DMN） | 26 研究の実測応答（T4/T5 の方向選択性など）を予測 | PyTorch（MIT）: [TuragaLab/flyvis](https://github.com/TuragaLab/flyvis) |
| Pugliese et al. 2025–26 (bioRxiv) | VNC（4 データセット） | 生物学的にあり得るパラメータを総当たりで振ったシミュレーション | **DNg100 → 6 脚のリズム**。抑制 1 + 興奮 2 の 3 ニューロン CPG を同定。DNb08 の予測を光遺伝学で確認 | [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.09.12.675944) |
| Li et al. 2026 (bioRxiv) | FlyWire 全脳 | 頭部固定ハエの自発活動カルシウムイメージングに全脳モデルをフィット | 抑制ハブ少数が自発活動の維持に必要十分、と予測 | [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.08.21.745055v1) |
| Chen & Xi 2025 (bioRxiv) | FlyWire v783 | LIF | 逃避時の摂食抑制が分散回路で実現される | [bioRxiv](https://www.biorxiv.org/content/10.64898/2025.12.14.694122v1) |
| Jin et al. 2026 (arXiv) | 全脳コネクトーム | 配線を**グラフ NN の構造**として使い、強化学習で歩行制御を学習 | 生物の動態ではなく「構造の事前分布」として使う路線 | [arXiv:2602.17997](https://arxiv.org/abs/2602.17997) |

**注意喚起の論文**：connectome-constrained モデルを課題誤差で選ぶと、実測のチューニングと一致しないことがある（[bioRxiv 2026.04.18.717873](https://www.biorxiv.org/content/10.64898/2026.04.18.717873v1)）。課題が解けることと、生物として正しいことは別物。

### 3.2 シミュレータ・ライブラリ

| 用途 | ライブラリ | メモ |
|---|---|---|
| データ取得 | `neuprint-python`、`caveclient`（BANC/FlyWire）、`fafbseg`、`navis`、`flybrains`（JRC2018 テンプレート変換）、R の natverse | MaleCNS は Feather の一括 DL が一番速い |
| 閲覧 | neuPrint、Codex（`?dataset=mc` で MaleCNS、`/banc`）、Clio、Neuroglancer、Male CNS Cell Type Explorer、NeuronBridge（EM ↔ GAL4 系統の対応） | |
| スパイキング | Brian2 / Brian2CUDA、GeNN、NEST GPU | 検証用の基準実装として |
| 微分可能 / GPU | PyTorch、JAX | 本プロジェクトの本命。多個体のバッチ化、パラメータ学習 |
| 多コンパートメント | Jaxley（*Nature Methods* 2025、JAX で微分可能）、NEURON | MaleCNS の SWC スケルトンがそのまま使える。特定回路の精密化に |

### 3.3 脳モデルの正しさ（まとめ）
- **検証済みの回路**：味覚 → 摂食、触角機械感覚 → グルーミング、運動視 (T4/T5)、歩行リズム (VNC)。
- **未検証・不明**：中心複合体のリングアトラクタ（方位）が全脳 LIF の中で自然に維持されるか、学習、内部状態、求愛の意思決定。
- Eon 自身が、全脳モデルの内部動態はリングアトラクタや CPG のような生物学的特徴で検証されていないと明記している。

---

## 4. 身体モデル

### 4.1 NeuroMechFly / FlyGym（EPFL Ramdya lab）
- **身体**：成虫メスの micro-CT から作成 [論文/公式]。脚は 1 本 7 自由度の位置制御（全 42 自由度）が標準。脚先の接着をスタンス/スイングで切り替えられる。
- **感覚**：
  - 六角格子の複眼（NMF v2 論文の設定で片眼 721 個眼。flyvis と同じ格子）
  - 触角・小顎肢での匂い強度（複雑なプルーム環境あり）
  - 関節角・力・接触力・体節ごとの接触力
- **FlyGym 2.x**
  - 2.0.0（2026-04-02）で完全に書き直された（後方互換なし）。
  - 2.1.0 で dm_control を捨て MuJoCo ネイティブの MjSpec へ移行した。追加機能：
    - MuJoCo Warp 3.9 による GPU 実行（約 300 倍速の主張）
    - **FlyMimic**（筋駆動の筋骨格モデルによる模倣学習、実験的）
    - **flybody モデルの取り込み**（実験的）
    - ブラウザビューア
  - 視覚は CPU のみ対応で、GPU 版の視覚は保留中 [公式 changelog]。
- **脳との接続**：脳 ↔ VNC の階層。下行 (DN) / 上行 (AN) のインターフェースを想定した設計。v2 論文では flyvis を載せて他個体を追従させた。
- リンク: [neuromechfly.org](https://neuromechfly.org/)、[NeLy-EPFL/flygym](https://github.com/NeLy-EPFL/flygym)

### 4.2 flybody（Janelia + Google DeepMind、Vaxenburg et al. *Nature* 2025）
- **身体**：成虫メスの共焦点画像（Congo Red で外骨格と腱を染色）から作成。剛体 67、関節 66、**102 自由度**、0.983 mg [論文]。
- **アクチュエータ**：翅はトルク、他は位置。腹部と跗節は腱で連動。脚先と上唇に接着。**筋は無い**。
- **翅と飛行**：準定常の流体力モデル（付加質量・粘性抗力・Kutta 揚力など）。強化学習（DMPO）の模倣学習で、実際の飛行・歩行軌跡を追従する。
  - 飛行の参照データは *D. hydei*。翅打ち周波数の変化は実物より小さい。
- **感覚**：
  - 視覚：32×32 グレースケール、150° 視野
  - 前庭：速度・加速度・ジャイロ・重力
  - 固有感覚：関節角・速度
  - いずれも理想化されたセンサ
- 論文自身がコネクトームとの統合を将来課題としている。Apache-2.0: [TuragaLab/flybody](https://github.com/TuragaLab/flybody)

### 4.3 筋骨格・翅ヒンジ（MN → 力の橋）
- **脚の筋骨格モデル**（arXiv:2509.06426、ICLR 2026 [二次]）：
  - 前脚の 15 筋腱ユニット（胸 7・基節 6・腿節 2。19 筋群中 12 をカバー）を Hill 型で表現し、OpenSim と MuJoCo に実装。
  - 中・後脚は形状再構成のみ。
  - 論文中の概数は「脚の約 19 筋を約 69 MN が制御」。ただし MN と筋の具体的な対応表は含まれない。
- **翅ヒンジ**（Melis, Siwanowicz & Dickinson, *Nature* 2024）：
  - **12 本の舵取り筋（各 1 MN）の活動 → 翅の運動** を CNN で予測するモデル。物理シミュレーションで自由飛行に似た機動を再現した。
  - MaleCNS には、この 12 本のうち 11 本（b1–3、i1–2、iii1、iii3、hg1–4）に対応する MN 型が左右 1 個ずつある [実測]。残る iii4 は名前付きの型が無く、系統名（MNwm35/36 など）に含まれている可能性がある。ほかに tp1–2、ps1–2、tpn もある。**MN 活動から体を直接動かせる見込みが最も高い部位**。

### 4.4 性の不一致
MaleCNS はオス、身体モデルは 2 つともメス由来。オスは小さく、前脚に性櫛があり、腹部の形も違う。当面はスケールを合わせるだけで許容する。求愛歌（翅振動）をやるなら翅まわりの差を確認する必要がある。

---

## 5. 感覚器 ↔ 感覚ニューロンの対応状況

| モダリティ | MaleCNS 側 [実測] | 身体シミュ側 | 変換モデル | 評価 |
|---|---|---|---|---|
| 視覚 | R1–6 は約 1/3 のみ、R7/R8 2,714、L1–L5 全カラム、六角座標 L879/R892 | NMF 六角格子複眼、flybody 32×32 カメラ | flyvis（DMN）。Zhao et al. 2025（*Nature*）の micro-CT による個眼 ↔ カラム 1 対 1 写像 | 構造 ◎、入力 △（R1–6 欠損を補う必要あり） |
| 嗅覚 | ORN 2,635、53 型 | NMF の触角・小顎肢センサ、匂いプルーム | 受容体 × 匂い物質の応答表（DoOR 等） | 配線 ◎、応答 ○、プルーム物理 △ |
| 味覚 | GRN 1,428（受容体型は一部 putative） | 跗節・口吻の接触 | Shiu 2024 の糖・苦味・水 GRN セット | ○ |
| 触覚（剛毛） | 2,206 | 体節ごとの接触力 | 剛毛の体表位置 ↔ ニューロンの対応が粗い | △ |
| 固有感覚 | 弦音器官 425、鐘状感覚子 426、毛板 113 | 関節角・速度・トルク | 弦音器官 (FeCO) の claw/hook/club 符号化モデルは文献にあり、実装が必要 | △ |
| JO（聴覚・風・重力） | 114 + 475 | 触角は受動。音場が無い | 要実装 | △（求愛歌に必須） |
| 平均棍・前庭 | 205 | flybody のジャイロ | 要実装 | △ |
| 温度・湿度 | 25 / 66 | 無し | 環境側も要実装 | × |
| 内受容（腹部など） | 1,143 ほか | 無し | 無し | × |

---

## 6. 運動ニューロン ↔ 身体の対応状況

| 部位 | MN [実測] | 筋同定 | 身体側の受け口 | ギャップ |
|---|---|---|---|---|
| 脚 ×6 | 381 | 名前付きで約半分（Ti/Tr flexor・extensor など） | NMF の 7 DoF/脚（位置制御）、前脚の筋骨格（実験的） | MN プール → 筋 → トルクの変換を作る必要あり |
| 翅 | 67 | パワー筋と舵取り筋を同定済み | flybody の翅トルクと流体力 | Melis 2024 のヒンジモデルで舵取り筋 → 翅運動をつなげる |
| 平均棍 | 16 | 一部 | 無し | |
| 首 | 44 | 系統名中心 | 頭部関節 | |
| 口吻・咽頭 | 67 | MN1–13（MN9 など） | 口吻関節（NMF）、上唇の接着（flybody） | 摂食はここで閉じられる |
| 触角 | 13 | | 触角関節 | |
| 腹部 | 214 | | flybody の腹部（腱で連動） | |
| 内分泌 | 94 | IPC、DH44 など | 無し | ホルモンを状態変数として扱う |

---

## 7. 既存の「コネクトーム × 身体／ゲーム」統合例

| 事例 | 脳 | 体・環境 | インターフェース | 学べる点 |
|---|---|---|---|---|
| **Eon Systems** embodied fly（2026-03） | Shiu LIF（FlyWire 脳）+ flyvis | NeuroMechFly（MuJoCo） | 感覚：GRN・ORN・JO に「イベント」を注入。運動：DNa01/02（旋回）、oDN1（前進）、MN9（摂食）、触角グルーミング DN → 下位コントローラ。**15 ms ごとに脳と体を同期**。すべて手設計 | VNC なし、DN は少数。視覚は「飾り」に近いと本人たちが認めている。逃避は脳内では起きるが体では未実装。コード非公開 [公式] |
| NeuroMechFly v2 論文 | flyvis + 手設計 | NMF | 旋回バイアスを DN 経由で渡す | |
| webgpu-fly | FlyWire LIF + MANC LIF（WebGPU） | flybody（MuJoCo を wasm で） | 不明 | ブラウザでも全脳 + VNC が回る [二次] |
| FlyBrain（Minecraft mod） | MaleCNS の飛行関連サブネット | Minecraft | 不明 | [二次] |
| Connectome Fly（Minecraft mod） | BANC | Minecraft | 不明 | [二次] |
| DOOMFLY | MaleCNS LIF | Doom | 画面 → 感覚ニューロン、被弾 → PPL101（DAN）に報酬刺激 | 本人も「本物の脳ではない」と認めている [二次] |
| Jin et al. 2026 | 全脳グラフを RL の方策として使う | 物理ハエ | RL で学習 | 動態ではなく構造の事前分布として使う |

**共通の近道**：どの事例も VNC と MN を通さず、「DN → 手設計の運動プリミティブ」で体を動かしている。MaleCNS は脳 → DN → VNC 前運動回路 → MN が同一個体でつながっているので、**MN レベルでの身体駆動**はまだ誰もやっていない領域になる。

---

## 8. 計算コストの見積り（実データから）

| 量 | 値 |
|---|---|
| ニューロン数 N | 166,700 |
| エッジ数 E | 25.58M（平均出次数 約 153） |
| CSR（fp32 重み + int32 index） | **約 0.20 GB / 匹**。全個体が同じ配線なら重みは共有でき、個体ごとに増えるのは状態ベクトル（N × 数 float ≒ 数 MB）だけ |

- **イベント駆動 LIF**：シナプスイベント/秒 ≈ N × 平均発火率 × 153 ≈ 2.6×10⁷ × r。r = 5 Hz で 1.3×10⁸ イベント/秒（シミュレーション内時間）。GPU なら数十〜百匹規模で実時間に近づく見込み（要ベンチマーク）。
- **クロック駆動のレートモデル**（dt = 1 ms）：SpMV が 2.6×10¹⁰ MAC/秒・匹。複数個体を SpMM にまとめれば、1 GPU で数匹〜十数匹は実時間に届く見込み。
- **律速は脳ではなく物理**：flybody は CPU 1 コアで飛行 10 ms に 0.42 秒かかる [論文]。多個体なら MuJoCo Warp（FlyGym 2.x）でのバッチ化が前提。
- 脳 ↔ 体の同期：Eon は 15 ms。逃避 (GF) や翅の制御（約 200 Hz）には粗すぎるので、1 ms 程度を目標にする。

---

## 9. プシュケー・チップへの含意

### 9.1 チップの抽象モデル（叩き台）
```
read(fly_i)  : 部位 S_read のニューロン集合の活動 → 符号化（スパイク列 / 発火率ベクトル / 低次元射影）
channel      : 送信強度 P、距離減衰 1/r^α、雑音、遅延、帯域（チャネル数）、指向性、
               個体ごとの役割（送信のみ / 受信のみ / 送受信 / 中継）
write(fly_j) : 受信信号 → 部位 S_write へ電流注入 / Poisson 駆動 / ゲイン変調
               （上書き・加算・ゲート付き）
```

### 9.2 読み書きの候補部位（MaleCNS のニューロン数は [実測]）

| 部位 | 数 | 何を表すか | 生理的な裏付け | 実験アイデア |
|---|---|---|---|---|
| EPG（中心複合体） | 50 | 方位（リングアトラクタのバンプ） | 強い（ただし全脳 LIF でのバンプ維持は未検証） | **方位共有 → 群れの整列**。蔵本モデル的な同期転移を、結合強度と距離減衰で相図にする |
| LPLC2 → GF (DNp01) | 185 → 2 | 迫ってくる物体 → 逃避 | 強い（※GF → TTMn はギャップ結合なので手で追加する） | **恐怖の伝染**。送信専用の「見張り役」個体を置いて群れの生存率を比べる |
| PPL1 DAN（PPL101 など）/ MBON | 各数個 | 価値・罰 | 強い（キノコ体学習） | **テレパシー条件づけ**。A の罰が B の記憶になる |
| pC1 / P1 系、求愛歌回路 | 156 + | 覚醒・求愛 | 中程度 | 本物の通信（求愛歌・cVA フェロモン）と人工チャネルの情報量を比べる |
| IPC（インスリン） | 16 | 栄養状態 | 中程度（ホルモン系は未モデル） | **空腹の共有**。採餌の分業が出るか |
| DN 群 | 1,314 | 運動意図 | 強い（一部の型） | **操り人形**。A の DN を B に書く。B の VNC が A の意図をどう解釈するか |
| L1–L3（視覚入力） | 各約 1,775 | 網膜像 | 強い | **視界の共有**。A の見たものを B の視葉に書き込む |

**どこに書き込むと何が起こるかの予測可能性**は、検証済みの回路（視覚・逃避・摂食・グルーミング・歩行 CPG）に近いほど高い。中心複合体・キノコ体・求愛系はモデルの正しさ自体が実験対象になる。むしろ「テレパシーが社会行動を生むか」を、モデルの検証（リングアトラクタが自律的に保たれるか等）と同時に進める題材として面白い。

---

## 10. ロードマップ案

| 段階 | 内容 | 合格基準 |
|---|---|---|
| **P0 データ基盤** | MaleCNS → 疎行列 + 型・部位インデックス（`tools/malecns_inventory.py` を発展させる）。BANC・FlyWire との型対応表 | 任意の型集合を名前で引ける |
| **P1 脳単体（GPU）** | PyTorch/JAX で全中枢 LIF（視葉は graded/rate との混成も可）。ギャップ結合の手動追加機構 | 検証スイート：① 糖 GRN → MN9（Shiu を MaleCNS で再現）、② JO/剛毛 → 触角グルーミング DN、③ DNg100 → 6 脚交互リズム（Pugliese）、④ LPLC2 → GF → TTMn、⑤ EPG バンプの保持（未知、結果がそのまま知見になる） |
| **P2 身体化 L1** | FlyGym 2.x と閉ループ。DN → 既存コントローラ（Eon 方式）。感覚：接触 → GRN/剛毛、視覚 → L1–L3、匂い → ORN | 採餌・グルーミング・逃避が閉ループで出る |
| **P2' 身体化 L2（本命の挑戦）** | **MN 活動 → 筋 → 関節トルク**。脚は筋名付き MN プール + 筋骨格モデル + 学習デコーダ。翅は舵取り筋 MN → Melis ヒンジモデル | DN を直接触らず、VNC が出す MN 活動だけで歩く・方向転換する |
| **P3 多個体 + チップ** | 共有重みで N 匹をバッチ実行。チップの API とチャネルモデル | §9.2 の実験を再現できる形で回せる |
| **P4 同定** | flyvis 方式で、型ごとのパラメータを行動・既知の生理に合わせて学習。チップ ON/OFF の比較を生物的に妥当な範囲で行う | 学習前後で検証スイートが壊れない |

---

## 11. 未確認事項（次に調べる）
- BANC のライセンスと、Nature 版の最終的な数値。
- FlyWire の CC BY-NC が本プロジェクトの公開形態に与える影響（MaleCNS 中心なら回避可能）。
- 脚 MN の筋同定（MANC/FANC 由来）の完全な対応表。MaleCNS の残り約 400 MN の扱い。
- NMF の複眼（721 個眼の格子）と MaleCNS の実カラム（約 880–890 個）の写像方法。
- 全脳イメージング（行動中）のうち、細胞型レベルで検証に使える公開データ。

---

## 参考文献・リンク

**コネクトーム**
- MaleCNS: [公式](https://male-cns.janelia.org/) / [Download](https://male-cns.janelia.org/download/) / [bioRxiv 2025.10.09.680999](https://www.biorxiv.org/content/10.1101/2025.10.09.680999v2) / [Cell 2026](https://www.cell.com/cell/fulltext/S0092-8674(26)00942-6) / [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12636603/) / [Cell Type Explorer](https://github.com/reiserlab/celltype-explorer-drosophila-male-cns) / [Google Research](https://sites.research.google/neural-mapping/)
- BANC: [bioRxiv 2025.07.31.667571](https://www.biorxiv.org/content/10.1101/2025.07.31.667571) / [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12324551/) / [FlyWire blog](https://blog.flywire.ai/?p=1823)
- MANC: [Janelia](https://www.janelia.org/project-team/flyem/manc-connectome) / [Cheong et al. eLife](https://elifesciences.org/articles/96084)
- データセット一覧: [Virtual Fly Brain EM data](https://flybrain-ndb.virtualflybrain.org/docs/data/em/)
- 伝達物質予測: [Eckstein et al. 2024 Cell (PMC)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11106717/)

**脳モデル**
- Shiu et al. 2024: [doi:10.1038/s41586-024-07763-9](https://doi.org/10.1038/s41586-024-07763-9) / [code](https://github.com/philshiu/Drosophila_brain_model) / [Eon fork](https://github.com/eonsystemspbc/fly-brain)
- flyvis: [Lappalainen et al. 2024](https://doi.org/10.1038/s41586-024-07939-3) / [code](https://github.com/TuragaLab/flyvis)
- VNC CPG: [Pugliese et al.](https://www.biorxiv.org/content/10.1101/2025.09.12.675944) / [The Transmitter](https://www.thetransmitter.org/systems-neuroscience/long-sought-walking-circuit-found-in-fruit-flies/)
- [Li et al. 2026 自発活動フィット](https://www.biorxiv.org/content/10.64898/2026.08.21.745055v1) / [Chen & Xi 2025](https://www.biorxiv.org/content/10.64898/2025.12.14.694122v1) / [再現性の注意](https://www.biorxiv.org/content/10.64898/2026.04.18.717873v1) / [Jin et al. 2026](https://arxiv.org/abs/2602.17997)
- Jaxley: [doi:10.1038/s41592-025-02895-w](https://doi.org/10.1038/s41592-025-02895-w)
- 総説: [State of Brain Emulation Report 2025](https://arxiv.org/abs/2510.15745) / [The embodied brain (Wang-Chen & Ramdya 2026)](https://arxiv.org/abs/2601.08056)

**身体**
- NeuroMechFly: [公式](https://neuromechfly.org/) / [changelog](https://neuromechfly.org/changelog/) / [flygym](https://github.com/NeLy-EPFL/flygym) / [NMF v2 Nature Methods](https://doi.org/10.1038/s41592-024-02497-y)
- flybody: [code](https://github.com/TuragaLab/flybody) / [Vaxenburg et al. 2025 (PMC)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12310536/)
- 筋骨格の脚: [arXiv:2509.06426](https://arxiv.org/abs/2509.06426)
- 翅ヒンジ: [Melis et al. 2024](https://doi.org/10.1038/s41586-024-07293-4)
- 眼の地図: [Zhao et al. 2025](https://doi.org/10.1038/s41586-025-09276-5)

**統合例**
- Eon: [How the Eon Team Produced a Virtual Embodied Fly](https://eon.systems/updates/embodied-brain-emulation) / [More flies are getting uploaded](https://eon.systems/updates/more-flies-are-getting-uploaded)
- [webgpu-fly](https://abgunaydin-webgpu-fly.static.hf.space/index.html) / [Connectome Fly (Modrinth)](https://modrinth.com/project/ZZxOe21T) / [FlyBrain (CurseForge)](https://www.curseforge.com/minecraft/mc-mods/flybrain) / [Doom 報道 (TechSpot)](https://www.techspot.com/news/113780-google-mapped-fly-nervous-system-developers-using-play.html)
