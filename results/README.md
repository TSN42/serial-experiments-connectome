# results/

実験スクリプトの数値出力。大きい中間ファイル（`.npy`, `.parquet`）とログは git 管理外で、スクリプトを再実行すると再生成される。

## structure/（`experiments/01_structure.py`）
| ファイル | 中身 |
|---|---|
| `summary.json` | 規模、伝達物質の内訳と出どころ、抑制入力割合、左右差の集計、結合行列の先頭モード、影響度プロファイルの相関 |
| `superclass_by_nt.csv` | スーパークラス × 伝達物質の細胞数 |
| `input_ei_by_superclass.csv` | スーパークラスごとの抑制入力割合・符号なし入力割合・入力シナプス数（中央値） |
| `lr_asymmetry.csv` | 型ごとの左右の 1 細胞あたり入出力シナプス数と log2 比（再構成の欠けの検出） |
| `traversal_sensory_to_effector.csv` | 感覚群 → 効果器群の到達ステップ（確率的走査、中央値）と到達率 |
| `influence_sensory_to_effector.csv` | 線形影響度（符号なし/付き × λ = 0.99/0.5、種 1 細胞あたり）の効果器群平均 |
| `effector_input_sources.csv` | 効果器群が受ける入力の、送り手スーパークラス別の割合 |

## lif/（`experiments/02_lif.py`, `04_feeding_dose.py`, `lif_figures.py`）
| ファイル | 中身 |
|---|---|
| `sweep.csv` | Shiu 型の重み（α = 0）で倍率を振った回路テストの全指標（縦持ち：w_scale, test, condition, metric, value） |
| `sweep_alpha0.5.csv` | 細胞の大きさ補正 α = 0.5 |
| `sweep_adapt1.csv` | スパイク頻度順応 b = 1 mV |
| `verdicts*.csv` | 上の各スイープの合否表（基準は `lif_figures.verdicts`） |
| `feeding_dose.csv` | 摂食回路の用量反応（糖の頻度、苦味の同時刺激、水） |
| `influence_vs_lif_auc.json` | 糖刺激時、線形影響度が LIF の活動細胞を見分ける AUC |
| `meta*.json` | 刺激・読み出しに使った細胞数、ギャップ結合の対 |
