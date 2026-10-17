# serial-experiments-connectome

ショウジョウバエのコネクトーム（主に **MaleCNS v1.0**）から全中枢モデルを組み、物理シミュレーション上の身体に載せる。最終的には、脳の一部を読み出して周囲へ送信し、受信した信号を別の部位へ書き込む「プシュケー・チップ」を載せた個体群を仮想空間に放つ。

## ドキュメント
- [docs/01_survey_landscape.md](docs/01_survey_landscape.md) — 調査報告 #01。データ・モデル・ライブラリと、それぞれの解剖学的な正しさ
- [docs/02_malecns_deep_dive.md](docs/02_malecns_deep_dive.md) — 調査報告 #02。MaleCNS の構造解析と全中枢 LIF の検証、公開後に出てきた道具、正しい扱い方、LIF か HH か

## コード
```
secon/            ライブラリ（data: 取得と読み込み / structure: 構造解析 / lif: 全中枢 LIF / groups: 感覚・効果器の定義）
experiments/      01_structure → 02_lif → 03_figures の順に実行。04_feeding_dose は摂食回路の用量反応
tools/            malecns_inventory.py（データの棚卸し）
tests/            LIF の単体テスト
results/          実験の数値出力（CSV/JSON）。大きい .npy/.parquet は git 管理外
```

## 動かし方
```bash
pip install -r requirements.txt
python -m pytest -q tests                       # LIF の単体テスト
python experiments/01_structure.py              # 構造解析（初回は MaleCNS を ~1.2 GB ダウンロード、約 15 分）
python experiments/02_lif.py                    # LIF 検証スイート（Shiu 型の重み、11 倍率 × 6 テスト、3 並列で約 15 分）
python experiments/02_lif.py --alpha 0.5 --scales 0.44 0.6 0.8 1.0 1.3 1.6 2.0 2.6   # 細胞の大きさ補正あり
python experiments/04_feeding_dose.py           # 摂食回路の用量反応
python experiments/03_figures.py                # docs/figures/ の図
```
データは `data/` 以下にキャッシュされる（git 管理外）。MaleCNS は CC-BY 4.0。
