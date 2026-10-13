# serial-experiments-connectome

ショウジョウバエのコネクトーム（主に **MaleCNS v1.0**）から全中枢モデルを組み、物理シミュレーション上の身体に載せる。最終的には、脳の一部を読み出して周囲へ送信し、受信した信号を別の部位へ書き込む「プシュケー・チップ」を載せた個体群を仮想空間に放つ。

## ドキュメント
- [docs/01_survey_landscape.md](docs/01_survey_landscape.md) — 調査報告 #01。データ・モデル・ライブラリと、それぞれの解剖学的な正しさ

## ツール
```bash
pip install -r requirements.txt
python tools/malecns_inventory.py            # MaleCNS v1.0 の注釈を棚卸し (~57 MB DL)
python tools/malecns_inventory.py --weights  # 全結合表も解析 (~1.1 GB DL)
```
データは `data/` 以下にキャッシュされる（git 管理外）。MaleCNS は CC-BY 4.0。
