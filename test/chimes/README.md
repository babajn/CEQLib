# EQCLib と CHIMES の比較

CHIMES との比較は外部コードを使う任意の検証です。 [EQCLib の基本テスト](../README.md#基本テスト)

| 比較 | 計算条件 | 作図ノートブック | 出力 |
|---|---|---|---|
| [無照射の単一セル](optically-thin/README.md) | 初期温度と組成をそろえた4条件の非平衡計算 | [PlotEQCLibChimes.ipynb](../../plot/PlotEQCLibChimes.ipynb) | `test/out/chimes/optically-thin/` |
| [FUV・遮蔽ありのグリッド](shielded/README.md) | 現行の平衡モードを中心に、固定温度の平衡、時間更新、等積冷却履歴を比較 | [PlotShieldedChimes.ipynb](../../plot/PlotShieldedChimes.ipynb) | `test/out/chimes/shielded/` |
| [低温ガス](shielded/README.md#低温ガスの現行モデル比較) | 同じ元素構成での H₂ 分率、熱収支、平衡温度、冷却履歴と更新時間 | [PlotH2Comparison.ipynb](../../plot/PlotH2Comparison.ipynb) | `test/out/chimes/shielded/{h2_current,thermal_balance_current}/` |

再計算には CHIMES のソース、Git LFS で取得する HDF5 形式の反応データ、HDF5 と倍精度の SUNDIALS が必要です。`chimes/` と `chimes-data/` を含むディレクトリを `CHIMES_DIR` で指定します。既定値は EQCLib と同じ階層にある `Chimes/` です。ビルドと再実行の手順は各比較の README にあります。ノートブックは表の出力先に生成された CSV を読み込みます。 [CHIMES 公式の取得手順](https://richings.bitbucket.io/chimes/user_guide/GettingStarted/download.html)

ダスト熱過程・収束・CIE 冷却表・H₂ 処方の開発時の比較も、当時の図と再計算手順を参照できます。 [開発時の検証](shielded/README.md#開発時の検証)
