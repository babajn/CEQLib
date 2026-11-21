# EQCLib と CHIMES の比較

| 比較 | 計算条件 | 作図ノートブック | 出力 |
|---|---|---|---|
| [無照射の単一セル](optically-thin/README.md) | 初期温度と組成をそろえた4条件の非平衡計算 | [PlotEQCLibChimes.ipynb](../../plot/PlotEQCLibChimes.ipynb) | `test/out/chimes/optically-thin/` |
| [FUV・遮蔽ありのグリッド](shielded/README.md) | 現行の平衡モードを中心に、固定温度の平衡、時間更新、等積冷却履歴を比較 | [PlotShieldedChimes.ipynb](../../plot/PlotShieldedChimes.ipynb) | `test/out/chimes/shielded/` |
| [低温ガス](shielded/README.md#低温ガスの現行モデル比較) | 同じ元素構成での H₂ 分率、熱収支、平衡温度、冷却履歴と更新時間 | [PlotH2Comparison.ipynb](../../plot/PlotH2Comparison.ipynb) | `test/out/chimes/shielded/{h2_current,thermal_balance_current}/` |

各比較では外部の CHIMES ソースと反応データを `CHIMES_DIR` で参照します。既定値は EQCLib と同じ階層にある `Chimes/` です。計算条件、ビルド、再実行方法はそれぞれの README にあります。生成した CSV と図の出力先は Git の管理対象外です。

ダスト熱過程・収束・CIE 冷却表・H₂ 処方の開発時の比較も、当時の図と再計算手順を参照できます。 [開発時の検証](shielded/README.md#開発時の検証)
