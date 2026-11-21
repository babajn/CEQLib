# EQCLib のテストと比較計算

コマンドは EQCLib のルートで実行します。初回は `make -C test notebook-data` で、基本テストと H II 領域、H₂・CO 化学、非平衡化学、Wolfire (2003) の再計算をまとめて実行できます。作図ノートブックは `test/out/` の計算結果を読み込みます。Cloudy と Wolfire (1995) の比較には保存済みの参照データもあります。 [ビルドと JupyterLab の手順](../README.md#テスト結果を-jupyter-で確認する)

| 検証・比較 | コードと参照データ | 実行コマンド | 作図ノートブック |
| --- | --- | --- | --- |
| 平衡・水素化学・熱化学モード | [unit/](unit/) | `make -C test check`、`make -C test cell-check` | |
| 熱化学更新の計算時間 | [benchmark/](benchmark/) | `make -C test benchmark` | |
| H II 領域 | [hii/](hii/) | `make -C test hii`、`make -C test hii-benchmark` | [平衡・冷却・反応率](../plot/PlotHIIRegion.ipynb) |
| Cloudy C23.01 との H II 領域比較 | [hii/cloudy/](hii/cloudy/) | [実行手順](../README.md#cloudy-との比較) | [H II 領域ノートブックの Cloudy 比較](../plot/PlotHIIRegion.ipynb) |
| H₂・CO 化学と反応設定 | [co/](co/) | `make -C test co`、`make -C test molecular-settings-benchmark`、[設定比較の計算ノートブック](co/RunMolecularSettingsBenchmark.ipynb) | [H₂・CO 化学](../plot/PlotCOChemistry.ipynb)、[設定比較](../plot/PlotMolecularSettings.ipynb) |
| 非平衡化学と高温ガスの冷却 | [noneq/](noneq/) | `make -C test noneq`、`make -C test noneq-benchmark`、`make -C test noneq-hot` | [非平衡化学](../plot/PlotNonEquilibrium.ipynb) |
| Wolfire et al. (1995) | [wolfire1995/](wolfire1995/) | [実行手順](../README.md#wolfire-et-al-1995-の検証計算と描画) | [熱平衡と熱収支](../plot/PlotWolfire1995.ipynb) |
| Wolfire et al. (2003) | [wolfire2003/](wolfire2003/) | `make -C test/wolfire2003 run` | [熱平衡曲線](../plot/PlotWolfire2003.ipynb) |
| CHIMES との比較 | [chimes/](chimes/README.md) | [無照射の単一セル](chimes/optically-thin/README.md)、[遠紫外線・遮蔽ありのグリッド](chimes/shielded/README.md) | [比較別のノートブック一覧](chimes/README.md) |

`common/` にはテスト共通の C ヘッダーがあります。`reference/EQCooling/` には冷却率や反応率の保存済み診断データがあります。H II 領域や Wolfire の比較用データは各テストのフォルダにあります。
