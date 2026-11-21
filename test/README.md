# EQCLib のテストと比較計算

コマンドは EQCLib のルートで実行します。初回は `make -C test notebook-data` で、基本テストと H II 領域、H₂・CO 化学、非平衡化学、Wolfire (2003) の再計算をまとめて実行できます。作図ノートブックは `test/out/` の計算結果を読み込みます。Cloudy と Wolfire (1995) の比較には保存済みの参照データもあります。 [ビルドと JupyterLab の手順](../README.md#テスト結果を-jupyter-で確認する) · [反応・冷却・加熱の物理モデル](../docs/physics.md)

## テスト一覧

| 検証・比較 | コードと参照データ | 実行コマンド | 作図ノートブック |
| --- | --- | --- | --- |
| [平衡・水素化学・熱化学モード](#基本テスト) | [unit/](unit/) | `make -C test check`、`make -C test cell-check` | |
| [熱化学更新の計算時間](#熱化学更新の並列ベンチマーク) | [benchmark/](benchmark/) | `make -C test benchmark` | |
| [H II 領域](#h-ii-領域の検証計算と描画) | [hii/](hii/) | `make -C test hii`、`make -C test hii-benchmark` | [平衡・冷却・反応率](../plot/PlotHIIRegion.ipynb) |
| [Cloudy C23.01 との H II 領域比較](#cloudy-との比較) | [hii/cloudy/](hii/cloudy/) | [実行手順](#cloudy-との比較) | [H II 領域ノートブックの Cloudy 比較](../plot/PlotHIIRegion.ipynb) |
| [H₂・CO 化学と反応設定](#h₂co-化学の検証計算と描画) | [co/](co/) | `make -C test co`、`make -C test molecular-settings-benchmark`、[設定比較の計算ノートブック](co/RunMolecularSettingsBenchmark.ipynb) | [H₂・CO 化学](../plot/PlotCOChemistry.ipynb)、[設定比較](../plot/PlotMolecularSettings.ipynb) |
| [非平衡化学と高温ガスの冷却](#非平衡化学の検証計算と描画) | [noneq/](noneq/) | `make -C test noneq`、`make -C test noneq-benchmark`、`make -C test noneq-hot` | [非平衡化学](../plot/PlotNonEquilibrium.ipynb) |
| [Wolfire et al. (1995)](#wolfire-et-al-1995-の検証計算と描画) | [wolfire1995/](wolfire1995/) | [実行手順](#wolfire-et-al-1995-の検証計算と描画) | [熱平衡と熱収支](../plot/PlotWolfire1995.ipynb) |
| [Wolfire et al. (2003)](#wolfire-et-al-2003-との熱平衡曲線の比較) | [wolfire2003/](wolfire2003/) | `make -C test/wolfire2003 run` | [熱平衡曲線](../plot/PlotWolfire2003.ipynb) |
| [CHIMES との比較](#chimes-との比較) | [chimes/](chimes/README.md) | [無照射の単一セル](chimes/optically-thin/README.md)、[遠紫外線・遮蔽ありのグリッド](chimes/shielded/README.md) | [比較別のノートブック一覧](chimes/README.md) |

`common/` にはテスト共通の C ヘッダーがあります。`reference/EQCooling/` には冷却率や反応率の保存済み診断データがあります。H II 領域や Wolfire の比較用データは各テストのフォルダにあります。

## 基本テスト

`check`、H II 領域、H₂・CO、非平衡などの光電離を含むテストでは、スペクトルのディレクトリを `PEGASE_DIR` に絶対パスで指定します。既定値は手元の `data/PEGASE2/` です。 [必要なスペクトルと設定例](../README.md#pegase2-のスペクトルを指定する)

```sh
make -C test check cell-check
```

`make -C test check` では、一定照射下での既存の化学平衡解への収束を検証します。照射切り替え後は、水素化学、炭素・ヘリウムの電離・再結合、温度、輝線冷却の応答を確認します。

`make -C test check` は、17,600 K で温度と組成を固定して電離光を切り替えたときに、酸素輝線と CIE の冷却率が連続することと、照射下で熱化学更新が進むことを確認します。 [熱化学のテスト](unit/HydrogenChemistryTest.c)

`make -C test cell-check` で平衡・非平衡の両モードの初期組成と極短ステップでの温度を検証できます。

## 熱化学更新の並列ベンチマーク

熱化学更新の並列性能は `OMP_NUM_THREADS=12 make -C test benchmark` で測定できます。このベンチマークは OpenMP を使い、原子ガスと分子ガスそれぞれの325点を、10年・1万年・100万年の時間刻みで3ステップずつ進めます。各ステップの収束数、実行時間、逐次と並列の比較を表示します。

## H II 領域の検証計算と描画

```sh
make -C test hii
```

このコマンドは `PEGASE_DIR` のスペクトルを読み込みます。ノートブックの再描画には、ここで得た CSV を使います。 [スペクトルの指定方法](../README.md#pegase2-のスペクトルを指定する)

[HIIRegionTest.c](hii/HIIRegionTest.c) は、ダストと分子を含まない単一セルを、太陽金属量の PEGASE2 電離スペクトルと case B 再結合で計算します。入射電離光の減衰、非電離 FUV、宇宙線、X 線は含めません。水素核密度 `n_H = 1, 10, 100, 1000 cm⁻³`、ガス金属量 `Z/Z☉ = 0, 0.1, 1`、電離パラメータ `U = NLyn/(n_H c)` の `log10(U) = −4〜−2` の108条件で、等積の熱平衡温度と各元素の電離分率を求めます。さらに `n_H = 100 cm⁻³`、`Z/Z☉ = 1`、`log10(U) = −3` で約3,160〜31,600 K の温度範囲を走査します。

- `test/out/hii/equilibrium.csv` と `thermal_curve.csv` に電離分率と過程別の冷却・加熱率（erg cm⁻³ s⁻¹）を出力します。元素を含まない場合、その電離分率は `nan` になります。
- `test/out/hii/reactions.csv` に H、He、C、O、Si の光電離、電子衝突電離、電子再結合の反応率（cm⁻³ s⁻¹）を出力します。各点で電離と再結合の釣り合いを確認します。
- Jupyter で [PlotHIIRegion.ipynb](../plot/PlotHIIRegion.ipynb) を実行すると、平衡温度・電離分率、冷却・加熱の内訳、反応率を表示します。保存済みの計算結果を描く場合は、最初のコードセルの `run_dir` を変更してください。描画には Python 3、NumPy、Matplotlib 3.7 以降が必要です。

PNG を生成するには `make -C test hii-plot` を実行してください。計算済みの CSV は [plot_hii_region.py](../plot/plot_hii_region.py) でも描画できます。PNG は `plot/out/hii/` に保存されます。

### H II 領域の計算時間

```sh
make -C test hii-benchmark
```

[HIIRegionBenchmark.c](hii/HIIRegionBenchmark.c) は、電離光なし、照射継続中、照射開始直後を同じ密度・元素組成・初期温度で比較します。照射開始直後の粒子には、照射前の化学組成と内部エネルギーを保持したまま、照射継続中と同じ `NLyn` を与えます。ダストと分子なしの108条件で、固定温度の化学平衡・冷却加熱率と、時間刻み10年・1万年・100万年の陰的エネルギー更新を測ります。

`test/out/hii/hii_benchmark.csv` に、各条件群を繰り返した7ブロックの単一スレッド CPU 時間（1呼び出しあたりの平均 µs）を出力します。測定区間はセルの更新処理です。同じ Xeon w5-3423 の1コアで採った測定例も CSV で参照できます。 [探索の変更前](hii/reference/hii_benchmark_before_xeon_w5_3423.csv)、[変更後](hii/reference/hii_benchmark_xeon_w5_3423.csv)

H II 領域で設定を切り替えた測定では、電荷交換を無効にすると更新が約 20〜35% 速くなりました。化学平衡の許容誤差の上限を `1e-5` に緩めた場合は、更新が約 10〜25% 速くなり、温度の変化は 2×10⁻⁵ 程度以下でした。

CO を無効にした `make -C test hii-benchmark` と同じ条件で、電荷交換を無効にして化学平衡の許容誤差を `1e-5` にすると、エネルギー更新1回あたりの時間は既定値の 0.6〜0.8 倍になりました。電荷交換を加える前の版と比べても 0.3〜0.6 倍です。測定値は単一スレッドの例で、計算機によって変わります。

### Cloudy との比較

[Cloudy C23.01](https://data.nublado.org/cloudy_releases/c23/) をビルドし、EQCLib のルートで次を実行します。

```sh
make -C test hii
python3 test/hii/cloudy/run_cloudy_hii.py --cloudy /path/to/c23.01/source/cloudy.exe \
  --pegase-dir /absolute/path/to/pegase-spectra
```

`--pegase-dir` の省略時は環境変数 `PEGASE_DIR` を参照し、設定がなければ手元の `data/PEGASE2/` を使います。

[run_cloudy_hii.py](hii/cloudy/run_cloudy_hii.py) は EQCLib の `equilibrium.csv` から密度、電離パラメータ、元素組成を読み取ります。PEGASE2 の連続スペクトルを EQCLib と同じ水素電離光子数で年齢ごとに規格化し、0〜30 Myr の18時点を等重みで平均します。Cloudy には case B、ダストと分子なし、厚さ 1 cm の1ゾーンを指定します。密度、入射電離光子束、光学的厚さ、熱収支を確認し、冷却の出力から O II と O III の正味の衝突線冷却を集計します。EQCLib の比較量には H II 領域で使われる温度重みを含めています。

各ゾーンの入力、ログ、集計 CSV は `test/out/hii/cloudy/` に出力します。[PlotHIIRegion.ipynb](../plot/PlotHIIRegion.ipynb) は再計算の出力を優先し、出力がなければ [Cloudy C23.01 の保存済み参照結果](hii/cloudy/cloudy_hii_c23.01.csv) を使います。参照結果で描画する場合は EQCLib の `make -C test hii` の出力があれば実行できます。Cloudy のデータ更新と計算方法は [C23 リリース論文](https://doi.org/10.22201/ia.01851101p.2023.59.02.12) と [C23.01 リリース論文](https://doi.org/10.3847/2515-5172/ad0e75) を参照してください。対応する PDF は [reference/](../reference/index.tsv) に保存しています。

H–O 電荷交換を無効にすると、Cloudy との比較における太陽金属量の平衡温度の差は、中央値で 1.4% から 5.3% に広がりました。

## H₂・CO 化学の検証計算と描画

```sh
make -C test co
```

[MolecularTest.c](co/MolecularTest.c) は、H II 領域・H₂・CO・ダスト・宇宙線・軟X線・H–O 電荷交換をすべて有効にした設定で化学平衡を解き、反応収支を確かめます。次の2種類の条件を計算します。

- 密度・温度・電離光子束の格子：`n_H = 1〜10⁵ cm⁻³`（0.25 dex 刻み）、`T = 30, 100, 300, 1000 K`、電離光子なしと `log10(U) = −5, −3` の252条件です。`G0 = 1.7`、全柱密度は `N_H = 2×10²¹ (n_H/100)^{1/2} cm⁻²` です。結果は `test/out/co/equilibrium.csv` に出力します。
- 柱密度走査では電離光子を入れず、`T = 30 K`、`n_H = 100, 1000 cm⁻³` と `G0 = 1.7, 17` の組で `N_H = 10¹⁹〜10²³ cm⁻²` の164点を調べます。各点は一つ前の解から始めます。結果は `test/out/co/column_scan.csv` に出力します。

各平衡状態について、H₂ と CO の生成・消滅、電荷交換を含む H と O の電離・再結合の反応量を、公開されている反応率関数から独立に計算し直します。どれかの相対差が `10⁻⁶` を超えるか、最終的な H₂ 量で CO を解き直したときに CO が `10⁻⁶` 以上変わると、試験は失敗します。CSV には存在比、各反応量（cm⁻³ s⁻¹）、相対差を出力します。

Jupyter で [PlotCOChemistry.ipynb](../plot/PlotCOChemistry.ipynb) を実行すると、反応収支、密度に対する分子の割合、柱密度に対する H/H₂ と C⁺/C/CO の遷移、生成と消滅の反応量を表示します。PNG を生成するには `make -C test co-plot` を実行してください。[plot_co_chemistry.py](../plot/plot_co_chemistry.py) が PNG を `plot/out/co/` に保存します。

### 柱密度走査での CO モデルの確認

柱密度走査では、`A_V ≈ 3〜10` で `CO/C > 0.99` となり、`A_V = 34` では `CO/C` が 0.05 になります。 [CO モデルの適用範囲](../docs/physics.md#co-モデルの適用範囲)

宇宙線による C の電離と He⁺ + CO → C⁺ + O + He を KC08e に加え、走査の反応率を使った試算では、深部の `CO/C` は 0.8〜0.98 に保たれますが、C⁺ を除く反応が遅い C⁺ + H₂ だけなので、`A_V = 3〜50` で C⁺ が炭素の 2〜20% 残ります。実際の暗黒雲の `C⁺/C ≈ 10⁻⁴` より大きく、[C II] 冷却も過大になります。

### H₂・CO・H–O 電荷交換の設定比較

```sh
make -C test molecular-settings-plot
```

「H₂ と CO」「H₂ のみ」「分子なし」のそれぞれを、H–O 電荷交換の有効・無効と組み合わせて比較します。全設定で H II モード、ダスト、宇宙線、軟 X 線、CIE 冷却を有効にします。太陽金属量で、水素核数密度は `10⁻²〜10⁴ cm⁻³`（0.5 dex 刻み）、初期温度は `30, 100, 300, 1000, 3000, 8000, 15000 K`、電離光子束は照射なしと `log10(U) = −5, −3` です。全柱密度は `N_H = 10²⁰ sqrt(n_H / cm⁻³) cm⁻²`、`G0 = 1.7`、`X0 = 1` です。各設定で同じ初期温度の化学平衡を作り、その組成から内部エネルギーを設定してから、10 年・1 万年・100 万年をそれぞれ独立に1ステップ更新します。更新時の許容誤差は化学平衡が `1e-5`、エネルギー方程式が `1e-4` です。

`make -C test molecular-settings-benchmark` は条件別の結果と収束ステータスを `test/out/co/molecular_settings.csv` に保存します。計算時間は単一スレッドの CPU 時間です。全設定で収束した共通条件について、5回の測定の平均時間の中央値と範囲を図に表示します。温度・内部エネルギー変化・化学組成の差は、比較する両設定で収束した条件を対応づけて集計します。データの図と集計表は `plot/out/co/settings/` に保存します。描画済みの時間刻みを変える場合は、リポジトリのルートで `python3 plot/plot_molecular_settings.py --dt 1e6` のように指定します。

1 万年刻みの単一スレッド測定では、H₂ と H–O 電荷交換を有効にしたまま CO を無効にすると、更新1回あたりの CPU 時間の中央値は 453.3 µs から 126.3 µs に短くなりました。

保存した比較では、全6設定で819条件中818条件が収束しました。初期温度 30 K、水素核密度 0.01 cm⁻³、`log10(U) = −5`、時間刻み100万年の1条件で平衡解が得られず、比較集計から外れました。

Jupyter で計算を実行して設定ごとの収束数と失敗した条件を確認できます。 [RunMolecularSettingsBenchmark.ipynb](co/RunMolecularSettingsBenchmark.ipynb)

描画用ノートブックは同じ CSV から計算時間と結果の差を表示します。結果の図に使う時間刻みは `dt_years` で切り替えられます。 [PlotMolecularSettings.ipynb](../plot/PlotMolecularSettings.ipynb)

## 非平衡化学の検証計算と描画

```sh
make -j1 -C test noneq-plot
```

[NonEquilibriumTest.c](noneq/NonEquilibriumTest.c) は、固定密度での照射切り替えと数値積分の精度を検証します。`make -C test noneq` で時間発展と冷却率の比較を、`make -j1 -C test noneq-benchmark` で計算時間の比較を実行できます。CSV は `test/out/noneq/` に保存されます。

- `evolution.csv` は、恒星電離光と遠紫外線の照射・消灯後の時間発展です。恒星電離光のセルは初期温度 8,000 K、`nH = 10 cm⁻³`、ダストなしで、水素の光電離率を `10⁻⁷ s⁻¹` として約16年照射し、消灯後100万年まで追います。遠紫外線のセルは初期温度 100 K、`nH = 100 cm⁻³`、`Nt = 2×10²⁰ cm⁻²`、初期の H₂ 水素核分率 0.8、ダストありで、`G0 = 1000` を30万年照射し、消灯後3,000万年まで追います。両セルの He/H、C/H、O/H 数比はそれぞれ `0.1`、`3×10⁻⁴`、`5×10⁻⁴` で、宇宙線と軟 X 線はオフです。切り替え直前と直後には、同じ `time_yr` に組成が同じ2行を記録します。`phase_time_yr` が照射開始または消灯からの経過時間です。`H2` 列は水素核分率 `2 n(H₂)/nH`、`H_`・`L_` 列は erg cm⁻³ s⁻¹ です。
- `recombination.csv` は、遮蔽と外部照射のない純水素を固定温度1,000 K、`nH = 10 cm⁻³` で再結合させ、case B の解析解 `x_HII = 1/(1 + t/t_rec)` と比較します。`t_rec = 1/(α_B nH)`、初期の `x_HII = 1` です。
- `convergence.csv` は H・He・C・O の混在セル270条件で、同じ時間幅の1ステップと8分割を比較します。`nH = 0.03, 1, 10³ cm⁻³`、初期温度 `30, 100, 8000, 10⁵, 10⁶ K`、時間幅 `100, 2×10⁴, 10⁶ yr`、照射なし・遠紫外線・恒星電離光と遠紫外線、H₂ のオン・オフを組み合わせます。宇宙線はオン、軟 X 線はオフです。温度・比内部エネルギーの相対差と、元素ごとに規格化した化学種分率の最大絶対差を出力します。
- `evolution_direct.csv`、`recombination_direct.csv`、`convergence_direct.csv` は直接評価で同じ計算をした結果です。接尾辞のない CSV は表引きの結果です。時間発展の488時点で温度と全冷却率の相対差、化学種・電子分率の絶対差を評価し、270条件の時間刻み比較でも両モードの到達状態を突き合わせます。
- `benchmark.csv` は、同じ初期セルから時間刻み `10, 10³, 10⁵ yr` でエネルギーを1回更新する場合の実時間を記録します。`equilibrium` は平衡、`nonequilibrium` は表引き、`nonequilibrium_direct` は非平衡の直接評価です。点灯時の中性セルと遠紫外線セルは上記と同条件です。消灯後の再結合セルは `nH = 10 cm⁻³`、初期温度 13,500 K、H⁺ 分率 0.98 で、He・C・O の中性・一価・二価の分率はそれぞれ 0.10、0.20、0.70 です。

冷却率の比較では、非平衡計算の各時点で温度・密度・元素存在量・柱密度・照射を固定して化学平衡を解きます。`evolution.csv` の `L_tot` は時間発展した組成、`eq_L_tot` と `eq_L_*` は平衡組成を共通の冷却式で評価した率です。`eq_native_L_tot` は同じ平衡組成を既存の平衡モードの式で評価した全冷却率です。既存式は [O II]・[O III] を光電離の寄与で重み付けし、[C II]・[O I] の背景放射や H₂ と He の衝突の扱いにも違いがあります。両モードで金属 CIE テーブルと CO はオフです。

[NonEquilibriumBenchmark.c](noneq/NonEquilibriumBenchmark.c) は各モード固有のエネルギー更新 API を使います。両方に精度引数 `10⁻⁶` を渡します。平衡側はエネルギー方程式の収束精度、非平衡側は時間積分の許容誤差に使われます。初期化とウォームアップを済ませた後、1スレッドで毎回同じ初期状態に戻して、更新時間を各条件・モードで5バッチ計ります。図の点は `elapsed_s / calls` の中央値、エラーバーは最小値から最大値までです。測定値は時間刻み、照射、初期組成、実行環境に依存します。

Jupyter の [PlotNonEquilibrium.ipynb](../plot/PlotNonEquilibrium.ipynb) では、照射への応答、平衡・非平衡の冷却率と計算時間、解析解と時間刻みの検証を図と表で確認できます。再計算時の非平衡モードには既定の温度係数の表引きを使います。既定の出力先で必要なデータが不足していればテスト計算を実行します。再計算するときは最初のコードセルの `rerun` を `True` にし、保存済みの別の結果を描くときは `data_dir` を変更してください。

[plot_noneq.py](../plot/plot_noneq.py) は CSV から時間発展図 `ionizing_pulse.png`・`fuv_pulse.png`、数値検証図 `validation.png`、平衡との冷却率比較図 `ionizing_cooling_comparison.png`・`fuv_cooling_comparison.png`、計算時間図 `benchmark.png` を `plot/out/noneq/` に生成します。開発時の表引きの検証には、直接評価との時系列比較図 `ionizing_temperature_lookup.png`・`fuv_temperature_lookup.png` と計算時間図 `temperature_lookup_benchmark.png` を使えます。切り替え時刻や対数軸の表示範囲外の率は CSV に記録します。計算済みのデータを PNG に描き直す場合は `python3 plot/plot_noneq.py --data-dir test/out/noneq` を使います。

### 高温ガスの冷却と消灯後の酸素輝線

```sh
make -j1 -C test noneq-hot
```

同じ固定密度のセルで、CIE 切り替えを導入する前の冷却式（`legacy`）、線形補間で CIE 冷却に切り替える実装（`cie_linear`）、Steffen 補間で CIE 冷却に切り替える現行実装（`cie_bridge`）を比較します。CSV は `test/out/noneq/` に出力します。

従来版の高温セルでは [O III] が全冷却の大きな割合を占めていました。CIE を使う更新版ではその寄与を除いたうえで高温の金属冷却を補うため、この比較では全冷却率が増え、セルがより早く冷えます。初期温度 10⁶ K のセルでは、初期の全冷却率が従来版の約 1.1×10⁻²¹ から現行版の約 1.5×10⁻²⁰ erg cm⁻³ s⁻¹ に増えます。

- `hot_cooling.csv` は、初期温度 10⁶ K と 10⁷ K、水素核密度 10 cm⁻³、太陽金属量の無照射セルの冷却履歴です。初期組成には H II・He III・C III・O III を使います。
- `hot_benchmark.csv` は、これらの高温セルと初期温度 20,000 K のセルを、10年・1,000年・10万年の時間刻みで1回更新するときの CPU 時間を記録します。3方式を同じコンパイル設定と積分の許容誤差で比較し、ウォームアップ後、各条件を単一スレッドで5バッチ測ります。
- `oxygen_after_switch.csv` は、8,000 K のセルを約16年照射し、消灯直前・直後とその後の酸素の電離分率と [O II]・[O III] 冷却率を記録します。

旧 CIE 表での計測では、初期温度 20,000 K の3条件で Steffen 補間が線形補間より1更新の CPU 時間を11〜23%短縮しました。高温の冷却履歴で 25,000 K を下回る時刻の差は0.16%以下でした。同じ時刻の温度を比べると、10⁶ K からの急冷区間では最大約15%の差が出ました。消灯試験の酸素輝線の履歴は3方式で一致しました。

ノートブックでは現行実装の高温冷却履歴と消灯後の酸素輝線を表示します。開発時の実装比較には上記の CSV と、`plot_noneq.py` が出力する `hot_cooling.png`、`hot_benchmark.png`、`oxygen_after_switch.png` を使えます。 [PlotNonEquilibrium.ipynb](../plot/PlotNonEquilibrium.ipynb)

## Wolfire et al. (1995) の検証計算と描画

```sh
make -C test start
mkdir -p test/out/wolfire
(cd test/out/wolfire && ../../test.out >wolfire.stdout 2>wolfire.log)
```

出力ファイルは `test/out/wolfire/out/` に作成されます。求解に失敗した点の詳細は `wolfire.log` に記録されます。失敗点があるとコマンドは非ゼロで終了し、成功した点は出力に残ります。ノートブックに保存した実行では3,300点中3,291点を描画し、柱密度 `10¹⁸ cm⁻²` の低密度側9点で平衡解が得られませんでした。

ノートブックの失敗点の表は、実行ディレクトリにある `rejected.tsv` を読み込みます。このファイルは保存した実行のログを集計したものです。同じディレクトリで再計算する場合は、以前の `rejected.tsv` を別の場所に保管し、今回の失敗点は新しい `wolfire.log` で確認してください。

Jupyter で [PlotWolfire1995.ipynb](../plot/PlotWolfire1995.ipynb) を実行すると、求解できたデータ点が表示されます。冷却・加熱率の図は水素核数密度に対する体積あたりの各成分と合計を表示します。既定は基準条件で、`rate_case` で計算条件を切り替えられます。保存済みの実行を描く場合は、最初のコードセルの `run_dir` を変更してください。旧形式のデータでは記録済みの成分が表示され、追加成分を得るには上の手順で再計算します。描画には Python 3、NumPy、Matplotlib 3.7 以降が必要です。

PNG を生成するには `make -C test wolfire-plot` を実行してください。`test/out/wolfire` の計算結果または同梱の [参照データ](wolfire1995/reference/) から `thermal_rates.png` などを生成し、`plot/out/wolfire/` に保存します。 [plot_wolfire1995.py](../plot/plot_wolfire1995.py)

## Wolfire et al. (2003) との熱平衡曲線の比較

```sh
make -C test/wolfire2003 run
```

太陽円の基準条件で化学・熱平衡を求め、冷相と暖相の安定枝および二相共存圧力を論文の Figure 8 と Table 3 に照合します。比較用の実行ファイルは GNU リンカの `--wrap` で論文の電離・加熱率に近似的にそろえたモデルと、EQCLib の既定の処方で計算した対照を出力します。計算結果は `test/out/wolfire2003/equilibrium.csv` に保存されます。ノートブックの実行には Jupyter、NumPy、pandas、SciPy、Matplotlib が必要です。 [PlotWolfire2003.ipynb](../plot/PlotWolfire2003.ipynb) · [Wolfire et al. (2003)](https://doi.org/10.1086/368016)

## CHIMES との比較

無照射の単一セルと、FUV・遮蔽ありの温度・密度グリッドで EQCLib と CHIMES を比較します。`test/chimes/` の比較ドライバーが、既定では EQCLib の隣にある `Chimes/chimes/` の C ライブラリを呼び出し、`Chimes/chimes-data/` の反応データを使って計算します。`plot/` の比較ノートブックは `test/out/chimes/` に出力した CSV から図と集計を作ります。 [CHIMES 比較の一覧と実行手順](chimes/README.md)

CHIMES の標準比較では銀河シミュレーション向けの設定を土台に、軟 X 線と局所電離光束をゼロにします。 [FUV・遮蔽ありのグリッドの比較条件](chimes/shielded/README.md) · [銀河シミュレーション向けの設定](../README.md#シミュレーションで速度を優先する設定)

遠紫外線・遮蔽ありの格子で、修正前後の EQCLib と CHIMES full の冷却率および等積冷却履歴を比較できます。 [開発時のダスト冷却検証](chimes/shielded/PlotDustCooling.ipynb)
