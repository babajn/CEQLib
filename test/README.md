# EQCLib のテストと比較計算

コマンドは EQCLib のルートで実行します。 [JupyterLab の起動手順](../README.md#テスト結果を-jupyter-で確認する) · [反応・冷却・加熱の物理モデル](../docs/physics.md)

## 基本テスト

インストール時の必須の動作確認です。通常のビルド環境、Python 3、同梱の入力データで実行できます。

```sh
make -C test check cell-check
```

`make -C test check` では、一定照射下での既存の化学平衡解への収束を検証します。照射切り替え後は、水素化学、炭素・ヘリウムの電離・再結合、温度、輝線冷却の応答を確認します。

`make -C test check` は、17,600 K で温度と組成を固定して電離光を切り替えたときに、酸素輝線と CIE の冷却率が連続することと、照射下で熱化学更新が進むことを確認します。 [熱化学のテスト](unit/HydrogenChemistryTest.c)

`make -C test cell-check` は、平衡・非平衡の両モードの初期組成と極短ステップでの温度を検証します。合成スペクトルを使い、UVB の既定値がオフであること、有効時に指定した赤方偏移と両側の遮蔽が結果に反映されることも確認します。

`make -C test check` には合成べき乗スペクトルを使う UVB テストも含まれます。別の数値積分との照合、片側と両側の規格化、柱密度に対する遮蔽とスペクトルの硬化、赤方偏移の補間と入力エラー、エネルギー上限、UVB だけでの電離を確認します。二次電離については、光電子のエネルギー収支、H I と He I が消える極限、電離度の上下限、固定温度での H・He の平衡解と長時間 BDF 解の一致を、無遮蔽と遮蔽ありで検証します。恒星光と UVB を併用する平衡計算では、中性と電離の初期値から同じ組成を得ることも確認します。FG20 が欠落して HM12 だけがある場合や、スペクトルが破損している場合は、再初期化がエラーを返して直前の率を保持することを確認します。 [UVB の単体テスト](unit/UVBackgroundTest.c)

`[expected test diagnostics]` に続くメッセージは、欠落・破損したスペクトルや不正な設定を検出するテストの出力です。各テストの `checks passed` とコマンドの終了コードで成否を確認できます。

### 軟 X 線の打ち切りの回帰テスト

通常の `check` では、電離率の単調な減衰、両端でのなめらかな接続、既存の加熱率を確認します。旧モデルで電子分率 `xe = 0.5` の段差にかかって失敗したセルについて、平衡・熱更新が成功することも検証します。次のコマンドは比較用の `EQC_LEGACY_XRAY_CUTOFF` ビルド `src/libEQCLibLegacyXray.a` を作り、従来の打ち切りでの戻り値と失敗時の状態保持を検証します。FG20 データは不要です。 [軟 X 線の電離モデル](../docs/physics.md#軟-x-線)

```sh
make -C test xray-legacy-check
```

## EQCLib 単体で実行できる検証と作図

`make -C test notebook-data` は、基本テストと H II 領域、H₂・CO 化学、非平衡化学、Wolfire (2003) の計算を順に実行します。ノートブックは `test/out/` の計算結果を読み込みます。Wolfire (1995) と (2003) は、EQCLib の結果を論文や同梱データと照合する検証です。

光電離を含むテストは、既定で `data/PEGASE2/` の同梱スペクトルを使います。別のディレクトリを使う場合は `PEGASE_DIR` に絶対パスを指定します。 [スペクトルの指定方法と設定例](../README.md#pegase2-のスペクトルを指定する)

| 検証・比較 | コードと参照データ | 実行コマンド | 作図ノートブック |
| --- | --- | --- | --- |
| [熱化学更新の計算時間](#熱化学更新の並列ベンチマーク) | [benchmark/](benchmark/) | `make -C test benchmark` | |
| [H II 領域](#h-ii-領域の検証計算と描画) | [hii/](hii/) | `make -C test hii`、`make -C test hii-benchmark` | [平衡・冷却・反応率](../plot/PlotHIIRegion.ipynb) |
| [H₂・CO 化学と反応設定](#h₂co-化学の検証計算と描画) | [co/](co/) | `make -C test co`、`make -C test molecular-settings-benchmark`、[設定比較の計算ノートブック](co/RunMolecularSettingsBenchmark.ipynb) | [H₂・CO 化学](../plot/PlotCOChemistry.ipynb)、[設定比較](../plot/PlotMolecularSettings.ipynb) |
| [非平衡化学と高温ガスの冷却](#非平衡化学の検証計算と描画) | [noneq/](noneq/) | `make -C test noneq`、`make -C test noneq-benchmark`、`make -C test noneq-hot` | [非平衡化学](../plot/PlotNonEquilibrium.ipynb) |
| [宇宙紫外線背景と遮蔽](#fg20-の遮蔽と計算時間) | [uvb/](uvb/) | `make -C test uvb-data`（FG20 スペクトルが必要） | [遮蔽・電離・円盤ガス](../plot/PlotUVBackground.ipynb) |
| [Wolfire et al. (1995)](#wolfire-et-al-1995-の検証計算と描画) | [wolfire1995/](wolfire1995/) | [実行手順](#wolfire-et-al-1995-の検証計算と描画) | [熱平衡と熱収支](../plot/PlotWolfire1995.ipynb) |
| [Wolfire et al. (2003)](#wolfire-et-al-2003-との熱平衡曲線の比較) | [wolfire2003/](wolfire2003/) | `make -C test/wolfire2003 run` | [熱平衡曲線](../plot/PlotWolfire2003.ipynb) |

`common/` にはテスト共通の C ヘッダーがあります。`reference/EQCooling/` には冷却率や反応率の保存済み診断データがあります。H II 領域や Wolfire の比較用データは各テストのフォルダにあります。

### 熱化学更新の並列ベンチマーク

熱化学更新の並列性能は `OMP_NUM_THREADS=12 make -C test benchmark` で測定できます。このベンチマークは OpenMP を使い、原子ガスと分子ガスそれぞれの325点を、10年・1万年・100万年の時間刻みで3ステップずつ進めます。各ステップの収束数、実行時間、逐次と並列の比較を表示します。

### H II 領域の検証計算と描画

```sh
make -C test hii
```

このコマンドは同梱スペクトルを読み込みます。`PEGASE_DIR` を指定すると照射スペクトルを切り替えられます。ノートブックの再描画には、ここで得た CSV を使います。 [スペクトルの指定方法](../README.md#pegase2-のスペクトルを指定する)

[HIIRegionTest.c](hii/HIIRegionTest.c) は、ダストと分子を含まない単一セルを、太陽金属量の PEGASE2 電離スペクトルと case B 再結合で計算します。入射電離光の減衰、非電離 FUV、宇宙線、X 線は含めません。水素核密度 `n_H = 1, 10, 100, 1000 cm⁻³`、ガス金属量 `Z/Z☉ = 0, 0.1, 1`、電離パラメータ `U = NLyn/(n_H c)` の `log10(U) = −4〜−2` の108条件で、等積の熱平衡温度と各元素の電離分率を求めます。さらに `n_H = 100 cm⁻³`、`Z/Z☉ = 1`、`log10(U) = −3` で約3,160〜31,600 K の温度範囲を走査します。

- `test/out/hii/equilibrium.csv` と `thermal_curve.csv` に電離分率と過程別の冷却・加熱率（erg cm⁻³ s⁻¹）を出力します。元素を含まない場合、その電離分率は `nan` になります。
- `test/out/hii/reactions.csv` に H、He、C、O、Si の光電離、電子衝突電離、電子再結合の反応率（cm⁻³ s⁻¹）を出力します。各点で電離と再結合の釣り合いを確認します。
- Jupyter で [PlotHIIRegion.ipynb](../plot/PlotHIIRegion.ipynb) を実行すると、平衡温度・電離分率、冷却・加熱の内訳、反応率を表示します。保存済みの計算結果を描く場合は、最初のコードセルの `run_dir` を変更してください。描画には Python 3、NumPy、Matplotlib 3.7 以降が必要です。

PNG を生成するには `make -C test hii-plot` を実行してください。計算済みの CSV は [plot_hii_region.py](../plot/plot_hii_region.py) でも描画できます。PNG は `plot/out/hii/` に保存されます。

#### H II 領域の計算時間

```sh
make -C test hii-benchmark
```

[HIIRegionBenchmark.c](hii/HIIRegionBenchmark.c) は、電離光なし、照射継続中、照射開始直後を同じ密度・元素組成・初期温度で比較します。照射開始直後の粒子には、照射前の化学組成と内部エネルギーを保持したまま、照射継続中と同じ `NLyn` を与えます。ダストと分子なしの108条件で、固定温度の化学平衡・冷却加熱率と、時間刻み10年・1万年・100万年の陰的エネルギー更新を測ります。

`test/out/hii/hii_benchmark.csv` に、各条件群を繰り返した7ブロックの単一スレッド CPU 時間（1呼び出しあたりの平均 µs）を出力します。測定区間はセルの更新処理です。同じ Xeon w5-3423 の1コアで採った測定例も CSV で参照できます。 [探索の変更前](hii/reference/hii_benchmark_before_xeon_w5_3423.csv)、[変更後](hii/reference/hii_benchmark_xeon_w5_3423.csv)

H II 領域で設定を切り替えた測定では、電荷交換を無効にすると更新が約 20〜35% 速くなりました。化学平衡の許容誤差の上限を `1e-5` に緩めた場合は、更新が約 10〜25% 速くなり、温度の変化は 2×10⁻⁵ 程度以下でした。

CO を無効にした `make -C test hii-benchmark` と同じ条件で、電荷交換を無効にして化学平衡の許容誤差を `1e-5` にすると、エネルギー更新1回あたりの時間は既定値の 0.6〜0.8 倍になりました。電荷交換を加える前の版と比べても 0.3〜0.6 倍です。測定値は単一スレッドの例で、計算機によって変わります。

### H₂・CO 化学の検証計算と描画

```sh
make -C test co
```

[MolecularTest.c](co/MolecularTest.c) は、H II 領域・H₂・CO・ダスト・宇宙線・軟X線・H–O 電荷交換をすべて有効にした設定で化学平衡を解き、反応収支を確かめます。次の2種類の条件を計算します。

- 密度・温度・電離光子束の格子：`n_H = 1〜10⁵ cm⁻³`（0.25 dex 刻み）、`T = 30, 100, 300, 1000 K`、電離光子なしと `log10(U) = −5, −3` の252条件です。`G0 = 1.7`、全柱密度は `N_H = 2×10²¹ (n_H/100)^{1/2} cm⁻²` です。結果は `test/out/co/equilibrium.csv` に出力します。
- 柱密度走査では電離光子を入れず、`T = 30 K`、`n_H = 100, 1000 cm⁻³` と `G0 = 1.7, 17` の組で `N_H = 10¹⁹〜10²³ cm⁻²` の164点を調べます。各点は一つ前の解から始めます。結果は `test/out/co/column_scan.csv` に出力します。

各平衡状態について、H₂ と CO の生成・消滅、電荷交換を含む H と O の電離・再結合の反応量を、公開されている反応率関数から独立に計算し直します。どれかの相対差が `10⁻⁶` を超えるか、最終的な H₂ 量で CO を解き直したときに CO が `10⁻⁶` 以上変わると、試験は失敗します。CSV には存在比、各反応量（cm⁻³ s⁻¹）、相対差を出力します。

Jupyter で [PlotCOChemistry.ipynb](../plot/PlotCOChemistry.ipynb) を実行すると、反応収支、密度に対する分子の割合、柱密度に対する H/H₂ と C⁺/C/CO の遷移、生成と消滅の反応量を表示します。PNG を生成するには `make -C test co-plot` を実行してください。[plot_co_chemistry.py](../plot/plot_co_chemistry.py) が PNG を `plot/out/co/` に保存します。

#### 柱密度走査での CO モデルの確認

柱密度走査では、`A_V ≈ 3〜10` で `CO/C > 0.99` となり、`A_V = 34` では `CO/C` が 0.05 になります。 [CO モデルの適用範囲](../docs/physics.md#co-モデルの適用範囲)

宇宙線による C の電離と He⁺ + CO → C⁺ + O + He を KC08e に加え、走査の反応率を使った試算では、深部の `CO/C` は 0.8〜0.98 に保たれますが、C⁺ を除く反応が遅い C⁺ + H₂ だけなので、`A_V = 3〜50` で C⁺ が炭素の 2〜20% 残ります。実際の暗黒雲の `C⁺/C ≈ 10⁻⁴` より大きく、[C II] 冷却も過大になります。

#### H₂・CO・H–O 電荷交換の設定比較

```sh
make -C test molecular-settings-plot
```

「H₂ と CO」「H₂ のみ」「分子なし」のそれぞれを、H–O 電荷交換の有効・無効と組み合わせて比較します。全設定で H II モード、ダスト、宇宙線、軟 X 線、CIE 冷却を有効にします。太陽金属量で、水素核数密度は `10⁻²〜10⁴ cm⁻³`（0.5 dex 刻み）、初期温度は `30, 100, 300, 1000, 3000, 8000, 15000 K`、電離光子束は照射なしと `log10(U) = −5, −3` です。全柱密度は `N_H = 10²⁰ sqrt(n_H / cm⁻³) cm⁻²`、`G0 = 1.7`、`X0 = 1` です。各設定で同じ初期温度の化学平衡を作り、その組成から内部エネルギーを設定してから、10 年・1 万年・100 万年をそれぞれ独立に1ステップ更新します。更新時の許容誤差は化学平衡が `1e-5`、エネルギー方程式が `1e-4` です。

`make -C test molecular-settings-benchmark` は条件別の結果と収束ステータスを `test/out/co/molecular_settings.csv` に保存します。計算時間は単一スレッドの CPU 時間です。全設定で収束した共通条件について、5回の測定の平均時間の中央値と範囲を図に表示します。温度・内部エネルギー変化・化学組成の差は、比較する両設定で収束した条件を対応づけて集計します。データの図と集計表は `plot/out/co/settings/` に保存します。描画済みの時間刻みを変える場合は、リポジトリのルートで `python3 plot/plot_molecular_settings.py --dt 1e6` のように指定します。

1 万年刻みの単一スレッド測定では、H₂ と H–O 電荷交換を有効にしたまま CO を無効にすると、更新1回あたりの CPU 時間の中央値は 453.3 µs から 126.3 µs に短くなりました。

保存した比較では、全6設定で819条件中818条件が収束しました。初期温度 30 K、水素核密度 0.01 cm⁻³、`log10(U) = −5`、時間刻み100万年の1条件で平衡解が得られず、比較集計から外れました。

Jupyter で計算を実行して設定ごとの収束数と失敗した条件を確認できます。 [RunMolecularSettingsBenchmark.ipynb](co/RunMolecularSettingsBenchmark.ipynb)

描画用ノートブックは同じ CSV から計算時間と結果の差を表示します。結果の図に使う時間刻みは `dt_years` で切り替えられます。 [PlotMolecularSettings.ipynb](../plot/PlotMolecularSettings.ipynb)

### 非平衡化学の検証計算と描画

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

#### 高温ガスの冷却と消灯後の酸素輝線

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

### Wolfire et al. (1995) の検証計算と描画

```sh
make -C test start
mkdir -p test/out/wolfire
(cd test/out/wolfire && ../../test.out >wolfire.stdout 2>wolfire.log)
```

全33条件の再計算は検証環境で約8分（496秒）でした。出力ファイルは `test/out/wolfire/out/` に作成され、求解に失敗した点は `wolfire.log` に記録されます。同梱の参照計算では3,300点中3,290点が収束し、柱密度 `10¹⁸ cm⁻²` の低密度側10点で求解に失敗しました。失敗点があると `test.out` は終了コード `1` を返し、収束した点は出力に残ります。

ノートブックは再計算時の `wolfire.log` から失敗点を抽出します。同梱の参照データを描画するときは、参照計算の失敗点を記録した `rejected.log` を使います。

Jupyter で [PlotWolfire1995.ipynb](../plot/PlotWolfire1995.ipynb) を実行すると、求解できたデータ点が表示されます。冷却・加熱率の図は水素核数密度に対する体積あたりの各成分と合計を表示します。既定は基準条件で、`rate_case` で計算条件を切り替えられます。保存済みの実行を描く場合は、最初のコードセルの `run_dir` を変更してください。旧形式のデータでは記録済みの成分が表示され、追加成分を得るには上の手順で再計算します。描画には Python 3、NumPy、Matplotlib 3.7 以降が必要です。

PNG を生成するには `make -C test wolfire-plot` を実行してください。`test/out/wolfire` の計算結果または同梱の [参照データ](wolfire1995/reference/) から `thermal_rates.png` などを生成し、`plot/out/wolfire/` に保存します。 [plot_wolfire1995.py](../plot/plot_wolfire1995.py)

### Wolfire et al. (2003) との熱平衡曲線の比較

```sh
make -C test/wolfire2003 run
```

太陽円の基準条件で化学・熱平衡を求め、冷相と暖相の安定枝および二相共存圧力を論文の Figure 8 と Table 3 に照合します。比較用の実行ファイルは GNU リンカの `--wrap` で論文の電離・加熱率に近似的にそろえたモデルと、EQCLib の既定の処方で計算した対照を出力します。計算結果は `test/out/wolfire2003/equilibrium.csv` に保存されます。ノートブックの実行には Jupyter、NumPy、pandas、SciPy、Matplotlib が必要です。 [PlotWolfire2003.ipynb](../plot/PlotWolfire2003.ipynb) · [Wolfire et al. (2003)](https://doi.org/10.1086/368016)

### FG20 の遮蔽と計算時間

FG20 の実スペクトルを使う任意の検証です。`fg20_cloudy.ascii` を置いたディレクトリを `UVB_DIR` に指定します。既定値はリポジトリルートの `data/FG20/` です。 [FG20 の取得方法](../README.md#fg20-のスペクトルを用意する) · [CHIMES・Cloudy との照合](#uvb-の外部データ外部コードとの照合)

#### 検証計算と作図

```sh
make -C test uvb-data UVB_DIR=/absolute/path/to/FG20
```

`uvb-data` は FG20 で EQCLib の検証計算を行い、作図データを `test/out/uvb/` に書き出します。[PlotUVBackground.ipynb](../plot/PlotUVBackground.ipynb) は、遮蔽による一次電離率と光電子の運動エネルギー供給率の減衰、固定温度の H・He 平衡、二次電離、円盤ガスでの UVB の有無による差を表示します。`run_dir` を変更すると保存済みの出力を読み込めます。

`disk-history.csv` には円盤格子144条件の各12更新について、温度、H II・電子分率、各更新に渡した片側の遮蔽柱密度、加熱率を記録します。ノートブックで履歴を描く密度と軟 X 線の有無を選べます。`make -C test uvb-plot` は、計算から `plot/out/uvb/` への PNG 出力までを行います。計算済みのデータから描画するには [plot_uvb.py](../plot/plot_uvb.py) を使います。

`cooling-curve.csv` は、密度 `nH = 0.01 cm⁻³`、上下それぞれの UVB 遮蔽柱密度 `10¹⁸ cm⁻²` を固定した温度走査です。`z = 0, 2` で `10²〜10⁶ K` を 0.05 dex 刻みで走査し、各温度で化学平衡を求めます。ノートブックと `cooling_curve.png` では、軟 X 線の設定をそろえ、UVB オン・オフの全冷却率 `L_tot` とその比を温度に対して描きます。

`z = 0, 1, 2, 3` について、無遮蔽と上下の柱密度が対称・非対称な場合の一次電離率を `rates.csv` に、無遮蔽・固定温度10⁴ Kの H・He 平衡組成を `equilibrium.csv` に保存します。`soft-xray.csv` は同じ組成での UVB と `X0 = 1` の軟 X 線の加熱率、`UVBackgroundRates.dat` は運動エネルギー100 eV以上の光電子が供給エネルギーに占める割合を記録します。

`secondary.csv` は、赤方偏移0と2、片側の柱密度 `10¹⁹〜10²¹ cm⁻²`、H II 分率 `10⁻³〜10⁻¹` で、熱、二次電離、励起などへの配分を記録します。全スペクトルと、光子エネルギーを 300 eV で打ち切った場合を比べます。

`disk-updates.csv` は、100 pc の局所柱密度と前ステップの組成を使い、1 Myr 刻みで12回の平衡熱更新をした結果です。旧モデルで失敗した希薄セル（`z = 0`、`nH = 0.01 cm⁻³`）が軟 X 線と UVB を併用しても全ステップ成功することを確認します。分子ガスでは、推奨する `xHI + 2*xH2` の有効柱密度と、比較用の `xHI` だけの指定を両側に与えます。

`disk-grid.csv` は、宇宙線を有効にした `z = 0, 2`、`nH = 10⁻³〜10 cm⁻³`（0.5 dex 刻み）での更新結果です。UVB と軟 X 線のオン・オフ、平衡・非平衡を組み合わせた144条件で、各12ステップの成功を確認します。温度・電子分率・加熱の内訳と全冷却率 `L_tot` を出力します。ノートブックと `disk_cooling.png` では、同じ密度で UVB オン・オフの全冷却率とその比を表示します。比較には12 Myr 後の温度と組成の違いも含まれます。FG20 だけの放射寄与を調べる基準は軟 X 線オフです。 [円盤セルの再現条件](uvb/UVBackgroundReference.c)

#### 金属吸収・光電離と実行時間の追加検証

```sh
make -C test uvb-metal uvb-metal-sensitivity uvb-metal-photoionization-disk
make -C test uvb-benchmark
```

`uvb-metal` は Python 3 の NumPy、SciPy、h5py が必要です。`Z/Z☉ = 0, 0.1, 1, 2`、片側の柱密度 `10¹⁹〜10²² cm⁻²` に加え、赤方偏移・金属量・上下の柱密度を同時に補間する点を、独立の Python 積分と比較します。総断面積を Morrison & McCammon (1983) の表2と照合し、金属のパワーを含むエネルギー収支と、公開 API の金属量ゼロでの互換性を検証します。`metal-comparison.csv`、`metal-cross-section-comparison.csv`、`metal-deep-comparison.csv`、`metal-deep-bands.csv`、`metal-chemistry-comparison.csv` に結果、`metal-timing.csv` に初期化と赤方偏移の設定時間、保持メモリを保存します。 [独立積分コード](uvb/metal_compare.py)

`uvb-metal-sensitivity` は太陽金属量、宇宙線・軟 X 線ありで、FG20 全域と光子上限 500 eV を比較します。`z = 0, 2`、`G0 = 0.1, 1`、`X0 = 0.1, 1`、`nH = 10⁻³〜10² cm⁻³`（0.5 dex 刻み）、平衡・非平衡の全352ケースを12ステップ進め、`disk-500ev.csv` に出力します。軟 X 線を両方の計算に含めるので、この比較は高エネルギー側の感度を測るものです。 [円盤試験](uvb/UVBackgroundReference.c)

太陽金属量、`nH = 1 cm⁻³`、`G0 = 0.1`、`X0 = 1` の円盤試験では、Wolfire の項を含めたまま FG20 を 500 eV で打ち切ると、平衡モードで12ステップ後の温度が `z = 0` で `116.96 → 116.67 K`、`z = 2` で `169.52 → 149.59 K` になりました。この差には FG20 の高エネルギー光子の寄与全体が含まれます。

`uvb-metal-photoionization-disk` は円盤格子に H–O 電荷交換と、UVB を含む酸素輝線の重みの比較を追加します。恒星光は `NLyn = 0`、各ケース12更新で、温度、電子分率、イオン分率と輝線冷却を `disk-metal-photoionization.csv` と `disk-metal-photoionization-oxygen-weight.csv` に保存します。後者は冷却モジュールの試験用ビルドです。 [円盤試験](uvb/UVBackgroundReference.c)

電荷交換オフの比較では、`z = 2`、`nH = 0.01 cm⁻³` の平衡円盤で `L_res` が97%下がり、温度が `15151 → 15830 K` になりました。

`uvb-benchmark` は既存の非平衡ベンチマークを使い、同じ初期条件で UVB を無効・有効にしたときの平衡・非平衡両モードのセル更新あたりの CPU 時間を測ります。表の初期化を除いた測定で、照射で異なる組成・温度へ収束する計算負荷も含みます。各モードの比較は `test/out/uvb/benchmark-summary.csv` に保存します。

`uvb_reference.out` と `uvb_asura_benchmark.out` は、オプション `--zero-metallicity` を付けるとセルの `Z` をゼロにします。既存の元素数密度やダスト量を保持したケースで、金属の追加前後を比較できます。

#### UVB 照射下の H–O 電荷交換

```sh
make -C test uvb-charge-transfer-disk uvb-charge-transfer-bdf
make -C test uvb-charge-transfer-benchmark BENCHMARK_CPU=10
```

`uvb-charge-transfer-disk` は太陽金属量、`z = 0, 2`、`nH = 10⁻³〜10² cm⁻³`（0.5 dex 刻み）、`G0 = 0.1, 1`、`X0 = 0.1, 1` で12回ずつ更新します。恒星光は `NLyn = 0` とし、電荷交換オフ、H II 領域モードを使わない電荷交換オン、H II 領域モードでの電荷交換オンを平衡・非平衡の両方で比較します。`disk-charge-transfer.csv` に各ケースの最終状態と [O I] 冷却、`disk-charge-transfer-budget.csv` に各更新での電荷残差と H–O 電荷交換の反応収支を保存します。 [円盤試験](uvb/UVBackgroundReference.c)

太陽金属量の円盤88条件を12ステップ進めて双方向解と比べると、最終温度の相対差は最大0.148%でした。このケースは `z = 2`、`nH = 1 cm⁻³`、`G0 = 1`、`X0 = 0.1` で、電荷交換による H の正味反応率 `F_CT` は通常の H 反応のフラックスの0.54%でした。最終状態で双方の分率が `10⁻⁴` 以上の H I・H II・O I・O II・O III の相対差は最大0.82%でした。微量な段階では相対差が大きくなります。全更新で、イオンの電荷の合計と試行電子密度の相対残差は指定した化学平衡の許容誤差 `10⁻⁵` 以下でした。 [電荷交換の式](../docs/physics.md#ho-電荷交換)

前ステップの組成から次の柱密度を計算する円盤試験では、遮蔽が進む途中に差が大きくなりました。`z = 2`、`nH = 0.0316 cm⁻³`、`G0 = X0 = 1` の第7ステップでは温度差が0.753%でした。同じ開始状態と柱密度を与えた1回の更新では、差は最大0.028%でした。`X0 = 0.1` の第8ステップでは H II 分率が一方向で0.285、双方向で0.345でした。

`uvb-charge-transfer-bdf` は温度固定で H・He・O の非平衡 BDF 積分を長時間進め、一方向と双方向の平衡解との比較を `uvb-charge-transfer-bdf.csv` に保存します。酸素なし、電荷交換オフ、H II 領域モード、通常の非平衡更新の出力は、ビット比較用の `uvb-charge-transfer-regression.csv` に16進の浮動小数点形式で保存します。 [固定温度試験](uvb/UVBChargeTransferReference.c)

`uvb-charge-transfer-benchmark` は、同じ ASURA の初期状態から FG20 を全域で使い、一方向・双方向の更新時間を `asura-charge-transfer.csv` と `asura-charge-transfer-bidirectional.csv` に出力します。比較用の初期状態は電荷交換オフで暖機します。 [計測の条件](uvb/ASURABenchmark.c)

#### ASURA 設定でのセル更新時間

```sh
make -C test uvb-asura-benchmark BENCHMARK_CPU=10
```

`uvb-asura-benchmark` は `taskset` で計算を1個の論理 CPU に固定します。`BENCHMARK_CPU` の既定値はプロセスに割り当てられた最初の論理 CPU です。ASURA の推奨設定を使い、`G0 = X0 = 1`、柱密度 `Nt = nH × 100 pc`、`z = 0, 2`、`nH = 10⁻³〜10² cm⁻³`（0.5 dex 刻み）を評価します。前ステップの `xHI + 2*xH2` を柱密度に掛け、時間刻み10万年で6回の平衡更新を行い、その組成と温度から平衡・非平衡をそれぞれ1回進める時間を測ります。化学平衡の許容誤差は `1e-5`、エネルギーの許容誤差は `1e-4` です。

UVB 無効（`off`）、光子上限110 eV（`cut110`）、FG20 全域（`full`）の132行を `test/out/uvb/asura-benchmark.csv` に出力します。110 eV 上限では H・He から生じる光電子は100 eV未満です。金属量が正の場合は、一部の金属から100 eV以上の電子も生じます。`cpu_us` は1セル更新あたりの CPU 時間で、5回の計測バッチの中央値です。平衡では別の計測用呼び出しで温度・電子密度の試行回数と二次電離の固定点反復回数を記録します。 [ベンチマークの条件](uvb/ASURABenchmark.c)

## 外部コードとの比較（任意）

比較計算に応じて、外部コードや外部の反応データを用意します。

| 比較 | 再計算に必要なもの | 作図ノートブック |
| --- | --- | --- |
| [Cloudy C23.01 との H II 領域比較](#cloudy-との比較) | Cloudy の実行ファイル | [H II 領域ノートブック](../plot/PlotHIIRegion.ipynb) |
| [CHIMES との比較](#chimes-との比較) | CHIMES のソース、反応データ、HDF5、SUNDIALS | [比較別のノートブック一覧](chimes/README.md) |
| [UVB の外部データ・外部コードとの照合](#uvb-の外部データ外部コードとの照合) | FG20 と CHIMES の反応データ、または Cloudy c25.00 のデータ・実行ファイル | |

### Cloudy との比較

[Cloudy C23.01](https://data.nublado.org/cloudy_releases/c23/) をビルドし、EQCLib のルートで次を実行します。

```sh
make -C test hii
python3 test/hii/cloudy/run_cloudy_hii.py --cloudy /path/to/c23.01/source/cloudy.exe
```

スペクトルは `--pegase-dir`、環境変数 `PEGASE_DIR`、同梱の `data/PEGASE2/` の順に参照します。

[run_cloudy_hii.py](hii/cloudy/run_cloudy_hii.py) は EQCLib の `equilibrium.csv` から密度、電離パラメータ、元素組成を読み取ります。PEGASE2 の連続スペクトルを EQCLib と同じ水素電離光子数で年齢ごとに規格化し、0〜30 Myr の18時点を等重みで平均します。Cloudy には case B、ダストと分子なし、厚さ 1 cm の1ゾーンを指定します。密度、入射電離光子束、光学的厚さ、熱収支を確認し、冷却の出力から O II と O III の正味の衝突線冷却を集計します。EQCLib の比較量には H II 領域で使われる温度重みを含めています。

各ゾーンの入力、ログ、集計 CSV は `test/out/hii/cloudy/` に出力します。[PlotHIIRegion.ipynb](../plot/PlotHIIRegion.ipynb) は再計算の出力を優先し、出力がなければ [Cloudy C23.01 の保存済み参照結果](hii/cloudy/cloudy_hii_c23.01.csv) を使います。参照結果で描画する場合は EQCLib の `make -C test hii` の出力があれば実行できます。Cloudy のデータ更新と計算方法は [C23 リリース論文](https://doi.org/10.22201/ia.01851101p.2023.59.02.12) と [C23.01 リリース論文](https://doi.org/10.3847/2515-5172/ad0e75) を参照してください。対応する PDF は [reference/](../reference/index.tsv) に保存しています。

H–O 電荷交換を無効にすると、Cloudy との比較における太陽金属量の平衡温度の差は、中央値で 1.4% から 5.3% に広がりました。

### CHIMES との比較

無照射の単一セルと、FUV・遮蔽ありの温度・密度グリッドで EQCLib と CHIMES を比較します。`test/chimes/` の比較ドライバーが、既定では EQCLib の隣にある `Chimes/chimes/` の C ライブラリを呼び出し、`Chimes/chimes-data/` の反応データを使って計算します。`plot/` の比較ノートブックは `test/out/chimes/` に出力した CSV から図と集計を作ります。 [CHIMES 比較の一覧と実行手順](chimes/README.md)

CHIMES の標準比較では銀河シミュレーション向けの設定を土台に、軟 X 線と局所電離光束をゼロにします。 [FUV・遮蔽ありのグリッドの比較条件](chimes/shielded/README.md) · [銀河シミュレーション向けの設定](../README.md#シミュレーションで速度を優先する設定)

遠紫外線・遮蔽ありの格子で、修正前後の EQCLib と CHIMES full の冷却率および等積冷却履歴を比較できます。 [開発時のダスト冷却検証](chimes/shielded/PlotDustCooling.ipynb)

### UVB の外部データ・外部コードとの照合

各ターゲットには FG20 のスペクトルが必要です。`uvb` と `uvb-metal-photoionization` は CHIMES の HDF5 データを読み、`uvb-metal-cloudy` は Cloudy c25.00 を実行します。Python 3 の NumPy、SciPy、h5py を用意します。FG20 を `data/FG20/` に置いた場合は、下の `export` 行を省略できます。

```sh
export UVB_DIR=/absolute/path/to/FG20
make -C test uvb CHIMES_DATA_DIR=/absolute/path/to/chimes-data
make -C test uvb-metal-photoionization CLOUDY_DIR=/absolute/path/to/c25.00 CHIMES_DATA_DIR=/absolute/path/to/chimes-data
make -C test uvb-metal-cloudy CLOUDY_DIR=/absolute/path/to/c25.00
```

`CHIMES_DATA_DIR` は `chimes_main_data.hdf5` と `FG20_cross_sections/` を含むディレクトリです。既定では EQCLib の隣にある `Chimes/chimes-data/` を使います。`CLOUDY_DIR` の既定では EQCLib の隣にある `c25.00/` を使います。 [CHIMES データの取得方法](chimes/README.md)

`uvb` は `uvb-data` を実行した上で、一次電離率・光電子パワー・二次電離を独立の平板積分と照合します。H・He の固定温度10⁴ Kでの平衡組成は、各イオンの定常比と電荷保存から求めた参照解と比較します。無遮蔽での一次電離率は CHIMES の表と照合し、結果を `comparison.csv` と `secondary-comparison.csv` などに保存します。 [独立積分と比較コード](uvb/compare.py)

`uvb-metal-photoionization` は C・O・Si の6つの率を、Cloudy c25 の `phfit.dat` と `mewe_nelectron.dat` に基づく独立積分で照合します。しきい値、遮蔽、補間点、金属の光電離を加えたエネルギー収支を検証し、`metal-photoionization-integrals.csv`、`metal-photoionization-auger.csv`、`metal-photoionization-chimes.csv` と `metal-photoionization-impact-estimate.csv` に結果を出力します。CHIMES の `photoion_euv` では C II、O I、O II、Si II を比較し、単一電離と対応する Auger チャンネルを併記します。C I と Si I は CHIMES の FUV グループに入り、その率は13.6 eV未満も含みます。 [独立積分コード](uvb/metal_photoionization_compare.py)

CHIMES の表生成コードでは、最初の内殻端より下でも、外殻の総断面積に C・O の 2s や Si の 3s を別途加えています。この領域の Verner et al. (1996) のフィットにはそれらの殻がすでに含まれます。EQCLib は Verner の `phfit2.f` と同じ殻の切り替えを使うため、生成済みの CHIMES 表との率には系統差が出ます。 [CHIMES の表生成コード](https://bitbucket.org/richings/chimes-tools/src/411fd88eb38ae05bd47adf598f1db5dc08fb4579/generate_cross_sections/generate_cross_sections.py) · [Verner の配布ページ](https://www.pa.uky.edu/~verner/photo.html)

`uvb-metal-cloudy` は `source/cloudy.exe` と `data/hm12_galaxy.ascii` を使い、HM12 と FG20 のそれぞれで同じ13.6 eV以上のスペクトルを EQCLib と Cloudy に与えます。`z = 0, 2`、`T = 10⁴ K`、`nH = 10⁻³〜1 cm⁻³` の単一ゾーンのイオン分率を `metal-photoionization-cloudy.csv` に出力します。Cloudy の IV 価以上の割合も保存します。 [単一ゾーン比較](uvb/metal_cloudy.py)
