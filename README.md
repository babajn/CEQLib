# EQCLib (Equilibrium Cooling)

星間ガスの化学平衡と冷却・加熱率を計算する C ライブラリです。固定温度での化学平衡と、化学平衡を仮定したエネルギー方程式の陰的更新に対応します。

## ビルドと検証

C99 対応の GCC と GNU Scientific Library（GSL）が必要です。以下の `GSL_DIR` は GSL のインストール先に合わせてください。ベンチマークは OpenMP を使用します。

```sh
make -C src GSL_DIR="$HOME/lib/gsl" libEQCLib.a
make -C test GSL_DIR="$HOME/lib/gsl" check
OMP_NUM_THREADS=12 make -C test GSL_DIR="$HOME/lib/gsl" benchmark
```

ベンチマークは、原子ガスと分子ガスそれぞれの325点を、10年・1万年・100万年の時間刻みで3ステップずつ進めます。各ステップの収束数、実行時間、逐次と並列の比較を表示します。

## シミュレーションコードからの利用

公開 API は [ChemicalEquilibriumCooling.h](src/ChemicalEquilibriumCooling.h) にあります。ホストコードを `libEQCLib.a`、`libgsl`、`libgslcblas`、`libm` とリンクします。
粒子の状態は [EQCDataStructure.h](src/EQCDataStructure.h) の `StructEQCData` で渡し、ステータスを返す API では `EQC_SUCCESS` で成功を判定します。

実行パラメータを設定し、並列計算を始める前に `EQCLibTryInitialization` を呼び出してください。H II 領域の計算には PEGASE2 の入力ディレクトリを指定できます。エネルギー更新では `EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating` の成功を確認してから、呼び出し側で `U += dt * Du` を適用します。

H II 領域の酸素冷却を使うには、初期化前に `EQCLibSetRunParamHIIRegion(true)` を呼び、各粒子の `NLyn` に局所の水素電離光子束を設定します。光電離ガスの CIE 金属冷却を [O II] と [O III] の冷却に徐々に置き換え、15,000〜25,000 K の間で CIE に戻します。各遷移の光学的に薄い放射率は [HIIregionLineCooling.h](src/HIIregionLineCooling.h) の `CoolingLinesOII` と `CoolingLinesOIII` で取得できます。

この設定では、O⁰ + H⁺ ⇄ O⁺ + H⁰ と O²⁺ + H⁰ → O⁺ + H⁺ の電荷交換も解きます。H の電離分率を1変数とし、O の電離平衡と H I・H II の反応平衡を、電荷交換を含めて同時に満たします。H₂ 生成を有効にしている場合は、H の電離分率を変えるたびに H₂ の生成と消滅の釣り合いも解き直すので、H₂、H、O の反応がすべて同時に釣り合います。

### シミュレーションで速度を優先する設定

分子雲内部を分解せずに銀河全体の熱進化を追う ASURA 計算では、**H₂ オン・CO オフ・H–O 電荷交換オン**を推奨します。ダスト遮蔽も有効にし、初期化前に次のように設定します。局所電離光子束 `NLyn` を渡す場合は `useLocalIonizingFlux` を `true` にします。

```c
EQCLibSetRunParamRecombinationCaseB(true);
EQCLibSetRunParamDustAttenuation(true);
EQCLibSetRunParamGrainSurfaceProcesses(true);
EQCLibSetRunParamCosmicRayIonization(true);
EQCLibSetRunParamSoftXrayIonization(true);
EQCLibSetRunParamCIEResonanceCooling(true);
EQCLibSetRunParamHIIRegion(useLocalIonizingFlux);
EQCLibSetRunParamH2Formation(true, 1.0);
EQCLibSetRunParamCOFormation(false, 1.0);
EQCLibSetRunParamChargeTransfer(true);
EQCLibSetRunParamChemistryAccuracy(1.0e-5);
```

陰的エネルギー更新の許容誤差は `1e-4` を目安にします。1 万年刻みの単一スレッド測定では、CO を無効にすると EQCLib の更新時間は 1 回あたり 385 µs から 128 µs に短くなりました。 [設定比較ノートブック](plot/PlotMolecularSettings.ipynb)

CO 輝線を予測する場合は CO を有効にします。遮蔽された低温ガスの熱収支が重要な場合は、代表的な粒子で CO の有無による冷却率の差を確認してください。 [CO モデルの適用範囲](#co-モデルの適用範囲)

### 各設定の詳細

| 設定 | 既定値 | 設定例 | 効果 |
| --- | --- | --- | --- |
| `EQCLibSetRunParamCOFormation` | `false` | `true, 1.0` | CO の化学平衡と輝線冷却を有効にします。第2引数は CO 生成率の倍率です。 |
| `EQCLibSetRunParamChargeTransfer` | `true` | `false` | H II 領域の更新が約 20〜35% 速くなります。物理モデルが変わり、H I が残る条件で O III が過剰になります。Cloudy との比較では、太陽金属量の平衡温度の差が中央値 1.4% から 5.3% に広がります。H II 領域の設定が無効なときは効果がありません。 |
| `EQCLibSetRunParamChemistryAccuracy` | `1e-8` | `1e-5` | 陰的なエネルギー更新の中で解く化学平衡の許容誤差の上限です。実際には `min(値, 0.1 × エネルギーの精度)` を使います。更新が約 10〜25% 速くなり、温度の変化は 2×10⁻⁵ 程度以下です。正味冷却率の API は内部で精度 1e-8 を指定するので、この設定にかかわらず 1e-9 以下で解きます。 |

CO 化学を使う場合は、`EQCLibSetRunParamH2Formation(true, 1.0)` と `EQCLibSetRunParamCOFormation(true, 1.0)` を設定します。CO を無効にするには `EQCLibSetRunParamCOFormation(false, 1.0)` を呼びます。無効にすると、更新後の CO 存在量と CO 輝線冷却はゼロになり、炭素・酸素の化学状態や冷却率にも影響します。[CO モデルの適用範囲](#co-モデルの適用範囲)

CO を無効にした `make hii-benchmark` と同じ条件で、電荷交換を無効にして化学平衡の許容誤差を `1e-5` にすると、エネルギー更新1回あたりの時間は既定値の 0.6〜0.8 倍になりました。電荷交換を加える前の版と比べても 0.3〜0.6 倍です。測定値は単一スレッドの例で、計算機によって変わります。

### H₂・CO・H–O 電荷交換の設定比較

```sh
make -C test molecular-settings-plot
```

「H₂ と CO」「H₂ のみ」「分子なし」のそれぞれを、H–O 電荷交換の有効・無効と組み合わせて比較します。全設定で H II モード、ダスト、宇宙線、軟 X 線、CIE 冷却を有効にします。太陽金属量で、水素核数密度は `10⁻²〜10⁴ cm⁻³`（0.5 dex 刻み）、初期温度は `30, 100, 300, 1000, 3000, 8000, 15000 K`、電離光子束は照射なしと `log10(U) = −5, −3` です。全柱密度は `N_H = 10²⁰ sqrt(n_H / cm⁻³) cm⁻²`、`G0 = 1.7`、`X0 = 1` です。各設定で同じ初期温度の化学平衡を作り、その組成から内部エネルギーを設定してから、10 年・1 万年・100 万年をそれぞれ独立に1ステップ更新します。更新時の許容誤差は化学平衡が `1e-5`、エネルギー方程式が `1e-4` です。

`make -C test molecular-settings-benchmark` は条件別の結果と収束ステータスを `test/out/co/molecular_settings.csv` に保存します。計算時間は単一スレッドの CPU 時間です。全設定で収束した共通条件について、5回の測定の平均時間の中央値と範囲を図に表示します。温度・内部エネルギー変化・化学組成の差は、比較する両設定で収束した条件を対応づけて集計します。データの図と集計表は `plot/out/co/settings/` に保存します。描画済みの時間刻みを変える場合は、リポジトリのルートで `python3 plot/plot_molecular_settings.py --dt 1e6` のように指定します。

Jupyter では [RunMolecularSettingsBenchmark.ipynb](test/RunMolecularSettingsBenchmark.ipynb) で計算を実行し、設定ごとの収束数と失敗した条件を確認できます。[PlotMolecularSettings.ipynb](plot/PlotMolecularSettings.ipynb) は同じ CSV を読み込み、計算時間と結果の差を表示します。結果の図に使う時間刻みは、ノートブックの `dt_years` で切り替えられます。

## H II 領域の検証計算と描画

```sh
make -C test GSL_DIR="$HOME/lib/gsl" hii
```

[HIIRegionTest.c](test/HIIRegionTest.c) は、ダストと分子を含まない単一セルを、太陽金属量の PEGASE2 電離スペクトルと case B 再結合で計算します。入射電離光の減衰、非電離 FUV、宇宙線、X 線は含めません。水素核密度 `n_H = 1, 10, 100, 1000 cm⁻³`、ガス金属量 `Z/Z☉ = 0, 0.1, 1`、電離パラメータ `U = NLyn/(n_H c)` の `log10(U) = −4〜−2` の108条件で、等積の熱平衡温度と各元素の電離分率を求めます。さらに `n_H = 100 cm⁻³`、`Z/Z☉ = 1`、`log10(U) = −3` で約3,160〜31,600 K の温度範囲を走査します。

- `test/out/hii/equilibrium.csv` と `thermal_curve.csv` に電離分率と過程別の冷却・加熱率（erg cm⁻³ s⁻¹）を出力します。元素を含まない場合、その電離分率は `nan` になります。
- `test/out/hii/reactions.csv` に H、He、C、O、Si の光電離、電子衝突電離、電子再結合の反応率（cm⁻³ s⁻¹）を出力します。各点で電離と再結合の釣り合いを確認します。
- Jupyter で [PlotHIIRegion.ipynb](plot/PlotHIIRegion.ipynb) を実行すると、平衡温度・電離分率、冷却・加熱の内訳、反応率を表示します。保存済みの計算結果を描く場合は、最初のコードセルの `run_dir` を変更してください。描画には Python 3、NumPy、Matplotlib 3.7 以降が必要です。

PNG を生成するには `make -C test GSL_DIR="$HOME/lib/gsl" hii-plot` を実行してください。計算済みの CSV は [plot_hii_region.py](plot/plot_hii_region.py) でも描画できます。PNG は `plot/out/hii/` に保存されます。

### H II 領域の計算時間

```sh
make -C test GSL_DIR="$HOME/lib/gsl" hii-benchmark
```

[HIIRegionBenchmark.c](test/HIIRegionBenchmark.c) は、電離光なし、照射継続中、照射開始直後を同じ密度・元素組成・初期温度で比較します。照射開始直後の粒子には、照射前の化学組成と内部エネルギーを保持したまま、照射継続中と同じ `NLyn` を与えます。ダストと分子なしの108条件で、固定温度の化学平衡・冷却加熱率と、時間刻み10年・1万年・100万年の陰的エネルギー更新を測ります。

`test/out/hii/hii_benchmark.csv` に、各条件群を繰り返した7ブロックの単一スレッド CPU 時間（1呼び出しあたりの平均 µs）を出力します。図はその中央値と最小〜最大を表示します。同じ Xeon w5-3423 の1コアで採った [探索の変更前](test/reference/hii_benchmark_before_xeon_w5_3423.csv) と [変更後](test/reference/hii_benchmark_xeon_w5_3423.csv) の比較も掲載しています。ローカルの CSV がない場合、変更後の測定例を使います。初期化、照射前後の初期状態の準備、ASURA 側の光子束の計算は計測範囲外です。 [PlotHIIRegion.ipynb](plot/PlotHIIRegion.ipynb)

### Cloudy との比較

[Cloudy C23.01](https://data.nublado.org/cloudy_releases/c23/) をビルドし、EQCLib のルートで次を実行します。

```sh
make -C test GSL_DIR="$HOME/lib/gsl" hii
python3 test/run_cloudy_hii.py --cloudy /path/to/c23.01/source/cloudy.exe
```

[run_cloudy_hii.py](test/run_cloudy_hii.py) は EQCLib の `equilibrium.csv` から密度、電離パラメータ、元素組成を読み取ります。PEGASE2 の連続スペクトルを EQCLib と同じ水素電離光子数で年齢ごとに規格化し、0〜30 Myr の18時点を等重みで平均します。Cloudy には case B、ダストと分子なし、厚さ 1 cm の1ゾーンを指定します。密度、入射電離光子束、光学的厚さ、熱収支を確認し、冷却の出力から O II と O III の正味の衝突線冷却を集計します。EQCLib の比較量には H II 領域で使われる温度重みを含めています。

各ゾーンの入力、ログ、集計 CSV は `test/out/hii/cloudy/` に出力します。[PlotHIIRegion.ipynb](plot/PlotHIIRegion.ipynb) は再計算の出力を優先し、出力がなければ [Cloudy C23.01 の保存済み参照結果](test/reference/cloudy_hii_c23.01.csv) を使います。参照結果で描画する場合は EQCLib の `make -C test hii` の出力があれば実行できます。Cloudy のデータ更新と計算方法は [C23 リリース論文](https://doi.org/10.22201/ia.01851101p.2023.59.02.12) と [C23.01 リリース論文](https://doi.org/10.3847/2515-5172/ad0e75) を参照してください。対応する PDF は [reference/](reference/index.tsv) に保存しています。

## H₂・CO 化学の検証計算と描画

```sh
make -C test GSL_DIR="$HOME/lib/gsl" co
```

[MolecularTest.c](test/MolecularTest.c) は、H II 領域・H₂・CO・ダスト・宇宙線・軟X線・H–O 電荷交換をすべて有効にした設定で化学平衡を解き、反応収支を確かめます。次の2種類の条件を計算します。

- 密度・温度・電離光子束の格子：`n_H = 1〜10⁵ cm⁻³`（0.25 dex 刻み）、`T = 30, 100, 300, 1000 K`、電離光子なしと `log10(U) = −5, −3` の252条件です。`G0 = 1.7`、全柱密度は `N_H = 2×10²¹ (n_H/100)^{1/2} cm⁻²` です。結果は `test/out/co/equilibrium.csv` に出力します。
- 柱密度走査：電離光子なし、`T = 30 K`、`n_H = 100, 1000 cm⁻³` と `G0 = 1.7, 17` の組で、`N_H = 10¹⁹〜10²³ cm⁻²` を変えます。各点は一つ前の解から始めます。結果は `test/out/co/column_scan.csv` に出力します。

各平衡状態について、H₂ と CO の生成・消滅、電荷交換を含む H と O の電離・再結合の反応量を、公開されている反応率関数から独立に計算し直します。どれかの相対差が `10⁻⁶` を超えるか、最終的な H₂ 量で CO を解き直したときに CO が `10⁻⁶` 以上変わると、試験は失敗します。CSV には存在比、各反応量（cm⁻³ s⁻¹）、相対差を出力します。

Jupyter で [PlotCOChemistry.ipynb](plot/PlotCOChemistry.ipynb) を実行すると、反応収支、密度に対する分子の割合、柱密度に対する H/H₂ と C⁺/C/CO の遷移、生成と消滅の反応量を表示します。PNG を生成するには `make -C test GSL_DIR="$HOME/lib/gsl" co-plot` を実行してください。[plot_co_chemistry.py](plot/plot_co_chemistry.py) が PNG を `plot/out/co/` に保存します。

### CO モデルの適用範囲

CO 化学は Keto & Caselli (2008) の平衡モデル（[Glover & Clark 2012](https://doi.org/10.1111/j.1365-2966.2011.20260.x) の KC08e）です。反応は C⁺ + H₂ → CO、CO の光解離（CO → C）、C の光電離（C → C⁺）の3つだけで、宇宙線による反応と C⁺ の再結合を含みません。平衡では `x_C/x_CO = Γ_CO/Γ_C` となり、ダスト減光の指数が CO（2.5）と C（3.0）で異なるため、この比は `A_V` とともに `e^{0.5 A_V}` で増えます。

- `A_V ≲ 20`：CO は KC08e の想定どおり振る舞います。柱密度走査では `A_V ≈ 3〜10` で `CO/C > 0.99` です。
- `A_V ≳ 25`：CO が減り、C I が主になります。走査では `CO/C` が `A_V = 34` で 0.05 です。実際の高密度コアでは CO が主成分として残るので、この領域の C・CO の存在比と冷却は信頼できません。
- 原子状炭素は全体に少なめです。Glover & Clark (2012) は、KC08 が C⁺ の再結合を含まないことを理由に、観測より小さい C/CO を報告しています。

Glover & Clark (2012) は、銀河規模の計算では KC08e のような簡単なモデルで十分であり、モデル間の差は分解能による不確かさより小さいと結論しています。ASURA などで使う場合は、粒子の `A_V` が 20 を超える割合を確認してください。

宇宙線による C の電離と He⁺ + CO → C⁺ + O + He を KC08e に加えるだけでは解決しません。走査の反応率を使った試算では、深部の `CO/C` は 0.8〜0.98 に保たれますが、C⁺ を除く反応が遅い C⁺ + H₂ だけなので、`A_V = 3〜50` で C⁺ が炭素の 2〜20% 残ります。実際の暗黒雲の `C⁺/C ≈ 10⁻⁴` より大きく、[C II] 冷却も過大になります。深部まで扱うには、C⁺ の再結合や CHx・OHx を含む Nelson & Langer (1999) 型のモデル（Glover & Clark 2012 の NL99）が必要です。

## Wolfire et al. (1995) の検証計算と描画

```sh
make -C test GSL_DIR="$HOME/lib/gsl" start
mkdir -p test/out/wolfire
(cd test/out/wolfire && ../../test.out >wolfire.stdout 2>wolfire.log)
```

出力ファイルは `test/out/wolfire/out/` に作成されます。求解に失敗した点の詳細は `wolfire.log` に記録され、失敗点があるとコマンドは非ゼロで終了します。

Jupyter で [PlotWolfire1995.ipynb](plot/PlotWolfire1995.ipynb) を実行すると、求解できたデータ点が表示されます。冷却・加熱率の図は水素核数密度に対する体積あたりの各成分と合計を表示します。既定は基準条件で、`rate_case` で計算条件を切り替えられます。保存済みの実行を描く場合は、最初のコードセルの `run_dir` を変更してください。旧形式のデータでは記録済みの成分が表示され、追加成分を得るには上の手順で再計算します。描画には Python 3、NumPy、Matplotlib 3.7 以降が必要です。

PNG を生成するには `make -C test wolfire-plot` を実行してください。[plot_wolfire1995.py](plot/plot_wolfire1995.py) が `test/out/wolfire` の結果を描画し、`thermal_rates.png` などを `plot/out/wolfire/` に保存します。2つのノートブックとスクリプトは、[plot_style.py](plot/plot_style.py) の共通の図のスタイルを使います。

## 参考文献

### 主要論文

| 用途 | 文献 |
| --- | --- |
| 中性原子ガスの検証、宇宙線・軟X線 | [Wolfire et al. (1995), *The Neutral Atomic Phases of the Interstellar Medium*, ApJ 443, 152–168](https://doi.org/10.1086/175510) |
| 原子の電離・再結合・冷却 | [Cen (1992)](https://doi.org/10.1086/191630)、[Voronov (1997)](https://doi.org/10.1006/adnd.1997.0732)、[Hui & Gnedin (1997)](https://arxiv.org/abs/astro-ph/9612232) |
| 化学反応係数と [C II]・[O I]・[Si II] の衝突率 | [Glover & Jappsen (2007)](https://arxiv.org/abs/0705.0182) |
| ダストの光電加熱と粒子表面での再結合 | [Bakes & Tielens (1994)](https://ui.adsabs.harvard.edu/abs/1994ApJ...427..822B/abstract)、[Weingartner & Draine (2001)](https://arxiv.org/abs/astro-ph/0105237) |
| H₂ の自己遮蔽と回転振動冷却 | [Draine & Bertoldi (1996)](https://arxiv.org/abs/astro-ph/9603032)、[Glover & Abel (2008)](https://arxiv.org/abs/0803.1768) |
| CO の化学反応・簡略化モデル | [Nelson & Langer (1997)](https://doi.org/10.1086/304167)、[Keto & Caselli (2008)](https://lweb.cfa.harvard.edu/~eketo/StarFormation/lowmass/keto_caselli_2008.pdf)、[Glover et al. (2010)](https://arxiv.org/abs/0907.4081)、[Glover & Clark (2012)](https://doi.org/10.1111/j.1365-2966.2011.20260.x) |
| CO の遮蔽表と回転振動冷却 | [Lee et al. (1996)](https://ui.adsabs.harvard.edu/abs/1996A%26A...311..690L/abstract)、[Hollenbach & McKee (1989)](https://ui.adsabs.harvard.edu/abs/1989ApJ...342..306H/abstract) |
| 高温ガスの衝突電離平衡（CIE）冷却表 | [Sutherland & Dopita (1993)](https://doi.org/10.1086/191823) |
| H II 領域の [O II]・[O III] 輝線 | O II の遷移確率は [Zeippen (1982)](https://adsabs.harvard.edu/abs/1982MNRAS.198..111Z)、有効衝突強度は [Kisielius et al. (2009)](https://doi.org/10.1111/j.1365-2966.2009.14989.x)。O III の遷移確率は [Froese Fischer & Tachiev (2004)](https://doi.org/10.1016/j.adt.2004.02.001) と [Storey & Zeippen (2000)](https://doi.org/10.1046/j.1365-8711.2000.03184.x)、有効衝突強度は [Storey et al. (2014)](https://doi.org/10.1093/mnras/stu777) |
| 原子データの選定と準位エネルギー | [Morisset et al. (2020), PyNeb](https://doi.org/10.3390/atoms8040066)、[NIST Atomic Spectra Database](https://doi.org/10.18434/T4W30F) |
| H II 領域の O–H 電荷交換 | O²⁺ + H⁰ は [Kingdon & Ferland (1996)](https://ui.adsabs.harvard.edu/abs/1996ApJS..106..205K/abstract)、O⁰ + H⁺ ⇄ O⁺ + H⁰ は [Millar et al. (1997)](https://ui.adsabs.harvard.edu/abs/1997A%26AS..121..139M/abstract) |
| Hα 放射率 | [Peters et al. (2017)](https://doi.org/10.1093/mnras/stw3216) |

上表の原典 PDF、Cloudy の比較用リリース論文、既存モデルの全参考文献は [EQCLib/reference/index.tsv](reference/index.tsv) に取得元とともに整理しています。

Nelson & Langer (1997) は出版社の PDF 直リンクが自動取得時に HTTP 403 を返すため、`reference/Nelson_Langer_1997.pdf` に公開リーダー経由の本文保持版を収録しています。組版画像が必要な場合は、索引に記載した DOI の出版社ページを参照してください。

### 教科書

- Draine (2011), *Physics of the Interstellar and Intergalactic Medium* — 星間ガスの電離・再結合、ダスト、輝線過程。 [Princeton University Press](https://press.princeton.edu/books/paperback/9780691122144/physics-of-the-interstellar-and-intergalactic-medium)
- Tielens (2005), *The Physics and Chemistry of the Interstellar Medium* — 分子化学、粒子表面反応、光解離領域。 [Cambridge University Press](https://www.cambridge.org/core/books/physics-and-chemistry-of-the-interstellar-medium/B71008A57D96A2875E851DAEA70B16AE)
