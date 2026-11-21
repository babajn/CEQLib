# CEQLib

星間ガスの化学平衡と冷却・加熱率を計算する C ライブラリです。固定温度での化学平衡と、化学平衡を仮定したエネルギー方程式の陰的更新に対応します。

## ビルドと検証

C99 対応の GCC と GNU Scientific Library（GSL）が必要です。以下の `GSL_DIR` は GSL のインストール先に合わせてください。ベンチマークは OpenMP を使用します。

```sh
make -C src GSL_DIR="$HOME/lib/gsl" libCEQLib.a
make -C test GSL_DIR="$HOME/lib/gsl" check
OMP_NUM_THREADS=12 make -C test GSL_DIR="$HOME/lib/gsl" benchmark
```

ベンチマークは、原子ガスと分子ガスそれぞれの325点を、10年・1万年・100万年の時間刻みで3ステップずつ進めます。各ステップの収束数、実行時間、逐次と並列の比較を表示します。

## シミュレーションコードからの利用

公開 API は [ChemicalEquilibriumCooling.h](src/ChemicalEquilibriumCooling.h) にあります。ホストコードを `libCEQLib.a`、`libgsl`、`libgslcblas`、`libm` とリンクします。

実行パラメータを設定し、並列計算を始める前に `CEQLibTryInitialization` を呼び出してください。H II 領域の計算には PEGASE2 の入力ディレクトリを指定できます。エネルギー更新では `CEQLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating` の成功を確認してから、呼び出し側で `U += dt * Du` を適用します。

H II 領域の酸素冷却を使うには、初期化前に `CEQLibiSetRunParamHIIRegion(true)` を呼び、各粒子の `NLyn` に局所の水素電離光子束を設定します。光電離ガスの CIE 金属冷却を [O II] と [O III] の冷却に徐々に置き換え、15,000〜25,000 K の間で CIE に戻します。各遷移の光学的に薄い放射率は [HIIregionLineCooling.h](src/HIIregionLineCooling.h) の `CoolingLinesOII` と `CoolingLinesOIII` で取得できます。

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

[Cloudy C23.01](https://data.nublado.org/cloudy_releases/c23/) をビルドし、CEQLib のルートで次を実行します。

```sh
make -C test GSL_DIR="$HOME/lib/gsl" hii
python3 test/run_cloudy_hii.py --cloudy /path/to/c23.01/source/cloudy.exe
```

[run_cloudy_hii.py](test/run_cloudy_hii.py) は CEQLib の `equilibrium.csv` から密度、電離パラメータ、元素組成を読み取ります。PEGASE2 の連続スペクトルを CEQLib と同じ水素電離光子数で年齢ごとに規格化し、0〜30 Myr の18時点を等重みで平均します。Cloudy には case B、ダストと分子なし、厚さ 1 cm の1ゾーンを指定します。密度、入射電離光子束、光学的厚さ、熱収支を確認し、冷却の出力から O II と O III の正味の衝突線冷却を集計します。CEQLib の比較量には H II 領域で使われる温度重みを含めています。

各ゾーンの入力、ログ、集計 CSV は `test/out/hii/cloudy/` に出力します。[PlotHIIRegion.ipynb](plot/PlotHIIRegion.ipynb) は再計算の出力を優先し、出力がなければ [Cloudy C23.01 の保存済み参照結果](test/reference/cloudy_hii_c23.01.csv) を使います。参照結果で描画する場合は CEQLib の `make -C test hii` の出力があれば実行できます。Cloudy のデータ更新と計算方法は [C23 リリース論文](https://doi.org/10.22201/ia.01851101p.2023.59.02.12) と [C23.01 リリース論文](https://doi.org/10.3847/2515-5172/ad0e75) を参照してください。対応する PDF は [reference/](reference/index.tsv) に保存しています。

## Wolfire et al. (1995) の検証計算と描画

```sh
make -C test GSL_DIR="$HOME/lib/gsl" start
mkdir -p test/out/wolfire
(cd test/out/wolfire && ../../test.out >wolfire.stdout 2>wolfire.log)
```

出力ファイルは `test/out/wolfire/out/` に作成されます。求解に失敗した点の詳細は `wolfire.log` に記録され、失敗点があるとコマンドは非ゼロで終了します。

Jupyter で [PlotWolfire1995.ipynb](plot/PlotWolfire1995.ipynb) を実行すると、求解できたデータ点が表示されます。保存済みの実行を描く場合は、最初のコードセルの `run_dir` を変更してください。描画には Python 3、NumPy、Matplotlib 3.7 以降が必要です。

PNG を生成するには `make -C test wolfire-plot` を実行してください。[plot_wolfire1995.py](plot/plot_wolfire1995.py) が `test/out/wolfire` の結果を描画し、PNG を `plot/out/wolfire/` に保存します。2つのノートブックとスクリプトは、[plot_style.py](plot/plot_style.py) の共通の図のスタイルを使います。

## 参考文献

### 主要論文

| 用途 | 文献 |
| --- | --- |
| 中性原子ガスの検証、宇宙線・軟X線 | [Wolfire et al. (1995), *The Neutral Atomic Phases of the Interstellar Medium*, ApJ 443, 152–168](https://doi.org/10.1086/175510) |
| 原子の電離・再結合・冷却 | [Cen (1992)](https://doi.org/10.1086/191630)、[Voronov (1997)](https://doi.org/10.1006/adnd.1997.0732)、[Hui & Gnedin (1997)](https://arxiv.org/abs/astro-ph/9612232) |
| 化学反応係数と [C II]・[O I]・[Si II] の衝突率 | [Glover & Jappsen (2007)](https://arxiv.org/abs/0705.0182) |
| ダストの光電加熱と粒子表面での再結合 | [Bakes & Tielens (1994)](https://ui.adsabs.harvard.edu/abs/1994ApJ...427..822B/abstract)、[Weingartner & Draine (2001)](https://arxiv.org/abs/astro-ph/0105237) |
| H₂ の自己遮蔽と回転振動冷却 | [Draine & Bertoldi (1996)](https://arxiv.org/abs/astro-ph/9603032)、[Glover & Abel (2008)](https://arxiv.org/abs/0803.1768) |
| CO の化学反応・簡略化モデル | [Glover et al. (2010)](https://arxiv.org/abs/0907.4081)、[Glover & Clark (2012)](https://doi.org/10.1111/j.1365-2966.2011.20260.x) |
| CO の遮蔽表と回転振動冷却 | [Lee et al. (1996)](https://ui.adsabs.harvard.edu/abs/1996A%26A...311..690L/abstract)、[Hollenbach & McKee (1989)](https://ui.adsabs.harvard.edu/abs/1989ApJ...342..306H/abstract) |
| 高温ガスの衝突電離平衡（CIE）冷却表 | [Sutherland & Dopita (1993)](https://doi.org/10.1086/191823) |
| H II 領域の [O II]・[O III] 輝線 | O II の遷移確率は [Zeippen (1982)](https://adsabs.harvard.edu/abs/1982MNRAS.198..111Z)、有効衝突強度は [Kisielius et al. (2009)](https://doi.org/10.1111/j.1365-2966.2009.14989.x)。O III の遷移確率は [Froese Fischer & Tachiev (2004)](https://doi.org/10.1016/j.adt.2004.02.001) と [Storey & Zeippen (2000)](https://doi.org/10.1046/j.1365-8711.2000.03184.x)、有効衝突強度は [Storey et al. (2014)](https://doi.org/10.1093/mnras/stu777) |
| 原子データの選定と準位エネルギー | [Morisset et al. (2020), PyNeb](https://doi.org/10.3390/atoms8040066)、[NIST Atomic Spectra Database](https://doi.org/10.18434/T4W30F) |
| Hα 放射率 | [Peters et al. (2017)](https://doi.org/10.1093/mnras/stw3216) |

上表の原典 PDF、Cloudy の比較用リリース論文、既存モデルの全参考文献は [CEQLib/reference/index.tsv](reference/index.tsv) に取得元とともに整理しています。

### 教科書

- Draine (2011), *Physics of the Interstellar and Intergalactic Medium* — 星間ガスの電離・再結合、ダスト、輝線過程。 [Princeton University Press](https://press.princeton.edu/books/paperback/9780691122144/physics-of-the-interstellar-and-intergalactic-medium)
- Tielens (2005), *The Physics and Chemistry of the Interstellar Medium* — 分子化学、粒子表面反応、光解離領域。 [Cambridge University Press](https://www.cambridge.org/core/books/physics-and-chemistry-of-the-interstellar-medium/B71008A57D96A2875E851DAEA70B16AE)
