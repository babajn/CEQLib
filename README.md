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

## Wolfire et al. (1995) の検証計算と描画

```sh
make -C test GSL_DIR="$HOME/lib/gsl" start
mkdir -p test/out/wolfire
(cd test/out/wolfire && ../../test.out >wolfire.stdout 2>wolfire.log)
```

出力ファイルは `test/out/wolfire/out/` に作成されます。求解に失敗した点の詳細は `wolfire.log` に記録され、失敗点があるとコマンドは非ゼロで終了します。

NumPy、pandas、Matplotlib、IPython が入った Jupyter 環境で [PlotWolfire1995.ipynb](plot/PlotWolfire1995.ipynb) を実行すると、求解できたデータ点が表示されます。保存済みの実行を描く場合は、最初のコードセルの `run_dir` を変更してください。

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

### 教科書

- Draine (2011), *Physics of the Interstellar and Intergalactic Medium* — 星間ガスの電離・再結合、ダスト、輝線過程。 [Princeton University Press](https://press.princeton.edu/books/paperback/9780691122144/physics-of-the-interstellar-and-intergalactic-medium)
- Tielens (2005), *The Physics and Chemistry of the Interstellar Medium* — 分子化学、粒子表面反応、光解離領域。 [Cambridge University Press](https://www.cambridge.org/core/books/physics-and-chemistry-of-the-interstellar-medium/B71008A57D96A2875E851DAEA70B16AE)
