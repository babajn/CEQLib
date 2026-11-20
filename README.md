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

## Wolfire の計算と描画

```sh
make -C test GSL_DIR="$HOME/lib/gsl" start
mkdir -p test/out/wolfire
(cd test/out/wolfire && ../../test.out >wolfire.stdout 2>wolfire.log)
```

出力ファイルは `test/out/wolfire/out/` に作成されます。求解に失敗した点の詳細は `wolfire.log` に記録され、失敗点があるとコマンドは非ゼロで終了します。

NumPy、pandas、Matplotlib、IPython が入った Jupyter 環境で [PlotWolfire1995.ipynb](plot/PlotWolfire1995.ipynb) を実行すると、求解できたデータ点が表示されます。保存済みの実行を描く場合は、最初のコードセルの `run_dir` を変更してください。
