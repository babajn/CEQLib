# EQCLib (Equilibrium Cooling)

星間ガスの化学平衡と冷却・加熱率を計算する C ライブラリです。固定温度での化学平衡、化学平衡を仮定したエネルギー方程式の陰的更新、固定密度での水素・ヘリウム・炭素・酸素の非平衡化学と熱の時間発展に対応します。宇宙紫外線背景（UVB）の光電離・加熱と、セルの両側からの遮蔽を扱えます。

EQCLib の開発者は Junichi Baba（馬場 淳一）です。

## 開発方針

EQCLib は化学平衡モードを既定とします。その時点の温度・密度・照射条件で化学平衡の組成を求め、その組成で評価した冷却・加熱率を使ってガスの内部エネルギーを時間発展させます。化学反応の時間積分にかかる計算負荷を抑え、銀河スケールの流体計算で粒子ごとに繰り返し使える速度と、ガスの熱進化を再現するのに十分な冷却精度の両立を目指します。

より多くの化学種と反応を時間発展させる CHIMES を比較対象とし、無照射の単一セルと、遠紫外線・遮蔽を含む温度・密度グリッドで検証します。冷却・加熱率と温度履歴を比べてモデル間の差がどの条件で大きくなるかを調べ、1セル更新あたりの CPU 時間を測って計算負荷を評価します。 [CHIMES との比較条件と結果](test/chimes/README.md)

## ディレクトリ構成

- `src/` に C ライブラリと単一セルの実行コード、`examples/` に設定例があります。
- `data/` に計算入力、`reference/` に文献索引があります。
- [物理モデル](docs/physics.md) に化学反応と冷却・加熱過程の率式、単位、適用条件を記載しています。
- `test/` に検証項目ごとのドライバーと参照データがあり、計算結果は `test/out/` に出力します。実行コマンドと作図ノートブックは [テスト一覧](test/README.md) から参照できます。
- `plot/` に作図スクリプトとノートブックがあります。

## ビルドと検証

ビルドには GNU Make、C99 対応の GCC、Steffen 補間（`gsl_interp_steffen`）に対応した GNU Scientific Library（GSL）の開発ファイルが必要です。テストと作図には Python 3 を使います。GSL のパスは `gsl-config` があれば自動で取得します。別の場所にインストールした GSL を指定する場合は、ビルドの前に `export GSL_DIR=/path/to/gsl` を実行してください。

ソースを取得し、EQCLib のリポジトリルートで静的ライブラリをビルドします。

```sh
make -C src libEQCLib.a
```

### インストール後の基本テスト（必須）

ビルド後は、リポジトリのルートで基本テストを実行します。GSL、Python 3、同梱の入力データで化学平衡・非平衡の計算と単一セルの実行を確認します。UVB の基本テストは実行時に生成する合成スペクトルを使います。 [基本テストの内容と成否の見方](test/README.md#基本テスト)

```sh
make -C test check cell-check
```

### テスト結果を Jupyter で確認する

`notebook-data` は基本テストも実行し、EQCLib 単体で H II 領域、H₂・CO 化学、非平衡化学、Wolfire (2003) の作図データを `test/out/` に生成します。H II 領域だけを確認する場合は `make -C test hii` で作図データを用意できます。

```sh
make -C test notebook-data
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r plot/requirements.txt
jupyter lab
```

JupyterLab で `plot/PlotHIIRegion.ipynb` などを開き、仮想環境の Python 3 カーネルで **Run → Run All Cells** を選びます。Wolfire (1995) の図と、H II 領域ノートブックの Cloudy 比較には同梱の参照結果を使えます。 [単体で実行できる検証とノートブックの一覧](test/README.md#eqclib-単体で実行できる検証と作図)

UVB の図を再計算する場合は、FG20 のスペクトルを用意し、`make -C test uvb-data` を実行してからノートブックを開きます。 [FG20 の取得方法](#fg20-のスペクトルを用意する) · [検証と作図の手順](test/README.md#fg20-の遮蔽と計算時間) · [PlotUVBackground.ipynb](plot/PlotUVBackground.ipynb)

### 外部コードとの比較（任意）

比較結果を再計算したい場合は、CHIMES のソースと反応データ、または Cloudy の実行ファイルを用意し、対象の比較手順を実行します。 [外部コードとの比較手順](test/README.md#外部コードとの比較任意)

### PEGASE.2 のスペクトルを指定する

`data/PEGASE2/` に、EQCLib の開発者が PEGASE.2 で生成した Kroupa IMF のスペクトル7本を同梱しています。光電離を含むテストは既定でこの7本を使います。PEGASE.2 のプログラムと生成に使った入力資料はローカルで管理します。

```text
spectra_Kroupa_Z0.0001.dat
spectra_Kroupa_Z0.0004.dat
spectra_Kroupa_Z0.0040.dat
spectra_Kroupa_Z0.0080.dat
spectra_Kroupa_Z0.0200.dat
spectra_Kroupa_Z0.0500.dat
spectra_Kroupa_Z0.1000.dat
```

別のスペクトルを使う場合は、同じ名前の7本を置いたディレクトリの絶対パスを `PEGASE_DIR` で指定します。

```sh
make -C test hii PEGASE_DIR="/absolute/path/to/pegase-spectra"
```

### FG20 のスペクトルを用意する

Faucher-Giguère (2020) の `fg20_102720_rescaled.zip` を取得し、展開した `fg20_cloudy.ascii` を `data/FG20/` または任意のローカルディレクトリに置きます。EQCLib はこのスペクトルから率を積分します。 [著者の配布ページ](https://galaxies.northwestern.edu/uvb-fg20/) · [公式アーカイブ](https://galaxies.northwestern.edu/files/2020/10/fg20_102720_rescaled.zip) · [Faucher-Giguère (2020)](https://arxiv.org/abs/1903.08657)

| ファイル | SHA-256 |
| --- | --- |
| `fg20_102720_rescaled.zip` | `587a27854de50c95d79901b514d9acee421119ca41ff483eaeea84b6d695d49d` |
| `fg20_cloudy.ascii` | `40026f69c90c1aa76038f93664e8f3dea47269cc61ce0a5bc7b0fb6aa056879d` |

## 設定ファイルで平衡・非平衡を切り替える

固定密度・固定照射の単一セルは、付属の `eqclib-cell` で実行できます。リポジトリのルートで次を実行してください。

```sh
make -C src eqclib-cell
cp examples/cell.conf cell.conf
./src/eqclib-cell cell.conf > cell.csv
```

[設定例](examples/cell.conf) の `chemistry_mode` で平衡（`equilibrium`）、非平衡（`nonequilibrium`）、H₂ のみ非平衡（`h2_nonequilibrium`）を選びます。項目の省略時は平衡を使います。ファイルは `名前 = 値` を1行ずつ記述し、`#` から行末まではコメントです。引数を省略すると下表の既定値で実行します。

| 項目 | 既定値 | 意味 |
| --- | ---: | --- |
| `nH` | 10 | 水素核数密度 [cm⁻³] |
| `T` | 8000 | 初期温度 [K] |
| `dt_years` | 1000 | 時間刻み [年] |
| `steps` | 10 | 更新回数 |
| `x_HII` | 0.5 | 非平衡での初期水素電離分率 `n(H⁺)/nH`。H₂ のみ非平衡では `1 - x_H2` を上限とする初期電荷探索の手掛かり |
| `x_H2` | 0 | 非平衡での初期 H₂ 水素核分率 `2 n(H₂)/nH` |
| `G0` | 0 | 遠紫外線強度 |
| `NLyn` | 0 | 水素電離光子束 [cm⁻² s⁻¹] |
| `Nt` | 0 | 水素柱密度 [cm⁻²] |
| `Zd` | 0 | EQCLib のダスト量パラメータ。銀河系並みのダスト量は `0.019` |
| `h2_cr_dissociation_rate` | 0 | H₂ 1分子あたりの宇宙線直接破壊率 [s⁻¹]。正の場合は宇宙線を有効にして `X0=1` を使う |
| `uvb` | `off` | `on` で宇宙紫外線背景を有効化 |
| `uvb_directory` | `../data/FG20` | `fg20_cloudy.ascii` のあるディレクトリ。相対パスは実行ディレクトリが基準 |
| `uvb_redshift` | `on` の場合は必須 | UVB の赤方偏移 |
| `uvb_column_top` | `on` の場合は必須 | 上側境界からセルまでの遮蔽柱密度 [cm⁻²] |
| `uvb_column_bottom` | `on` の場合は必須 | 下側境界からセルまでの遮蔽柱密度 [cm⁻²] |

平衡モードは指定された `T` で化学平衡を解き、その組成から初期の内部エネルギーを計算します。非平衡モードは `x_HII` と `x_H2` を初期値に使い、He・C・O は中性から始めます。H₂ のみ非平衡モードは初期分率 `x_H2` を固定して原子・イオンを平衡化します。He/H、C/H、O/H の数比は 0.1、`3×10⁻⁴`、`5×10⁻⁴` です。H₂ の反応は `Zd` が正の場合、非平衡で初期 H₂ がある場合、H₂ のみ非平衡モードで有効になります。

`NLyn` が正の場合は、第2引数にスペクトルのディレクトリを渡します。同梱データを使う場合は、リポジトリのルートから `./src/eqclib-cell cell.conf data/PEGASE2` と実行します。標準出力の CSV には初期状態と各ステップ後の温度、組成、加熱・冷却率 [erg cm⁻³ s⁻¹] が入ります。 [必要なスペクトル](#pegase2-のスペクトルを指定する)

`chemistry_mode` は `eqclib-cell` の起動時に読み込みます。この実行例では、`nonequilibrium` を選ぶと高温側の金属 CIE 冷却を有効にします。 [単一セルの CIE 設定](docs/physics.md#元素別-cie-冷却と酸素輝線の切り替え)

### 単一セルで UVB を有効にする

リポジトリのルートから実行する場合は、`cell.conf` の `uvb = off` を `on` に変更し、次のように設定します。FG20 の赤方偏移は `0 ≤ z ≤ 10` です。柱密度 `0` はその側が無遮蔽であることを表します。指定した赤方偏移と柱密度は、初期状態と全ステップで共通です。 [FG20 の取得方法](#fg20-のスペクトルを用意する) · [遮蔽モデル](docs/physics.md#宇宙紫外線背景)

```ini
uvb = on
uvb_directory = data/FG20
uvb_redshift = 2
uvb_column_top = 1e18
uvb_column_bottom = 1e18
```

`./src/eqclib-cell cell.conf > cell.csv` で実行します。恒星光も使う場合は `NLyn` を正にし、第2引数に PEGASE.2 のディレクトリを渡します。

## シミュレーションコードからの利用

公開 API は [ChemicalEquilibriumCooling.h](src/ChemicalEquilibriumCooling.h) にあります。ホストコードを `libEQCLib.a`、`libgsl`、`libgslcblas`、`libm` とリンクします。
粒子の状態は [EQCDataStructure.h](src/EQCDataStructure.h) の `StructEQCData` で渡し、ステータスを返す API では `EQC_SUCCESS` で成功を判定します。

### 共通 API での設定と時間更新

ASURA などのホストコードでは、起動時に設定を読み込み、粒子ごとの冷却更新を `EQCLibTryAdvanceCell` に接続します。この接続を一度行えば、以降はホストの設定ファイルの1項目でモードを選べます。このリポジトリの単一セル実行例 `eqclib-cell` も同じ更新 API を使います。ASURA への接続はホスト側で行います。

```ini
chemistry_mode = nonequilibrium
```

設定値は `equilibrium`（`EQC_CHEMISTRY_EQUILIBRIUM`）、`nonequilibrium`（`EQC_CHEMISTRY_NONEQUILIBRIUM`）、`h2_nonequilibrium`（`EQC_CHEMISTRY_H2_NONEQUILIBRIUM`）です。項目を省略すると従来の平衡モードを使います。`EQCLibTryLoadChemistryMode` は `名前 = 値` 形式で `chemistry_mode` を読み取り、ホストの他の設定項目や `#` 以降のコメントは読み飛ばします。設定ファイルを使わない場合はパスに `NULL` を渡せます。パスを指定してファイルを開けない場合や、`chemistry_mode` の値が不正または重複している場合はエラーを返します。ホストに設定パーサーがある場合は、前後の空白を除いた値を `EQCLibTryParseChemistryMode` に渡せます。その場合も項目を省略できるように、モード変数を `EQC_CHEMISTRY_EQUILIBRIUM` で初期化します。

`h2_nonequilibrium` は `U` と `nmol[ipH2]` をステップ間で引き継ぐ試験モードです。初回は正の温度の初期推定値を `T` に入れるか、状態方程式の値を設定します。H₂ 生成を有効にし、CO 生成を無効、初期 CO 存在量をゼロにします。所定の温度で初期状態を作るときは、`EQCLibTryCalculateFixedH2CoolingHeatingRates` で H₂ を保持したまま原子・イオン組成と加熱・冷却率を求めます。ホスト側で温度制限などにより `U` を補正した場合は、`EQCLibTryNetCoolingRateFixedH2` で補正後の `U` に整合する温度と率を求めます。 [H₂ の時間更新と適用条件](docs/physics.md#h₂-のみを非平衡で進めるモデル) · [試験モードの検証](test/README.md#h₂-のみ非平衡と宇宙線解離のテスト)

ライブラリの実行パラメータを設定し、並列計算を始める前に `EQCLibTryInitialization` を呼びます。H II 領域を有効にした場合は、第2引数に PEGASE.2 のスペクトルがあるディレクトリを渡します。ホスト側でモードを保持し、各粒子の更新に渡します。 [必要なスペクトル](#pegase2-のスペクトルを指定する)

```c
#include "ChemicalEquilibriumCooling.h"

/* 実行パラメータを設定した後、起動時に一度呼ぶ。 */
static EQCLibStatus InitializeCooling(const char *configPath, const char *pegaseDirectory,
                                      EQCLibChemistryMode *mode)
{
    EQCLibStatus status = EQCLibTryLoadChemistryMode(configPath, mode);
    if (status != EQC_SUCCESS) return status;
    return EQCLibTryInitialization(0, pegaseDirectory);
}

/* 初期化の成功を確認した後、粒子ごとに呼ぶ。cell->dt の単位は秒。 */
static EQCLibStatus AdvanceCooling(StructEQCData *cell, double accuracy,
                                   EQCLibChemistryMode mode)
{
    EQCLibStatus status = EQCLibTryAdvanceCell(cell, accuracy, mode);
    if (status != EQC_SUCCESS) return status;
    /* 更新後の cell->U、cell->T、化学種の数密度を粒子に反映する。 */
    return EQC_SUCCESS;
}
```

共通 API は各モードで `U`、温度、化学種などを更新し、`Du` にステップ平均の比内部エネルギー変化率を格納します。`U` の時間更新はこの呼び出しで完結します。`dt = 0` の場合とエラー時は入力状態を保持します。

同じ粒子データでモードを切り替える場合は、CO 生成を既定の無効設定にして `nmol[ipCO]` をゼロにします。`h2_nonequilibrium` では Si の平衡組成を求めます。`nonequilibrium` での Si の入力と扱いは下記の節に示します。

平衡モードで指定温度から新規に計算を始める場合は、その温度で `EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates` を呼び、成功を確認してから `EQCLibConvertTtoU` で平衡組成に対応する `U` を求めます。非平衡モードでは初期組成に対応する `U` を設定し、次のステップや再開計算に更新後の `U` と組成を引き継ぎます。モードは起動時に決めます。非平衡組成から平衡計算を再開する場合は、指定温度なら組成をそろえて `U` を再計算し、`U` を引き継ぐ場合は再平衡化で温度が変わり得ることを考慮してください。 [初期化の実装例](src/CellRunner.c)

従来の平衡エネルギー API `EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating` を直接使うホストコードでは、引き続き成功を確認してから呼び出し側で `U += dt * Du` を適用します。

### 宇宙紫外線背景と遮蔽

UVB は既定で無効です。すべての化学モードで有効にできます。UVB の追加で `StructEQCData` の ABI が変わるため、ライブラリとホストコードを同じヘッダーで再ビルドしてください。 [FG20 の取得方法](#fg20-のスペクトルを用意する) · [物理モデルと適用範囲](docs/physics.md#宇宙紫外線背景)

1. 初期化前に `EQCLibSetRunParamUVBackground(true, uvbDirectory)` を呼びます。ディレクトリに `NULL` を渡すと、実行ディレクトリからの相対パス `../data/FG20` を使います。光子エネルギーの上限を設ける場合は、初期化前に `EQCLibSetRunParamUVBackgroundMaximumEnergy(energyEV)` を呼びます。既定値の `0` ではスペクトルの上端まで積分します。
2. `EQCLibTryInitialization` の成功を確認します。初期化直後は表の最小赤方偏移を使います。赤方偏移を設定・更新するには、並列のセル更新に入る前に `EQCLibTrySetUVBackgroundRedshift(z)` を呼び、戻り値を確認します。背景の表は並列更新中には読み取り専用です。
3. 各セルの `NHIuvb[0]`、`NHIuvb[1]` に、上下それぞれの境界からセルまでの有限かつ非負の H I 柱密度 [cm⁻²] を渡します。無遮蔽の場合は両方に `0` を設定します。`cell->Z` にはダスト中の金属も含む有限かつ非負の全金属質量分率を渡します。太陽組成は `Z = 0.0134` です。
4. 通常どおり `EQCLibTryAdvanceCell` で更新します。UVB の加熱率は `H_UVB` に格納され、`H_tot` に含まれます。恒星光の光電離加熱は `H_PI` で確認できます。 [公開 API](src/ChemicalEquilibriumCooling.h) · [セルの入出力](src/EQCDataStructure.h)

片側の一次電離率と配分前の光電子パワーを診断するには `EQCLibTryGetUVBackgroundHemisphereRates`（`Z = 0`）または `EQCLibTryGetUVBackgroundHemisphereRatesWithMetallicity` を使います。後者は金属に由来する光電子パワーも返します。 [率の定義と単位](docs/physics.md#宇宙紫外線背景)

ASURA で方向別の柱密度を局所値で近似する場合は、`CalcNcol` で求めた全水素核柱密度 `localColumn` に、更新前の H I と H₂ に含まれる水素核の分率を掛けます。両側に同じ有効柱密度を渡す例です。 [遮蔽近似の適用範囲](docs/physics.md#uvb-の適用範囲) · [円盤ガスの検証](test/README.md#fg20-の遮蔽と計算時間)

```c
const double nH = cell->n[ipHYDROGEN];
const double xHI = cell->natom[ipHYDROGEN] / nH;
const double xH2 = cell->nmol[ipH2] / nH;
const double effectiveColumn = localColumn * (xHI + 2.0 * xH2);
cell->NHIuvb[0] = effectiveColumn;
cell->NHIuvb[1] = effectiveColumn;
EQCLibStatus status = EQCLibTryAdvanceCell(cell, accuracy, mode);
```

### H II 領域の平衡冷却

H II 領域の酸素冷却を使うには、初期化前に `EQCLibSetRunParamHIIRegion(true)` を呼び、各粒子の `NLyn` に局所の水素電離光子束 [cm⁻² s⁻¹] を設定します。各遷移の光学的に薄い放射率は `CoolingLinesOII` と `CoolingLinesOIII` で取得できます。 [輝線の API](src/HIIregionLineCooling.h) · [酸素輝線と CIE 冷却](docs/physics.md#金属輝線と高温-cie-冷却) · [H–O 電荷交換](docs/physics.md#ho-電荷交換)

### 固定温度・密度での水素の非平衡計算

`EQCLibTryAdvanceHydrogenChemistry(&state, 1e-8)` は、単一セルの水素化学を `state.dt` 秒進めます。`T` と `n[ipHYDROGEN]` に温度（K）と水素核数密度を、`nion[ipHYDROGEN]` と `nmol[ipH2]` に初期の H⁺・H₂ 数密度（cm⁻³）を設定してください。水素以外の気相化学種はゼロにします。H₂ を扱う場合は `EQCLibSetRunParamH2Formation(true, 1.0)` を呼び、粒子表面での生成には `Zd` を設定します。電離光 `NLyn` を使う場合は、初期化前に H II 領域の設定を有効にします。

各呼び出しは温度・密度・外部照射を固定して積分し、H₂ の自己遮蔽をその時点の H₂ 存在量から計算します。`EQC_SUCCESS` が返ると、H I・H⁺・H₂・電子の数密度が更新されます。照射の切り替えは `NLyn` や `G0` を変更して次のステップを呼び出します。 [H₂ の反応と遮蔽](docs/physics.md#h₂-化学と熱過程)

### 固定密度での非平衡化学と熱の連成

`EQCLibTryAdvanceHydrogenThermochemistry(&state, 1e-8)` は、水素・ヘリウム・炭素・酸素の化学反応と比内部エネルギー `U`（erg/g）を連成します。正の水素核数密度、`U`、気相質量密度 `Rho`（g/cm³）を設定してください。He・C・O については、対応する `n` と `nion` に核数密度と一価イオン数密度を、`nHeIII`、`nCIII`、`nOIII` に二価イオン数密度を設定します（すべて cm⁻³）。炭素は C I〜C III を扱います。中性原子と電子の数密度は元素・電荷の保存から、初期温度は `U` と初期組成から求めます。密度と外部照射は各ステップ中に固定します。

Si を含む場合は `n[ipSILICON]`、`nion[ipSILICON]`、`nSiIII` に Si 核、Si II、Si III の数密度（cm⁻³）を設定します。Si II と Si III の数密度は非負で、合計が核数密度以下になるようにします。両イオンの数密度をステップ中も固定し、Si I は核数密度から両イオンの数密度を差し引いて求めます。この組成を電子密度、状態方程式、[Si II] 35 µm 冷却に反映します。

非平衡の熱化学更新では CO 生成を無効にして `nmol[ipCO] = 0` とします。CO 生成が有効、または入力の CO 数密度がゼロ以外の状態で正の時間刻みを指定すると、共通 API と熱化学 API は `EQC_INVALID_INPUT` を返し、入力状態を保持します。CO の化学と冷却を扱う場合は平衡モードを使います。

高温の金属 CIE 冷却を使うには `EQCLibSetRunParamCIEResonanceCooling(true)` を呼び、粒子の `Z` に物理的な全金属質量分率を設定します。ライブラリでは既定で無効です。C・O・Si の気相核数密度とダスト量パラメータ `Zd` は個別に設定します。 [CIE 冷却と金属量の基準](docs/physics.md#元素別-cie-冷却と酸素輝線の切り替え)

成功時は `U`、温度、組成、状態方程式の諸量、ステップ末の瞬時冷却・加熱率と診断量が更新されます。水素のみの化学更新と熱化学更新は、いずれもエラー時に入力状態を保持します。非平衡更新の係数表は既定で有効です。`EQCLibSetRunParamNonEquilibriumTemperatureTables(false)` で直接評価に切り替えられます。 [非平衡での温度係数表](docs/physics.md#非平衡での温度係数表) · [軟 X 線と積分の注意点](docs/physics.md#軟-x-線)

### シミュレーションで速度を優先する設定

ダスト過程は `Zd / 0.019` を銀河系に対する相対ダスト量として使います。ホストコードが物理的なダスト対ガス質量比を持つ場合は、その比をホスト側の銀河系の基準値で割り、`0.019` を掛けて `Zd` に設定します。`Z` には物理的な全金属質量分率を渡します。 [ダストと気相元素の扱い](docs/physics.md#ダスト量と気相元素量)

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

陰的エネルギー更新の許容誤差は `1e-4` を目安にします。 [設定比較の測定結果](test/README.md#h₂coho-電荷交換の設定比較) · [ASURA での平衡冷却の実測](test/README.md#asura-での平衡冷却の実測)

この設定例は、アプリケーションが初期化前に指定する銀河計算向けの推奨設定です。C API の初期値では H₂ と CIE は無効、化学平衡許容誤差は `1e-8` です。

`useLocalIonizingFlux = false` の平衡モードでも、UVB と H–O 電荷交換を有効にした酸素を含むセルでは、H の平衡組成を使って O の電離に電荷交換を反映します。H II 領域モードでは H と O を双方向に解きます。非平衡では `EQCLibSetRunParamChargeTransfer` が電荷交換を直接制御します。 [UVB での H–O 電荷交換](docs/physics.md#ho-電荷交換)

FG20 と軟 X 線は個別に設定できます。両方を有効にすると銀河系外軟 X 線の寄与が一部重複します。 [UVB の適用範囲](docs/physics.md#uvb-の適用範囲) · [併用時の再現テスト](test/README.md#金属吸収光電離と実行時間の追加検証)

CO 輝線を予測する場合は CO を有効にします。遮蔽された低温ガスの熱収支が重要な場合は、代表的な粒子で CO の有無による冷却率の差を確認してください。 [CO モデルの適用範囲](docs/physics.md#co-モデルの適用範囲)

### 各設定の詳細

| 設定 | 既定値 | 設定例 | 効果 |
| --- | --- | --- | --- |
| `EQCLibSetRunParamCOFormation` | `false` | `true, 1.0` | 平衡モードで CO の化学平衡と輝線冷却を有効にします。第2引数は CO 生成率の倍率です。 |
| `EQCLibSetRunParamChargeTransfer` | `true` | `false` | 平衡では H II 領域モードで双方向に結合し、H II 領域モードの外では UVB が有効で O があるセルに一方向で適用します。非平衡では照射の有無にかかわらず使われます。 [H–O 電荷交換](docs/physics.md#ho-電荷交換) |
| `EQCLibSetRunParamChemistryAccuracy` | `1e-8` | `1e-5` | 陰的なエネルギー更新の中で解く化学平衡の許容誤差の上限です。実際には `min(値, 0.1 × エネルギーの精度)` を使います。正味冷却率の API は内部で精度 1e-8 を指定するので、この設定にかかわらず 1e-9 以下で解きます。 |

CO 化学を使う場合は、`EQCLibSetRunParamH2Formation(true, 1.0)` と `EQCLibSetRunParamCOFormation(true, 1.0)` を設定します。CO を無効にするには `EQCLibSetRunParamCOFormation(false, 1.0)` を呼びます。無効にすると、更新後の CO 存在量と CO 輝線冷却はゼロになり、炭素・酸素の化学状態や冷却率にも影響します。 [CO モデルの適用範囲](docs/physics.md#co-モデルの適用範囲)

## 参考文献

各過程の原典は物理モデルの該当節から参照できます。論文と取得元は文献索引にまとめています。 [物理モデル](docs/physics.md) · [文献索引](reference/index.tsv)
