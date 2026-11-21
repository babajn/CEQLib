# EQCLib (Equilibrium Cooling)

星間ガスの化学平衡と冷却・加熱率を計算する C ライブラリです。固定温度での化学平衡、化学平衡を仮定したエネルギー方程式の陰的更新、固定密度での水素・ヘリウム・炭素・酸素の非平衡化学と熱の時間発展に対応します。

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

ソースを取得し、EQCLib のリポジトリルートで次を実行すると、静的ライブラリをビルドして基本テストを実行できます。

```sh
make -C src libEQCLib.a
make -C test check cell-check
```

### テスト結果を Jupyter で確認する

リポジトリのルートで次を実行します。`notebook-data` はライブラリのビルドと基本テストも行い、H II 領域、H₂・CO 化学、非平衡化学、Wolfire (2003) の作図データを `test/out/` に出力します。計算時間はターゲットを順番に実行して測ります。

```sh
make -C test notebook-data
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r plot/requirements.txt
jupyter lab
```

JupyterLab で `plot/PlotHIIRegion.ipynb` などを開き、仮想環境の Python 3 カーネルで **Run → Run All Cells** を選ぶと、保存済みの計算結果から図と集計表を更新できます。Wolfire (1995) と Cloudy の比較図には同梱の参照データも使えます。各比較を再計算するときは実行手順を参照してください。CHIMES の再計算には外部の CHIMES ライブラリと反応データが必要です。 [テストとノートブックの一覧](test/README.md) · [Wolfire (1995) の実行手順](test/README.md#wolfire-et-al-1995-の検証計算と描画) · [Cloudy の実行手順](test/README.md#cloudy-との比較) · [CHIMES の実行手順](test/chimes/README.md)

## 設定ファイルで平衡・非平衡を切り替える

固定密度・固定照射の単一セルは、付属の `eqclib-cell` で実行できます。リポジトリのルートで次を実行してください。

```sh
make -C src eqclib-cell
cp examples/cell.conf cell.conf
./src/eqclib-cell cell.conf > cell.csv
```

[設定例](examples/cell.conf) の `chemistry_mode` を `nonequilibrium` に変更すると非平衡、`equilibrium` にすると平衡です。項目の省略時は平衡を使います。ファイルは `名前 = 値` を1行ずつ記述し、`#` から行末まではコメントです。引数を省略すると下表の既定値で実行します。

| 項目 | 既定値 | 意味 |
| --- | ---: | --- |
| `nH` | 10 | 水素核数密度 [cm⁻³] |
| `T` | 8000 | 初期温度 [K] |
| `dt_years` | 1000 | 時間刻み [年] |
| `steps` | 10 | 更新回数 |
| `x_HII` | 0.5 | 非平衡での初期水素電離分率 `n(H⁺)/nH` |
| `x_H2` | 0 | 非平衡での初期 H₂ 水素核分率 `2 n(H₂)/nH` |
| `G0` | 0 | 遠紫外線強度 |
| `NLyn` | 0 | 水素電離光子束 [cm⁻² s⁻¹] |
| `Nt` | 0 | 水素柱密度 [cm⁻²] |
| `Zd` | 0 | ダスト対ガス質量比 |

平衡モードは指定された `T` で化学平衡を解き、その組成から初期の内部エネルギーを計算します。非平衡モードは `x_HII` と `x_H2` を初期値に使い、He・C・O は中性から始めます。He/H、C/H、O/H の数比は 0.1、`3×10⁻⁴`、`5×10⁻⁴` です。H₂ 生成は `Zd` が正の場合と、非平衡で初期 H₂ がある場合に有効になります。

`NLyn` が正の場合は、第2引数に PEGASE2 のディレクトリを指定します。リポジトリのルートからなら `./src/eqclib-cell cell.conf data/PEGASE2` です。標準出力の CSV には初期状態と各ステップ後の温度、組成、加熱・冷却率 [erg cm⁻³ s⁻¹] が入ります。

`chemistry_mode` は `eqclib-cell` の起動時に読み込みます。この実行例では、非平衡モードを選ぶと高温側の金属 CIE 冷却を有効にします。平衡モードでは無効です。 [単一セルの CIE 設定](docs/physics.md#元素別-cie-冷却と酸素輝線の切り替え)

## シミュレーションコードからの利用

公開 API は [ChemicalEquilibriumCooling.h](src/ChemicalEquilibriumCooling.h) にあります。ホストコードを `libEQCLib.a`、`libgsl`、`libgslcblas`、`libm` とリンクします。
粒子の状態は [EQCDataStructure.h](src/EQCDataStructure.h) の `StructEQCData` で渡し、ステータスを返す API では `EQC_SUCCESS` で成功を判定します。

### 共通 API での設定と時間更新

ASURA などのホストコードでは、起動時に設定を読み込み、粒子ごとの冷却更新を `EQCLibTryAdvanceCell` に接続します。この接続を一度行えば、以降はホストの設定ファイルの1項目でモードを選べます。このリポジトリの単一セル実行例 `eqclib-cell` も同じ更新 API を使います。ASURA への接続はホスト側で行います。

```ini
chemistry_mode = nonequilibrium
```

設定値は `equilibrium`（`EQC_CHEMISTRY_EQUILIBRIUM`）と `nonequilibrium`（`EQC_CHEMISTRY_NONEQUILIBRIUM`）です。項目を省略すると従来の平衡モードを使います。`EQCLibTryLoadChemistryMode` は `名前 = 値` 形式で `chemistry_mode` を読み取り、ホストの他の設定項目や `#` 以降のコメントは読み飛ばします。設定ファイルを使わない場合はパスに `NULL` を渡せます。パスを指定してファイルを開けない場合や、`chemistry_mode` の値が不正または重複している場合はエラーを返します。ホストに設定パーサーがある場合は、前後の空白を除いた値を `EQCLibTryParseChemistryMode` に渡せます。その場合も項目を省略できるように、モード変数を `EQC_CHEMISTRY_EQUILIBRIUM` で初期化します。

ライブラリの実行パラメータを設定し、並列計算を始める前に `EQCLibTryInitialization` を呼びます。H II 領域の計算には PEGASE2 の入力ディレクトリを指定できます。ホスト側でモードを保持し、各粒子の更新に渡します。

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

共通 API は両モードで `U`、温度、化学種などを更新し、`Du` にステップ平均の比内部エネルギー変化率を格納します。`U` の時間更新はこの呼び出しで完結します。`dt = 0` の場合とエラー時は入力状態を保持します。

同じ粒子データでモードを切り替える場合は、CO 生成を既定の無効設定にして `nmol[ipCO]` をゼロにします。Si を含む粒子は両モードで使えます。非平衡モードでの Si の入力と扱いは下記の節に示します。

平衡モードで指定温度から新規に計算を始める場合は、その温度で `EQCLibTryCalculateChemicalAbundanceCoolingHeatingRates` を呼び、成功を確認してから `EQCLibConvertTtoU` で平衡組成に対応する `U` を求めます。非平衡モードでは初期組成に対応する `U` を設定し、次のステップや再開計算に更新後の `U` と組成を引き継ぎます。モードは起動時に決めます。非平衡組成から平衡計算を再開する場合は、指定温度なら組成をそろえて `U` を再計算し、`U` を引き継ぐ場合は再平衡化で温度が変わり得ることを考慮してください。 [初期化の実装例](src/CellRunner.c)

従来の平衡エネルギー API `EQCLibTryEnergyEquationSolverChemcalEquilibriumCoolingHeating` を直接使うホストコードでは、引き続き成功を確認してから呼び出し側で `U += dt * Du` を適用します。

### H II 領域の平衡冷却

H II 領域の酸素冷却を使うには、初期化前に `EQCLibSetRunParamHIIRegion(true)` を呼び、各粒子の `NLyn` に局所の水素電離光子束 [cm⁻² s⁻¹] を設定します。各遷移の光学的に薄い放射率は `CoolingLinesOII` と `CoolingLinesOIII` で取得できます。 [輝線の API](src/HIIregionLineCooling.h) · [酸素輝線と CIE 冷却](docs/physics.md#金属輝線と高温-cie-冷却) · [H–O 電荷交換](docs/physics.md#ho-電荷交換)

### 固定温度・密度での水素の非平衡計算

`EQCLibTryAdvanceHydrogenChemistry(&state, 1e-8)` は、単一セルの水素化学を `state.dt` 秒進めます。`T` と `n[ipHYDROGEN]` に温度（K）と水素核数密度を、`nion[ipHYDROGEN]` と `nmol[ipH2]` に初期の H⁺・H₂ 数密度（cm⁻³）を設定してください。水素以外の気相化学種はゼロにします。H₂ を扱う場合は `EQCLibSetRunParamH2Formation(true, 1.0)` を呼び、粒子表面での生成には `Zd` を設定します。電離光 `NLyn` を使う場合は、初期化前に H II 領域の設定を有効にします。

各呼び出しは温度・密度・外部照射を固定して積分し、H₂ の自己遮蔽をその時点の H₂ 存在量から計算します。`EQC_SUCCESS` が返ると、H I・H⁺・H₂・電子の数密度が更新されます。照射の切り替えは `NLyn` や `G0` を変更して次のステップを呼び出します。 [H₂ の反応と遮蔽](docs/physics.md#h₂-化学と熱過程)

### 固定密度での非平衡化学と熱の連成

`EQCLibTryAdvanceHydrogenThermochemistry(&state, 1e-8)` は、水素・ヘリウム・炭素・酸素の化学反応と比内部エネルギー `U`（erg/g）を連成します。正の水素核数密度、`U`、気相質量密度 `Rho`（g/cm³）を設定してください。He・C・O については、対応する `n` と `nion` に核数密度と一価イオン数密度を、`nHeIII`、`nCIII`、`nOIII` に二価イオン数密度を設定します（すべて cm⁻³）。炭素は C I〜C III を扱います。中性原子と電子の数密度は元素・電荷の保存から、初期温度は `U` と初期組成から求めます。密度と外部照射は各ステップ中に固定します。

Si を含む場合は `n[ipSILICON]`、`nion[ipSILICON]`、`nSiIII` に Si 核、Si II、Si III の数密度（cm⁻³）を設定します。Si II と Si III の数密度は非負で、合計が核数密度以下になるようにします。両イオンの数密度をステップ中も固定し、Si I は核数密度から両イオンの数密度を差し引いて求めます。この組成を電子密度、状態方程式、[Si II] 35 µm 冷却に反映します。

非平衡の熱化学更新では CO 生成を無効にして `nmol[ipCO] = 0` とします。CO 生成が有効、または入力の CO 数密度がゼロ以外の状態で正の時間刻みを指定すると、共通 API と熱化学 API は `EQC_INVALID_INPUT` を返し、入力状態を保持します。CO の化学と冷却を扱う場合は平衡モードを使います。

高温の金属 CIE 冷却を使うには `EQCLibSetRunParamCIEResonanceCooling(true)` を呼び、粒子の `Z` に物理的な全金属質量分率を設定します。ライブラリでは既定で無効です。C・O・Si の数密度とダスト対ガス質量比 `Zd` は個別に設定します。 [CIE 冷却と金属量の基準](docs/physics.md#元素別-cie-冷却と酸素輝線の切り替え)

成功時は `U`、温度、組成、状態方程式の諸量、ステップ末の瞬時冷却・加熱率と診断量が更新されます。水素のみの化学更新と熱化学更新は、いずれもエラー時に入力状態を保持します。非平衡更新の係数表は既定で有効です。`EQCLibSetRunParamNonEquilibriumTemperatureTables(false)` で直接評価に切り替えられます。 [非平衡での温度係数表](docs/physics.md#非平衡での温度係数表) · [軟 X 線と積分の注意点](docs/physics.md#軟-x-線)

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

陰的エネルギー更新の許容誤差は `1e-4` を目安にします。 [設定比較の測定結果](test/README.md#h₂coho-電荷交換の設定比較)

この設定例は、アプリケーションが初期化前に指定する銀河計算向けの推奨設定です。C API の初期値では H₂ と CIE は無効、化学平衡許容誤差は `1e-8` です。

CO 輝線を予測する場合は CO を有効にします。遮蔽された低温ガスの熱収支が重要な場合は、代表的な粒子で CO の有無による冷却率の差を確認してください。 [CO モデルの適用範囲](docs/physics.md#co-モデルの適用範囲)

### 各設定の詳細

| 設定 | 既定値 | 設定例 | 効果 |
| --- | --- | --- | --- |
| `EQCLibSetRunParamCOFormation` | `false` | `true, 1.0` | 平衡モードで CO の化学平衡と輝線冷却を有効にします。第2引数は CO 生成率の倍率です。 |
| `EQCLibSetRunParamChargeTransfer` | `true` | `false` | 平衡では H II 領域を有効にしたときに使われます。非平衡では照射の有無にかかわらず使われます。 [H–O 電荷交換](docs/physics.md#ho-電荷交換) |
| `EQCLibSetRunParamChemistryAccuracy` | `1e-8` | `1e-5` | 陰的なエネルギー更新の中で解く化学平衡の許容誤差の上限です。実際には `min(値, 0.1 × エネルギーの精度)` を使います。正味冷却率の API は内部で精度 1e-8 を指定するので、この設定にかかわらず 1e-9 以下で解きます。 |

CO 化学を使う場合は、`EQCLibSetRunParamH2Formation(true, 1.0)` と `EQCLibSetRunParamCOFormation(true, 1.0)` を設定します。CO を無効にするには `EQCLibSetRunParamCOFormation(false, 1.0)` を呼びます。無効にすると、更新後の CO 存在量と CO 輝線冷却はゼロになり、炭素・酸素の化学状態や冷却率にも影響します。 [CO モデルの適用範囲](docs/physics.md#co-モデルの適用範囲)

## 参考文献

各過程の原典は物理モデルの該当節から参照できます。論文と取得元は文献索引にまとめています。 [物理モデル](docs/physics.md) · [文献索引](reference/index.tsv)
