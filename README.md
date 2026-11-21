# EQCLib (Equilibrium Cooling)

星間ガスの化学平衡と冷却・加熱率を計算する C ライブラリです。固定温度での化学平衡、化学平衡を仮定したエネルギー方程式の陰的更新、固定密度での水素・ヘリウム・炭素・酸素の非平衡化学と熱の時間発展に対応します。

## 開発方針

EQCLib は化学平衡モードを既定とします。その時点の温度・密度・照射条件で化学平衡の組成を求め、その組成で評価した冷却・加熱率を使ってガスの内部エネルギーを時間発展させます。化学反応の時間積分にかかる計算負荷を抑え、銀河スケールの流体計算で粒子ごとに繰り返し使える速度と、ガスの熱進化を再現するのに十分な冷却精度の両立を目指します。

より多くの化学種と反応を時間発展させる CHIMES を比較対象とし、無照射の単一セルと、遠紫外線・遮蔽を含む温度・密度グリッドで検証します。冷却・加熱率と温度履歴を比べてモデル間の差がどの条件で大きくなるかを調べ、1セル更新あたりの CPU 時間を測って計算負荷を評価します。 [CHIMES との比較条件と結果](test/chimes/README.md)

## ディレクトリ構成

- `src/` に C ライブラリと単一セルの実行コード、`examples/` に設定例があります。
- `data/` に計算入力、`reference/` に文献索引があります。
- `test/` に検証項目ごとのドライバーと参照データがあり、計算結果は `test/out/` に出力します。実行コマンドと作図ノートブックは [テスト一覧](test/README.md) から参照できます。
- `plot/` に作図スクリプトとノートブックがあります。

## ビルドと検証

C99 対応の GCC と Steffen 補間（`gsl_interp_steffen`）に対応した GNU Scientific Library（GSL）が必要です。以下の `GSL_DIR` は GSL のインストール先に合わせてください。`make -C test benchmark` は OpenMP を使用します。

```sh
make -C src GSL_DIR="$HOME/lib/gsl" libEQCLib.a
make -C test GSL_DIR="$HOME/lib/gsl" check
OMP_NUM_THREADS=12 make -C test GSL_DIR="$HOME/lib/gsl" benchmark
```

ベンチマークは、原子ガスと分子ガスそれぞれの325点を、10年・1万年・100万年の時間刻みで3ステップずつ進めます。各ステップの収束数、実行時間、逐次と並列の比較を表示します。

### 検証ノートブックの再計算

リポジトリのルートで次を実行すると、H₂・CO、H II 領域、非平衡化学の検証データと計算時間を `test/out/` に出力します。計算時間はターゲットを順番に実行して測ります。

```sh
make -C test GSL_DIR="$HOME/lib/gsl" co hii noneq
make -j1 -C test GSL_DIR="$HOME/lib/gsl" molecular-settings-benchmark hii-benchmark noneq-benchmark noneq-hot
```

Wolfire の検証と Cloudy の比較は個別の手順で実行します。計算後、Jupyter で各ノートブックのセルを先頭から実行し、保存すると図と集計表が更新されます。 [Wolfire (1995) の実行手順](#wolfire-et-al-1995-の検証計算と描画)、[Wolfire (2003) の実行手順](#wolfire-et-al-2003-との熱平衡曲線の比較)、[Cloudy の実行手順](#cloudy-との比較)

## 設定ファイルで平衡・非平衡を切り替える

固定密度・固定照射の単一セルは、付属の `eqclib-cell` で実行できます。リポジトリのルートで次を実行してください。

```sh
make -C src GSL_DIR="$HOME/lib/gsl" eqclib-cell
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

`chemistry_mode` は `eqclib-cell` の起動時に読み込みます。非平衡モードでは高温側の金属 CIE 冷却を有効にし、セルの C と O の数密度を使います。追跡しない N・Ne・Mg・S・Ca・Fe の分は、O/H を Asplund et al. (2009) の太陽値と比べて `Z = 0.0134 × (O/H) / 10⁻³·³¹` とした金属量で補います。平衡モードでは金属 CIE 冷却を無効にします。`make -C test cell-check` で両モードの初期組成と極短ステップでの温度を検証できます。 [非平衡化学と熱の連成](#固定密度での非平衡化学と熱の連成)

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

H II 領域の酸素冷却を使うには、初期化前に `EQCLibSetRunParamHIIRegion(true)` を呼び、各粒子の `NLyn` に局所の水素電離光子束を設定します。光電離ガスの CIE 金属冷却を [O II] と [O III] の冷却に徐々に置き換え、15,000〜25,000 K の間で CIE に戻します。各遷移の光学的に薄い放射率は [HIIregionLineCooling.h](src/HIIregionLineCooling.h) の `CoolingLinesOII` と `CoolingLinesOIII` で取得できます。

この設定では、O⁰ + H⁺ ⇄ O⁺ + H⁰ と O²⁺ + H⁰ → O⁺ + H⁺ の電荷交換も解きます。H の電離分率を1変数とし、O の電離平衡と H I・H II の反応平衡を、電荷交換を含めて同時に満たします。H₂ 生成を有効にしている場合は、H の電離分率を変えるたびに H₂ の生成と消滅の釣り合いも解き直すので、H₂、H、O の反応がすべて同時に釣り合います。

### 固定温度・密度での水素の非平衡計算

`EQCLibTryAdvanceHydrogenChemistry(&state, 1e-8)` は、単一セルの水素化学を `state.dt` 秒進めます。`T` と `n[ipHYDROGEN]` に温度（K）と水素核数密度を、`nion[ipHYDROGEN]` と `nmol[ipH2]` に初期の H⁺・H₂ 数密度（cm⁻³）を設定してください。水素以外の気相化学種はゼロにします。H₂ を扱う場合は `EQCLibSetRunParamH2Formation(true, 1.0)` を呼び、粒子表面での生成には `Zd` を設定します。電離光 `NLyn` を使う場合は、初期化前に H II 領域の設定を有効にします。

各呼び出しは温度・密度・外部照射を固定して積分し、H₂ の自己遮蔽をその時点の H₂ 存在量から計算します。`EQC_SUCCESS` が返ると、H I・H⁺・H₂・電子の数密度が更新されます。照射の切り替えは `NLyn` や `G0` を変更して次のステップを呼び出します。`make -C test check` では照射切り替え後の応答と、一定照射下での既存の化学平衡解への収束を検証します。

ダスト上の H₂ 生成には Glover & Jappsen (2007) の反応60を使います。体積生成率は `k₆₀ × n_H × n(H I)` で、`k₆₀` に掛かる係数は `3.0×10⁻¹⁸ √(T/K) cm³ s⁻¹` です。粘着・熱脱離の効率、ダスト対ガス比の太陽値に対する比、`EQCLibSetRunParamH2Formation` のクランピング係数 `Cp` が掛かります。 [Glover & Jappsen (2007), Table 3](https://arxiv.org/pdf/0705.0182#page=47)

H₂ の光解離には、Draine 型スペクトルを仮定した遮蔽前の率 `3.3×10⁻¹¹ G0 s⁻¹` と、断面積 `2.0×10⁻²¹ cm²/H` をダスト量でスケールした減光を使います。自己遮蔽は Draine & Bertoldi (1996) の式37を柱密度ゼロで1になるよう規格化します。熱ドップラー幅には H₂ の質量を使い、乱流ドップラー幅 3 km/s と二乗和で合成します。平衡更新では、入力の H₂ 量を初期値とし、遮蔽式の微分を使って安定な生成・消失収支を探します。 [Glover & Mac Low (2007), §2.2.1](https://arxiv.org/pdf/astro-ph/0605120#page=13) · [Draine & Bertoldi (1996), 式37](https://arxiv.org/pdf/astro-ph/9603032#page=12)

### 固定密度での非平衡化学と熱の連成

`EQCLibTryAdvanceHydrogenThermochemistry(&state, 1e-8)` は、水素・ヘリウム・炭素・酸素の化学反応と比内部エネルギー `U`（erg/g）を連成します。正の水素核数密度、`U`、気相質量密度 `Rho`（g/cm³）を設定してください。He・C・O については、対応する `n` と `nion` に核数密度と一価イオン数密度を、`nHeIII`、`nCIII`、`nOIII` に二価イオン数密度を設定します（すべて cm⁻³）。炭素は C I〜C III を扱います。中性原子と電子の数密度は元素・電荷の保存から、初期温度は `U` と初期組成から求めます。密度と外部照射は各ステップ中に固定します。

Si を含む場合は `n[ipSILICON]`、`nion[ipSILICON]`、`nSiIII` に Si 核、Si II、Si III の数密度（cm⁻³）を設定します。Si II と Si III の数密度は非負で、合計が核数密度以下になるようにします。両イオンの数密度をステップ中も固定し、Si I は核数密度から両イオンの数密度を差し引いて求めます。この組成を電子密度、状態方程式、[Si II] 35 µm 冷却に反映します。Si の光電離は起こらないものとして扱い、Si I・Si II の光電離加熱は加熱率に含めません。

非平衡の熱化学更新では CO 生成を無効にして `nmol[ipCO] = 0` とします。CO 生成が有効、または入力の CO 数密度がゼロ以外の状態で正の時間刻みを指定すると、共通 API と熱化学 API は `EQC_INVALID_INPUT` を返し、入力状態を保持します。CO の化学と冷却を扱う場合は平衡モードを使います。

積分中は、その時点の組成で状態方程式と冷却・加熱率を評価し、`dU/dt = (H_tot - L_tot) / Rho` で熱を進めます。He I・He II の光電離加熱と、He の励起・電離・再結合などの冷却には現在の組成と電子密度を使います。H₂ の回転振動冷却で衝突相手になる中性ヘリウムにも He I の現在値を使います。既定で H–O 電荷交換の3反応を含めます。消灯後も、残っている C⁺ から [C II] 158 µm、O⁺・O²⁺ から [O II]・[O III] 冷却をその時点の温度・組成で評価します。

高温側の金属冷却を補うには `EQCLibSetRunParamCIEResonanceCooling(true)` を呼び、粒子の全金属の質量分率 `Z` に正の値を設定します。`Z` は Asplund et al. (2009) の太陽値 `Z☉ = 0.0134` を基準とする全金属量として扱います。ダストの `Zd` は従来どおり `0.019` を基準とします。ライブラリの CIE 設定の既定値は `false` です。CIE を無効にした場合や `Z = 0` の場合、`L_res` はゼロです。非平衡モードではその場合、全温度で [O II] と [O III] の輝線冷却をそのまま使います。

ホストが物理的な全金属質量分率を持っている場合は、その値を `Z` に渡します。旧基準 `0.019` で太陽比を表している場合は、同じ太陽比になるよう `Z = 0.0134 × (Z_old / 0.019)` を渡します。たとえば旧基準の `Z_old = 0.019` は `Z = 0.0134` に対応します。C・O・Si の数密度とダスト対ガス質量比 `Zd` は従来どおり個別に設定します。

CIE が有効で `Z > 0` の場合、Gnat & Ferland (2012) の元素別係数 `Λ_X(T)` を使います。C・O・Si は、電子数密度と `n[]` の各元素の全数密度を掛けて足し合わせ、他の冷却過程と同じ元素組成を使います。EQCLib が追跡しない N・Ne・Mg・S・Ca・Fe は、Asplund et al. (2009) の太陽組成（ダストへの取り込みなし）で重み付けした係数の和に `n_e n_H Z / 0.0134` を掛けて加えます。この補完項は 10⁴ K 以下でゼロ、10⁵ K 以上で全量とし、その間は `s = log10(T / 10⁴ K)` として `s²(3 − 2s)` で滑らかに変えます。10⁴ K 以下の冷却率は C・O・Si だけの場合と同じです。H・He の冷却はそれぞれの冷却過程で計算します。 [元素別 CIE 係数の出典](https://doi.org/10.1088/0067-0049/199/1/20)

非平衡モードで CIE が有効かつ `Z > 0` のとき、酸素輝線冷却の重みは 15,000 K 以下で1、25,000 K 以上で0です。その間は `s = (25,000 K - T) / 10,000 K` として `s²(3 - 2s)` で滑らかに変えます。CIE 冷却率には、酸素輝線の重みを1から引いた値を掛けて `L_res` に入れます。消灯直後の温度と組成はそのままなので、輝線冷却も連続し、その後は再結合と温度の変化に応じて変わります。

補間範囲は 10⁴〜10⁸ K です。平衡モードでは対数温度に対して係数の対数を線形補間し、非平衡モードでは各区間で単調な Steffen 補間を使います。Steffen 補間では格子点でも係数と傾きが連続し、積分器が細かいステップを取り直す回数を抑えます。10⁴ K 未満では端点の係数に `exp(100 × (log10(T/K) − 4))` を掛け、10⁸ K より上は 10⁹ K まで `√T` で外挿します。 [CIE の実装](src/CoolingCIE.c)

非平衡モードの温度による重みは切り替え区間内の光電離ガスにも適用します。高温の金属冷却には衝突電離平衡の係数を使い、酸素の時間発展は O I〜O III を追います。

組み込みの表は `python3 data/GnatFerland2012/generate_cie_table.py > src/CoolingTableGnatFerland.h` で再生成できます。初期化時に `Check=1` を指定すると、実行ディレクトリに `EQCooling/` が作成されます。その中の `CIECooling.dat` に対数温度、C・O・Si の各係数、6元素の太陽組成での補完係数（重みを掛ける前）が出力されます。現行設定の冷却率と冷却履歴は、遮蔽ありのグリッドで CHIMES と比較しています。 [CHIMES 比較ノートブック](plot/PlotShieldedChimes.ipynb)

`make -C test check` は、17,600 K で温度と組成を固定して電離光を切り替えたときに、酸素輝線と CIE の冷却率が連続することと、照射下で熱化学更新が進むことを確認します。 [熱化学のテスト](test/unit/HydrogenChemistryTest.c)

外部照射のない冷却ガスでも、再結合の遅れにより非平衡冷却率が CIE の値より低くなる温度域があります。 [Gnat & Sternberg (2007)](https://arxiv.org/abs/astro-ph/0608181)

C I の光電離には、`NLyn` が有効な場合は恒星スペクトルを使い、それ以外では遠紫外線 `G0` を使います。両方の照射を止める場合は `NLyn` と `G0` をゼロにします。[C II] の正味冷却では、宇宙マイクロ波背景の平均比強度を放射エネルギー密度に換算して吸収・誘導放射を評価します。成功時は API 内で `U`、温度、組成、状態方程式の諸量と、ステップ末の瞬時冷却・加熱率を更新します。`make -C test check` で照射切り替え後の炭素・ヘリウムの電離・再結合、温度、輝線冷却の応答を検証できます。

積分には GSL の陰的 BDF 法を使います。光電離係数と H₂ のダスト遮蔽は更新の開始時に一度求め、積分途中では温度や組成に依存する率を更新します。初期イオンがなく、電離源もゼロと判定できるセルでは、積分中も電子密度を厳密にゼロに保ちます。H₂ 生成がオフの場合は H₂ を積分変数から外します。Hα と H I 21 cm の診断量はステップ末の状態で計算します。

両 API ともエラー時は入力状態を保持します。軟 X 線を有効にすると、既存の軟 X 線電離率が電子分率 0.5 で不連続なため、この付近で積分が停止して `EQC_NO_CONVERGENCE` を返す場合があります。

非平衡更新では既定で温度係数を表引きします。`EQCLibSetRunParamNonEquilibriumTemperatureTables(false)` で直接評価を選べます。初期化時に 1〜10⁸ K の対数温度格子（1 dex あたり1,024区間）で係数を計算します。補間には隣接4点を使い、式の切り替え点をまたぐ区間、格子外、ゼロへのアンダーフロー付近、区間内の検査点で相対差が `5×10⁻⁹` を超える区間は直接評価します。再結合の case A と case B はそれぞれの係数を用意します。

対象は H・He・C・O の電子衝突電離と再結合、H–O 電荷交換、H₂ の衝突解離、H・He と H₂ の冷却係数、H₂ 加熱の臨界密度、[C II]・[O I] の衝突係数、[O II]・[O III] の励起に使うボルツマン因子です。電子密度や各化学種の存在量による重み、H₂ の自己遮蔽、準位分布などは積分途中の状態から求めます。

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

陰的エネルギー更新の許容誤差は `1e-4` を目安にします。1 万年刻みの単一スレッド測定では、H₂ と H–O 電荷交換を有効にしたまま CO を無効にすると、更新1回あたりの CPU 時間の中央値は 453.3 µs から 126.3 µs に短くなりました。 [設定比較ノートブック](plot/PlotMolecularSettings.ipynb)

この設定例は、アプリケーションが初期化前に指定する銀河計算向けの推奨設定です。C API の初期値では H₂ と CIE は無効、化学平衡許容誤差は `1e-8` です。CHIMES の標準比較では上の設定を土台に、軟 X 線と局所電離光束をゼロにします。 [比較条件](test/chimes/shielded/README.md)

CO 輝線を予測する場合は CO を有効にします。遮蔽された低温ガスの熱収支が重要な場合は、代表的な粒子で CO の有無による冷却率の差を確認してください。 [CO モデルの適用範囲](#co-モデルの適用範囲)

### 各設定の詳細

| 設定 | 既定値 | 設定例 | 効果 |
| --- | --- | --- | --- |
| `EQCLibSetRunParamCOFormation` | `false` | `true, 1.0` | 平衡モードで CO の化学平衡と輝線冷却を有効にします。第2引数は CO 生成率の倍率です。 |
| `EQCLibSetRunParamChargeTransfer` | `true` | `false` | H II 領域の更新が約 20〜35% 速くなります。物理モデルが変わり、H I が残る条件で O III が過剰になります。Cloudy との比較では、太陽金属量の平衡温度の差が中央値 1.4% から 5.3% に広がります。H II 領域の設定が無効なときは効果がありません。 |
| `EQCLibSetRunParamChemistryAccuracy` | `1e-8` | `1e-5` | 陰的なエネルギー更新の中で解く化学平衡の許容誤差の上限です。実際には `min(値, 0.1 × エネルギーの精度)` を使います。更新が約 10〜25% 速くなり、温度の変化は 2×10⁻⁵ 程度以下です。正味冷却率の API は内部で精度 1e-8 を指定するので、この設定にかかわらず 1e-9 以下で解きます。 |

H₂ の化学平衡は Newton 法で解き、化学平衡の許容誤差が `1e-12` 以下の場合や Newton 法が収束しない場合は区間探索を使います。

CO 化学を使う場合は、`EQCLibSetRunParamH2Formation(true, 1.0)` と `EQCLibSetRunParamCOFormation(true, 1.0)` を設定します。CO を無効にするには `EQCLibSetRunParamCOFormation(false, 1.0)` を呼びます。無効にすると、更新後の CO 存在量と CO 輝線冷却はゼロになり、炭素・酸素の化学状態や冷却率にも影響します。[CO モデルの適用範囲](#co-モデルの適用範囲)

CO を無効にした `make hii-benchmark` と同じ条件で、電荷交換を無効にして化学平衡の許容誤差を `1e-5` にすると、エネルギー更新1回あたりの時間は既定値の 0.6〜0.8 倍になりました。電荷交換を加える前の版と比べても 0.3〜0.6 倍です。測定値は単一スレッドの例で、計算機によって変わります。

### H₂・CO・H–O 電荷交換の設定比較

```sh
make -C test molecular-settings-plot
```

「H₂ と CO」「H₂ のみ」「分子なし」のそれぞれを、H–O 電荷交換の有効・無効と組み合わせて比較します。全設定で H II モード、ダスト、宇宙線、軟 X 線、CIE 冷却を有効にします。太陽金属量で、水素核数密度は `10⁻²〜10⁴ cm⁻³`（0.5 dex 刻み）、初期温度は `30, 100, 300, 1000, 3000, 8000, 15000 K`、電離光子束は照射なしと `log10(U) = −5, −3` です。全柱密度は `N_H = 10²⁰ sqrt(n_H / cm⁻³) cm⁻²`、`G0 = 1.7`、`X0 = 1` です。各設定で同じ初期温度の化学平衡を作り、その組成から内部エネルギーを設定してから、10 年・1 万年・100 万年をそれぞれ独立に1ステップ更新します。更新時の許容誤差は化学平衡が `1e-5`、エネルギー方程式が `1e-4` です。

`make -C test molecular-settings-benchmark` は条件別の結果と収束ステータスを `test/out/co/molecular_settings.csv` に保存します。計算時間は単一スレッドの CPU 時間です。全設定で収束した共通条件について、5回の測定の平均時間の中央値と範囲を図に表示します。温度・内部エネルギー変化・化学組成の差は、比較する両設定で収束した条件を対応づけて集計します。データの図と集計表は `plot/out/co/settings/` に保存します。描画済みの時間刻みを変える場合は、リポジトリのルートで `python3 plot/plot_molecular_settings.py --dt 1e6` のように指定します。

保存した比較では、全6設定で819条件中818条件が収束しました。初期温度 30 K、水素核密度 0.01 cm⁻³、`log10(U) = −5`、時間刻み100万年の1条件で平衡解が得られず、比較集計から外れました。

Jupyter で計算を実行して設定ごとの収束数と失敗した条件を確認できます。 [RunMolecularSettingsBenchmark.ipynb](test/co/RunMolecularSettingsBenchmark.ipynb)

描画用ノートブックは同じ CSV から計算時間と結果の差を表示します。結果の図に使う時間刻みは `dt_years` で切り替えられます。 [PlotMolecularSettings.ipynb](plot/PlotMolecularSettings.ipynb)

## H II 領域の検証計算と描画

```sh
make -C test GSL_DIR="$HOME/lib/gsl" hii
```

[HIIRegionTest.c](test/hii/HIIRegionTest.c) は、ダストと分子を含まない単一セルを、太陽金属量の PEGASE2 電離スペクトルと case B 再結合で計算します。入射電離光の減衰、非電離 FUV、宇宙線、X 線は含めません。水素核密度 `n_H = 1, 10, 100, 1000 cm⁻³`、ガス金属量 `Z/Z☉ = 0, 0.1, 1`、電離パラメータ `U = NLyn/(n_H c)` の `log10(U) = −4〜−2` の108条件で、等積の熱平衡温度と各元素の電離分率を求めます。さらに `n_H = 100 cm⁻³`、`Z/Z☉ = 1`、`log10(U) = −3` で約3,160〜31,600 K の温度範囲を走査します。

- `test/out/hii/equilibrium.csv` と `thermal_curve.csv` に電離分率と過程別の冷却・加熱率（erg cm⁻³ s⁻¹）を出力します。元素を含まない場合、その電離分率は `nan` になります。
- `test/out/hii/reactions.csv` に H、He、C、O、Si の光電離、電子衝突電離、電子再結合の反応率（cm⁻³ s⁻¹）を出力します。各点で電離と再結合の釣り合いを確認します。
- Jupyter で [PlotHIIRegion.ipynb](plot/PlotHIIRegion.ipynb) を実行すると、平衡温度・電離分率、冷却・加熱の内訳、反応率を表示します。保存済みの計算結果を描く場合は、最初のコードセルの `run_dir` を変更してください。描画には Python 3、NumPy、Matplotlib 3.7 以降が必要です。

PNG を生成するには `make -C test GSL_DIR="$HOME/lib/gsl" hii-plot` を実行してください。計算済みの CSV は [plot_hii_region.py](plot/plot_hii_region.py) でも描画できます。PNG は `plot/out/hii/` に保存されます。

### H II 領域の計算時間

```sh
make -C test GSL_DIR="$HOME/lib/gsl" hii-benchmark
```

[HIIRegionBenchmark.c](test/hii/HIIRegionBenchmark.c) は、電離光なし、照射継続中、照射開始直後を同じ密度・元素組成・初期温度で比較します。照射開始直後の粒子には、照射前の化学組成と内部エネルギーを保持したまま、照射継続中と同じ `NLyn` を与えます。ダストと分子なしの108条件で、固定温度の化学平衡・冷却加熱率と、時間刻み10年・1万年・100万年の陰的エネルギー更新を測ります。

`test/out/hii/hii_benchmark.csv` に、各条件群を繰り返した7ブロックの単一スレッド CPU 時間（1呼び出しあたりの平均 µs）を出力します。測定区間はセルの更新処理です。同じ Xeon w5-3423 の1コアで採った測定例も CSV で参照できます。 [探索の変更前](test/hii/reference/hii_benchmark_before_xeon_w5_3423.csv)、[変更後](test/hii/reference/hii_benchmark_xeon_w5_3423.csv)

### Cloudy との比較

[Cloudy C23.01](https://data.nublado.org/cloudy_releases/c23/) をビルドし、EQCLib のルートで次を実行します。

```sh
make -C test GSL_DIR="$HOME/lib/gsl" hii
python3 test/hii/cloudy/run_cloudy_hii.py --cloudy /path/to/c23.01/source/cloudy.exe
```

[run_cloudy_hii.py](test/hii/cloudy/run_cloudy_hii.py) は EQCLib の `equilibrium.csv` から密度、電離パラメータ、元素組成を読み取ります。PEGASE2 の連続スペクトルを EQCLib と同じ水素電離光子数で年齢ごとに規格化し、0〜30 Myr の18時点を等重みで平均します。Cloudy には case B、ダストと分子なし、厚さ 1 cm の1ゾーンを指定します。密度、入射電離光子束、光学的厚さ、熱収支を確認し、冷却の出力から O II と O III の正味の衝突線冷却を集計します。EQCLib の比較量には H II 領域で使われる温度重みを含めています。

各ゾーンの入力、ログ、集計 CSV は `test/out/hii/cloudy/` に出力します。[PlotHIIRegion.ipynb](plot/PlotHIIRegion.ipynb) は再計算の出力を優先し、出力がなければ [Cloudy C23.01 の保存済み参照結果](test/hii/cloudy/cloudy_hii_c23.01.csv) を使います。参照結果で描画する場合は EQCLib の `make -C test hii` の出力があれば実行できます。Cloudy のデータ更新と計算方法は [C23 リリース論文](https://doi.org/10.22201/ia.01851101p.2023.59.02.12) と [C23.01 リリース論文](https://doi.org/10.3847/2515-5172/ad0e75) を参照してください。対応する PDF は [reference/](reference/index.tsv) に保存しています。

## H₂・CO 化学の検証計算と描画

```sh
make -C test GSL_DIR="$HOME/lib/gsl" co
```

[MolecularTest.c](test/co/MolecularTest.c) は、H II 領域・H₂・CO・ダスト・宇宙線・軟X線・H–O 電荷交換をすべて有効にした設定で化学平衡を解き、反応収支を確かめます。次の2種類の条件を計算します。

- 密度・温度・電離光子束の格子：`n_H = 1〜10⁵ cm⁻³`（0.25 dex 刻み）、`T = 30, 100, 300, 1000 K`、電離光子なしと `log10(U) = −5, −3` の252条件です。`G0 = 1.7`、全柱密度は `N_H = 2×10²¹ (n_H/100)^{1/2} cm⁻²` です。結果は `test/out/co/equilibrium.csv` に出力します。
- 柱密度走査では電離光子を入れず、`T = 30 K`、`n_H = 100, 1000 cm⁻³` と `G0 = 1.7, 17` の組で `N_H = 10¹⁹〜10²³ cm⁻²` の164点を調べます。各点は一つ前の解から始めます。結果は `test/out/co/column_scan.csv` に出力します。

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

出力ファイルは `test/out/wolfire/out/` に作成されます。求解に失敗した点の詳細は `wolfire.log` に記録されます。失敗点があるとコマンドは非ゼロで終了し、成功した点は出力に残ります。ノートブックに保存した実行では3,300点中3,291点を描画し、柱密度 `10¹⁸ cm⁻²` の低密度側9点で平衡解が得られませんでした。

ノートブックの失敗点の表は、実行ディレクトリにある `rejected.tsv` を読み込みます。このファイルは保存した実行のログを集計したものです。同じディレクトリで再計算する場合は、以前の `rejected.tsv` を別の場所に保管し、今回の失敗点は新しい `wolfire.log` で確認してください。

Jupyter で [PlotWolfire1995.ipynb](plot/PlotWolfire1995.ipynb) を実行すると、求解できたデータ点が表示されます。冷却・加熱率の図は水素核数密度に対する体積あたりの各成分と合計を表示します。既定は基準条件で、`rate_case` で計算条件を切り替えられます。保存済みの実行を描く場合は、最初のコードセルの `run_dir` を変更してください。旧形式のデータでは記録済みの成分が表示され、追加成分を得るには上の手順で再計算します。描画には Python 3、NumPy、Matplotlib 3.7 以降が必要です。

PNG を生成するには `make -C test wolfire-plot` を実行してください。`test/out/wolfire` の計算結果または同梱の [参照データ](test/wolfire1995/reference/) から `thermal_rates.png` などを生成し、`plot/out/wolfire/` に保存します。 [plot_wolfire1995.py](plot/plot_wolfire1995.py)

## Wolfire et al. (2003) との熱平衡曲線の比較

```sh
make -C test/wolfire2003 GSL_DIR="$HOME/lib/gsl" run
```

太陽円の基準条件で化学・熱平衡を求め、冷相と暖相の安定枝および二相共存圧力を論文の Figure 8 と Table 3 に照合します。比較用の実行ファイルは GNU リンカの `--wrap` で論文の電離・加熱率に近似的にそろえたモデルと、EQCLib の既定の処方で計算した対照を出力します。計算結果は `test/out/wolfire2003/equilibrium.csv` に保存されます。ノートブックの実行には Jupyter、NumPy、pandas、SciPy、Matplotlib が必要です。 [PlotWolfire2003.ipynb](plot/PlotWolfire2003.ipynb) · [Wolfire et al. (2003)](https://doi.org/10.1086/368016)

## 非平衡化学の検証計算と描画

```sh
make -j1 -C test noneq-plot
```

[NonEquilibriumTest.c](test/noneq/NonEquilibriumTest.c) は、固定密度での照射切り替えと数値積分の精度を検証します。`make -C test noneq` で時間発展と冷却率の比較を、`make -j1 -C test noneq-benchmark` で計算時間の比較を実行できます。CSV は `test/out/noneq/` に保存されます。

- `evolution.csv` は、恒星電離光と遠紫外線の照射・消灯後の時間発展です。恒星電離光のセルは初期温度 8,000 K、`nH = 10 cm⁻³`、ダストなしで、水素の光電離率を `10⁻⁷ s⁻¹` として約16年照射し、消灯後100万年まで追います。遠紫外線のセルは初期温度 100 K、`nH = 100 cm⁻³`、`Nt = 2×10²⁰ cm⁻²`、初期の H₂ 水素核分率 0.8、ダストありで、`G0 = 1000` を30万年照射し、消灯後3,000万年まで追います。両セルの He/H、C/H、O/H 数比はそれぞれ `0.1`、`3×10⁻⁴`、`5×10⁻⁴` で、宇宙線と軟 X 線はオフです。切り替え直前と直後には、同じ `time_yr` に組成が同じ2行を記録します。`phase_time_yr` が照射開始または消灯からの経過時間です。`H2` 列は水素核分率 `2 n(H₂)/nH`、`H_`・`L_` 列は erg cm⁻³ s⁻¹ です。
- `recombination.csv` は、遮蔽と外部照射のない純水素を固定温度1,000 K、`nH = 10 cm⁻³` で再結合させ、case B の解析解 `x_HII = 1/(1 + t/t_rec)` と比較します。`t_rec = 1/(α_B nH)`、初期の `x_HII = 1` です。
- `convergence.csv` は H・He・C・O の混在セル270条件で、同じ時間幅の1ステップと8分割を比較します。`nH = 0.03, 1, 10³ cm⁻³`、初期温度 `30, 100, 8000, 10⁵, 10⁶ K`、時間幅 `100, 2×10⁴, 10⁶ yr`、照射なし・遠紫外線・恒星電離光と遠紫外線、H₂ のオン・オフを組み合わせます。宇宙線はオン、軟 X 線はオフです。温度・比内部エネルギーの相対差と、元素ごとに規格化した化学種分率の最大絶対差を出力します。
- `evolution_direct.csv`、`recombination_direct.csv`、`convergence_direct.csv` は直接評価で同じ計算をした結果です。接尾辞のない CSV は表引きの結果です。時間発展の488時点で温度と全冷却率の相対差、化学種・電子分率の絶対差を評価し、270条件の時間刻み比較でも両モードの到達状態を突き合わせます。
- `benchmark.csv` は、同じ初期セルから時間刻み `10, 10³, 10⁵ yr` でエネルギーを1回更新する場合の実時間を記録します。`equilibrium` は平衡、`nonequilibrium` は表引き、`nonequilibrium_direct` は非平衡の直接評価です。点灯時の中性セルと遠紫外線セルは上記と同条件です。消灯後の再結合セルは `nH = 10 cm⁻³`、初期温度 13,500 K、H⁺ 分率 0.98 で、He・C・O の中性・一価・二価の分率はそれぞれ 0.10、0.20、0.70 です。

冷却率の比較では、非平衡計算の各時点で温度・密度・元素存在量・柱密度・照射を固定して化学平衡を解きます。`evolution.csv` の `L_tot` は時間発展した組成、`eq_L_tot` と `eq_L_*` は平衡組成を共通の冷却式で評価した率です。`eq_native_L_tot` は同じ平衡組成を既存の平衡モードの式で評価した全冷却率です。既存式は [O II]・[O III] を光電離の寄与で重み付けし、[C II]・[O I] の背景放射や H₂ と He の衝突の扱いにも違いがあります。両モードで金属 CIE テーブルと CO はオフです。

[NonEquilibriumBenchmark.c](test/noneq/NonEquilibriumBenchmark.c) は各モード固有のエネルギー更新 API を使います。両方に精度引数 `10⁻⁶` を渡します。平衡側はエネルギー方程式の収束精度、非平衡側は時間積分の許容誤差に使われます。初期化とウォームアップを済ませた後、1スレッドで毎回同じ初期状態に戻して、更新時間を各条件・モードで5バッチ計ります。図の点は `elapsed_s / calls` の中央値、エラーバーは最小値から最大値までです。測定値は時間刻み、照射、初期組成、実行環境に依存します。

Jupyter の [PlotNonEquilibrium.ipynb](plot/PlotNonEquilibrium.ipynb) では、照射への応答、平衡・非平衡の冷却率と計算時間、解析解と時間刻みの検証を図と表で確認できます。再計算時の非平衡モードには既定の温度係数の表引きを使います。既定の出力先で必要なデータが不足していればテスト計算を実行します。再計算するときは最初のコードセルの `rerun` を `True` にし、保存済みの別の結果を描くときは `data_dir` を変更してください。

[plot_noneq.py](plot/plot_noneq.py) は CSV から時間発展図 `ionizing_pulse.png`・`fuv_pulse.png`、数値検証図 `validation.png`、平衡との冷却率比較図 `ionizing_cooling_comparison.png`・`fuv_cooling_comparison.png`、計算時間図 `benchmark.png` を `plot/out/noneq/` に生成します。開発時の表引きの検証には、直接評価との時系列比較図 `ionizing_temperature_lookup.png`・`fuv_temperature_lookup.png` と計算時間図 `temperature_lookup_benchmark.png` を使えます。切り替え時刻や対数軸の表示範囲外の率は CSV に記録します。計算済みのデータを PNG に描き直す場合は `python3 plot/plot_noneq.py --data-dir test/out/noneq` を使います。

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

ノートブックでは現行実装の高温冷却履歴と消灯後の酸素輝線を表示します。開発時の実装比較には上記の CSV と、`plot_noneq.py` が出力する `hot_cooling.png`、`hot_benchmark.png`、`oxygen_after_switch.png` を使えます。 [PlotNonEquilibrium.ipynb](plot/PlotNonEquilibrium.ipynb)

### CHIMES との比較

無照射の単一セルと、FUV・遮蔽ありの温度・密度グリッドで EQCLib と CHIMES を比較します。`test/chimes/` の比較ドライバーが、既定では EQCLib の隣にある `Chimes/chimes/` の C ライブラリを呼び出し、`Chimes/chimes-data/` の反応データを使って計算します。`plot/` の比較ノートブックは `test/out/chimes/` に出力した CSV から図と集計を作ります。 [CHIMES 比較の一覧と実行手順](test/chimes/README.md)

### 高温でのダスト再結合冷却と光電加熱

ダスト再結合冷却 `L_DustRec` と光電加熱 `H_PE` には、ガス温度が `10⁴ K` 以下で Bakes & Tielens (1994) の式をそのまま使います。`10⁴〜10⁵ K` では対数温度上で重みを滑らかに 1 から 0 に下げ、`10⁵ K` 以上では両項をゼロにします。光電加熱の式が検証された温度域を保ち、2つの熱収支項の高温外挿を抑えるための近似です。 [Bakes & Tielens (1994)](https://ui.adsabs.harvard.edu/abs/1994ApJ...427..822B/abstract)

遠紫外線・遮蔽ありの格子で、修正前後の EQCLib と CHIMES full の冷却率および等積冷却履歴を比較できます。 [開発時のダスト冷却検証](test/chimes/shielded/PlotDustCooling.ipynb)

## 参考文献

### 主要論文

| 用途 | 文献 |
| --- | --- |
| 中性原子ガスの検証、宇宙線・軟X線 | [Wolfire et al. (1995), *The Neutral Atomic Phases of the Interstellar Medium*, ApJ 443, 152–168](https://doi.org/10.1086/175510) |
| 原子の電離・再結合・冷却 | [Cen (1992)](https://doi.org/10.1086/191630)、[Voronov (1997)](https://doi.org/10.1006/adnd.1997.0732)、[Hui & Gnedin (1997)](https://arxiv.org/abs/astro-ph/9612232) |
| 化学反応係数と [C II]・[O I]・[Si II] の衝突率 | [Glover & Jappsen (2007)](https://arxiv.org/abs/0705.0182) |
| ダストの光電加熱と粒子表面での再結合 | [Bakes & Tielens (1994)](https://ui.adsabs.harvard.edu/abs/1994ApJ...427..822B/abstract)、[Weingartner & Draine (2001)](https://arxiv.org/abs/astro-ph/0105237) |
| H₂ の自己遮蔽と回転振動冷却 | [Draine & Bertoldi (1996)](https://arxiv.org/abs/astro-ph/9603032)、[Glover & Abel (2008)](https://arxiv.org/abs/0803.1768) |
| CO の化学反応・簡略化モデル | [Nelson & Langer (1997)](https://doi.org/10.1086/304167)、[Nelson & Langer (1999)](https://doi.org/10.1086/307823)、[Keto & Caselli (2008)](https://lweb.cfa.harvard.edu/~eketo/StarFormation/lowmass/keto_caselli_2008.pdf)、[Glover et al. (2010)](https://arxiv.org/abs/0907.4081)、[Glover & Clark (2012)](https://doi.org/10.1111/j.1365-2966.2011.20260.x) |
| CO の遮蔽表と回転振動冷却 | [Lee et al. (1996)](https://ui.adsabs.harvard.edu/abs/1996A%26A...311..690L/abstract)、[Hollenbach & McKee (1989)](https://ui.adsabs.harvard.edu/abs/1989ApJ...342..306H/abstract) |
| 高温ガスの元素別衝突電離平衡（CIE）冷却表 | [Gnat & Ferland (2012)](https://doi.org/10.1088/0067-0049/199/1/20) |
| CIE で補う6元素の太陽組成と基準金属量 | [Asplund et al. (2009)](https://doi.org/10.1146/annurev.astro.46.060407.145222) |
| 非平衡 CIE 表の単調補間 | [Steffen (1990), “A simple method for monotonic interpolation in one dimension”](https://adsabs.harvard.edu/pdf/1990A%26A...239..443S) |
| H II 領域の [O II]・[O III] 輝線 | O II の遷移確率は [Zeippen (1982)](https://adsabs.harvard.edu/abs/1982MNRAS.198..111Z)、有効衝突強度は [Kisielius et al. (2009)](https://doi.org/10.1111/j.1365-2966.2009.14989.x)。O III の遷移確率は [Froese Fischer & Tachiev (2004)](https://doi.org/10.1016/j.adt.2004.02.001) と [Storey & Zeippen (2000)](https://doi.org/10.1046/j.1365-8711.2000.03184.x)、有効衝突強度は [Storey et al. (2014)](https://doi.org/10.1093/mnras/stu777) |
| 原子データの選定と準位エネルギー | [Morisset et al. (2020), PyNeb](https://doi.org/10.3390/atoms8040066)、[NIST Atomic Spectra Database](https://doi.org/10.18434/T4W30F) |
| H II 領域の O–H 電荷交換 | O²⁺ + H⁰ は [Kingdon & Ferland (1996)](https://ui.adsabs.harvard.edu/abs/1996ApJS..106..205K/abstract)、O⁰ + H⁺ ⇄ O⁺ + H⁰ は [Millar et al. (1997)](https://ui.adsabs.harvard.edu/abs/1997A%26AS..121..139M/abstract) |
| Hα 放射率 | [Peters et al. (2017)](https://doi.org/10.1093/mnras/stw3216) |

上表の原典 PDF、Cloudy の比較用リリース論文、既存モデルの全参考文献は [EQCLib/reference/index.tsv](reference/index.tsv) に取得元とともに整理しています。

### 教科書

- Draine (2011), *Physics of the Interstellar and Intergalactic Medium* — 星間ガスの電離・再結合、ダスト、輝線過程。 [Princeton University Press](https://press.princeton.edu/books/paperback/9780691122144/physics-of-the-interstellar-and-intergalactic-medium)
- Tielens (2005), *The Physics and Chemistry of the Interstellar Medium* — 分子化学、粒子表面反応、光解離領域。 [Cambridge University Press](https://www.cambridge.org/core/books/physics-and-chemistry-of-the-interstellar-medium/B71008A57D96A2875E851DAEA70B16AE)
