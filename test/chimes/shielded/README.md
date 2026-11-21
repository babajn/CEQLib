# FUV・遮蔽ありのグリッド比較

太陽金属量、`G0 = 1.7` Habing、水素の一次宇宙線電離率 `1.8e-17 s⁻¹`、遮蔽柱密度 `10²⁰ √(nH / cm⁻³) cm⁻²` の条件で比較します。CHIMES の `matched` は H・He・C・O・Si、`full` は全元素のネットワークです。どちらも B87 の EUV と Auger-EUV の断面積をゼロにします。固定温度の化学平衡には、B87 の EUV を含む `matched-euv` と `full-euv` の対照も含みます。

EQCLib の主系列は化学平衡モードで、H₂・ダスト遮蔽・ダスト熱過程・宇宙線・H–O 電荷交換・高温 CIE 冷却を有効にし、CO を無効にします。エネルギー更新の許容誤差は `1e-4`、化学平衡の許容誤差の上限は `1e-5` です。標準比較では軟 X 線と水素の光電離を無効にします。EQCLib は `NLyn = 0` とし、化学平衡で H–O 電荷交換を利用するために H II 領域のフラグを有効にします。EQCLib の CO 有効と非平衡モードも同じ現行版で比較します。 [銀河シミュレーション向けの設定](../../../README.md#シミュレーションで速度を優先する設定)

温度・密度格子での平衡、1セル更新の CPU 時間と到達温度、等積冷却履歴を比較します。銀河計算に近い時間刻みとして、粒子質量 10² と 10⁴ M☉ の SPH 粒子のクーラン条件 `dt = 0.3 h / c_s`（`h = (m/ρ)^{1/3}`、上限 1 Myr）での更新も測ります。時間更新では各コードが自身の平衡組成から出発します。CPU 時間は成功条件を照合し、合計時間の比と、条件ごとの比の中央値を集計します。

## 再実行

EQCLib のリポジトリルートから実行します。

```sh
make -C test/chimes/shielded run
```

個別のターゲットは `equilibrium`、`timing`、`courant`、`cooling` です。従来の `make -C test/chimes run` からも同じ比較を実行できます。`CHIMES_DIR` にはビルド済みの倍精度 `chimes/libchimes.so` と `chimes-data/` が必要です。GSL は `GSL_DIR`、HDF5 と SUNDIALS は `CHIMES_LOCAL`（既定値 `$CHIMES_DIR/.local`）で指定します。計測に使う CPU は `BENCH_CPU`（既定値 `2`）で指定します。

## 作図と出力

保存済みの CSV から図と集計を再生成できます。ノートブックで全セルを実行するか、PNG と Markdown を出力する場合は `make -C test/chimes/shielded plot` を実行します。作図には Python 3、NumPy、Matplotlib が必要です。 [PlotShieldedChimes.ipynb](../../../plot/PlotShieldedChimes.ipynb)

CSV、PNG、集計は `test/out/chimes/shielded/` に保存します。別の計算出力を使う場合は Make の `OUT` とノートブックの `data_dir` を指定できます。以前の EUV ありの計測データは、出力ディレクトリ内の `with_euv_run1/` にあります。

## 低温ガスの現行モデル比較

密度 3、10、30、100 cm⁻³ で、H₂ 分率、冷却率、熱平衡温度、5,000 K からの冷却履歴、クーラン刻みの更新時間を比較します。CHIMES matched と EQCLib で H・He・C・O・Si の数比、入射 FUV 強度、柱密度、宇宙線電離率、ダスト温度 10 K、乱流ドップラー幅 3 km/s をそろえます。EQCLib の CIE は同じ元素構成に制限します。同一温度・指定組成での [C II] 冷却の照合も含みます。 [PlotH2Comparison.ipynb](../../../plot/PlotH2Comparison.ipynb)

```sh
bash test/chimes/shielded/run_h2_prescription_comparison.sh --current-only
bash test/chimes/shielded/run_thermal_balance_comparison.sh --current-only
```

出力先は `test/out/chimes/shielded/h2_current/` と `test/out/chimes/shielded/thermal_balance_current/` です。

## 開発時の検証

### 高温のダスト冷却・加熱の比較

EQCLib の修正前後と CHIMES full で、固定温度の冷却・加熱率と等積冷却履歴を比較します。EQCLib の修正前は `bc72f90` のダスト実装を使い、両ビルドで同じ旧 CIE 表を使います。EQCLib と CHIMES は同じ密度・温度・照射条件から始め、化学組成はそれぞれの固定温度平衡から求めます。 [作図ノートブック](PlotDustCooling.ipynb)

```sh
bash test/chimes/shielded/run_dust_comparison.sh
```

CSV は `test/out/chimes/shielded/dust/` に出力します。ノートブックはこのディレクトリの CSV から結果を再描画します。別の出力先を使う場合はスクリプトの引数とノートブックの `data_dir` に同じパスを指定します。修正前後で `10⁴ K` 以下の固定温度の冷却・加熱率と、`10⁴ K` から始めた冷却履歴の出力を照合します。

### 化学平衡モードの収束検証

陰的エネルギー更新は、通常の化学平衡許容誤差と機械精度に近い許容誤差の両方で `EQC_NO_CONVERGENCE` になった場合、通常の値から1桁ずつ厳しくした中間の化学平衡許容誤差で再試行します。熱収支の残差には指定された許容誤差を使います。

```sh
bash test/chimes/shielded/run_convergence_comparison.sh
```

ダスト処方と旧 CIE 表を両ビルドでそろえ、修正前の `ThermalEquilibrium.c`（`bc72f90`）と現在の実装を比較します。固定温度格子、6条件の冷却履歴、固定刻みとクーラン刻みの計測結果を `test/out/chimes/shielded/convergence/` に出力します。速度は1コアに固定し、修正前後の実行順を交互に入れ替えて7回測ります。`BENCH_CPU` で使用コアを指定できます。 [収束と速度の比較ノートブック](PlotDustCooling.ipynb)

### CIE 冷却表の元素組成の比較

ダストの高温処方と収束修正を共通にし、Sutherland & Dopita の旧表、Gnat & Ferland の C・O・Si のみの表、6元素を補った現行表を比較します。旧表と C・O・Si のみの版には当時の `Z = 0.019` を、現行版には Asplund et al. (2009) の基準 `Z = 0.0134` を渡します。`Zd` と追跡元素の数比は全ビルドで共通です。

CHIMES `matched` と `full` は、EQCLib と同じ H・He・C・O・Si の数比を使います。`full` が追加する6元素には Wiersma et al. (2009) の組成を使います。両ネットワークとも EUV を切り、各コードの固定温度平衡から等積冷却を始め、6条件を10億年まで追います。

```sh
bash test/chimes/shielded/run_cie_comparison.sh
```

CSV は `test/out/chimes/shielded/cie_supplement/` に出力します。`before`、`after`、`supplemented` はそれぞれ旧表、C・O・Si のみの表、6元素を補った現行表です。速度は `BENCH_CPU` の1コアで EQCLib 3版の実行順を入れ替えて7回測り、CHIMES の各ネットワークは同じコアで1回測ります。ノートブックに冷却率と冷却履歴の図、および低温側と速度の集計があります。 [比較ノートブック](PlotDustCooling.ipynb)

### 低温ガスの H₂ 分率と冷却

水素核密度 3、10、30、100 cm⁻³、温度 30〜10⁴ K で、H₂ の化学平衡と冷却を比較します。両コードの H・He・C・O・Si の元素数比、入射 FUV 強度、照射の柱密度、宇宙線電離率をそろえます。CHIMES は EUV を除いた B87 の `matched`、EQCLib は追跡元素の CIE のみで CO をオフにします。CHIMES の分子反応には、これらの元素を含む OH・CO などの経路があります。入射光は同じ Habing 強度に規格化し、遮蔽前の H₂ 光解離率を各コードから測ります。最初の原因調査は当時の EQCLib の H₂ 処方と CHIMES の乱流ドップラー幅 7.1 km/s を使います。

```sh
bash test/chimes/shielded/run_h2_comparison.sh
```

CSV は `test/out/chimes/shielded/h2/` に出力します。固定温度では各コードの平衡組成に加え、同じ H₂ 分率と衝突相手の組成を与えた状態でも各項を測ります。生成・消失は H₂ 分子数 cm⁻³ s⁻¹、熱の寄与は erg cm⁻³ s⁻¹ です。等積履歴は 5,000 K における各コードの化学平衡から始め、同じ時刻列で 10 億年まで進めます。生成だけの感度試験では、EQCLib のダスト上の H₂ 生成係数を `Cp=0.25` で変更して再計算します。 [H₂ の比較ノートブック](PlotH2Development.ipynb)

### 文献に沿った H₂ 処方の再比較

旧版を基準に、ダスト上の生成係数の訂正、遮蔽前の光解離率を B87 に合わせる対照、Draine & Bertoldi の式37と H₂ の熱ドップラー幅を使う自己遮蔽を一つずつ適用します。生成と自己遮蔽を更新した現行版は Draine 型スペクトルの遮蔽前光解離率を使います。この版の光解離率だけを B87 に合わせた対照も計算します。CHIMES の乱流ドップラー幅は、現行 EQCLib と同じ 3 km/s です。

```sh
bash test/chimes/shielded/run_h2_prescription_comparison.sh
```

CSV は `test/out/chimes/shielded/h2_prescription/` に出力します。固定温度の平衡、熱平衡温度、5,000 K からの等積履歴を比較します。速度は 10² と 10⁴ 太陽質量のクーラン刻みで、EQCLib の各版の平衡モードを1コアで実行順を入れ替えて7回測ります。CHIMES matched は同じコアで1回計測します。 [H₂ の比較ノートブック](PlotH2Development.ipynb)

### 低温ガスの熱平衡を決める加熱・冷却の内訳

文献に沿った H₂ 処方の比較と同じ密度・元素量・照射・遮蔽条件を使います。固定温度で各コードの化学平衡を求め、ダスト光電加熱、宇宙線加熱、H₂ の熱過程、[C II]・[O I] などの冷却を記録します。89、129、211 K は水素核密度 10 cm⁻³ の3方式の平衡温度付近を直接比較する点です。

```sh
bash test/chimes/shielded/run_thermal_balance_comparison.sh
```

CSV は `test/out/chimes/shielded/thermal_balance/` に出力します。`native_*.csv` は各方式の平衡組成を使った熱収支です。`shared_*.csv` では両コードに同じ H₂、水素イオン、電子、C⁺、O I、Si⁺ の分率を与えます。`eq_legacy_h2` ではさらに H₂ と H I の配分だけを旧処方に戻します。`cool_cii_e_only` は電子衝突だけを使った [C II] の診断値です。指定組成は各項を照合するためのもので、熱平衡温度は各コードの化学平衡組成で得た熱収支から求めます。 [H₂ の比較ノートブック](PlotH2Development.ipynb)
