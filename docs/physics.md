# EQCLib の物理モデル

通常のビルドで化学組成と熱収支に使う反応・過程を記載します。式中の係数、温度の切り替え、成分の加算は EQCLib の実装に従います。ホストコードでの初期化と設定、各モデルの再現計算は別の手順から確認できます。 [シミュレーションコードからの利用](../README.md#シミュレーションコードからの利用) · [検証計算と結果](../test/README.md) · [文献索引](../reference/index.tsv)

- [状態量と時間発展](#状態量と時間発展) · [原子・イオンの反応](#原子イオンの反応)
- [放射場、宇宙線、軟 X 線](#放射場宇宙線軟-x-線) · [H₂ 化学と熱過程](#h₂-化学と熱過程) · [CO 化学と冷却](#co-化学と冷却)
- [原子ガスの冷却](#原子ガスの冷却) · [金属輝線と高温 CIE 冷却](#金属輝線と高温-cie-冷却) · [ダスト再結合冷却と光電加熱](#高温でのダスト再結合冷却と光電加熱)
- [診断輝線と熱収支の集計](#診断輝線と熱収支の集計) · [非平衡での温度係数表](#非平衡での温度係数表) · [参考文献](#参考文献)

## 状態量と時間発展

数密度 `n(X)` は cm⁻³、温度 `T` とダスト温度 `Td` は K、柱密度 `N(X)` は cm⁻²です。フィット式に代入する温度は K 単位での数値です。水素核密度は `nH = n(H I) + n(H II) + 2 n(H₂)`、電子分率は `xe = ne/nH` です。二体反応の体積反応率は `k(T) n(A) n(B)` で、`k` は cm³ s⁻¹です。一体反応は `Γ n(A)` で、`Γ` は s⁻¹です。冷却 `L` と加熱 `H` の各成分と合計は erg cm⁻³ s⁻¹です。

| モード | 化学種と計算方法 |
| --- | --- |
| 化学平衡 | 温度・密度・照射を与え、H I・H II、He I–III、C I–III、O I–III、Si I–III の定常反応収支と電荷保存を解きます。H₂ と CO は個別の設定で有効になります。CO に入った炭素と酸素は原子・イオンの総量から差し引きます。 |
| 非平衡 | H II と有効時の H₂、He・C・O の一価・二価イオンを時間積分し、中性成分を元素保存から求めます。Si I–III の組成は入力から固定し、電荷、状態方程式、[Si II] 35 µm 冷却に反映します。CO を含む状態は平衡モードで扱います。 |
| 水素のみの非平衡 API | `EQCLibTryAdvanceHydrogenChemistry` は固定温度で H I・H II と有効時の H₂ を積分します。入力する気相元素は水素のみです。 |

電子密度は一価イオン密度の和と二価イオン密度の2倍の和です。温度と比内部エネルギー `U` [erg g⁻¹] は `U = f k_B T / (μ m_p)` で変換します。`μ` は陽子質量を単位とする平均分子量、`f = (3 n_mono + 5 n_mol)/(2 n_tot)` です。`n_mono` は原子・イオン・電子の粒子数密度、`n_mol` は分子の粒子数密度です。 [状態方程式と電荷](../src/ChemicalEquilibriumCooling.c#L497) · [化学平衡](../src/ChemicalEquilibriumSolver.c#L91) · [非平衡の状態更新](../src/HydrogenChemistry.c#L64)

熱の更新は固定密度で `dU/dt = (H_tot − L_tot)/ρ` です。`ρ` は `Rho` に渡す質量密度 [g cm⁻³] です。平衡モードでは更新先の温度で組成と熱収支を再評価し、陰的な差分方程式を解きます。非平衡モードでは組成と熱収支を GSL の陰的 BDF 法で連成します。照射、密度、入力の柱密度は1ステップ中固定し、組成に比例する分子柱密度は反応の進行とともに変わります。 [平衡モードのエネルギー更新](../src/ThermalEquilibrium.c#L24) · [非平衡モードの時間積分](../src/HydrogenChemistry.c#L164)

## 原子・イオンの反応

隣接する電離段階間の正味反応率は `F_j = n(X^j) I_j − n(X^(j+1)) R_(j+1)` [cm⁻³ s⁻¹] です。通常の `I_j` [s⁻¹] は電子衝突項 `ne k_ion` と光電離項の和で、下表の段階に宇宙線・軟 X 線の項を加えます。`R_(j+1)` [s⁻¹] は気相項 `ne α` と、一価イオンに対する粒子表面の項 `nH α_gr` の和です。H–O 電荷交換は別に加算します。平衡では正味反応率をゼロにし、非平衡では各段階への流入と流出の差を積分します。 [電離・再結合の係数と段階](../src/EquilibriumAbundance.c#L628) · [Si の平衡](../src/ChemicalEquilibriumSolver.c#L76) · [段階間のフラックス](../src/HydrogenChemistry.c#L116)

| 元素 | 電子衝突・恒星光による電離 | 宇宙線 | 軟 X 線 | 粒子表面での再結合 |
| --- | --- | --- | --- | --- |
| H | H I → H II | H I | H I | H II → H I |
| He | He I → II → III | He I | He I | He II → He I |
| C | C I → II → III。C I には遠紫外線の分岐もあります | C I | | C II → C I |
| O | O I → II → III | O I | | O II → O I |
| Si | 平衡モードで Si I → II → III | 平衡モードで Si I・II | | 平衡モードで Si II → Si I |

CO が有効な低温側では炭素に簡略化した分子ネットワークを使います。電離光があるときの遊離炭素の扱いは [CO モデル](#co-化学と冷却) に記載します。

### 電子衝突電離と気相再結合

H I・He I・He II の衝突電離係数は `k_ion = A sqrt(T) exp(−T₀/T)/(1 + sqrt(10⁻⁵ T))` です。`(A, T₀)` は順に `(5.85×10⁻¹¹, 157807 K)`、`(2.38×10⁻¹¹, 285335 K)`、`(5.68×10⁻¹², 631515 K)` です。C I・II、O I・II、Si I・II には `k_ion = A(1 + P sqrt(u)) u^K exp(−u)/(X + u)` を使います。`u = ΔE/(k_B T)` です。どちらの式も `k_ion` は cm³ s⁻¹です。 [Cen (1992)][cen] · [Voronov (1997) の係数表][voronov] · [衝突電離の実装](../src/CollisionalGasPhaseChemistry.c#L69)

| イオン | `ΔE` [eV] | `A` [cm³ s⁻¹] | `P` | `X` | `K` |
| --- | ---: | ---: | ---: | ---: | ---: |
| C I | 11.3 | 6.85×10⁻⁸ | 0 | 0.193 | 0.25 |
| C II | 24.4 | 1.86×10⁻⁸ | 1 | 0.286 | 0.24 |
| O I | 13.6 | 3.59×10⁻⁸ | 0 | 0.073 | 0.34 |
| O II | 35.1 | 1.39×10⁻⁸ | 1 | 0.212 | 0.22 |
| Si I | 8.2 | 1.88×10⁻⁷ | 1 | 0.376 | 0.25 |
| Si II | 16.4 | 6.43×10⁻⁸ | 1 | 0.632 | 0.20 |

| 上位イオン | 気相再結合係数 `α` [cm³ s⁻¹] |
| --- | --- |
| H II、He III | `α = A λ^p [1 + (aλ)^b]^(−c)`、`λ = 2 T₀/T`。H II では case A に `(A,p,a,b,c) = (1.269×10⁻¹³, 1.503, 1.91571, 0.470, 1.923)`、case B に `(2.753×10⁻¹⁴, 1.500, 0.36496, 0.407, 2.242)` を使います。He III は `T₀ = 631515 K` と各 case の `2A` を使います。 [Hui & Gnedin (1997)][hui] · [実装](../src/CollisionalGasPhaseChemistry.c#L267) |
| He II | `λ = 2 × 285335/T` とし、case A は `3.0×10⁻¹⁴ λ^0.654`、case B は `1.26×10⁻¹⁴ λ^0.750` です。組成の He II → He I にこの係数を使い、He I の冷却には二電子再結合冷却も加えます。 [Hui & Gnedin (1997)][hui] · [実装](../src/CollisionalGasPhaseChemistry.c#L285) |
| C II、O II | 順に `4.40×10⁻¹² (0.0033333 T)^−0.61`、`3.40×10⁻¹² (0.0033333 T)^−0.63` です。 [Millar et al. (1997)][millar] · [実装](../src/CollisionalGasPhaseChemistry.c#L337) |
| C III、O III | 下に示す放射再結合と二電子再結合のフィットを使います。O III は Nahar (1999) の酸素データに基づきます。C III の出典候補は Nahar & Pradhan (1997) で、フィット元は未同定です。 [酸素データ][nahar-o] · [炭素データ][nahar-c] · [文献索引](../reference/index.tsv) · [実装](../src/CollisionalGasPhaseChemistry.c#L350) |
| Si II | `α = A(T/300)^b` です。`(A,b)` は `T ≤ 2000`、`2000 < T ≤ 10⁴`、`T > 10⁴` の順に `(7.50×10⁻¹², −0.55)`、`(4.86×10⁻¹², −0.32)`、`(9.08×10⁻¹⁴, 0.818)` です。 [Glover & Jappsen (2007) の反応31][gj07] · [実装](../src/CollisionalGasPhaseChemistry.c#L388) |
| Si III | `α = 1.75×10⁻¹² (T/10⁴)^−0.6346 + T_e^−1.5 [2.2552×10⁻¹¹ exp(−2.76/T_e) + 5.6058×10⁻⁹ exp(−10.13/T_e)]`。ここで `T_e` は `k_B T` を eV で表した数値です。 [Glover & Jappsen (2007) の反応47][gj07] · [実装](../src/CollisionalGasPhaseChemistry.c#L395) |

```text
α_CIII = 2.0×10⁻¹¹ (T/300)^−0.61
       + 2.6×10⁻⁶ T^−1.5 exp(−5150/T) [1 + 0.5 exp(−1000/T)]
       + 8.5×10⁻³ T^−1.7 exp(−125000/T) [1 + 4.5 exp(−55000/T)]
       + 1.9×10⁻¹² exp(−10⁻⁶ T)
α_OIII = 5.5×10⁻¹⁰ T^−0.585
       + 4.8×10⁻³ T^−1.55 exp(−205000/T) [1 + 0.55 exp(−150000/T)]
```

H と He の case A は再結合光子が逃げる近似、case B は基底準位への再結合で生じた電離光子が局所で再吸収される近似です。`EQCLibSetRunParamRecombinationCaseB` で選びます。 [Hui & Gnedin (1997)][hui] · [係数の切り替え](../src/EquilibriumAbundance.c#L571)

### H–O 電荷交換

| 反応 | 係数 `k` [cm³ s⁻¹] | 出典 |
| --- | --- | --- |
| O I + H II → O II + H I | `7.0×10⁻¹⁰ exp(−232/T)` | [Millar et al. (1997)][millar] |
| O II + H I → O I + H II | `7.7×10⁻¹⁰` | [Millar et al. (1997)][millar] |
| O III + H I → O II + H II | `1.04×10⁻⁹ t₄^0.27 [1 + 2.02 exp(−5.92 t₄)]`。ここで `t₄ = T/10⁴` は 0.01〜10 に制限します | [Kingdon & Ferland (1996)][kingdon] |

各係数に反応物2種の密度を掛けて、H と O の状態に対応する符号で反映します。平衡モードでは H II 領域と電荷交換の両方を有効にすると、H の電離分率を変えながら O の電離平衡と H の収支を解きます。H₂ 生成が有効なら、各試行で H₂ の生成・消滅の釣り合いも解き直します。非平衡の熱化学更新では、照射の有無にかかわらず電荷交換を有効にすると3反応を積分します。 [反応係数](../src/CollisionalGasPhaseChemistry.c#L422) · [平衡側の有効条件](../src/ChemicalEquilibriumSolver.c#L42) · [非平衡での電荷交換](../src/HydrogenChemistry.c#L135) · [電荷交換の有無による検証結果](../test/README.md#h₂coho-電荷交換の設定比較)

### 粒子表面での再結合

`α_gr = 10⁻¹⁴ D C₀ / [1 + C₁(ψ^C₂ + C₃ T^C₄ ψ^(C₂−C₅−C₆ ln T))]` [cm³ s⁻¹] を使います。`ψ = G0 sqrt(T)/ne` は指定された入射 `G0` から計算し、`D = Zd/0.019` です。H II、He II、C II、Si II はイオンごとの係数を使い、O II の `α_gr` は H II の0.25倍です。粒子表面過程を有効にし `Zd > 0` のときに加え、`G0 = 0` では `α_gr = 10⁻¹⁴ D C₀` です。`G0 > 0` かつ `ne = 0` のときはゼロです。 [Weingartner & Draine (2001)][wd01] · [Glover & Jappsen (2007)][gj07] · [粒子表面再結合の係数と実装](../src/GrainSurfaceChemistry.c#L48)

## 放射場、宇宙線、軟 X 線

### 恒星光による光電離と加熱

H II 領域を有効にしたセルでは、入力の H 電離光子束 `NLyn` [photons cm⁻² s⁻¹] に規格化した恒星スペクトルを使います。光子スペクトルを `φ(E)`、元素・電離段階ごとのしきい値を `E_i`、光電離断面積を `σ_i(E)` [cm²] とすると、実装のイオン1個あたりの率は次の形です。`Q_H = ∫_(13.60 eV)^∞ φ(E) dE` で規格化した係数を年齢ごとに求め、それらを等重みで平均します。

```text
Γ_i = NLyn ⟨∫_(E_i)^∞ φ(E) σ_i(E) dE / Q_H⟩               [s⁻¹]
h_i = NLyn ⟨∫_(E_i)^∞ φ(E) σ_i(E) (E − E_i) dE / Q_H⟩     [erg s⁻¹]
H_PI = Σ_i n(i) h_i                                      [erg cm⁻³ s⁻¹]
```

積分には PEGASE2 の Kroupa IMF、0〜30 Myr のスペクトルと Verner の断面積を使います。実行時には恒星スペクトルの金属量 `Z = 0.0200` の係数を使い、セルの金属量は吸収する各化学種の密度に反映します。断面積で重み付けするのは H I、He I・II、C I・II、O I・II と、平衡モードの Si I・II です。伝播や入射電離光の減衰はホスト側で `NLyn` に反映する扱いです。 [PEGASE2][pegase] · [Kroupa (2001)][kroupa] · [Verner et al. (1996)][verner] · [光電離・加熱率](../src/PhotochemicalGasPhaseChemistry.c#L110) · [スペクトル積分](../src/PhotochemicalGasPhaseChemistry.c#L430) · [同梱スペクトルと指定方法](../README.md#pegase2-のスペクトルを指定する)

### 遠紫外線とダスト遮蔽

`G0` は Habing 場を単位とする無次元の入射遠紫外線強度、`Nt` は水素核柱密度です。ダスト減光が有効なとき、可視減光は `A_V = 5.35×10⁻²² D Nt` [mag]、一般の減光係数は `exp(−γ A_V)` です。H₂ 光解離には別のダスト減光 `exp(−2.0×10⁻²¹ D Nt)` を使います。ダスト減光を無効にすると `A_V = 0` となり、ダストによる減光係数はいずれも1です。分子の線遮蔽はそのときも柱密度から評価します。 [Glover et al. (2010)][g10] · [Glover & Mac Low (2007)][gml] · [ダスト減光の実装](../src/MoleculeChemistry.c#L151)

遊離炭素を原子ネットワークで解くとき、C I には H II 領域が有効で `NLyn > 0` なら恒星光の `Γ_CI`、それ以外では `3.1×10⁻¹⁰ G0 exp(−3 A_V)` [s⁻¹] を使います。CO モデル内部の C I 光電離には同じ遠紫外線の式を使います。 `H_PI` に足すのは恒星スペクトルで求めた光電離加熱です。 [Glover et al. (2010) の反応171][g10] · [C I の光電離の分岐](../src/EquilibriumAbundance.c#L513) · [光電離加熱の集計](../src/PhotochemicalGasPhaseChemistry.c#L160)

### 宇宙線

宇宙線を有効にすると、H I、He I、C I、O I、Si I、Si II の電離率に `X0 ζ_i` [s⁻¹] を加えます。`X0` は宇宙線と軟 X 線に共通の無次元の強度倍率です。基準 `ζ_HI = 1.8×10⁻¹⁷ s⁻¹` に対して順に `1, 1.1, 3.8, 5.7, 6.5, 2.5` 倍です。Si の電離は平衡モードで評価します。

加熱は `H_CR = nH(1 − xe) X0 ζ_HI E_heat(xe)` で、`0 ≤ xe < 1` の範囲で評価します。`E_heat = (35 eV) [f₁ + 1/(56.4/f₃ + 1/(1−f₁))]` は一次電子1個あたりの加熱エネルギーで、率の計算では erg に換算します。`f₁ = 0.9971 [1 − (1 − xe^0.2663)^1.3163]` です。`f₃` は `xe = 10⁻⁴, 10⁻³, 10⁻², 10⁻¹` の順に `3.25, 7.20, 16.2, 31.0` とし、区分線形で補間・外挿します。 [Le Teuff et al. (2000)][lte] · [Glover et al. (2010)][g10] · [Glover & Jappsen (2007)][gj07] · [Wolfire et al. (1995) の式4と付録A][wolfire95] · [宇宙線の実装](../src/CosmicRay.c)

### 軟 X 線

軟 X 線を有効にすると、H I と He I にそれぞれ `0.5 X0 ζ_X` [s⁻¹] の電離率を加えます。`ζ_X = 10^F` は `0 < xe ≤ 0.5` で評価します。加熱は水素核1個あたりの `Q_X = 10^G` [erg s⁻¹] から `H_XR = nH X0 Q_X` とし、`xe > 0` で評価します。フィットの変数は `p₁ = max(−0.3, log10(10⁻¹⁸ Nt))` と `p₂ = log10(xe)` です。電離率が `xe = 0.5` をまたいで不連続なため、非平衡更新がそこで `EQC_NO_CONVERGENCE` になる場合があります。 [Wolfire et al. (1995) の式A6–A10][wolfire95] · [電離と加熱のフィット](../src/Xray.c#L35)

```text
 F = f₄ (−15.6 − 1.10 p₁ + 0.0913 p₁²) + 0.87 f₅ exp[−((p₁−0.41)/0.84)²]
 G = f₆ [−26.5 − 0.920 p₁ + 0.0589 p₁² + 0.96 exp[−((p₁−0.38)/0.87)²]]
f₄ = 1.06 + 0.0408 p₂ + 0.00651 p₂²
f₅ = 1.90 + 0.678 p₂ + 0.113 p₂²
f₆ = 0.990 − 0.00274 p₂ + 0.00113 p₂²
```

## H₂ 化学と熱過程

H₂ 生成を有効にすると、粒子表面での生成、遠紫外線による光解離、電子・H I・H₂ による衝突解離を組成に反映します。分子数の収支は `d n(H₂)/dt = R_form − n(H₂)(Γ_diss + Γ_col)` [cm⁻³ s⁻¹] です。平衡モードではゼロ点を求め、非平衡モードでは水素核分率 `2 n(H₂)/nH` を積分します。 [平衡の H₂ 反応収支](../src/EquilibriumAbundance.c#L209) · [非平衡の分子水素](../src/HydrogenChemistry.c#L222)

| 反応 | 率式と適用条件 |
| --- | --- |
| ダスト上での H₂ 生成 | `R_form = k_form nH n(H I)` [cm⁻³ s⁻¹]。 `k_form = 3.0×10⁻¹⁸ C_p D sqrt(T) S(T,Td) / [1 + 10⁴ exp(−600/Td)]` [cm³ s⁻¹]、`S = 1/[1 + 0.04 sqrt(T+Td) + 0.002 T + 8×10⁻⁶ T²]` です。`C_p` は生成率の倍率で、`Td` に正の値を渡さない場合は通常ビルドの 20 K を使います。 [Glover & Jappsen (2007) の反応60][gj07] · [生成率](../src/MoleculeChemistry.c#L229) · [ダスト温度](../src/GrainSurfaceChemistry.c#L207) |
| 遠紫外線による光解離 | `Γ_diss = 3.3×10⁻¹¹ G0 exp(−2.0×10⁻²¹ D Nt) f_sh` [s⁻¹]。Draine 型スペクトルのフィットに入力の `G0` を掛けます。式はダスト減光が有効な場合です。自己遮蔽にはセルの組成から求める `N(H₂) = Nt n(H₂)/nH` を使います。 [Glover & Mac Low (2007) の式29・31][gml] · [光解離率](../src/MoleculeChemistry.c#L278) |
| 衝突解離 | `Γ_col = ne k_e + n(H I) k_H + n(H₂) k_H2` [s⁻¹]。係数は `k_e = 3.73×10⁻⁹ T^0.1121 exp(−99430/T)`、`k_H = 6.67×10⁻¹² sqrt(T) exp(−1−63590/T)`、`k_H2 = 5.996×10⁻³⁰ T^4.1881 exp(−54657.4/T) (1+6.761×10⁻⁶ T)^−5.6881` [cm³ s⁻¹] です。 [Glover & Jappsen (2007) の反応8–10][gj07] · [衝突解離係数](../src/MoleculeChemistry.c#L305) |

自己遮蔽には Draine & Bertoldi の式37を `N(H₂) = 0` で1になるよう規格化して使います。`x = N(H₂)/(5×10¹⁴ cm⁻²)`、`b₅ = b/(1 km s⁻¹)` とすると、`f_sh = F(x,b₅)/F(0,b₅)` です。

```text
F(x,b₅) = 0.965/(1 + x/b₅)² + 0.035 exp(−8.5×10⁻⁴ sqrt(1+x))/sqrt(1+x)
b² = b_turb² + 2 k_B T/m(H₂)
```

乱流ドップラー幅 `b_turb` の初期値は 3 km s⁻¹です。平衡計算は初期の H₂ 分率を手掛かりに安定な生成・消失収支を探し、許容誤差が `1e-12` 以下の場合や Newton 法が収束しない場合は区間探索を使います。 [Draine & Bertoldi (1996) の式37][db96] · [遮蔽とドップラー幅](../src/MoleculeChemistry.c#L243) · [平衡の根の選択](../src/EquilibriumAbundance.c#L391)

| 熱過程 | 体積率 [erg cm⁻³ s⁻¹] |
| --- | --- |
| 回転振動冷却 | `L_rv = n(H₂) Σ_j n(j) Λ_j(T)`。低密度極限の `Λ_j` [erg cm³ s⁻¹] は H I・H₂・He・H II・電子の Table 8 の多項式です。H I・H₂・He は 6,000 K、H II・電子は 10,000 K の値を上限側に延長します。10 K 以下では H₂ 衝突項だけを使います。非平衡では He I の密度、平衡では He の全密度を係数に掛けます。 [Glover & Abel (2008)][ga08] · [温度範囲と衝突係数](../src/MoleculeChemistry.c#L400) · [集計](../src/ChemicalEquilibriumCooling.c#L426) |
| 衝突解離冷却 | `L_diss = (4.48 eV) n(H₂) Γ_col`。1分子の結合エネルギーを熱から引きます。 [Hollenbach & McKee (1979)][hm79] · [実装](../src/MoleculeChemistry.c#L707) |
| 光解離加熱 | `H_diss = 6.4×10⁻¹³ n(H₂) Γ_diss`。係数は1解離あたりの加熱 [erg] です。 [Glover & Mac Low (2007)][gml] · [実装](../src/MoleculeChemistry.c#L729) |
| 紫外線ポンピング加熱 | `H_pump = 8.5 × 3.2×10⁻¹² n(H₂) Γ_diss f_c`。`8.5` は励起と解離の比、`f_c = nH/(nH+n_crit)` は衝突脱励起の効率です。 [Glover & Mac Low (2007)][gml] · [Draine & Bertoldi (1996)][db96] · [実装](../src/MoleculeChemistry.c#L729) |
| 生成加熱 | `H_form = 7.2×10⁻¹² R_form f_c`。係数は1分子あたりのエネルギー [erg] です。 [Glover & Mac Low (2007)][gml] · [実装](../src/MoleculeChemistry.c#L729) |

`n_crit` [cm⁻³] は `2 / [(n(H I)/nH)/n_crit,HI + (n(H₂)/nH)/n_crit,H2]` です。`u = log10(T/10⁴)` とし、H I には `10^(4.00 − 0.416u − 0.327u²)`、H₂ には `10^(4.845 − 1.3u + 1.62u²)` を使います。熱収支には `L_H2 = L_rv + L_diss` と `H_H2 = H_diss + H_pump + H_form` を加えます。 [Lepp & Shull (1983)][ls83] · [Shapiro & Kang (1987)][sk87] · [臨界密度](../src/MoleculeChemistry.c#L319)

## CO 化学と冷却

CO は平衡モードで生成を有効にした場合に扱います。炭素の分子ネットワークは `T ≤ 10⁴ K` で使い、それより高温では原子・イオンのネットワークを解きます。低温側では C⁺ + H₂ の放射会合を起点とする CO 生成を1反応で近似し、CO の光解離を C I に戻し、C I の光電離を C⁺ に戻します。 [Keto & Caselli (2008)][kc08] · [Glover & Clark (2012) の KC08e][gc12] · [分子炭素の平衡](../src/EquilibriumAbundance.c#L1240)

| 過程 | 率式 |
| --- | --- |
| 有効 CO 生成 | `R_CO = 4.0×10⁻¹⁶ C_p,CO (T/300)^−0.2 n(C II) n(H₂)` [cm⁻³ s⁻¹]。`C_p,CO` は生成率の倍率です。 [Glover et al. (2010) の反応148][g10] · [実装](../src/EquilibriumAbundance.c#L536) |
| CO 光解離 | `Γ_CO = 2.0×10⁻¹⁰ G0 exp(−2.5 A_V) f_CO(N(CO)) f_H2→CO(N(H₂))` [s⁻¹]。CO 自己遮蔽と H₂ による遮蔽を掛けます。柱密度は `N(X) = Nt n(X)/nH` です。 [Lee et al. (1996) の表10・11][lee96] · [Glover & Clark (2012)][gc12] · [実装](../src/MoleculeChemistry.c#L786) |
| C I 光電離 | `Γ_CI = 3.1×10⁻¹⁰ G0 exp(−3.0 A_V)` [s⁻¹]。 [Glover et al. (2010)][g10] · [実装](../src/EquilibriumAbundance.c#L541) |

遮蔽表は柱密度と遮蔽係数の対数で線形補間します。H₂ による CO 遮蔽は表の高柱密度端から傾きを連続に延長し、Visser et al. の値を使ってその先の減光を決めます。CO の自己遮蔽は表の上端を超えると `10⁻⁷` です。

生成と光反応が進む定常状態では、分子ネットワーク内で `R_CO = Γ_CO n(CO) = Γ_CI n(C I)` を解きます。ここで `R_CO` に使う C II は分子ネットワーク内の量です。恒星電離光が有効なときは、CO 量を求めた後に、遊離炭素を原子・イオンのネットワークで再配分します。 [Visser et al. (2009) の表5][visser] · [遮蔽表と外挿](../src/MoleculeChemistry.c#L96) · [遊離炭素の扱い](../src/EquilibriumAbundance.c#L1325)

CO の回転冷却は `L_rot = n(CO) 4(k_B T)² A₀/[E₀(1+r+1.5 sqrt(r))]` [erg cm⁻³ s⁻¹] です。`A₀ = 9.7×10⁻⁸ s⁻¹`、`E₀/k_B = 2.76 K`、`r = n_crit/nH`、`n_crit = 3.3×10⁶ (T/10³)^0.75 cm⁻³` です。振動冷却は `L_vib = n(CO) E₁₀ [n(H I) γ_H + n(H₂) γ_H2]` で、`E₁₀/k_B = 3080 K` と次の係数を使います。`L_CO = L_rot + L_vib` です。 [McKee et al. (1982)][mckee82] · [Hollenbach & McKee (1979)][hm79] · [Hollenbach & McKee (1989)][hm89] · [CO 冷却の実装](../src/MoleculeChemistry.c#L938)

```text
γ_H  = 3.0×10⁻¹² sqrt(T) exp[−(2000/T)^3.43 − 3080/T]     [cm³ s⁻¹]
γ_H2 = 4.3×10⁻¹⁴ T exp[−(3.14×10⁵/T)^0.333 − 3080/T]    [cm³ s⁻¹]
```

### CO モデルの適用範囲

KC08e の3反応では、電離光を含まない定常状態の比 `n(C I)/n(CO) = Γ_CO/Γ_CI` です。ダスト減光の指数が CO で2.5、C I で3.0のため、分子の線遮蔽を固定すると比は `exp(0.5 A_V)` で増えます。この近似は C⁺ の再結合と宇宙線による CO の化学反応を省略します。

- `A_V ≲ 20` では CO は KC08e の想定どおり振る舞います。
- `A_V ≳ 25` ではモデル上 CO が減り、C I が主になります。実際の高密度コアでは CO が主成分として残るため、この領域の C・CO 存在比と冷却は適用範囲外です。
- Glover & Clark は、KC08 が C⁺ の再結合を省略したため、原子状炭素が観測より少なくなると報告しています。

銀河規模の計算で使う場合は、粒子の `A_V` が20を超える割合を確認してください。深部を扱うには、C⁺ の再結合や CHx・OHx を含む Nelson & Langer (1999) 型のようなネットワークが必要です。 [Glover & Clark (2012)][gc12] · [Nelson & Langer (1999)][nl99] · [柱密度走査と試算の結果](../test/README.md#柱密度走査での-co-モデルの確認)

## 原子ガスの冷却

H と He の電子衝突励起・電離、気相再結合と自由–自由放射を評価します。`L_HI`、`L_HeI`、`L_HeII` は、それぞれ励起・電離と隣の上位イオンからの再結合を足したものです。下表の式は体積冷却率 [erg cm⁻³ s⁻¹] です。 [集計](../src/ChemicalEquilibriumCooling.c#L366)

| 過程 | 率式と出典 |
| --- | --- |
| 自由–自由放射 `L_ff` | `1.42×10⁻²⁷ × 1.5 sqrt(T) ne [n(H II)+n(He II)+4 n(He III)]`。Gaunt 因子は1.5です。 [Cen (1992) の式16][cen] · [実装](../src/Bremsstrahlung.c#L7) |
| 衝突電離冷却 | H I・He I・He II に `ne n(i) A sqrt(T) exp(−T₀/T)/(1 + sqrt(10⁻⁵ T))` を使います。`T₀` は衝突電離の表と同じで、`A` は順に `1.27×10⁻²¹, 9.38×10⁻²², 4.95×10⁻²²` です。 [Cen (1992)][cen] · [実装](../src/CollisionalGasPhaseChemistry.c#L502) |
| 衝突励起冷却 | H I は `ne n(H I) 7.5×10⁻¹⁹ exp(−118348/T)/d`、He II は `ne n(He II) 5.54×10⁻¹⁷ T^−0.397 exp(−473638/T)/d` です。He I の三重項には、実装で `ne² n(He II) 9.10×10⁻²⁷ T^−0.1687 exp(−13179/T)/d` を使います。`d = 1 + sqrt(10⁻⁵ T)` です。通常の二体項の冷却係数は erg cm³ s⁻¹、He I の三重項の係数は erg cm⁶ s⁻¹です。 [Cen (1992)][cen] · [係数と密度因子](../src/CollisionalGasPhaseChemistry.c#L516) |
| 放射再結合冷却 | `ne n(i) η_i(T)`。H II と He III の `η` [erg cm³ s⁻¹] は下に示すフィット、He II は `η = k_B T α_HeII(T)` です。 [Hui & Gnedin (1997)][hui] · [実装](../src/CollisionalGasPhaseChemistry.c#L535) |
| He II の二電子再結合冷却 | `ne n(He II) × 0.75 k_B T₀ × α_DR`。`T₀ = 631515 K`、`α_DR = 1.90×10⁻³ T^−1.5 exp(−0.75 T₀/T) [1+0.3 exp(−0.15 T₀/T)]` [cm³ s⁻¹] です。 [Cen (1992) の式14][cen] · [実装](../src/CollisionalGasPhaseChemistry.c#L566) |

H II と He III の `λ = 2 T₀/T` には気相再結合と同じ `T₀` を使います。以下は H II の係数で、He III では数値係数を2倍します。

```text
η_A = 1.778×10⁻²⁹ T λ^1.965 [1 + (1.84843 λ)^0.502]^−2.697
η_B = 3.435×10⁻³⁰ T λ^1.970 [1 + (0.44444 λ)^0.376]^−3.720
```

## 金属輝線と高温 CIE 冷却

### 微細構造輝線

衝突脱励起率は `C_ul = Σ_j n(j) q_ul,j` [s⁻¹] です。`q_ul,j` [cm³ s⁻¹] は衝突相手ごとの係数で、励起率には詳細釣り合い `q_lu,j = (g_u/g_l) exp(−ΔE/k_B T) q_ul,j` を使います。準位数密度は衝突、自然放射、背景放射による誘導遷移の統計平衡から求めます。各線の正味冷却率は `L_ul = ΔE [(A_ul + B_ul J) n_u − B_lu J n_l]` [erg cm⁻³ s⁻¹] です。`J` は各ルーチンが Einstein B に渡す量です。背景放射は 2.735 K の宇宙マイクロ波背景を使い、非平衡の [C II] では平均比強度 `J_ν` を `4π J_ν/c` に換算して渡します。平衡の [C II]、[O I]、[Si II] にはコードにある従来の背景量を使います。 [Glover & Jappsen (2007) の表6][gj07] · [準位平衡](../src/StatisticalEquilibrium.c) · [背景放射と [C II] の実装](../src/FineStructureLines.c#L605)

| 熱収支に含む線 | 準位と衝突相手 |
| --- | --- |
| [C II] 158 µm | 2準位。電子、H I、オルト・パラ H₂。 [係数と実装](../src/FineStructureLines.c#L631) |
| [O I] 63・146 µm | 3準位。H I、H II、オルト・パラ H₂。44 µm 遷移も準位平衡に入ります。 [実装](../src/FineStructureLines.c#L713) |
| [Si II] 35 µm | 2準位。電子と H I。非平衡モードでは入力の固定 Si II 密度を使います。 [実装](../src/FineStructureLines.c#L869) |

微細構造輝線の H₂ 衝突項では、オルトとパラの比を 3:1 とします。熱収支の集計では [C I] の項はゼロで、[O I] 44 µm の放射は準位平衡への寄与までです。 [衝突相手の比率](../src/FineStructureLines.c#L18) · [冷却率の集計](../src/ChemicalEquilibriumCooling.c#L383)

### H II 領域の酸素輝線

O II と O III はそれぞれ5準位の統計平衡を解きます。電子脱励起係数は `q_ul = 8.629×10⁻⁶ Υ_ul/(g_u sqrt(T))` [cm³ s⁻¹] です。励起側には詳細釣り合いのボルツマン因子を掛け、衝突強度 `Υ` は対数温度上で線形補間し、範囲外では端点値を使います。自然放射を加えて準位占有率 `p_u` を求め、光学的に薄い放射率 `n(ion) Σ_(u>l) p_u A_ul ΔE_ul` [erg cm⁻³ s⁻¹] を計算します。 [準位・遷移データ](../src/HIIregionOxygenData.h) · [酸素輝線の実装](../src/HIIregionLineCooling.c#L73)

O II の遷移確率には Zeippen (1982)、O III には Froese Fischer & Tachiev (2004) と Storey & Zeippen (2000) を使います。有効衝突強度は O II が Kisielius et al. (2009)、O III が Storey et al. (2014) です。準位エネルギーには NIST の値を使い、使用した PyNeb のデータ表をヘッダーに記載しています。 [Zeippen (1982)][oii-data] · [Froese Fischer & Tachiev (2004)][fischer] · [Storey & Zeippen (2000)][storey-zeippen] · [Kisielius et al. (2009)][kisielius] · [Storey et al. (2014)][storey] · [NIST Atomic Spectra Database][nist] · [PyNeb][pyneb]

### 元素別 CIE 冷却と酸素輝線の切り替え

衝突電離平衡（CIE）では、光学的に薄い低密度ガスに対する温度別の金属冷却係数を使います。C・O・Si はセルの**全元素数密度**、補完の N・Ne・Mg・S・Ca・Fe は太陽組成と `Z` を使います。係数の単位は erg cm³ s⁻¹です。補完係数は太陽組成の水素核1個あたりに換算済みです。 [Gnat & Ferland (2012) の元素別 CIE 係数][gf12] · [Asplund et al. (2009) の太陽組成][asplund] · [テーブル生成](../data/GnatFerland2012/generate_cie_table.py) · [CIE 冷却の実装](../src/CoolingCIE.c#L72)

```text
L_CIE = ne [ Σ_(X=C,O,Si) n(X) Λ_X(T) + nH (Z/0.0134) w_extra(T) Λ_extra,solar(T) ]
```

`w_extra` は 10⁴ K 以下で0、10⁵ K 以上で1です。間は `s = log10(T/10⁴)` の `S(s) = s²(3−2s)` を使います。物理的な全金属質量分率を `Z` に渡します。旧基準 0.019 での太陽比は `Z = 0.0134 (Z_old/0.019)` に換算します。`Zd` はダストの比として別に設定します。単一セルの実行例 `eqclib-cell` では `Z = 0.0134 (O/H)/10⁻³·³¹` とし、非平衡モードで CIE を有効にします。 [補完項と基準金属量](../src/CoolingCIE.c#L61) · [定数](../src/Constants.h#L100) · [単一セルの金属量](../src/CellRunner.c#L155)

平衡モードは `log10 T` に対して `log10 Λ` を線形補間し、非平衡モードは Steffen 補間を使います。表の範囲は 10⁴〜10⁸ K です。10⁴ K 未満では冷却率に `exp[100(log10(T/K)−4)]` を掛け、10⁸〜10⁹ K では上端の係数を `sqrt(T/10⁸ K)` 倍します。CIE を有効にしている場合の上限は 10⁹ K です。`EQCLibSetRunParamCIEResonanceCooling(true)` と `Z > 0` で CIE 項を使います。 [Steffen (1990)][steffen] · [補間と範囲の処理](../src/CoolingCIE.c#L72)

酸素の輝線和に `w_O`、CIE 冷却に `1 − w_O` を掛けます。温度の重みは `w_T = S(clamp((25000−T)/10000, 0, 1))` です。

| モード | 酸素輝線の重みと適用条件 |
| --- | --- |
| 平衡 | `w_O = w_T × [n(H II)/nH] × Γ_π,H/(Γ_π,H + Γ_other,H)`。`Γ_other,H` は電子衝突・宇宙線・軟 X 線の H I 電離率の和です。H II 領域が有効で `NLyn > 0`、かつ光電離率が正のときに使います。 |
| 非平衡 | CIE 有効かつ `Z > 0` なら `w_O = w_T` です。照射の有無にかかわらず現在の O II・O III 密度を使います。CIE を無効にした場合、または `Z = 0` の場合は `w_O = 1` です。 |

`L_OII` と `L_OIII` は重みを掛けた明示的な酸素輝線、`L_res` は重みを掛けた CIE 冷却です。高温の CIE はより高い電離段階を含む冷却の近似で、時間積分する酸素の状態は O I–III です。 [平衡と非平衡の重み](../src/ChemicalEquilibriumCooling.c#L298) · [冷却率の切り替え](../src/ChemicalEquilibriumCooling.c#L387) · [非平衡ガスと CIE の比較][gs07] · [高温冷却の検証](../test/README.md#高温ガスの冷却と消灯後の酸素輝線)

組み込み表はリポジトリのルートで `python3 data/GnatFerland2012/generate_cie_table.py > src/CoolingTableGnatFerland.h` として生成します。最初の初期化で `Check=1` を渡すと、実行ディレクトリの `EQCooling/CIECooling.dat` に `log10(T/K)`、C・O・Si の係数、補完係数が出力されます。補完係数は温度の重みを掛ける前の値です。 [表の生成スクリプト](../data/GnatFerland2012/generate_cie_table.py) · [初期化と診断](../src/CoolingCIE.c#L25)

## 高温でのダスト再結合冷却と光電加熱

ダストに対する二つの熱過程は Bakes & Tielens の温度フィットを使います。減光が有効なら `G = G0 exp(−2.5 A_V)`、無効なら `G = G0` とし、`ψ = G sqrt(T)/(ne φ_PAH)`、`φ_PAH = 1` とします。次の率の単位は erg cm⁻³ s⁻¹です。

```text
H_PE = 1.0×10⁻²⁴ D nH G w_d ε
ε = 0.0487/(1 + 0.004 ψ^0.73) + 0.0365 (T/10⁴)^0.7/(1 + 0.0002 ψ)
L_DustRec = 3.49×10⁻³⁰ D ne nH φ_PAH T^0.944 ψ^β w_d
β = 0.735 T^−0.068
```

`w_d` は 10⁴ K 以下で1、10⁵ K 以上で0です。間は `w_d = 1 − S(log10(T/10⁴))` で滑らかに変え、高温への外挿を抑えます。`L_DustRec` は粒子表面過程が有効な場合に評価します。`H_PE` は入射遠紫外線がある場合に評価します。どちらも電子やダストがない場合はゼロです。この温度重みは熱過程に使い、粒子表面の化学再結合率には前述の式を使います。 [Bakes & Tielens (1994) の式42–44][bt94] · [Wolfire et al. (1995) の式11][wolfire95] · [熱過程の実装](../src/GrainSurfaceChemistry.c#L125) · [開発時のダスト冷却比較](../test/chimes/shielded/README.md#高温のダスト冷却加熱の比較)

## 診断輝線と熱収支の集計

`H_tot` は光電離、光電加熱、宇宙線、軟 X 線、H₂ の和です。`L_tot` は自由–自由放射、H・He の冷却、金属輝線、CIE、ダスト再結合、H₂、CO を足します。Hα と H I 21 cm は、組成・温度から計算する診断出力です。 [冷却・加熱率の集計](../src/ChemicalEquilibriumCooling.c#L333)

| 診断量 | 式と単位 |
| --- | --- |
| Hα `L_Halpha` | `4π × 2.82×10⁻²⁶ ne n(H II) t₄^(−0.942−0.031 ln t₄) + 1.30×10⁻¹⁷ ne n(H I) exp(−12.1 eV/k_B T) G₁₃(T)/sqrt(T)` [erg cm⁻³ s⁻¹]。`t₄ = T/10⁴` で、`G₁₃` は衝突励起に用いる区分多項式です。 [Peters et al. (2017) の式4・5][peters] · [実装](../src/HIIregionLineCooling.c#L175) |
| H I 21 cm `L_HI21cm` | `4πj = (3/4) hν A₁₀ n(H I)` [erg cm⁻³ s⁻¹]。 `ν = 1.4204×10⁹ Hz`、`A₁₀ = 2.8843×10⁻¹⁵ s⁻¹` です。 [Draine (2011)][draine] · [実装](../src/HyperFineStructureLines.c#L7) |
| H I 21 cm のスピン温度と吸収量 | `T_s = (T_R + y_c T + y_α T_α)/(1+y_c+y_α)` [K]。`T_R = 3.77 K`、`T_α = T` とし、衝突と Lyα の結合を評価します。`k_HI21cm = 3hc² A₁₀ n(H I)/(32π k_B ν T_s)` の次元は周波数積分した吸収量 [Hz cm⁻¹] です。周波数ごとの係数は、積分を1に規格化した線形状関数 [Hz⁻¹] を掛けて得ます。 [Field (1958)][field] · [Draine (2011)][draine] · [実装](../src/HyperFineStructureLines.c#L26) |

Hα と H I 21 cm は `L_tot` と別のフィールドに格納します。そのほか、ガス–ダスト間の熱輸送や H₂O 冷却などの補助関数は、上記の熱収支ルーチンとは別に実装されています。 [診断量の出力](../src/ChemicalEquilibriumCooling.c#L445) · [ガス–ダスト間の熱輸送](../src/GrainSurfaceChemistry.c#L178) · [H₂O の補助関数](../src/MoleculeChemistry.c#L1047)

## 非平衡での温度係数表

非平衡更新は既定で 1〜10⁸ K の対数温度格子（1 dex あたり1,024区間）を使い、隣接4点で補間します。温度区分の境界を含む区間、格子外、ゼロ近傍、区間内の検査点で相対差が `5×10⁻⁹` を超える区間では率式を直接評価します。`EQCLibSetRunParamNonEquilibriumTemperatureTables(false)` で直接評価を選べます。

対象は H・He・C・O の電子衝突電離と再結合、H–O 電荷交換、H₂ の衝突解離、H・He と H₂ の冷却係数、H₂ 加熱の臨界密度、[C II]・[O I] の衝突係数、[O II]・[O III] のボルツマン因子です。組成と電子密度、分子の自己遮蔽、準位分布は積分中の状態から計算します。 [温度係数表の実装](../src/TemperatureTable.c) · [非平衡化学の時間発展と表引きの検証](../test/README.md#非平衡化学の検証計算と描画)

## 参考文献

各節の出典に対応する原典 PDF、Cloudy の比較用リリース論文とモデル全体の文献は [EQCLib/reference/index.tsv](../reference/index.tsv) に取得元とともに整理しています。

### 教科書

- Draine (2011), *Physics of the Interstellar and Intergalactic Medium* — 星間ガスの電離・再結合、ダスト、輝線過程。 [Princeton University Press](https://press.princeton.edu/books/paperback/9780691122144/physics-of-the-interstellar-and-intergalactic-medium)
- Tielens (2005), *The Physics and Chemistry of the Interstellar Medium* — 分子化学、粒子表面反応、光解離領域。 [Cambridge University Press](https://www.cambridge.org/core/books/physics-and-chemistry-of-the-interstellar-medium/B71008A57D96A2875E851DAEA70B16AE)

[asplund]: https://doi.org/10.1146/annurev.astro.46.060407.145222
[bt94]: https://adsabs.harvard.edu/pdf/1994ApJ...427..822B
[cen]: https://doi.org/10.1086/191630
[db96]: https://arxiv.org/abs/astro-ph/9603032
[draine]: https://www.astro.princeton.edu/~draine/book/
[field]: https://www.astro.rug.nl/~saleem/courses/EoRCourse/Field1958.pdf
[fischer]: https://doi.org/10.1016/j.adt.2004.02.001
[ga08]: https://arxiv.org/abs/0803.1768
[gc12]: https://doi.org/10.1111/j.1365-2966.2011.20260.x
[gf12]: https://arxiv.org/abs/1111.6980
[gj07]: https://arxiv.org/abs/0705.0182
[gml]: https://arxiv.org/abs/astro-ph/0605120
[gs07]: https://arxiv.org/abs/astro-ph/0608181
[g10]: https://arxiv.org/abs/0907.4081
[hm79]: https://adsabs.harvard.edu/pdf/1979ApJS...41..555H
[hm89]: https://ui.adsabs.harvard.edu/abs/1989ApJ...342..306H/abstract
[hui]: https://arxiv.org/abs/astro-ph/9612232
[kc08]: https://lweb.cfa.harvard.edu/~eketo/StarFormation/lowmass/keto_caselli_2008.pdf
[kingdon]: https://ui.adsabs.harvard.edu/abs/1996ApJS..106..205K/abstract
[kisielius]: https://doi.org/10.1111/j.1365-2966.2009.14989.x
[kroupa]: https://arxiv.org/abs/astro-ph/0009005
[lee96]: https://adsabs.harvard.edu/pdf/1996A%26A...311..690L
[ls83]: https://adsabs.harvard.edu/abs/1983ApJ...270..578L
[lte]: https://doi.org/10.1051/aas:2000265
[mckee82]: https://adsabs.harvard.edu/abs/1982ApJ...259..647M
[millar]: https://ui.adsabs.harvard.edu/abs/1997A%26AS..121..139M/abstract
[nahar-c]: https://www.astronomy.ohio-state.edu/nahar.1/papers/apjs-pirc-cn.pdf
[nahar-o]: https://www.astronomy.ohio-state.edu/nahar.1/papers/apjs-pirc-oions.pdf
[nist]: https://doi.org/10.18434/T4W30F
[nl99]: https://doi.org/10.1086/307823
[oii-data]: https://adsabs.harvard.edu/abs/1982MNRAS.198..111Z
[pegase]: https://arxiv.org/abs/astro-ph/9912179
[peters]: https://arxiv.org/abs/1610.06569
[pyneb]: https://doi.org/10.3390/atoms8040066
[sk87]: https://adsabs.harvard.edu/abs/1987ApJ...318...32S
[steffen]: https://adsabs.harvard.edu/pdf/1990A%26A...239..443S
[storey]: https://doi.org/10.1093/mnras/stu777
[storey-zeippen]: https://doi.org/10.1046/j.1365-8711.2000.03184.x
[verner]: https://arxiv.org/abs/astro-ph/9601009
[visser]: https://arxiv.org/abs/0906.3699
[voronov]: https://www.pa.uky.edu/~verner/col.html
[wd01]: https://arxiv.org/abs/astro-ph/0105237
[wolfire95]: https://public.websites.umich.edu/~monnier/Courses/2003Fall_AY403_Pages/wolfire95.pdf
