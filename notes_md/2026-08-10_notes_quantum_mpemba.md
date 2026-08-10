# 量子ムペンバ効果 精読ノート — 3論文の記法対応

作成日: 2026-08-10
対象レポート: `notes_tex/report_quantum_mpemba/report_quantum_mpemba_ja.tex`

## 0. 対象論文と参照した一次資料

| 略称 | 書誌 | 参照した一次資料 |
|---|---|---|
| **LR** | Lu & Raz, *PNAS* **114**, 5083 (2017) | arXiv:1609.05271 absページ + ar5iv本文 |
| **CLL** | Carollo, Lasanta & Lesanovsky, *PRL* **127**, 060401 (2021) | arXiv:2103.05020 absページ + ar5iv本文 |
| **WSW** | Wang, Su & Wang, arXiv:2410.06669 | `references/arXiv-2410.06669v3/` のソース配布（`main.tex`, `bib.bib`） |

WSWは2026-08-10時点でjournal-refなし（v3: 2025-11-26）。提出前に再確認すること。

---

## 1. 記法対応表（本ノートの中核）

レポートでは統一記号として**状態 $\varrho$**、**生成子 $\mathcal{G}$** を導入した。
定義は `notes_tex/report_quantum_mpemba/macros.tex` に集約してある。

### 1.1 スペクトル分解まわり

| 概念 | LR原論文 | CLL原論文 | 統一記法 | 備考 |
|---|---|---|---|---|
| 状態 | $\vec p(t)$ | $\rho_t$ | $\varrho(t)$ | 確率ベクトル / 密度行列 |
| 生成子 | $R$（$R_{ij}$） | $\mathcal{L}$ | $\mathcal{G}$ | LRは行列、CLLは超演算子 |
| 運動方程式 | $\dot p_i=\sum_j R_{ij}p_j$ | $\dot\rho_t=\mathcal{L}[\rho_t]$ | $\dot\varrho=\mathcal{G}\varrho$ | |
| 定常状態 | $\vec\pi(T_b)$ | $r_1\;(=\rho_{\rm ss})$ | $\varrho_{\rm ss}$ | **CLLは$r_1$と書く**（$\pi$を使わない） |
| 右固有モード | $\vec v_k$ | $r_k$ | $r_k$ | CLLに合わせた |
| 左固有モード | （明示せず） | $\ell_k$ | $\ell_k$ | LRは左固有ベクトルを陽に書かない |
| 双対規格化 | — | ${\rm Tr}(\ell_k r_h)=\delta_{kh}$ | $\langle \ell_k,r_h\rangle=\delta_{kh}$ | |
| 重なり係数 | $a_k$ | ${\rm Tr}(\ell_k\rho_0)$ | $a_k=\langle\ell_k,\varrho(0)\rangle$ | **議論全体の中心量** |
| 固有値の順序 | $\lambda_1>\lambda_2\ge\cdots$ | $\lambda_k$（$\lambda_1=0$） | $\lambda_1=0>{\rm Re}\lambda_2\ge\cdots$ | LRは実数、CLLは複素 |
| 緩和時間 | — | $\tau=1/\lvert\lambda_2\rvert$ | 同左 | |

**注意した点**

- LRの$\lambda_1>\lambda_2\ge\cdots$は$\lambda_1=0$で残りが負。したがって$\lambda_2>\lambda_3$は「$\lvert\lambda_2\rvert<\lvert\lambda_3\rvert$」＝第2モードが最も遅い、の意。符号の向きを取り違えやすい。
- LRは左固有ベクトルを陽に導入せず、$\vec p(0)=\vec\pi+\sum_k a_k\vec v_k$ の展開係数として$a_k$を定義する。CLLの${\rm Tr}(\ell_k\rho_0)$と同一だが、書き方が違うだけ。統一記法では$\langle\ell_k,\cdot\rangle$に揃えた。

### 1.2 距離・判定条件

| 概念 | LR | CLL | 統一記法での扱い |
|---|---|---|---|
| 距離 | 明示的な距離関数は主役でない | **Hilbert–Schmidt距離** $\mathcal{E}_t(\rho,\rho_{\rm ss})=\big[{\rm Tr}(e^{t\mathcal{L}}[\rho]-\rho_{\rm ss})^2\big]^{1/2}$ | CLLの定義をそのまま採用 |
| 判定条件 | $\lambda_2>\lambda_3$ かつ $\lvert a_2^c\rvert>\lvert a_2^h\rvert$ | ${\rm Tr}(\ell_2 U\rho_0U^\dagger)=0$ | 前者は受動的、後者は能動的 |
| 加速率 | — | $e^{t(\lvert{\rm Re}\lambda_3\rvert-\lvert{\rm Re}\lambda_2\rvert)}$ | |

### 1.3 WSW（統一枠組みの外側）

WSWには生成子が存在しないため、上表と対応する量がない。対応表の最下段は空欄になる。

| 概念 | WSW原論文 | 統一記法 |
|---|---|---|
| 系 / 浴 | Majorana $\chi_i$（$N$個）/ $\psi_i$（$N^2$個） | 同左 |
| 反交換関係 | $\{\chi_i,\chi_j\}=\delta_{ij}$ | 同左（$2\delta_{ij}$ではない） |
| SYK結合 | $J$、分散 $\overline{J^2}=3!J^2/N^3$ | $J$ |
| 系–浴結合 | $V$、分散 $\overline{V^2}=n!V^2/N^{2n}$ | $V$ |
| 参照スケール | $J_0$、**$J_0=2J$** | $J_0$ |
| 結合の次数 | $n$（本文の図はすべて $n=3$） | $n$ |
| 実効逆温度 | $\beta(t)=\dfrac{2\,{\rm Im}G_K(\omega,t)}{\omega\,{\rm Im}G_R(\omega,t)}\Big\vert_{\omega\to0}$ | $\beta_{\rm eff}(t)$ |
| 分布関数 | $f_{\rm eff}(\omega,t)$, $f_{\rm eq}(\omega)$ | 同左 |
| 浴からの距離 | $D_{\rm bath}(t)=\int d\omega\,\lvert f_{\rm eff}-f_{\rm bath}\rvert$ | $D_{\rm bath}(t)$ |
| 同温平衡からの距離 | $D(t)=\int d\omega\,\lvert f_{\rm eff}(\omega,\beta_{\rm eff})-f_{\rm eq}(\omega,\beta_{\rm eff})\rvert$ | $D_{\rm eq}(t)$ |

**$J_0=2J$ は図を読むときに必須**。図の横軸 $J_0t$、結合 $V/J_0$ はすべてこの規約。

---

## 2. 各論文の主張（一次資料で確認済み）

### LR
- Markov緩和で、定常への接近の遅さは最緩和モードとの重なり $\lvert a_2\rvert$ が支配。
- 十分条件: $\lambda_2>\lambda_3$ かつ $\lvert a_2^c\rvert>\lvert a_2^h\rvert$ → 高温側が速く冷える。
- 逆ムペンバ効果（加熱）: $\lvert a_2^h\rvert>\lvert a_2^c\rvert$。
- 例示模型: **非対称二重井戸上の1次元Fokker–Planck** と **三状態系2種**。
  → レポートの最小模型に三状態系を選んだのは原論文と直結している。

### CLL
- $e^{t\mathcal{L}}[\rho_0]=r_1+\sum_{k\ge2}e^{t\lambda_k}{\rm Tr}(\ell_k\rho_0)r_k$。
- **緩和開始直前にユニタリ $U$ を施して ${\rm Tr}(\ell_2U\rho_0U^\dagger)=0$ を作る**制御プロトコル。
- 適用例: **Dicke模型**、全対全相互作用スピン系（一次転移近傍の準安定領域の回避）。

### WSW
- 弱結合: $\beta_{\rm eff}(t)\simeq\beta_f+\alpha e^{-\Gamma t}$（Eberlein 2017）で単調、交差なし。
- 強結合: $V$ が閾値を超えると振動＋MPC。さらに強いと**負の実効温度**。閾値は初期条件に依存しない。
- **Lindblad SYKではMPCが出ない**。無限温度でもGreen関数が $\propto J$ の有限幅を持つのに対し、Lindbladはこれをdelta関数に潰すため。
- 二浴の温度差を大きくすると閾値はむしろ**増大**（可積分系の既知の結果と逆）。

---

## 3. 照合で判明した注意点

### 3.1 自分の初期理解が誤っていた点

| 項目 | 誤 | 正 |
|---|---|---|
| CLLの距離 | トレース距離 | **Hilbert–Schmidt距離** |
| CLLの機構 | $a_2=0$となる初期状態を受動的に選ぶ | **ユニタリ変換で能動的に$a_2$を消す制御プロトコル** |
| CLLの適用例 | 量子ドット等 | **Dicke模型、全対全スピン系** |

### 3.2 WSWの `bib.bib` の記載誤り

- `maldacena2016bound` が `pages={1--17}` になっている。正しくは論文番号 **106**（DOI: `10.1007/JHEP08(2016)106`）。
  → WSWのbibを転記すると誤りが伝播する。**Crossrefで照合すること**。
- LRはarXiv題名（"Anomalous cooling and heating—the Mpemba effect and its inverse"）とPNAS題名が異なる。引用はPNAS題名を用いた。

### 3.3 概念上の注意（レポートの論旨に関わる）

**CLLは厳密には「温度」の話ではない。**
CLLの$\rho_0$は必ずしも異なる温度の熱平衡状態ではなく、ユニタリで回した任意の初期状態である。
ムペンバ効果との関係は*着想の類似*として述べられている（abstractも "inspired by" という表現）。
したがって「LRとCLLをまとめて初期条件機構」と括る際は、

- 共通するのは **$a_2$が緩和速度を支配するという数学的構造**
- 異なるのは **$a_2$を何で動かすか**（LR: 初期温度 / CLL: ユニタリ）

という切り分けを明示する必要がある。レポート本文ではCLLを「制御プロトコル」と呼んで区別してあるが、この非対称性は Discussion で触れる価値がある。

**$a_2$は状態に対して定義され、$\beta_{\rm eff}$は温度計の読み。**
この非対称性がレポートのOpen Questionの出発点。LR/CLLが実効温度の定義依存性の問題を免れているのは、$a_2$が状態そのものの量だから。

---

## 4. 数値計算の記録

### 4.1 三状態Markov模型（実施済み）

コード: `code/mpemba_three_state.py`（模型・スペクトル分解）, `code/run_three_state_markov.py`（ドライバ）
出力: `code/out/2026-08-10_154155_three_state_markov/`

**採用パラメータ**: $E=(0,1,5)$, $\Gamma=(0.3,2.0,0.001)$, $\beta_{\rm bath}=0.5$

パラメータ空間のグリッド走査で以下を順に課して選定した。

1. $\lambda_2>\lambda_3$（Lu–Razの第一条件）
2. $a_2(\beta_0)$ が $(0,\beta_{\rm bath})$ 内に**第二の零点**をもつ（強ムペンバ）
3. $0.05\le|\lambda_2|\le5$ — **これが重要**。この制約を課さないと $|\lambda_2|\sim10^{-3}$ の準安定解が「ギャップ比が発散するので最良」として上位を占めてしまう。緩和がほぼ起きない病的な解なので除外が必要
4. $3\le\lambda_3/\lambda_2\le25$

**結果**

| 量 | 値 |
|---|---|
| $\lambda_2,\ \lambda_3$ | $-0.600$, $-7.575$（比 12.6） |
| $\beta_0^*$（第二の零点） | $0.1246$ |
| $\beta_0^h=0.100$ | $a_2=-0.0110$, $d(0)=0.388$, $t_{\rm conv}=3.93$ |
| $\beta_0^c=0.294$ | $a_2=+0.0330$, $d(0)=0.172$, $t_{\rm conv}=5.76$ |
| $\beta_0=\beta_0^*$ | $a_2=0$, $d(0)=0.358$, $t_{\rm conv}=0.78$ |
| 距離曲線の交差時刻 | $t=0.305$ |

$t_{\rm conv}$ は $\|p(t)-\pi\|_1<10^{-3}$ に到達する時刻。
高温側は初期距離が**遠い**のに先に到達する＝ムペンバ効果。$\beta_0^*$ では $\lambda_3$ 支配で 5 倍速い＝強ムペンバ効果。

**実装上の注意**

- 詳細釣り合いより $R$ は対称行列 $S=D^{-1}RD$（$D={\rm diag}(\sqrt\pi)$）と相似。`eigh` を使うと固有値が実数で得られ、左右固有ベクトルの双正規化も自動で満たされる（$v_k=Du_k$, $\ell_k=D^{-1}u_k$）。一般の `eig` を使う理由はない
- $a_2$ の零点探索は $a_2\propto \ell_2\cdot e^{-\beta_0 E}$（$Z>0$）を使えば規格化を省けて行列積1回になる。素朴なループ実装から約4000倍速くなり、走査が現実的になった

### 4.2 共鳴準位模型（実施済み）

コード: `code/resonant_level.py`, `code/run_resonant_level.py`
出力: `code/out/2026-08-10_161046_resonant_level/`

**設定**: $H=\varepsilon_d d^\dagger d+\sum_k\varepsilon_k c_k^\dagger c_k+\sum_k(V_k d^\dagger c_k+{\rm h.c.})$、
Lorentz型混成 $\Gamma(\omega)=\Gamma_0\Lambda^2/(\omega^2+\Lambda^2)$、$\varepsilon_d=1$, $\Gamma_0=1$, $\beta_{\rm bath}=0.25$, $\beta_0\in\{0.02,0.20\}$

**厳密解の構造（実装の要）**

$\Delta^R(\omega)=\Gamma_0\Lambda/(\omega+i\Lambda)$ なので $G^R$ は有理関数になり、極は2つだけ。
したがって伝播振幅 $u(t)=\sum_j c_j e^{-iz_jt}$ は**厳密に2指数の和**、記憶核 $A(t,\omega)=\int_0^t u(t-s)e^{-i\omega s}ds$ も解析的。
微分方程式を数値積分する必要はまったくない。Cini公式と組み合わせて2時間相関を組み立てる。

**結果**

| 予言 | 結果 |
|---|---|
| (i) 広帯域($\Lambda=20$)で単調・交差なし | **確認**。$\beta_{\rm eff}\to0.24997$（$\beta_{\rm bath}=0.25$と4桁一致）、窓依存性$<10^{-4}$ |
| (ii) 狭帯域($\Lambda=1$)＋強結合で振動 | **確認**。負の実効温度の領域も出現 |
| (iii) 閾値が初期温度に依存しない | **不成立**。閾値を定義できない |

(iii)の内訳: 軌跡の隔たり最大 $0.16$ に対し、FDTフィット窓（$2.0$–$4.0$）を変えたときの系統的な幅は中央値 $0.043$、最大 $1.32$。
`separation_exceeds_spread_anywhere: False`、つまり**交差はどの時刻でも系統誤差に埋もれる**。

一方、状態量である $D_{\rm bath}(t)$（スペクトル重み付き）は**両方の帯域幅で単調減少**。
Wangら Fig.3(b) と同じ「状態は単調・温度計は振動」という構造を、相互作用も無秩序もカオスもない自由フェルミオン系で再現した。

**実装で踏んだ落とし穴（再発防止）**

1. **周波数グリッドの範囲**。裾が $\Gamma_0\Lambda^2/(\pi\omega^4)$ なので打ち切り誤差は $\sim W^{-3}$。最初カットオフを `60×bandwidth` にしたら $\Lambda=40$ で範囲2400・間隔0.8となり、$e^{-i\omega t}$（$t\le60$）を解像できず完全性が0.83破れた。**範囲不足は雑音ではなく滑らかで完全に誤った答えを出す**ので気づきにくい
2. **$\Lambda=\infty$ は使えない**。$\Gamma(\omega)$ が減衰しないので裾積分がどのグリッドでも収束しない。Markov極限は「大きな有限 $\Lambda$」で表現する。コード側で例外を投げるようにした
3. **メモリ**。$A$ は $N_t\times N_\omega$。自己診断を $t_{\max}=400$, $\Lambda=40$ で書いたら16 GBになりプロセスが落ちた（出力なしで終了するので原因が分かりにくい）。ガードを入れた
4. **FDTの符号規約を推測しない**。半Fourier変換の結果、比は $+\tanh(\beta\omega/2)$ で $\beta=2\times$傾き。当初 $-2\times$ としていて符号も倍率も誤っていた。平衡での自己診断（$\beta_{\rm eff}\to\beta_{\rm bath}$）で確定させるのが正しい手順
5. **相対時間窓 `rel_max` を全時刻で揃える**。半Fourier変換を $t'_{\max}$ で打ち切ると周波数方向に周期 $2\pi/t'_{\max}$ のリンギングが出る。時刻ごとに窓が違うと系統誤差も違い、**$\beta_{\rm eff}$ の見かけの時間依存性がその人工物に支配される**
6. **線形フィットではなく $\tanh$ フィット**。原点を通る直線で近似すると幅 $W$ の窓で $\beta$ を約 $(\beta W)^2/20$ だけ過小評価する。パーセント級のバイアスで、物理と誤認しうる

**本質的な制約（物理側）**

$\beta_{\rm eff}$ の抽出には $t'\lesssim2T$ が必要な一方、リンギング周期 $\pi/T$ が熱的スケール $1/\beta$ より十分小さいこと、すなわち $T\gg\pi\beta$ が要る。
過渡現象は $t\sim1/\Gamma_0$ で終わるので、$\Gamma_0\beta\gtrsim1$ だと過渡領域で実効温度が定義できない。
$\beta_{\rm bath}=0.25$（高温浴）にしたのはこの制約を回避するため。
**この制約はWangらのSYK計算にも同様に効くはずだが、原論文では議論されていない。**

---

## 5. 未確認・要フォロー

- [ ] WSWの掲載状況（提出前に再確認）
- [ ] WSWのSM（`SM_remove_doc_class.tex`）は未読。二浴の場合とLindblad SYKの導出、別定義の実効温度による確認がここにある
- [ ] LR / CLL の本文PDFは未取得（ar5iv経由で式のみ確認）。数値パラメータを引用する場合は要確認
- [ ] CLLがDicke模型で使った具体的な$\mathcal{L}$のパラメータ（最小模型の設定に流用できる可能性）
