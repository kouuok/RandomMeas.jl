# 元論文の解説 — この研究が使っている部分、広げた部分、食い違って見える部分

[README.md](README.md) に戻る / 論文全体の入門: [crm_shadow_kaisetsu.pdf](crm_shadow_kaisetsu.pdf)(CRMシャドウ入門)

**元論文**: B. Vermersch, A. Rath, B. Sundar, C. Branciard, J. Preskill, A. Elben, "Enhanced estimation of quantum properties with common randomized measurements", *PRX Quantum* **5**, 010352 (2024). [doi:10.1103/PRXQuantum.5.010352](https://doi.org/10.1103/PRXQuantum.5.010352) / [arXiv:2304.12292](https://arxiv.org/abs/2304.12292)

この文書は論文全体の紹介ではない(それは上の入門にある)。**この研究が論文のどの式・どの主張を使っているか、どこを広げたか、どこが食い違って見えるか**に絞って対応を付ける。

**式番号は2つの版を併記する。** arXiv 版(v1、2023年4月24日投稿)と PRX Quantum 版(2024年3月28日出版)では、付録の並びと式番号が違い、PRX 版には実験データの応用例が1つ加わっている(0章)。本文では「式(25) [PRX: (C7)]」のように、先に arXiv 版、角括弧の中に PRX 版の番号を書く。本文の式(1)〜(6)は両方の版で同じ番号である。

---

## 0. arXiv 版と PRX Quantum 版の対応

| 内容 | arXiv 版(v1) | PRX Quantum 版 |
|---|---|---|
| 古典シャドウの定義 | 式(1) | 式(1)(II.A 節) |
| CRM シャドウの定義 | 式(2)、(3) | 式(2)、(3)(II.B 節) |
| 単一パウリ文字列の分散の上界 | 式(4) | 式(4)(II.C 節) |
| 多コピー観測量の推定量と上界 | 式(5)、(6) | 式(5)、(6)(II.D 節) |
| 実験データによるエンタングルメントの検出( $p_3$ -PPT、イオントラップ10量子ビット) | **なし** | **III.A 節**(Fig. 1(b)) |
| フォン・ノイマン・エントロピーの多項式近似 | 例1 | III.B 節、式(7) |
| 忠実度の推定と prior の反復改善 | 例2 | III.C 節 |
| 別の実験(companion)から prior を作る | 付録 D | **付録 A**(式(A1)) |
| 分散の一般式(全分散の分解) | 付録 A.1、式(15)〜(17) | **付録 B.2**、式(B8)〜(B10) |
| バッチシャドウの分散の展開 | 式(10) | 式(B3) |
| 単一パウリ文字列の**厳密な**分散 | 付録 B.1、**式(25)** | **付録 C.1、式(C7)** |
| 一般の多コピー観測量の上界 | 付録 B.2 | 付録 C.2、式(C16) |
| シフト演算子の $O_A^{(1)}$ | 付録 C | 付録 E |
| エントロピーの多項式近似の係数 | 付録 E | 付録 D |

付録の並びは、arXiv 版の A, B, C, D, E が PRX 版の B, C, E, A, D に対応する。

---

## この文書の要点

1. **論文の中心は CRM シャドウ $\hat\rho_\sigma^{(r)}=\hat\rho^{(r)}-\sigma^{(r)}+\sigma$ (式(2))である。** 実験の影から「prior $\sigma$ を同じ基底で測ったときの影」を引き、 $\sigma$ そのものを足し戻す。 $\sigma$ が何であっても不偏である。
2. **単一パウリ文字列の分散は、論文の付録に厳密な等式として書かれている**(arXiv 版の付録 B.1 の式(25)、PRX 版の付録 C.1 の式(C7))。 この研究の「利得法則」 $G=\frac{(3^{\lvert A\rvert}-1)x^2+v_s}{(3^{\lvert A\rvert}-1)\Delta^2+v_s}$ は、その式を「標準のシャドウ( $\sigma$ の項を落とす)」と「CRM」で割った比そのもので、**新しい法則ではない**。この研究の文書は、以前は上界の式(4)に触れているだけで厳密な式(25) [PRX: (C7)] への帰属を書いていなかったが、2026年10月に書き足した(5章)。損得の境目 $\lvert\Delta\rvert\le\lvert x\rvert$ も論文の本文(II.C 節)に書かれている。
3. **論文は、CRM の効果を「測定の設定の数 $N_U$ が有限であることによる誤差」を減らすものと位置づけ、較正に時間がかかる実験で特に重要だと述べている。** この研究の利得 $G$ も、同じ $n_u$ ・同じ $n_m$ での分散比で、つまり**設定の数を何分の1にできるか**を測っている。1基底1ショット( $n_m=1$ )では天井が約 $1/(1-x^2)$ まで下がる(6章)。
4. **論文が忠実度を prior の良さの指標に使うのは、多コピー観測量と純粋な prior の文脈である。** 多コピー観測量の上界(式(6))は Hilbert–Schmidt 距離で書かれ、純粋な prior なら忠実度 $1/2$ 以上で標準より良くなる。この研究の「局所観測量では忠実度ではなく $\Delta$ が効く」は、その文脈の外(単一コピーの局所観測量、台に縮約すると混合になる prior)の話で、論文と矛盾しない(7章)。
5. **論文のエントロピーの数値例(arXiv 版の例1、PRX 版の III.B 節)の prior は、厳密な MPS を小さい結合次元に切り詰めたものである。** この研究の「切断 MPS」と同じ作り方で、この研究はそれが $\chi=4$ で電荷のゆらぎを捨てることを見つけた(8章)。
6. **論文が扱っていないもの**: 対称性を破った prior(UHF)とその群平均、観測量ごとの違い(電荷とスピン)、和の観測量での項ごとの誤差、距離依存、読み出し誤差に対する prior 側の補正。この研究の寄与はこれらにある(9章)。
7. **論文の展望**(重点サンプリングや適応的な手法との組み合わせ、物質の相の探査への応用、prior を反復で改善する手順)は、この研究の次の方針にそのままつながる(10章)。

---

## 1. 論文の設定と、この研究の記号との対応

### 1.1 ランダム測定のプロトコル

(PRX 版の II.A 節) $N$ 量子ビットの状態 $\rho$ に局所ランダムユニタリー $U=\bigotimes_{i=1}^NU_i$ を掛けて計算基底で測る。各 $U_i$ は

```math
U_i\in\Bigl\{\ \mathbb 1_2,\ \ \tfrac1{\sqrt2}\begin{pmatrix}1&1\\1&-1\end{pmatrix},\ \ \tfrac1{\sqrt2}\begin{pmatrix}1&-i\\1&+i\end{pmatrix}\ \Bigr\}
```

から一様に選ぶ(それぞれ $Z$ 、 $X$ 、 $Y$ の基底で測ることに当たる)。ユニタリーを $N_U$ 通り選び、それぞれで $N_M$ 回ずつ測る。 $r$ 番目のユニタリーでの古典的な影(classical shadow)は

```math
\hat\rho^{(r)}=\sum_{\mathbf s}\widehat P_\rho(\mathbf s\,|\,U^{(r)})\ \mathcal M^{-1}\bigl(U^{(r)\dagger}\lvert\mathbf s\rangle\langle\mathbf s\rvert U^{(r)}\bigr)\qquad(\text{式(1)})
```

で、 $\widehat P_\rho(\mathbf s\,|\,U)$ は $N_M$ 回の測定から数えた出力 $\mathbf s$ の頻度、 $\mathcal M^{-1}$ は測定チャネルの逆である。

### 1.2 記号の対応

| 論文 | この研究 | 意味 |
|---|---|---|
| $N_U$ | $n_u$ | 測定の設定(ランダムユニタリー)の数 |
| $N_M$ | $n_m$ | 1つの設定で測るショット数 |
| $N_A$ | $\lvert A\rvert$ | 観測量の台の量子ビット数 |
| $O=\gamma=\bigotimes_i\gamma_i$ ( $\gamma_i\in\{\mathbb 1,X,Y,Z\}$ ) | $P$ | パウリ文字列 |
| ${\rm Tr}(\rho\gamma)$ | $x=\langle P\rangle_\rho$ | 真の値 |
| ${\rm Tr}(O\sigma)$ | $y=\langle P\rangle_\sigma$ | prior の値 |
| ${\rm Tr}(O(\rho-\sigma))$ | $\Delta=x-y$ | prior の外れ |

---

## 2. CRM シャドウの定義と不偏性

### 2.1 定義

```math
\hat\rho_\sigma^{(r)}=\hat\rho^{(r)}-\sigma^{(r)}+\sigma\qquad(\text{式(2)})
```

```math
\sigma^{(r)}=\sum_{\mathbf s}P_\sigma(\mathbf s\,|\,U^{(r)})\ \mathcal M^{-1}\bigl(U^{(r)\dagger}\lvert\mathbf s\rangle\langle\mathbf s\rvert U^{(r)}\bigr)\qquad(\text{式(3)})
```

- $\hat\rho^{(r)}$ : 実験で得た影(式(1))。
- $\sigma^{(r)}$ : **同じユニタリー $U^{(r)}$ で prior $\sigma$ を測ったときの影**。ただし $P_\sigma(\mathbf s\,|\,U)$ は頻度ではなく、**古典計算で厳密に求めた出力の確率**である(ハットが付いていない)。
- $\sigma$ : prior そのもの。

**実験と同じくじ(同じユニタリー)を prior にも引かせて、その結果を差し引く**のが「共通の乱数(common random numbers)」という名前の由来である。

### 2.2 不偏性

ユニタリーについて平均すると $\mathbb E[\hat\rho^{(r)}]=\rho$ 、 $\mathbb E[\sigma^{(r)}]=\sigma$ なので

```math
\mathbb E\bigl[\hat\rho_\sigma^{(r)}\bigr]=\rho-\sigma+\sigma=\rho
```

で、**$\sigma$ をどう選んでも推定は偏らない**。prior が悪くても分散が増えるだけで、答えはずれない。

### 2.3 この研究の実装との対応 — prior 側は標本にしない

$\sigma^{(r)}$ は厳密な確率 $P_\sigma$ で作る。これを prior のサンプリングで代用すると、prior 側にもショット雑音が乗って分散が増える。この研究の指針「**prior 側は絶対に標本推定しない**」([SUMMARY.md](SUMMARY.md) 実務上の指針3)は、この定義をそのまま守るということである。

単一パウリ文字列 $P$ (台 $A$ )では、1つの設定での推定値は次のように書ける。基底が $P$ の各因子と一致したとき( $P$ の台のすべての量子ビットで基底が合ったとき、確率 $3^{-\lvert A\rvert}$ )だけ値を持ち、

```math
{\rm Tr}\bigl(P\hat\rho^{(r)}\bigr)=3^{\lvert A\rvert}\,\mathbb 1_{\rm 一致}\,\overline{m},\qquad
{\rm Tr}\bigl(P\sigma^{(r)}\bigr)=3^{\lvert A\rvert}\,\mathbb 1_{\rm 一致}\,y
```

( $\overline m$ は $N_M$ 回の測定での $\pm1$ の出力の積の平均)。したがって CRM の推定値は

```math
{\rm Tr}\bigl(P\hat\rho_\sigma^{(r)}\bigr)=3^{\lvert A\rvert}\,\mathbb 1_{\rm 一致}\,(\overline m-y)+y
```

である。この研究の [README_gainlaw.md](README_gainlaw.md) の模擬はこの式をそのまま使っている。

---

## 3. 分散の一般形 — 「基底のくじ」と「ショット雑音」

付録の「General variance formula」(arXiv 版の付録 A.1、PRX 版の付録 B)は、1つの設定の分散を2つに分ける:

```math
\mathbb V_1=\mathbb V_U\bigl[f_{\rho,\sigma}(U)\bigr]+\frac{\mathbb E_U\bigl[g_\rho(U)\bigr]}{N_M}\qquad(\text{式(15) [PRX: (B8)]})
```

```math
f_{\rho,\sigma}(U)=\sum_{\mathbf s}\bigl(P_\rho(\mathbf s|U)-P_\sigma(\mathbf s|U)\bigr)\bigl[\mathcal M^{-1}(O^{(1)})\bigr](U,\mathbf s)\qquad(\text{式(16) [PRX: (B9)]})
```

```math
g_\rho(U)=\sum_{\mathbf s}P_\rho(\mathbf s|U)\bigl([\mathcal M^{-1}(O^{(1)})](U,\mathbf s)\bigr)^2-\Bigl(\sum_{\mathbf s}P_\rho(\mathbf s|U)[\mathcal M^{-1}(O^{(1)})](U,\mathbf s)\Bigr)^2\qquad(\text{式(17) [PRX: (B10)]})
```

- **第1項**: どのユニタリー(基底)が当たるかのくじによるゆらぎ。 $f$ は「 $\rho$ と $\sigma$ の出力分布の差」だけで決まるので、**prior が良ければ小さくなる**。
- **第2項**: 1つの基底の中でのショット雑音。 $g$ は $\rho$ だけで決まり、**prior では減らせない**。

これはこの研究の [README_eigenstate.md](README_eigenstate.md) 1.6 の「全分散の法則」による分解(基底のくじの雑音 $Kx^2$ と、減らせない床 $v_s$ )と同じものである。

**単一パウリ文字列での具体的な計算。** パウリ文字列 $P$ (台 $A$ )について、 $[\mathcal M^{-1}(P)](U,\mathbf s)$ は「基底が台の上で全部一致したとき $3^{\lvert A\rvert}\times(\text{台の上の測定値の積})$ 、それ以外は 0」である(2.3 節)。基底が一致したとき、台の上の測定値の積は $\rho$ のもとで平均 $x={\rm Tr}(P\rho)$ 、 $\sigma$ のもとで平均 $y={\rm Tr}(P\sigma)$ の $\pm1$ の確率変数である。したがって

- $f(U)=3^{\lvert A\rvert}\mathbb 1_{\rm 一致}\,(x-y)=3^{\lvert A\rvert}\mathbb 1_{\rm 一致}\,\Delta$ 。一致する確率は $3^{-\lvert A\rvert}$ なので、 $\mathbb E_U[f]=\Delta$ 、 $\mathbb E_U[f^2]=3^{-\lvert A\rvert}\cdot3^{2\lvert A\rvert}\Delta^2=3^{\lvert A\rvert}\Delta^2$ 。
- $g(U)$ は1ショットの値の分散で、一致したときは $3^{2\lvert A\rvert}\times(\pm1\text{ の分散})=3^{2\lvert A\rvert}(1-x^2)$ 、一致しないときは 0。

これらから

```math
\mathbb V_U[f]=\mathbb E_U[f^2]-\mathbb E_U[f]^2=3^{\lvert A\rvert}\Delta^2-\Delta^2=(3^{\lvert A\rvert}-1)\Delta^2,\qquad
\mathbb E_U[g]=3^{-\lvert A\rvert}\cdot3^{2\lvert A\rvert}(1-x^2)=3^{\lvert A\rvert}(1-x^2)
```

となり、次の式(25) [PRX: (C7)] がそのまま出る(PRX 版では式(C2)〜(C6)がこの計算に当たる)。

---

## 4. 単一パウリ文字列の分散 — 式(25) [PRX: (C7)] と式(4)

付録の「Estimating single-copy Pauli observables」(arXiv 版の付録 B.1、PRX 版の付録 C.1)は、 $O=\gamma$ がパウリ文字列のとき、推定量 $\hat O=\frac1{N_U}\sum_r{\rm Tr}(O\hat\rho_\sigma^{(r)})$ の分散を**厳密に**与えている:

```math
\mathbb V(\hat O)=\frac1{N_U}\Bigl((3^{N_A}-1)\,{\rm Tr}\bigl(O(\rho-\sigma)\bigr)^2+\frac{3^{N_A}\bigl(1-{\rm Tr}(\rho\gamma)^2\bigr)}{N_M}\Bigr)
\le\frac{3^{N_A}}{N_U}\Bigl({\rm Tr}\bigl[O(\rho-\sigma)\bigr]^2+\frac1{N_M}\Bigr)\qquad(\text{式(25) [PRX: (C7)]、右辺が本文の式(4)})
```

本文は式(4)の上界の形だけを示し、そのすぐ後で次のように述べている(要約):

- **標準のシャドウでは、 $\sigma$ を 0 に置き換えれば同じ式が成り立つ**(classical shadows の元論文 Huang–Kueng–Preskill の定理2と整合する)。
- **$N_M$ の値によらず、 $\lvert{\rm Tr}[O(\rho-\sigma)]\rvert\le\lvert{\rm Tr}(O\rho)\rvert$ なら CRM の分散は標準のシャドウより小さい**(両方の版の本文に明記)。この研究の損得の境目 $\lvert\Delta\rvert\le\lvert x\rvert$ ( $\varepsilon\le1$ 、 $0\le y/x\le2$ )はこの条件そのものである。
- 統計誤差は、有限の $N_U$ と、1設定あたり有限の $N_M$ の両方から来る。**CRM が減らすのは有限の $N_U$ による分散で、これは較正に時間がかかる実験(イオントラップや超伝導量子ビット)で特に重要である。そこでは $N_U$ が限られる一方、 $N_M$ は大きく取れる**(PRX 版の II.C 節)。

---

## 5. この研究の「利得法則」は、式(25) [PRX: (C7)] の比である

式(25) [PRX: (C7)] を、 $\sigma$ の項を落とした標準のシャドウ( $\Delta\to x$ )と CRM で割ると、 $1/N_U$ が約分されて

```math
G=\frac{\mathbb V_{\rm 標準}}{\mathbb V_{\rm CRM}}=\frac{(3^{\lvert A\rvert}-1)x^2+v_s}{(3^{\lvert A\rvert}-1)\Delta^2+v_s},\qquad v_s=\frac{3^{\lvert A\rvert}(1-x^2)}{n_m}
```

となる。これがこの研究の [SUMMARY.md](SUMMARY.md) 「中心にある式」そのものである。

### 5.1 どこまでが論文の結果で、どこからがこの研究の結果か

| 内容 | 出どころ |
|---|---|
| 単一パウリ文字列の CRM と標準のシャドウの厳密な分散 | **論文**(式(25) [PRX: (C7)] と、本文の式(4)の直後の「 $\sigma$ を 0 に置き換える」) |
| 利得をその比 $G$ として書くこと、 $n_u$ が約分されること | 論文の式から直ちに出る(比を取っただけ) |
| 天井 $G_{\max}=1+Kx^2/v_s$ 、 $G=G_{\max}/(1+R)$ 、 $R=\varepsilon^2(G_{\max}-1)$ の形 | この研究(式(25) [PRX: (C7)] の書き換え) |
| 損得の境目 $\lvert\Delta\rvert\le\lvert x\rvert$ (この研究の書き方では $0<y/x<2$ ) | **論文**(本文 II.C 節) |
| 固有状態に近い prior が得をする条件 $\lvert x\rvert>1/2$ | この研究([README_eigenstate.md](README_eigenstate.md)、論文の条件からの代数) |
| 和の観測量の厳密な分散(項の積の期待値と $3^{\lvert A_k\cap A_l\rvert}$ の重み) | この研究([README_observables.md](README_observables.md) 5.1。論文の付録(arXiv 版 B.2、PRX 版 C.2)は一般の多コピー観測量の**上界**を与えている) |
| 実際の測定過程を模擬して式と突き合わせること(288点) | この研究([README_gainlaw.md](README_gainlaw.md)) |
| Hubbard 模型の具体的な prior(UHF、UHF-sym、PHF、切断・変分 MPS)での $\Delta$ の系統的な評価と、その物理的な説明 | この研究 |

### 5.2 この研究の文書の書き方(2026年10月に直した)

この研究の文書は、利得の式を「利得法則」と呼び、SUMMARY の主要結論1に「利得法則は約5.9桁で成り立つ」と書いていたが、**論文の厳密な式(25) [PRX: (C7)] への帰属が書かれていなかった**。式(25)は厳密な等式なので、実測との一致は「法則が物理として成り立つ」ことではなく、**この研究の実装(量子ビットの割り当て、JW 変換、推定量)が論文の式どおりに動いている**ことの確認である。確認としての価値はあるが、主張としては「元論文の式(25)(の比)を、この研究の実装と Hubbard 模型の系で確認した」と書くのが正確である。2026年10月に、SUMMARY・ROADMAP・論文原稿・各文書の書き方をそのように直した。この研究の独自の寄与は 5.1 の表の下半分にある。

---

## 6. $N_M$ の役割 — 利得は「設定の数」の削減を測っている

式(25) [PRX: (C7)] で、prior が減らせるのは第1項(基底のくじ)だけで、第2項(ショット雑音)は減らない。したがって、prior が完全( $\Delta=0$ )でも

```math
G\le G_{\max}=1+\frac{(3^{\lvert A\rvert}-1)\,x^2}{3^{\lvert A\rvert}(1-x^2)}\,n_m
```

で、**天井は $n_m$ にほぼ比例する**。論文の例も $N_M=150$ (実験データ、PRX 版のみ)、 $N_M=1000$ (エントロピー)、 $N_M=10^5$ (忠実度)と、1設定あたりのショットを多く取っている。

この研究の見出しの利得は $n_m=100$ での値で、 $n_m$ を変えると次のようになる( $U=12$ ):

| 観測量 | $n_m=1$ | 10 | 100 | 1000 |
|---|---|---|---|---|
| オンサイト ZZ(UHF prior、1次元 $L=12$ 周期境界) | **6.5** | 55 | 464 | 1879 |
| ブロックの電荷の偶奇 $P_{16}$ (台32、UHF prior、1次元 $L=64$ ) | **7.7** | 67 | 591 | 2862 |

**1基底1ショットの通常の classical shadow では、どれほど良い prior でも利得は約 $1/(1-x^2)$ (台が大きい極限)で頭打ちになる。** したがって $G=464$ は「測定の設定(回路)の数を約 464 分の1にできる」という意味で、「総ショット数を 464 分の1にできる」という意味ではない。これは論文自身の位置づけ(4章の2つめの点)と同じであり、この研究の見出しの数値にも併記すべきである。基底の切り替えにかかるコストを入れた最適な $n_m$ は [SUMMARY.md](SUMMARY.md) 実務上の指針4にある。

---

## 7. 多コピー観測量と忠実度 — 論文が忠実度を使う文脈

### 7.1 多コピー観測量の推定と上界

純度やエンタングルメントの指標のような多コピー観測量 ${\rm Tr}[O\rho^{\otimes n}]$ は、異なる設定の影を組み合わせた U 統計量で推定する:

```math
\hat O=\frac{(N_U-n)!}{N_U!}\sum_{t_1\ne\cdots\ne t_n}{\rm Tr}\Bigl[O\bigl(\hat\rho_\sigma^{[t_1]}\otimes\cdots\otimes\hat\rho_\sigma^{[t_n]}\bigr)\Bigr]\qquad(\text{式(5)})
```

その分散の上界は

```math
\mathbb V[\hat O]\le\frac{n^2\lVert O_A^{(1)}\rVert_2^2}{N_U}\Bigl(3^{N_A}\lVert\rho_A-\sigma_A\rVert_2^2+\frac{2^{N_A}}{N_M}\Bigr)+\mathcal O(1/N_U^2)\qquad(\text{式(6)})
```

で、**prior の良さは Hilbert–Schmidt 距離 $\lVert\rho_A-\sigma_A\rVert_2$ で入る**( $O_A^{(1)}$ は $\rho$ に依存する部分系 A の上の演算子、 $n$ はコピー数)。

### 7.2 忠実度 1/2 という条件

$\rho$ も $\sigma=\lvert\phi\rangle\langle\phi\rvert$ も純粋なら、 ${\rm Tr}\rho^2={\rm Tr}\sigma^2=1$ 、 ${\rm Tr}(\rho\sigma)=F$ ( $F$ は忠実度)なので

```math
\lVert\rho-\sigma\rVert_2^2={\rm Tr}\rho^2-2\,{\rm Tr}(\rho\sigma)+{\rm Tr}\sigma^2=2(1-F)
```

である。標準のシャドウは $\sigma$ の項を落とした( $\sigma\to0$ )ものなので、対応する量は $\lVert\rho\rVert_2^2=1$ である。式(6)の第1項が標準より小さくなる条件 $2(1-F)<1$ は、 $F>1/2$ である。論文の忠実度の例(arXiv 版の例2、PRX 版の III.C 節)の反復手順は、この「忠実度 $1/2$ 以上」を prior の合格ラインにしている(8章)。PRX 版の III.C 節は「式(6)を見れば、 $\mathcal F_\phi\ge1/2$ なら CRM の方が分散が小さい」と明記している。

### 7.3 この研究との関係 — 矛盾しない

この研究の中心的な主張の1つは「**局所観測量の利得は、prior の大域的な忠実度ではなく、観測量ごとの外れ $\Delta$ で決まる**」ことである([README_prior_survey.md](README_prior_survey.md)、[README_2m.md](README_2m.md))。これは論文と矛盾しない。

- 論文が忠実度を使うのは、**大域的な(または多コピーの)観測量と、純粋な prior** の文脈である。そこでは忠実度が式(6)の Hilbert–Schmidt 距離を正しく押さえる。
- この研究が扱うのは**単一コピーの局所観測量**で、効くのは式(25) [PRX: (C7)] の $\Delta$ だけである。
- さらに、prior を観測量の台に縮約すると一般に**混合状態**になり、純粋な prior で成り立っていた「忠実度が距離を押さえる」関係が崩れる([README_2m.md](README_2m.md) 準備4b〜4d)。

---

## 8. 論文の応用例と、この研究との対応

### 8.0 実験データによるエンタングルメントの検出(PRX 版の III.A 節のみ。arXiv 版にはない)

- **データ**: イオントラップ量子シミュレータ(10量子ビット)のクエンチ実験で測られたランダム測定データ( $N_U=500$ 、 $N_M=150$ )。
- **prior**: 実験のユニタリー時間発展を数値シミュレーションした純粋状態 $\sigma=\lvert\psi\rangle\langle\psi\rvert$ 。実験の状態はデコヒーレンスで混合状態になっているので、prior は近似にすぎない。
- **推定する量**: 2つの部分系の間のエンタングルメントを判定する $p_3$ -PPT 条件 $p_2^2/p_3>1$ ( $p_n$ は部分転置した密度行列のモーメント、多コピー観測量)。
- **結果**: 1標準偏差でエンタングルメントを確認するのに、CRM なら $N_U=53$ 、標準のシャドウなら $N_U=299$ で、**必要な設定数が約6分の1**になった(データ取得の時間が数時間短くなる)。

**この研究との対応**: 論文の中で唯一、実機のデータを使った例である。prior は実験の状態そのものではなく、ノイズのないシミュレーションでもよいことを示している。この研究は擬似実験(DMRG で作った参照状態)だけで、実機のデータや読み出し誤差は扱っていない(9章)。

### 8.1 フォン・ノイマン・エントロピーの多項式近似(arXiv 版の例1、PRX 版の III.B 節)

- **状態**: 臨界イジング鎖 $H=-\sum_i(Z_iZ_{i+1}+X_i)$ ( $N=16$ 量子ビット)の基底状態 $\lvert G\rangle$ を、結合次元 $\chi_G\approx40$ の MPS で表したもの。
- **prior**: $\lvert G\rangle$ を**結合次元 $\chi=1,2,3$ に切り詰めた MPS**。
- **推定する量**: 半分の部分系( $N_A=8$ )のモーメント $p_n={\rm Tr}[\rho_A^n]$ と、そこから作るエントロピーの多項式近似( $n_{\max}=3$ で $f_3(x)=\frac{137}{60}x-4x^2+\frac74x^3$ )。
- **設定**: $N_M=1000$ に固定して $N_U$ を変える。 $\chi=2,3$ で誤差が標準のシャドウより大きく減る。

**この研究との対応**: この研究の「切断 MPS」( $\chi_p$ )は、例1と同じ「厳密な状態を切り詰める」作り方である。この研究はさらに、

- 切り詰めた MPS が $\chi=4$ で**電荷のゆらぎを捨てる**(中央の結合の Schmidt 分解で、残る4状態がすべてスピンの成分になる)ため、オンサイトの電荷量では UHF に負けること、
- 2次元・周期境界・弱結合では $\chi_p=4,8$ でも損をする系があること、
- 厳密な状態なしで作れる**変分 MPS** は逆に電荷のゆらぎを正しく持つが、スピンの向きを選ぶこと

を示した([README_phf_mps.md](README_phf_mps.md))。例1のような切断 MPS の prior は「厳密な状態を知っている」ことが前提で、実際の実験の prior としては変分 MPS や平均場を考える必要がある。

### 8.2 忠実度の推定と、prior を反復で良くする手順(arXiv 版の例2、PRX 版の III.C 節)

- **状態**: ランダム量子回路で作った $N=30$ 量子ビットの状態(脱分極雑音 $p=10^{-3}$ 、 $10^{-4}$ を含む)。 $N_U=15$ 、 $N_M=10^5$ 。
- **推定する量**: 目標の純粋状態 $\lvert\psi\rangle$ との忠実度 $\mathcal F_\psi=\langle\psi\rvert\rho\lvert\psi\rangle$ ( $O=\lvert\psi\rangle\langle\psi\rvert$ )。
- **反復手順**:
  1. prior $\sigma=\lvert\phi\rangle\langle\phi\rvert$ から始め、CRM シャドウか標準のシャドウで $\mathcal F_\phi$ を推定する。
  2. $\mathcal F_\phi$ が閾値 $\mathsf F\ge1/2$ を下回ったら、新しい prior を作る(より多くの古典計算を使ってよい)。
  3. 1に戻る。十分に高い忠実度の prior が見つかったら、それを使って任意の多コピー観測量を高精度に推定する。

**この研究との対応**: 「測定データを見ながら prior を選び直す」という発想は、この研究の [README_2m.md](README_2m.md) の「同じデータから各 prior の $\hat\Delta$ を推定して選ぶ」手順につながっている(この研究では、局所観測量の選定には忠実度ではなく $\hat\Delta$ を使うのが正しいことを示した)。今後の方針「データから prior を作る」(10章)は、この反復手順を局所観測量向けに作り直すものである。

### 8.3 別の実験から prior を作る(arXiv 版の付録 D、PRX 版の付録 A)

prior $\sigma$ を古典計算で作る代わりに、**$\rho$ に近い状態を作る別の実験(companion experiment)** の classical shadow から作ることもできる。その実験が $\rho$ を正確に再現していなくても、より多くのユニタリーで測っておけば役に立つ、という議論である。この研究では扱っていない。

---

## 9. 論文が扱っていないもの — この研究の寄与

| テーマ | 論文 | この研究 |
|---|---|---|
| prior の種類 | 切り詰めた MPS、目標の純粋状態、別の実験 | UHF、UHF-sym(群平均)、RHF、PHF、切断 MPS、変分 MPS、変分 MPS の群平均 |
| 対称性を破った prior | 扱っていない | UHF はスピンの向きを選ぶので、スピンの量で損をする。群平均(UHF-sym)で損が消える([README_uhfsym.md](README_uhfsym.md)) |
| 観測量による違い | 扱っていない | 電荷の量は強結合で $\varepsilon\propto(t/U)^2$ で良くなり、スピンの量は飽和する([README_prior_survey.md](README_prior_survey.md) 結果2、[README_observables.md](README_observables.md)) |
| 和の観測量 | 一般の多コピー観測量の上界(arXiv 版の付録 B.2、PRX 版の付録 C.2) | 厳密な分散と、項ごとの誤差が効くこと(二重占有を4項の和で測ると UHF の利得は 1.8、群平均で 123) |
| 距離依存 | 扱っていない | 長距離秩序を持つ prior は遠くの相関で必ず損をし、群平均しても直らない([README_eigenstate.md](README_eigenstate.md)、[README_varsym_distance.md](README_varsym_distance.md)) |
| 量子的な相関を持つ prior の対称性の回復 | 扱っていない | 変分 MPS を群平均すると、電荷とスピンの両方で損をしない([README_varsym.md](README_varsym.md)) |
| 系の大きさと大域的な忠実度 | 例は16〜30量子ビット | 大域的な忠実度が指数的に落ちても局所観測量の利得は残る(128量子ビットまで) |
| 測定の誤差 | 忠実度の例は状態の用意に脱分極雑音を入れている。実験データの例(PRX 版)は実機のデコヒーレンスを含む | 読み出し誤差に対して prior を補正しないと、台の大きい観測量で CRM が標準より悪くなりうる([README_readout_noise.md](README_readout_noise.md)) |

---

## 10. 論文の展望と、この研究の次の方針

論文の結論は、今後の課題として次を挙げている。

- **重点サンプリング(importance sampling)や適応的な手法**との組み合わせ
- 変分量子アルゴリズムの勾配の推定や、**物質の相の探査**への応用
- 補助系を使った測定からの、より良い後処理

これらは、この研究の次の方針と直接つながる。

| 論文の展望 | この研究の方針 |
|---|---|
| 重点サンプリング・適応的な手法 | prior の情報で基底の選び方を偏らせる方法(locally biased classical shadows)と CRM の組み合わせ。特に、1基底1ショットでの天井 $1/(1-x^2)$ を超えられるか |
| 物質の相の探査 | Hubbard 模型のドープ系(ストライプ)で、平均場が壊れる領域での prior の選び方 |
| prior を反復で良くする(忠実度の例) | データの半分で prior を合わせ込み、残りで推定する(偏りを生まない分割) |
| (論文にない) | 読み出し誤差・状態の用意の誤差を prior 側に取り込む |

---

## 11. この文書の限界

- 論文の内容は arXiv 版(v1)の HTML と PRX Quantum 版の PDF の両方で確認した。式番号は両方を併記した(0章)。
- 数値例の細部(図の読み取りによる改善の倍率など)は要約にとどめ、論文の図の数値は引用していない。
- 6章の $n_m$ 依存の表と、9章最後の行の読み出し誤差の数値は、この研究の側で式(25) [PRX: (C7)] から計算したもので、論文には載っていない。
