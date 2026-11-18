# 元論文の詳細解説 — 共通ランダム測定(CRM)による量子状態の性質の推定の改善

[README.md](README.md) に戻る / この研究との対応だけを知りたいとき: [README_crm_paper.md](README_crm_paper.md) / 研究全体の見取り図: [README_overview.md](README_overview.md)

**元論文**: B. Vermersch, A. Rath, B. Sundar, C. Branciard, J. Preskill, A. Elben, "Enhanced Estimation of Quantum Properties with Common Randomized Measurements", *PRX Quantum* **5**, 010352 (2024). [doi:10.1103/PRXQuantum.5.010352](https://doi.org/10.1103/PRXQuantum.5.010352) / [arXiv:2304.12292](https://arxiv.org/abs/2304.12292)

この文書は、元論文の本文と付録の全体を、式の導出の途中を省かずに解説する。[README_crm_paper.md](README_crm_paper.md) が「この研究が論文のどこを使ったか」に絞っているのに対し、こちらは論文そのものを最初から最後まで読むための資料である。

- **読んだ範囲**: PRX Quantum 版の全14ページ(本文、付録 A〜E、参考文献)と、arXiv 版(v1)の本文と付録の全文。
- **式番号の書き方**: 他の文書と同じく「式(25) [PRX: (C7)]」のように、先に arXiv 版、角括弧の中に PRX 版の番号を書く。本文の式(1)〜(7)は両方の版で同じ番号である。節の番号は PRX 版に従う。
- **論文にない説明の区別**: 論文が結果だけを書いている箇所の途中計算、具体的な数値例、読むときの注意は、この文書で補ったものである。補った箇所には「(補足)」と書く。
- **前提知識**: 行列の計算(積、トレース、固有値)と、確率の平均・分散を仮定する。量子力学の記法は 0.5 節でまとめて説明する。もっと手前から読みたいときは、[README_basics.md](README_basics.md) の1〜4章が、同じ内容を式を減らして説明している。
- **読む範囲の目安**: この研究(Hubbard 模型への適用)が使っているのは、5章(パウリ観測量の分散)までである。6章以降は、純度やエントロピーのような「状態の多項式で書ける量」への拡張と応用例で、数学的に最も重いのは6章である。初めて読むときは、0〜5章と12章を読み、6〜10章は必要になってから読めばよい。

---

## 0. 論文の全体像

### 0.1 一言でいうと

量子コンピュータの状態 $\rho$ をランダムな基底で測ったデータがあるとする。 $\rho$ に近い状態 $\sigma$ を古典計算機が持っていれば、**「同じランダムな基底で $\sigma$ を測ったらどうなるか」を古典計算で求めて実験の結果から差し引く**ことで、推定を偏らせずに統計誤差だけを減らせる。これが共通ランダム測定(common randomized measurements、CRM)である。

### 0.2 論文の主な結果

| 結果 | 場所 | 内容 |
|---|---|---|
| CRM シャドウの定義 | 式(2)、(3) | $\hat\rho_\sigma^{(r)}=\hat\rho^{(r)}-\sigma^{(r)}+\sigma$ 。 $\sigma$ が何であっても不偏 |
| パウリ観測量の分散 | 式(4)、式(25) [PRX: (C7)] | 分散は $\mathrm{Tr}[O(\rho-\sigma)]^2$ に比例する項と、ショット雑音の項の和 |
| 多コピー観測量の分散の上界 | 式(6)、式(34) [PRX: (C16)] | prior の良さは Hilbert–Schmidt 距離 $\lVert\rho_A-\sigma_A\rVert_2$ で入る |
| 実験データでの実証(PRX 版のみ) | III.A 節 | イオントラップ10量子ビット。エンタングルメントの検出に必要な設定数が約6分の1 |
| エントロピーの推定(数値) | III.B 節 | 臨界イジング鎖。結合次元 2〜3 の MPS を prior にして誤差が大きく減る |
| 忠実度の推定(数値) | III.C 節 | 30量子ビットのランダム回路。prior を反復で選び直す手順 |
| 別の実験から prior を作る | 付録 D [PRX: 付録 A] | $\sigma$ の古典的な記述がなくても、 $\sigma$ を実現する実験のデータで代用できる |

### 0.3 論文の構成(PRX 版)と、この文書の章の対応

| PRX 版 | arXiv 版 | 内容 | この文書 |
|---|---|---|---|
| I 節 | 序論 | 動機、共通乱数法 | 1章、2章 |
| II.A 節 | 「Randomized measurements & classical shadows」 | 古典シャドウの復習 | 3章 |
| II.B 節 | 「Common randomized measurements」 | CRM シャドウ | 4章 |
| II.C 節 | 「Estimation of Pauli observables」 | パウリ観測量、式(4) | 5章 |
| II.D 節 | 「Estimation of MCOs with CRM shadows」 | 多コピー観測量、式(5)、(6) | 6章 |
| III.A 節 | なし | 実験データでのエンタングルメント検出 | 7章 |
| III.B 節 | 例1 | フォン・ノイマン・エントロピー | 8章 |
| III.C 節 | 例2 | 忠実度の推定 | 9章 |
| IV 節 | 結論 | 結論と展望 | 11章 |
| 付録 A | 付録 D | 別の実験から prior を作る | 10章 |
| 付録 B | 付録 A | 分散の一般式 | 5.2節、6.3節 |
| 付録 C | 付録 B | 局所パウリ測定での分散 | 5.3節、6.4節 |
| 付録 D | 付録 E | エントロピーの多項式近似の係数 | 8.2節 |
| 付録 E | 付録 C | シフト演算子の $O_A^{(1)}$ | 6.6節 |

### 0.4 記号

| 記号 | 意味 |
|---|---|
| $N$ | 量子ビットの数 |
| $\rho$ | 実験で用意された状態(密度行列)。知りたい対象 |
| $\sigma$ | $\rho$ の古典的に表現できる近似。prior と呼ぶ |
| $U=\bigotimes_iU_i$ | 測定の前に掛けるランダムなユニタリー。1つ選ぶことを「設定」と呼ぶ |
| $N_U$ | 設定の数 |
| $N_M$ | 1つの設定で繰り返す測定(ショット)の数 |
| $\mathbf s=(s_1,\dots,s_N)$ | 測定で得るビット列。 $s_i\in\lbrace0,1\rbrace$ |
| $P_\rho(\mathbf s\mid U)$ | 設定 $U$ でビット列 $\mathbf s$ が出る確率 $\langle\mathbf s\rvert U\rho U^\dagger\lvert\mathbf s\rangle$ |
| $\widehat P_\rho(\mathbf s\mid U)$ | 同じ確率を $N_M$ 回の測定の頻度で推定したもの(ハット付きは実験からの推定値) |
| $\mathcal M^{-1}$ | 測定チャネルの逆(3.2節) |
| $O$ | 観測量。 $n$ コピーに作用するものを多コピー観測量(MCO)と呼ぶ |
| $A$ 、 $N_A$ | 観測量の台(恒等演算子でない量子ビットの集合)と、その量子ビット数 |
| $\mathbb E_U$ 、 $\mathbb E_{QM}$ | ユニタリーの選び方についての平均と、設定を固定したときの測定結果についての平均 |
| $\mathbb V$ | 分散 |
| $\lVert X\rVert_2$ | Hilbert–Schmidt ノルム $\sqrt{\mathrm{Tr}(X^2)}$ ( $X$ はエルミート) |

### 0.5 読む前に — この文書で使う量子力学の記法(補足)

論文は量子情報の記法を説明なしに使う。ここで必要なものをまとめる。すでに知っている読者は飛ばしてよい。

**状態と期待値。**

- $\lvert\psi\rangle$ は縦ベクトル、 $\langle\psi\rvert$ はその複素共役を取って横に倒したもの、 $\langle\phi\vert\psi\rangle$ は内積である。 $\lvert\psi\rangle\langle\psi\rvert$ は「縦 × 横」で行列になる。
- **密度行列** $\rho$ は、状態を行列で表したものである。1つの状態ベクトルで書ける状態(**純粋状態**)なら $\rho=\lvert\psi\rangle\langle\psi\rvert$ 。いくつかの状態を確率 $p_k$ で混ぜたもの(**混合状態**)なら $\rho=\sum_kp_k\lvert\psi_k\rangle\langle\psi_k\rvert$ 。どちらも、固有値が 0 以上(半正定値)で、トレースが 1 のエルミート行列である。
- 観測量 $O$ (エルミート行列)の**期待値**は $\mathrm{Tr}(O\rho)$ である。 $\mathrm{Tr}$ はトレース(対角成分の和)で、 $\mathrm{Tr}(AB)=\mathrm{Tr}(BA)$ が成り立つ。
- **エルミート**とは、転置して複素共役を取ると元に戻る行列( $X^\dagger=X$ )のことである。固有値が実数になる。

**測定。**

- **計算基底**は、各量子ビットが 0 か 1 に決まった状態 $\lvert\mathbf s\rangle=\lvert s_1,\dots,s_N\rangle$ の集まりである。計算基底で測ると、ビット列 $\mathbf s$ が確率 $\langle\mathbf s\rvert\rho\lvert\mathbf s\rangle$ ( $\rho$ の対角成分)で出る。量子コンピュータが直接行えるのはこの測定である。
- **ユニタリー** $U$ は、 $U^\dagger U=\mathbb 1$ を満たす行列で、量子状態に対する回転に当たる。状態 $\rho$ に $U$ を掛けると $U\rho U^\dagger$ になる。測る前に $U$ を掛けることで、計算基底以外の基底で測ることができる。
- $\mathbb 1$ は恒等演算子(単位行列)、 $\mathbb 1_2$ は1量子ビットの $2\times2$ の単位行列である。

**複数の量子ビット。**

- $\otimes$ は**テンソル積**である。 $A\otimes B$ は「1番目の量子ビットに $A$ 、2番目に $B$ を作用させる」演算子を表す。 $\bigotimes_iU_i$ は、各量子ビットに別々の $U_i$ を作用させることを表す。 $\rho^{\otimes n}$ は、同じ状態 $\rho$ を $n$ 個並べたもの( $n$ コピー)である。
- **部分トレース**は、量子ビットの集合 $A$ の外( $\bar A$ と書く)を見なかったことにする操作で、 $\mathrm{Tr}_{\bar A}$ と書く。 $\rho$ にこの操作をした結果を $\rho_A$ と書き、**縮約密度行列**と呼ぶ。 $A$ の上だけに作用する観測量の期待値は、 $\rho_A$ だけで決まる。
- **パウリ文字列**は、各量子ビットに $\mathbb 1$ 、 $X$ 、 $Y$ 、 $Z$ のどれかを割り当てたテンソル積である。 $\mathbb 1$ 以外が割り当てられた量子ビットの集合が**台**である。

**その他の記号。**

- $\delta_{a,b}$ は、 $a=b$ のとき 1、そうでないとき 0 を取る記号(クロネッカーのデルタ)である。
- **ブロッホベクトル**は、1量子ビットの状態を $\rho=\frac12(\mathbb 1+r_XX+r_YY+r_ZZ)$ と書いたときの3つの実数 $(r_X,r_Y,r_Z)$ である。 $r_P$ はパウリ行列 $P$ の期待値に等しい。1量子ビットの状態は、この3つの数で完全に決まる。
- **チャネル**は、行列を行列に移す線形の写像のことである。この文書では「測定して結果を記録する」という操作全体を1つの写像と見るときに使う。
- **不偏**とは、推定値の平均が知りたい値に一致することである。

---

## 1. 動機(I 節)

### 1.1 ランダム測定で何を測るか

ランダム測定は、状態 $\rho$ にランダムなユニタリーを掛けてから計算基底で測ることを繰り返し、そのデータから $\rho$ の性質を推定する方法である。古典シャドウ(Huang–Kueng–Preskill, 2020)はその代表で、データから $\rho$ の「不偏な推定値」を作り、それを使って多数の観測量をまとめて推定する。

推定の対象は2種類ある。

- **線形の観測量** $\mathrm{Tr}(O\rho)$ 。少数の量子ビットに作用する観測量の期待値。変分量子固有値ソルバーで局所ハミルトニアンの期待値を測る場合がこれに当たる。
- **多コピー観測量(MCO)** $\mathrm{Tr}(O\rho^{\otimes n})$ ( $n\ge1$ )。 $\rho$ の $n$ 次の多項式で書ける量である。純度 $\mathrm{Tr}(\rho^2)$ 、Rényi エントロピー、部分転置のモーメント(混合状態のエンタングルメントの判定に使う)、量子 Fisher 情報量の下界などがこの形に書ける。

(補足)純度が2コピーの観測量であることは次のように確かめられる。2つのコピーを入れ替えるスワップ演算子 $\mathbb S\lvert i\rangle\lvert j\rangle=\lvert j\rangle\lvert i\rangle$ を使うと

```math
\mathrm{Tr}\bigl[\mathbb S(\rho\otimes\rho)\bigr]=\sum_{i,j}\langle i\rvert\langle j\rvert\,\mathbb S\,(\rho\otimes\rho)\lvert i\rangle\lvert j\rangle=\sum_{i,j}\langle j\rvert\rho\lvert i\rangle\langle i\rvert\rho\lvert j\rangle=\mathrm{Tr}(\rho^2)
```

となる。2番目の等号では、 $\mathbb S$ がエルミートなので左のブラに作用させて $\langle i\rvert\langle j\rvert\mathbb S=\langle j\rvert\langle i\rvert$ とした。

(補足)ここに出てくる量を一言ずつ説明する。

- **変分量子固有値ソルバー**は、量子コンピュータの上で試しの状態を作り、エネルギーの期待値を測り、それが下がるように状態を調節していく計算法である。エネルギーは局所的なパウリ文字列の和なので、その期待値を何度も測ることになる。
- **純度** $\mathrm{Tr}(\rho^2)$ は、状態がどれだけ「混ざっていないか」を表す。純粋状態で 1、混ざるほど小さくなる。
- **エントロピー**は、状態がどれだけ混ざっているかを表す量である。系の一部分 $A$ だけを見た状態 $\rho_A$ のエントロピーは、 $A$ と残りの部分の間の**エンタングルメント**(量子的な絡み合い)の強さを表す。全体が純粋状態でも、部分が残りと絡み合っていれば $\rho_A$ は混合状態になるからである。
- これらは $\rho$ の1次式ではなく、 $\rho$ を2回以上掛けた形をしている。「 $\rho$ のコピーを $n$ 個並べた $\rho^{\otimes n}$ の上の観測量」と見なせるので、多コピー観測量と呼ぶ。実際に $n$ 個の状態を同時に用意する必要はなく、1個ずつ測ったデータを後処理で組み合わせて推定する(6章)。

### 1.2 問題 — 必要な測定の数

ランダム測定の統計誤差を一定以下に抑えるのに必要な測定の数は、多コピー観測量では(部分)系の大きさについて指数的に増えるのが普通である。線形の観測量については誤差を減らす方法が多く提案されてきたが(基底の選び方を偏らせる方法、脱乱択化、観測量のグループ分けなど)、多コピー観測量に使える一般的な方法が求められていた。

論文の提案は、**$\rho$ についての近似的な知識を使って、線形・多コピーの両方で統計誤差を減らす枠組み**である。近似 $\sigma$ は、実験の理論的なモデル化(例えばテンソルネットワークによる近似的な古典シミュレーション)からでも、別の実験のデータからでも作れる。

---

## 2. 共通乱数法 — 論文の考え方の元(I 節)

CRM の元になっているのは、モンテカルロ法で使われる共通乱数法(common random numbers)である。

確率変数 $X$ の期待値 $\mathbb E[X]$ を、標本 $X_i$ の平均で推定したい。統計誤差の大きさは分散 $\mathbb V[X]$ で決まる。ここで、 $X$ と強く相関していて、**期待値 $\mathbb E[Y]$ が厳密にわかっている**別の確率変数 $Y$ が手元にあるとする。このとき $X$ の代わりに

```math
X'=X-Y+\mathbb E[Y]
```

を、同じ乱数から作った組 $(X_i,Y_i)$ について平均する。

- **偏らない**: $\mathbb E[X']=\mathbb E[X]-\mathbb E[Y]+\mathbb E[Y]=\mathbb E[X]$ 。
- **分散**: 定数 $\mathbb E[Y]$ は分散に効かないので

```math
\mathbb V[X']=\mathbb V[X-Y]=\mathbb V[X]+\mathbb V[Y]-2\,\mathrm{cov}(X,Y),\qquad\mathrm{cov}(X,Y)=\mathbb E[XY]-\mathbb E[X]\,\mathbb E[Y]
```

- **得をする条件**: $\mathbb V[X']<\mathbb V[X]$ は $\mathrm{cov}(X,Y)>\mathbb V[Y]/2$ と同じである(PRX 版の文献[41]の注)。 $Y$ が $X$ と十分に強く正に相関していれば分散が減り、 $Y=X$ なら分散は 0 になる。相関が弱いと $\mathbb V[Y]$ が加わる分だけ**かえって分散が増える**。

(補足)身近な例で言い直す。300 人の今年のテストの平均点を、くじで選んだ 10 人の採点だけから見積もるとする。 $X$ は選んだ 10 人の今年の平均点である。去年の点が全員分わかっていれば、 $Y$ として「同じ 10 人の去年の平均点」を取れる。 $\mathbb E[Y]$ は全員の去年の平均点で、厳密にわかっている。たまたま点の高い人が多く選ばれたという偶然は $X$ にも $Y$ にも同じように現れるので、 $X-Y$ では打ち消される。「同じ 10 人」を使うことが、「共通の乱数」に当たる。

(補足)得をする条件を数で見る。 $\mathbb V[X]=1$ 、 $\mathbb V[Y]=1$ とし、相関係数( $\mathrm{cov}$ を2つの標準偏差の積で割ったもの)を $c$ とすると、 $\mathbb V[X']=2-2c$ である。 $c=0.9$ なら 0.2(5分の1)、 $c=0.5$ なら 1(変わらない)、 $c=0$ なら 2(2倍に悪化)になる。

CRM では、 $X$ が「実験のデータから作った推定値」、 $Y$ が「同じランダムなユニタリーで $\sigma$ を測ったときの推定値(古典計算)」、 $\mathbb E[Y]$ が「 $\sigma$ での厳密な値」に当たる。 $X$ と $Y$ が共有する乱数は、ランダムなユニタリーの選び方である。

---

## 3. 古典シャドウの復習(II.A 節)

### 3.1 測定の手順とデータ

1. ランダムなユニタリー $U$ を選んで $\rho$ に掛ける。
2. 回転した状態 $U\rho U^\dagger$ を計算基底 $\lvert\mathbf s\rangle=\lvert s_1,\dots,s_N\rangle$ で測り、ビット列を得る。
3. 同じ $U$ で 2 を $N_M$ 回繰り返す。
4. $U$ を選び直して 1〜3 を $N_U$ 回繰り返す。

データは $N_UN_M$ 個のビット列 $\mathbf s^{(r,b)}$ ( $r=1,\dots,N_U$ は設定の番号、 $b=1,\dots,N_M$ はショットの番号)である。

論文の枠組みはどんなランダム測定にも使えるが、例ではパウリ測定を使う。 $U=\bigotimes_{i=1}^NU_i$ で、各 $U_i$ を次の3つから一様に選ぶ。

```math
U_i\in\mathcal U=\Bigl\{\ \mathbb 1_2,\ \ \tfrac1{\sqrt2}\begin{pmatrix}1&1\\1&-1\end{pmatrix},\ \ \tfrac1{\sqrt2}\begin{pmatrix}1&-i\\1&+i\end{pmatrix}\ \Bigr\},\qquad U_i^\dagger ZU_i=Z,\ X,\ Y
```

(補足)「 $U_i$ を掛けてから $Z$ を測る」ことは、「掛ける前の状態で $U_i^\dagger ZU_i$ を測る」ことと同じである。 $\langle\psi\rvert U_i^\dagger ZU_i\lvert\psi\rangle$ が、 $U_i\lvert\psi\rangle$ での $Z$ の期待値だからである。したがって上の3つは、それぞれ $Z$ 、 $X$ 、 $Y$ の基底で測ることに当たり、各量子ビットで測る基底を確率 $1/3$ ずつで選んでいることになる。

### 3.2 測定チャネルとその逆

この節の目的は、「測定の結果から、元の状態 $\rho$ の推定値をどう作るか」の公式を得ることである。考え方は次のとおり。ランダムな基底で測って結果を記録する、という操作を平均すると、元の状態 $\rho$ を決まった仕方でぼかしたもの $\mathcal M(\rho)$ になる。ぼかし方がわかれば、その逆の操作 $\mathcal M^{-1}$ を記録に施すことで、平均すれば $\rho$ に戻る推定値が作れる。

(補足)論文は逆チャネルの形 $\mathcal M^{-1}(\bigotimes_iO_i)=\bigotimes_i(3O_i-\mathrm{Tr}(O_i)\mathbb 1_2)$ を結果だけ引用している。ここでは1量子ビットで導く。

1量子ビットの状態を $\rho=\frac12(\mathbb 1+r_XX+r_YY+r_ZZ)$ と書く( $r_P=\mathrm{Tr}(P\rho)$ )。基底 $P\in\lbrace X,Y,Z\rbrace$ で測ると、結果 $\pm1$ が確率 $(1\pm r_P)/2$ で出て、測定後の状態は $P$ の固有状態 $\frac12(\mathbb 1\pm P)$ になる。「測定して、結果に応じた固有状態を記録する」操作を結果について平均すると

```math
\sum_{\pm}\frac{1\pm r_P}2\cdot\frac{\mathbb 1\pm P}2=\frac14\bigl[(1+r_P)(\mathbb 1+P)+(1-r_P)(\mathbb 1-P)\bigr]=\frac12\bigl(\mathbb 1+r_PP\bigr)
```

である。さらに基底の選び方( $X,Y,Z$ を確率 $1/3$ ずつ)について平均すると

```math
\mathcal M(\rho)=\frac13\sum_{P}\frac12\bigl(\mathbb 1+r_PP\bigr)=\frac12\Bigl(\mathbb 1+\frac{r_XX+r_YY+r_ZZ}3\Bigr)=\frac13\bigl(\rho+\mathbb 1\bigr)
```

となる。これが測定チャネル $\mathcal M$ で、 $\rho$ のブロッホベクトルを3分の1に縮める写像である。トレースが 1 でない行列にも線形に広げると $\mathcal M(X)=\frac13\bigl(X+\mathrm{Tr}(X)\mathbb 1\bigr)$ で、その逆は

```math
\mathcal M^{-1}(X)=3X-\mathrm{Tr}(X)\,\mathbb 1
```

である。確かめると、 $\mathcal M^{-1}(\mathcal M(X))=\bigl(X+\mathrm{Tr}(X)\mathbb 1\bigr)-\frac13\bigl(\mathrm{Tr}(X)+2\mathrm{Tr}(X)\bigr)\mathbb 1=X$ となる。 $N$ 量子ビットでは各量子ビットの基底を独立に選ぶので、積の形の演算子に対して

```math
\mathcal M^{-1}\Bigl(\bigotimes_iO_i\Bigr)=\bigotimes_i\bigl(3O_i-\mathrm{Tr}(O_i)\,\mathbb 1_2\bigr)\qquad(\text{式(19) [PRX: (C1)]})
```

で、積の形でない演算子には線形性で広げる。

この逆チャネルには、後で使う性質が2つある。

- **エルミート性を保つ**: $X$ がエルミートなら $\mathcal M^{-1}(X)$ もエルミートである。
- **自己共役である**: $\mathrm{Tr}[A\mathcal M^{-1}(B)]=\mathrm{Tr}[\mathcal M^{-1}(A)B]$ 。1量子ビットでは両辺とも $3\mathrm{Tr}(AB)-\mathrm{Tr}(A)\mathrm{Tr}(B)$ になる。

### 3.3 古典シャドウの定義と不偏性

$r$ 番目の設定のデータから作る古典シャドウは

```math
\hat\rho^{(r)}=\sum_{\mathbf s}\widehat P_\rho(\mathbf s\mid U^{(r)})\ \mathcal M^{-1}\bigl(U^{(r)\dagger}\lvert\mathbf s\rangle\langle\mathbf s\rvert U^{(r)}\bigr)\qquad(\text{式(1)})
```

```math
\widehat P_\rho(\mathbf s\mid U^{(r)})=\frac1{N_M}\sum_{b=1}^{N_M}\delta_{\mathbf s,\mathbf s^{(r,b)}}
```

である。 $\widehat P_\rho$ は、設定 $U^{(r)}$ の $N_M$ 回の測定でビット列 $\mathbf s$ が出た割合である。 $U^{(r)\dagger}\lvert\mathbf s\rangle\langle\mathbf s\rvert U^{(r)}$ は、「結果 $\mathbf s$ が出た」ことに対応する、回転する前の状態である。

**不偏性**: 設定を固定して測定結果について平均すると、頻度は確率になる( $\mathbb E_{QM}[\widehat P_\rho]=P_\rho$ )。続けて設定について平均すると、3.2節の測定チャネルの定義そのものが現れる。

```math
\mathbb E[\hat\rho^{(r)}]=\mathbb E_U\mathbb E_{QM}[\hat\rho^{(r)}]=\mathcal M^{-1}\Bigl(\underbrace{\mathbb E_U\sum_{\mathbf s}P_\rho(\mathbf s\mid U)\,U^\dagger\lvert\mathbf s\rangle\langle\mathbf s\rvert U}_{=\mathcal M(\rho)}\Bigr)=\rho
```

1つの $\hat\rho^{(r)}$ は状態ではない(負の固有値を持つ)が、平均すれば $\rho$ になる。

(補足)1量子ビット、1ショットの例。 $Z$ 基底を引いて結果が $s=0$ だったとする。記録される状態は $\lvert0\rangle\langle0\rvert$ で、対角成分が $(1,0)$ の行列である。3.2節の公式 $\mathcal M^{-1}(X)=3X-\mathrm{Tr}(X)\mathbb 1$ を施すと、対角成分が $(3-1,0-1)=(2,-1)$ の行列になる。固有値に $-1$ があるので、これは量子状態ではない。それでも、基底の選び方と測定結果について平均すれば、正しい $\rho$ になる。古典シャドウは、1個ずつは荒れているが平均すれば正しい「状態の推定値」である。観測量の推定値は $\mathrm{Tr}(O\hat\rho^{(r)})$ を $r$ について平均して作る。

### 3.4 パウリ文字列の推定値(補足)

ここが、この研究で使う「当たり・外れ」の描像の出どころである。結論を先に言うと、パウリ文字列の推定値は「選んだ基底が台の上で全部当たれば、測定値の積を $3^{N_A}$ 倍したもの。1つでも外れれば 0」になる。以下はその導出である。

後の章のために、パウリ文字列 $\gamma=\bigotimes_i\gamma_i$ ( $\gamma_i\in\lbrace\mathbb 1,X,Y,Z\rbrace$ 、台 $A$ )の推定値を具体的に書く。

- $\mathcal M^{-1}$ が自己共役なので、 $\mathrm{Tr}[\gamma\mathcal M^{-1}(X)]=\mathrm{Tr}[\mathcal M^{-1}(\gamma)X]$ である。
- パウリ行列はトレースが 0 なので $\mathcal M_i^{-1}(\gamma_i)=3\gamma_i$ 、恒等演算子は $\mathcal M_i^{-1}(\mathbb 1)=3\cdot\mathbb 1-2\cdot\mathbb 1=\mathbb 1$ である。したがって $\mathcal M^{-1}(\gamma)=3^{N_A}\gamma$ 。
- 残るのは $\gamma(U,\mathbf s)\equiv\mathrm{Tr}[\gamma U^\dagger\lvert\mathbf s\rangle\langle\mathbf s\rvert U]=\prod_i\langle s_i\rvert U_i\gamma_iU_i^\dagger\lvert s_i\rangle$ である。台の上の量子ビット $i$ について、 $U_i$ が $\gamma_i$ の基底に当たっていれば( $U_i^\dagger ZU_i=\gamma_i$ )、 $U_i\gamma_iU_i^\dagger=Z$ で因子は $(-1)^{s_i}$ になる。外れていれば $U_i\gamma_iU_i^\dagger$ は $Z$ 以外のパウリ行列(符号を除く)で、対角成分が 0 なので因子は 0 になる。

まとめると、台の上の基底を全部当てるユニタリーを $V_\gamma$ と書いて

```math
\gamma(U,\mathbf s)=\delta_{U_A,V_\gamma}\,(-1)^{\sum_{i\in A}s_i}\qquad(\text{式(20) [PRX: (C2)]})
```

である( $U_A=\bigotimes_{i\in A}U_i$ )。つまり1つの設定での推定値は

```math
\mathrm{Tr}\bigl(\gamma\,\hat\rho^{(r)}\bigr)=3^{N_A}\,\delta_{U_A^{(r)},V_\gamma}\ \overline m^{(r)},\qquad\overline m^{(r)}=\frac1{N_M}\sum_b(-1)^{\sum_{i\in A}s_i^{(r,b)}}
```

で、**基底が台の上で全部当たったとき(確率 $3^{-N_A}$ )だけ、測定値の積の平均 $\overline m$ を $3^{N_A}$ 倍した値になり、外れたときは 0 になる**。当たる確率と倍率が打ち消し合うので平均は $\mathrm{Tr}(\gamma\rho)$ になるが、「大きな値か 0 か」のくじなので分散が大きい。CRM が減らすのは、このくじによる分散である。

---

## 4. CRM シャドウ(II.B 節)

### 4.1 定義

$\rho$ の近似 $\sigma$ が古典計算機の上にあるとする。論文は $\sigma$ がエルミートであることだけを仮定し、半正定値であることもトレースが 1 であることも要求しない。つまり $\sigma$ は本物の量子状態でなくてもよく、「各観測量の値を何と予言するか」を決める行列でありさえすればよい。そのため $\sigma$ を「擬状態(pseudostate)」と呼んでいる。CRM シャドウは

```math
\hat\rho_\sigma^{(r)}=\hat\rho^{(r)}-\sigma^{(r)}+\sigma\qquad(\text{式(2)})
```

```math
\sigma^{(r)}=\sum_{\mathbf s}P_\sigma(\mathbf s\mid U^{(r)})\ \mathcal M^{-1}\bigl(U^{(r)\dagger}\lvert\mathbf s\rangle\langle\mathbf s\rvert U^{(r)}\bigr),\qquad P_\sigma(\mathbf s\mid U^{(r)})=\langle\mathbf s\rvert U^{(r)}\sigma U^{(r)\dagger}\lvert\mathbf s\rangle\qquad(\text{式(3)})
```

である。3つの項の意味は次のとおり。

| 項 | 何か | 何から作るか |
|---|---|---|
| $\hat\rho^{(r)}$ | 実験の古典シャドウ(式(1)) | 実験のデータ |
| $\sigma^{(r)}$ | **実験と同じユニタリー $U^{(r)}$ で $\sigma$ を測ったとしたときのシャドウ** | 古典計算。確率 $P_\sigma$ は頻度ではなく**厳密な値**(ハットが付いていない) |
| $\sigma$ | prior そのもの | 古典計算 |

2章の共通乱数法と見比べると、 $X=\hat\rho^{(r)}$ 、 $Y=\sigma^{(r)}$ 、 $\mathbb E[Y]=\sigma$ である。

### 4.2 不偏性

$\sigma^{(r)}$ は測定結果を含まないので、平均はユニタリーについてだけ取ればよい。3.3節と同じ計算で

```math
\mathbb E[\sigma^{(r)}]=\mathbb E_U[\sigma^{(r)}]=\mathcal M^{-1}\bigl(\mathcal M(\sigma)\bigr)=\sigma
```

である。(補足)この計算は $\mathcal M$ と $\mathcal M^{-1}$ が線形であることしか使っていないので、 $\sigma$ が状態でなくても成り立つ。したがって

```math
\mathbb E\bigl[\hat\rho_\sigma^{(r)}\bigr]=\rho-\sigma+\sigma=\rho
```

で、**$\sigma$ をどう選んでも CRM シャドウは $\rho$ の不偏な推定値である**。prior が悪いときに起こるのは、答えがずれることではなく、分散が増えることである。

### 4.3 データの取り方は変わらない

論文が強調している点である。

- 実験で行うのは普通のランダム測定で、 $\sigma$ に依存しない。 $\sigma$ は**後処理の段階でだけ**使う。
- したがって、**実験が終わった後で $\sigma$ を選べる**。より正確な理論モデルが後から手に入れば、同じデータに使える。複数の $\sigma$ を同じデータで試すこともできる(9章の反復手順)。
- これは、測る基底の選び方を prior に合わせて変える方法(重点サンプリング)との大きな違いである。

### 4.4 なぜ分散が減るか(直観)

ショット数が十分に多い( $N_M\gg1$ )とすると、頻度は確率に近づき、 $\hat\rho^{(r)}$ は「厳密な確率 $P_\rho$ で作ったシャドウ」 $\rho^{(r)}$ に近づく。 $\rho\approx\sigma$ なら、 $\rho^{(r)}$ と $\sigma^{(r)}$ は同じユニタリーから同じ式で作られているので、ほぼ同じように揺らぐ。差 $\hat\rho^{(r)}-\sigma^{(r)}$ ではその揺らぎが打ち消し合い、行列要素の分散は $\hat\rho^{(r)}$ だけのときより小さくなる。5章と6章は、この直観を分散の式として厳密にしたものである。

### 4.5 後処理のコスト

$\sigma^{(r)}$ を作る分だけ、後処理は標準のシャドウより重くなる。論文は、テンソルネットワークのような表現を使えば時間・メモリの両面で効率的に実行できると述べている。また、 $\sigma$ の理論的な記述がない場合に、 $\sigma$ を実現する別の実験のデータで代用する方法を付録に示している(10章)。

---

## 5. パウリ観測量の推定(II.C 節、付録 B.2 と C.1 [arXiv: A.2 と B.1])

この章が、この研究で使う利得の式の出どころである。流れを先に述べる。

1. 推定値のばらつきを、「どの設定を引いたか」による部分と「設定を固定した中でのショットの揺らぎ」による部分に分ける(5.2節)。
2. パウリ文字列の場合に、それぞれを具体的に計算する(5.3節)。前者は prior の外れの2乗に比例し、後者は prior によらない。
3. 結果を読む(5.4節)。

[README_basics.md](README_basics.md) の 3.4 節と 4.3 節に、同じ結果を1つのパウリ文字列に限って短く導いたものがある。先にそちらを読むと、以下の一般的な記号の意味がつかみやすい。

### 5.1 推定量と結果

パウリ文字列 $O=\gamma$ (台 $A$ 、 $N_A$ 量子ビット)の期待値 $\mathrm{Tr}(O\rho)$ を

```math
\hat O=\frac1{N_U}\sum_{r=1}^{N_U}\mathrm{Tr}\bigl(O\,\hat\rho_\sigma^{(r)}\bigr)
```

で推定する。論文の結果は

```math
\mathbb V[\hat O]=\frac1{N_U}\Bigl((3^{N_A}-1)\,\mathrm{Tr}\bigl[O(\rho-\sigma)\bigr]^2+\frac{3^{N_A}\bigl(1-\mathrm{Tr}(\rho O)^2\bigr)}{N_M}\Bigr)\le\frac{3^{N_A}}{N_U}\Bigl(\mathrm{Tr}\bigl[O(\rho-\sigma)\bigr]^2+\frac1{N_M}\Bigr)
```

で、左の等式が付録の式(25) [PRX: (C7)]、右の上界が本文の式(4)である。以下、これを順に導く。

(補足)この研究の記号との対応は、 $N_A=\lvert A\rvert$ 、 $\mathrm{Tr}(\rho O)=x$ 、 $\mathrm{Tr}[O(\rho-\sigma)]=\Delta$ 、 $N_U=n_u$ 、 $N_M=n_m$ である。括弧の中の第1項が「くじの揺らぎ」 $(3^{\lvert A\rvert}-1)\Delta^2$ 、第2項が「ショットの揺らぎ」 $v_s$ に当たる。

### 5.2 分散を2つに分ける(付録 B.2 [arXiv: A.2])

設定の番号 $(r)$ を省く。観測量を一般のエルミート演算子 $O^{(1)}$ とし(単一コピーでは $O^{(1)}=O$ 。多コピーでの意味は6章)、短く

```math
w(U,\mathbf s)\equiv\bigl[\mathcal M^{-1}(O^{(1)})\bigr](U,\mathbf s)=\mathrm{Tr}\bigl[\mathcal M^{-1}(O^{(1)})\,U^\dagger\lvert\mathbf s\rangle\langle\mathbf s\rvert U\bigr]
```

と書く。これは「設定 $U$ で結果 $\mathbf s$ が出たときに、1ショットが推定値に寄与する値」である。 $\mathcal M^{-1}$ が自己共役であることを使うと、式(2)、(3)から1つの設定での推定値は

```math
\mathrm{Tr}\bigl(O^{(1)}\hat\rho_\sigma\bigr)=\sum_{\mathbf s}\bigl(\widehat P_\rho(\mathbf s\mid U)-P_\sigma(\mathbf s\mid U)\bigr)\,w(U,\mathbf s)+\mathrm{Tr}(O^{(1)}\sigma)\qquad(\text{式(12) [PRX: (B5)]})
```

となる。最後の項は定数なので分散には効かない。この推定値には、乱数の出どころが2つある。**設定 $U$ のくじ**と、**設定を固定したときの測定結果の揺らぎ**(頻度 $\widehat P_\rho$ に入る)である。

**測定結果についての平均。** 頻度 $\widehat P_\rho(\mathbf s\mid U)$ は、独立な $N_M$ 回の測定で $\mathbf s$ が出た割合である。2つの頻度の積の平均は

```math
\mathbb E_{QM}\bigl[\widehat P_\rho(\mathbf s\mid U)\,\widehat P_\rho(\mathbf s'\mid U)\bigr]=P_\rho(\mathbf s\mid U)P_\rho(\mathbf s'\mid U)+\frac{\delta_{\mathbf s,\mathbf s'}P_\rho(\mathbf s\mid U)-P_\rho(\mathbf s\mid U)P_\rho(\mathbf s'\mid U)}{N_M}\qquad(\text{式(14) [PRX: (B7)]})
```

である。(補足)導出: $\widehat P_\rho(\mathbf s)\widehat P_\rho(\mathbf s')=N_M^{-2}\sum_{b,b'}\delta_{\mathbf s,\mathbf s^{(b)}}\delta_{\mathbf s',\mathbf s^{(b')}}$ を、 $b\ne b'$ の項と $b=b'$ の項に分ける。 $b\ne b'$ の $N_M(N_M-1)$ 個の項は独立なので平均が $P_\rho(\mathbf s)P_\rho(\mathbf s')$ 、 $b=b'$ の $N_M$ 個の項は同じショットなので平均が $\delta_{\mathbf s,\mathbf s'}P_\rho(\mathbf s)$ である。合わせると $(1-1/N_M)P_\rho(\mathbf s)P_\rho(\mathbf s')+\delta_{\mathbf s,\mathbf s'}P_\rho(\mathbf s)/N_M$ で、式(14)になる。

**結果。** 式(14)を分散の定義に代入して整理すると、1つの設定の分散は

```math
\mathbb V_1\equiv\mathbb V\bigl[\mathrm{Tr}(O^{(1)}\hat\rho_\sigma)\bigr]=\mathbb V_U\bigl[f_{\rho,\sigma}(U)\bigr]+\frac{\mathbb E_U\bigl[g_\rho(U)\bigr]}{N_M}\qquad(\text{式(15) [PRX: (B8)]})
```

```math
f_{\rho,\sigma}(U)=\sum_{\mathbf s}\bigl(P_\rho(\mathbf s\mid U)-P_\sigma(\mathbf s\mid U)\bigr)\,w(U,\mathbf s)\qquad(\text{式(16) [PRX: (B9)]})
```

```math
g_\rho(U)=\sum_{\mathbf s}P_\rho(\mathbf s\mid U)\,w(U,\mathbf s)^2-\Bigl(\sum_{\mathbf s}P_\rho(\mathbf s\mid U)\,w(U,\mathbf s)\Bigr)^2\qquad(\text{式(17) [PRX: (B10)]})
```

となる。

(補足)同じ結果は、全分散の法則「分散 = (条件付き平均の分散) + (条件付き分散の平均)」からすぐに読み取れる。

- 設定 $U$ を固定したときの推定値の平均は $f_{\rho,\sigma}(U)+\mathrm{Tr}(O^{(1)}\sigma)$ である(頻度を確率に置き換えたもの)。その $U$ についての分散が第1項 $\mathbb V_U[f]$ 。
- 設定 $U$ を固定したとき、推定値は「値 $w(U,\mathbf s)$ を確率 $P_\rho(\mathbf s\mid U)$ で取る確率変数」の $N_M$ 回の平均である。1回分の分散が $g_\rho(U)$ なので、 $N_M$ 回の平均の分散は $g_\rho(U)/N_M$ 。その $U$ についての平均が第2項。

2つの項の性質は対照的である。

| | 第1項 $\mathbb V_U[f_{\rho,\sigma}]$ | 第2項 $\mathbb E_U[g_\rho]/N_M$ |
|---|---|---|
| 出どころ | 設定のくじ(どのユニタリーを引いたか) | 1つの設定の中のショット雑音 |
| $N_M\to\infty$ で | 残る | 消える |
| prior への依存 | $\rho$ と $\sigma$ の出力分布の差で決まる。 $\sigma=\rho$ なら $f=0$ で消える | $\rho$ だけで決まる。prior では減らせない |

$N_U$ 個の設定は独立なので、単一コピーの観測量では $\mathbb V[\hat O]=\mathbb V_1/N_U$ である。

### 5.3 パウリ文字列での計算(付録 C.1 [arXiv: B.1])

$O^{(1)}=O=\gamma$ とする。3.4節より $\mathcal M^{-1}(\gamma)=3^{N_A}\gamma$ で、 $w(U,\mathbf s)=3^{N_A}\gamma(U,\mathbf s)=3^{N_A}\delta_{U_A,V_\gamma}(-1)^{\sum_{i\in A}s_i}$ である。

**第1項。** 基底が外れた設定では $f=0$ 。当たった設定では、 $\sum_{\mathbf s}P_\rho(\mathbf s\mid U)(-1)^{\sum_{i\in A}s_i}$ は「 $\gamma$ の基底で測った結果の積の期待値」なので $\mathrm{Tr}(\rho\gamma)$ に等しい。 $\sigma$ についても同じなので

```math
f_{\rho,\sigma}(U)=3^{N_A}\,\delta_{U_A,V_\gamma}\,\mathrm{Tr}\bigl((\rho-\sigma)\gamma\bigr)\qquad(\text{式(21) [PRX: (C3)]})
```

である。 $\delta_{U_A,V_\gamma}$ は確率 $3^{-N_A}$ で 1、それ以外で 0 を取る。 $\delta^2=\delta$ なので、その分散は

```math
\mathbb V_U[\delta_{U_A,V_\gamma}]=\mathbb E_U[\delta^2]-\mathbb E_U[\delta]^2=\frac1{3^{N_A}}-\frac1{9^{N_A}}
```

で、

```math
\mathbb V_U[f_{\rho,\sigma}(U)]=9^{N_A}\,\mathrm{Tr}\bigl((\rho-\sigma)\gamma\bigr)^2\Bigl(\frac1{3^{N_A}}-\frac1{9^{N_A}}\Bigr)=(3^{N_A}-1)\,\mathrm{Tr}\bigl((\rho-\sigma)\gamma\bigr)^2\qquad(\text{式(22) [PRX: (C4)]})
```

となる。

**第2項。** $\gamma(U,\mathbf s)^2=\delta_{U_A,V_\gamma}$ ( $\pm1$ の2乗は 1)なので、 $\sum_{\mathbf s}P_\rho w^2=9^{N_A}\delta_{U_A,V_\gamma}$ 、 $\bigl(\sum_{\mathbf s}P_\rho w\bigr)^2=9^{N_A}\delta_{U_A,V_\gamma}\mathrm{Tr}(\rho\gamma)^2$ である。したがって

```math
g_\rho(U)=9^{N_A}\,\delta_{U_A,V_\gamma}\bigl(1-\mathrm{Tr}(\rho\gamma)^2\bigr),\qquad\mathbb E_U[g_\rho(U)]=3^{N_A}\bigl(1-\mathrm{Tr}(\rho\gamma)^2\bigr)\qquad(\text{式(23)、(24) [PRX: (C5)、(C6)]})
```

となる。 $1-\mathrm{Tr}(\rho\gamma)^2$ は、平均 $\mathrm{Tr}(\rho\gamma)$ の $\pm1$ の確率変数の分散である。

**まとめ。** 2つを $\mathbb V[\hat O]=\mathbb V_1/N_U$ に入れると、5.1節の等式(式(25) [PRX: (C7)])になる。上界(式(4))は、 $3^{N_A}-1\le3^{N_A}$ と $1-\mathrm{Tr}(\rho\gamma)^2\le1$ を使って得られる。

### 5.4 結果の読み方

論文が式(4)の後で述べていることを整理する。

1. **標準のシャドウは $\sigma=0$ の場合である。** 式(2)、(3)で $\sigma=0$ とすると $\hat\rho_\sigma^{(r)}=\hat\rho^{(r)}$ になるので、標準のシャドウの分散は同じ式で $\mathrm{Tr}[O(\rho-\sigma)]$ を $\mathrm{Tr}(O\rho)$ に置き換えたものである。この結果は、多ショットのシャドウ推定についての既存の結果(Zhou–Liu, 2022 の定理2)と整合する。
2. **誤差の出どころは2つある。** 設定の数 $N_U$ が有限であること(第1項)と、1設定あたりのショット数 $N_M$ が有限であること(第2項)である。**CRM が減らすのは前者だけ**である。
3. **得をする条件。** $N_M$ の値によらず、 $\lvert\mathrm{Tr}[O(\rho-\sigma)]\rvert\le\lvert\mathrm{Tr}(O\rho)\rvert$ なら CRM の分散は標準のシャドウ以下である。(補足)等式(式(25) [PRX: (C7)])で比べると、第2項は共通で第1項だけが違うので、この条件は必要十分である。prior の予言 $\mathrm{Tr}(O\sigma)$ が、真の値と同じ符号で、真の値の2倍を超えなければよい。
4. **どんな実験で役に立つか。** 設定を切り替えるたびに較正の時間がかかる実験(イオントラップ、超伝導量子ビット)では、 $N_U$ は限られるが $N_M$ は大きく取れる。そこでは第1項が誤差の主な原因になるので、CRM の効果が大きい。

(補足)具体的な数値で見る。 $N_A=2$ 、 $\mathrm{Tr}(O\rho)=0.9$ 、 $\mathrm{Tr}(O\sigma)=0.8$ 、 $N_M=100$ とする。1設定あたりの分散は

| | 第1項(設定のくじ) | 第2項(ショット雑音) | 合計 |
|---|---|---|---|
| 標準のシャドウ | $8\times0.9^2=6.48$ | $9\times0.19/100=0.0171$ | 6.497 |
| CRM | $8\times0.1^2=0.08$ | 0.0171 | 0.0971 |

で、分散は約67分の1になる。同じ精度に必要な設定の数 $N_U$ が約67分の1になるということである。prior が完全( $\mathrm{Tr}(O\sigma)=0.9$ )でも第2項の 0.0171 は残るので、改善は $6.497/0.0171\approx380$ 倍で頭打ちになる。この頭打ちは $N_M$ を増やすと上がる。

---

## 6. 多コピー観測量の推定(II.D 節、付録 B と C.2 [arXiv: A と B.2])

(補足)この章は、純度やエントロピーのような、 $\rho$ を2回以上掛けた形の量を扱う。式が長いが、結論は次の2点である。

- 多コピー観測量の誤差の主な部分は、ある単一コピーの観測量 $O^{(1)}$ の推定の誤差に帰着する(6.3節)。したがって5章の「くじの揺らぎとショットの揺らぎ」の分解がそのまま使える。
- prior の良さは、台に縮約した状態どうしの Hilbert–Schmidt 距離 $\lVert\rho_A-\sigma_A\rVert_2$ として入る(6.4節)。これは、差の行列 $\rho_A-\sigma_A$ の全成分の2乗の和の平方根で、2つの行列が成分ごとにどれだけ違うかを表す。

Hubbard 模型へのこの研究の適用は単一コピーのパウリ文字列だけを使うので、この章は飛ばしても先を読める。

### 6.1 U 統計量 — なぜ異なる設定のシャドウを組み合わせるか

$n$ コピーの観測量 $\mathrm{Tr}(O\rho^{\otimes n})$ を推定するには、 $\rho^{\otimes n}$ の不偏な推定値が要る。独立な(異なる設定の)シャドウ $n$ 個のテンソル積は

```math
\mathbb E\bigl[\hat\rho_\sigma^{(r_1)}\otimes\cdots\otimes\hat\rho_\sigma^{(r_n)}\bigr]=\rho^{\otimes n}\qquad(r_1,\dots,r_n\ \text{がすべて異なるとき})
```

を満たす。独立な確率変数の積の平均は、平均の積だからである。同じ設定のシャドウを2回使うと( $r_1=r_2$ )、 $\mathbb E[\hat\rho\otimes\hat\rho]\ne\rho\otimes\rho$ となって偏る。そこで、**すべて異なる添字の組について平均する**。この形の推定量を U 統計量と呼ぶ。

(補足)数の例で確かめる。平均 $\mu$ 、分散 $v$ の独立な測定値 $x_1,x_2,\dots$ から $\mu^2$ を推定したいとする。同じ測定値を2回使った $x_1^2$ の平均は $\mu^2+v$ で、分散の分だけ大きく偏る。別々の測定値の積 $x_1x_2$ の平均は $\mu\times\mu=\mu^2$ で、偏らない。純度 $\mathrm{Tr}(\rho^2)$ のような $\rho$ の2次の量を推定するときに、異なる設定のシャドウどうしを掛けるのは同じ理由である。

### 6.2 バッチシャドウ(付録 B.1 [arXiv: A.1])

異なる添字の組は約 $N_U^n$ 通りあり、全部を計算するのは重い。論文は「バッチシャドウ」(Rath ら, 2023)を使う。 $N_U$ 個の CRM シャドウを $m$ 個のグループ( $m\ge n$ 、 $m$ は $N_U$ を割り切るとする)に分け、グループごとに平均する。

```math
\hat\rho_\sigma^{[t]}=\frac m{N_U}\sum_{r=1+(t-1)N_U/m}^{tN_U/m}\hat\rho_\sigma^{(r)},\qquad t=1,\dots,m\qquad(\text{式(8) [PRX: (B1)]})
```

(補足)論文の式(8) [PRX: (B1)]は前の因子が $N_U/m$ と印刷されているが、本文の説明「 $N_U/m$ 個のシャドウを平均する」のとおり、平均の因子 $m/N_U$ と読む。

推定量は、バッチシャドウの U 統計量である。

```math
\hat O=\frac{(m-n)!}{m!}\sum_{t_1\ne\cdots\ne t_n}\mathrm{Tr}\Bigl[O\bigl(\hat\rho_\sigma^{[t_1]}\otimes\cdots\otimes\hat\rho_\sigma^{[t_n]}\bigr)\Bigr]\qquad(\text{式(5)})
```

- 和は、 $1$ から $m$ までの互いに異なる $n$ 個の添字の並びすべてについて取る。並びの数は $m!/(m-n)!$ で、前の因子はそれで割って平均にしている。
- 異なるバッチは別々のシャドウから作るので独立で、それぞれ $\mathbb E[\hat\rho_\sigma^{[t]}]=\rho$ を満たす。したがって $\mathbb E[\hat O]=\mathrm{Tr}(O\rho^{\otimes n})$ で、**不偏**である。
- 計算する項の数は $\mathcal O(m^n)$ で、 $m$ を小さく取れば後処理が軽い。

3つの特別な場合がある。

| 場合 | 何になるか |
|---|---|
| $\sigma=0$ | 標準のバッチシャドウによる推定 |
| $\sigma=0$ かつ $m=N_U$ | 元々の古典シャドウの多コピー推定量(Huang ら 2020、Elben ら 2020) |
| $n=1$ | バッチに分ける意味がなくなり、5章の推定量と同じになる( $m$ に依存しない) |

### 6.3 分散の展開(付録 B.1 [arXiv: A.1])

バッチシャドウの推定量の分散は Rath ら(2023)が標準のシャドウについて計算している。その導出は $\mathbb E[\hat\rho^{(r)}]=\rho$ しか使っていないので、CRM シャドウにもそのまま当てはまる。

```math
\mathbb V[\hat O]=\sum_{\ell=1}^n\frac{\binom n\ell^2}{\binom m\ell}\Bigl(\frac m{N_U}\Bigr)^\ell\Bigl[\sum_{k=1}^\ell\binom\ell k(-1)^{\ell-k}\,\mathbb V_k\Bigr]=\frac{n^2}{N_U}\mathbb V_1+\frac{n^2(n-1)^2\frac m{m-1}}{2N_U^2}\bigl(\mathbb V_2-2\mathbb V_1\bigr)+\mathcal O\Bigl(\frac1{N_U^3}\Bigr)\qquad(\text{式(10) [PRX: (B3)]})
```

```math
\mathbb V_k=\mathbb V\Bigl[\mathrm{Tr}\Bigl(O^{(k)}\bigotimes_{r=1}^k\hat\rho_\sigma^{(r)}\Bigr)\Bigr],\qquad O^{(k)}=\mathrm{Tr}_{k+1,\dots,n}\Bigl[O_{\rm sym}\bigl(\mathbb 1^{\otimes k}\otimes\rho^{\otimes(n-k)}\bigr)\Bigr]\qquad(\text{式(11) [PRX: (B4)]})
```

記号の意味:

- $O_{\rm sym}=\frac1{n!}\sum_{\pi\in\mathcal S_n}W_\pi^\dagger OW_\pi$ は、 $O$ をコピーの入れ替えについて対称にしたもの( $W_\pi$ はコピーを置換 $\pi$ で入れ替える演算子、 $\mathcal S_n$ は $n$ 次の対称群)。 $\rho^{\otimes n}$ はコピーの入れ替えで変わらないので、 $O$ を $O_{\rm sym}$ に置き換えても期待値は同じである。
- $O^{(k)}$ は、 $n$ 個のコピーのうち $n-k$ 個を真の状態 $\rho$ で埋めて部分トレースを取り、残りの $k$ 個のコピーに作用する演算子にしたもの。 $\rho$ に依存する。
- $\mathbb V_k$ は、 $k$ 個のコピーだけをシャドウにし、残りを $\rho$ にしたときの分散。

(補足)右辺の展開を確かめる。 $\ell=1$ の項は $\frac{n^2}m\cdot\frac m{N_U}\cdot\mathbb V_1=\frac{n^2}{N_U}\mathbb V_1$ 。 $\ell=2$ の項は $\binom n2^2=\frac{n^2(n-1)^2}4$ 、 $\binom m2=\frac{m(m-1)}2$ 、角括弧の中が $\mathbb V_2-2\mathbb V_1$ なので、 $\frac{n^2(n-1)^2}{4}\cdot\frac2{m(m-1)}\cdot\frac{m^2}{N_U^2}(\mathbb V_2-2\mathbb V_1)$ で、式の第2項に一致する。

(補足)**主要項 $n^2\mathbb V_1/N_U$ の意味。** バッチシャドウを $\hat\rho_\sigma^{[t]}=\rho+\delta_t$ と書く( $\delta_t$ は平均 0 の揺らぎ)。テンソル積を $\delta$ の1次まで展開すると、 $n$ 個のコピーのうち1つだけが $\delta$ で残りが $\rho$ の項が $n$ 個出る。 $O$ を対称にしてあるので、どのコピーに $\delta$ があっても寄与は同じ形 $\mathrm{Tr}(O^{(1)}\delta)$ になる。すべてのバッチについて平均すると、推定量の揺らぎは1次で

```math
\hat O-\mathrm{Tr}(O\rho^{\otimes n})\approx n\,\mathrm{Tr}\Bigl(O^{(1)}\Bigl[\frac1{N_U}\sum_{r=1}^{N_U}\hat\rho_\sigma^{(r)}-\rho\Bigr]\Bigr)
```

となる。右辺は「単一コピーの観測量 $O^{(1)}$ を全部のシャドウで推定したときの揺らぎ」の $n$ 倍なので、分散は $n^2\mathbb V_1/N_U$ である。**多コピー観測量の主要な誤差は、 $\rho$ に依存する単一コピーの観測量 $O^{(1)}$ の推定の誤差に帰着する。** そのため5.2節の分解(式(15) [PRX: (B8)])がそのまま使えて

```math
\mathbb V[\hat O]=\frac{n^2}{N_U}\Bigl(\mathbb V_U\bigl[f_{\rho,\sigma}(U)\bigr]+\frac{\mathbb E_U\bigl[g_\rho(U)\bigr]}{N_M}\Bigr)+\mathcal O\Bigl(\frac1{N_U^2}\Bigr)\qquad(\text{式(18) [PRX: (B11)]})
```

となる。

論文の注意:

- 主要項は**バッチの数 $m$ に依存しない**。そこで後処理が最も軽い $m=n$ を選ぶのが普通である。
- $\mathcal O(1/N_U^2)$ の項を具体的に抑えるのは難しいが、 $N_M\to\infty$ かつ $\sigma\to\rho$ では 0 になることは保証される。このとき $\hat\rho^{(r)}\to\sigma^{(r)}$ で、CRM シャドウ $\hat\rho_\sigma^{(r)}\to\sigma$ が揺らがない定数になるからである。
- $n=1$ では $\mathcal O(1/N_U^2)$ の項は最初からない。

**台について。** 多コピーの演算子 $O$ の台 $A=\mathrm{supp}(O)$ を、「少なくとも1つのコピーで $O$ が恒等演算子でない量子ビットの集合」と定める。 $O^{(k)}$ の台は $A$ に含まれるので、 $O^{(k)}=O_A^{(k)}\otimes\mathbb 1_{\bar A}^{\otimes k}$ と書ける( $\bar A$ は $A$ の補集合)。

### 6.4 局所パウリ測定での上界の導出(付録 C.2 [arXiv: B.2])

式(18) [PRX: (B11)]の2つの項を、パウリ測定の場合に上から抑える。 $O^{(1)}=O_A^{(1)}\otimes\mathbb 1_{\bar A}$ とする。

**手順1: パウリ文字列で展開する。** $A$ の上の $4^{N_A}$ 個のパウリ文字列 $\gamma_A$ は $\mathrm{Tr}(\gamma_A\gamma_A')=2^{N_A}\delta_{\gamma_A,\gamma_A'}$ を満たす直交基底なので

```math
O_A^{(1)}=\frac1{2^{N_A}}\sum_{\gamma_A}c_\gamma\,\gamma_A,\qquad c_\gamma\equiv\mathrm{Tr}\bigl(O_A^{(1)}\gamma_A\bigr),\qquad\lVert O_A^{(1)}\rVert_2^2=\frac1{2^{N_A}}\sum_{\gamma_A}c_\gamma^2
```

と書ける( $c_\gamma$ はこの文書での略記)。最後の等式は、展開を2乗してトレースを取り、直交性を使えば出る。 $\gamma_A$ の台を $\Gamma\subseteq A$ 、その量子ビット数を $N_\Gamma$ と書くと、3.4節より $\mathcal M^{-1}(\gamma_A)=3^{N_\Gamma}\gamma_A$ なので

```math
w(U,\mathbf s)=\frac1{2^{N_A}}\sum_{\gamma_A}3^{N_\Gamma}\,c_\gamma\,\gamma_A(U_A,\mathbf s_A)\qquad(\text{式(26) [PRX: (C8)]})
```

である。 $w$ は $A$ の上のユニタリー $U_A$ とビット列 $\mathbf s_A$ にしか依存しない。

**手順2: $f$ は $A$ に縮約した状態だけで決まる。** $w$ が $\mathbf s_A$ にしか依存しないので、 $\bar A$ のビットについての和は確率にだけ掛かり、 $\sum_{\mathbf s_{\bar A}}P_\rho(\mathbf s\mid U)=P_{\rho_A}(\mathbf s_A\mid U_A)$ となる( $\rho_A=\mathrm{Tr}_{\bar A}\rho$ は縮約密度行列)。したがって

```math
f_{\rho,\sigma}(U)=\sum_{\mathbf s_A}\bigl(P_{\rho_A}(\mathbf s_A\mid U_A)-P_{\sigma_A}(\mathbf s_A\mid U_A)\bigr)\,w(U_A,\mathbf s_A)\qquad(\text{式(27) [PRX: (C9)]})
```

**手順3: Cauchy–Schwarz の不等式。** $\bigl(\sum_sa_sb_s\bigr)^2\le\bigl(\sum_sa_s^2\bigr)\bigl(\sum_sb_s^2\bigr)$ を、 $a$ を確率の差、 $b$ を $w$ として使う(式(28) [PRX: (C10)])。

**手順4: 1つめの因子。** $X=U_A(\rho_A-\sigma_A)U_A^\dagger$ と置くと、確率の差は対角成分 $\langle\mathbf s_A\rvert X\lvert\mathbf s_A\rangle$ である。 $X$ がエルミートなので $\langle\mathbf s\rvert X^2\lvert\mathbf s\rangle=\sum_{\mathbf s'}\lvert\langle\mathbf s\rvert X\lvert\mathbf s'\rangle\rvert^2\ge\langle\mathbf s\rvert X\lvert\mathbf s\rangle^2$ で、

```math
\sum_{\mathbf s_A}\bigl(P_{\rho_A}(\mathbf s_A\mid U_A)-P_{\sigma_A}(\mathbf s_A\mid U_A)\bigr)^2\le\sum_{\mathbf s_A}\langle\mathbf s_A\rvert X^2\lvert\mathbf s_A\rangle=\mathrm{Tr}(X^2)=\lVert\rho_A-\sigma_A\rVert_2^2\qquad(\text{式(29) [PRX: (C11)]})
```

となる。ここで prior の良さが Hilbert–Schmidt 距離として入る。

**手順5: 2つめの因子。** $w^2$ を展開すると、2つのパウリ文字列 $\gamma_A$ 、 $\gamma_A'$ の積の和 $\sum_{\mathbf s_A}\gamma_A(U_A,\mathbf s_A)\gamma_A'(U_A,\mathbf s_A)$ が現れる。式(20) [PRX: (C2)]より、それぞれは「基底が当たっていれば台の上のビットの符号の積、外れていれば 0」である。符号の積の和は量子ビットごとの積に分かれ、量子ビット $i$ が片方の台にだけ入っていると $\sum_{s_i}(-1)^{s_i}=0$ で消える。両方に入っているか、どちらにも入っていなければ 2 になる。したがって台が同じ( $\Gamma=\Gamma'$ )ときだけ残り、同じ台で同じ $U_A$ に当たっている2つのパウリ文字列は同じものなので

```math
\sum_{\mathbf s_A}\gamma_A(U_A,\mathbf s_A)\,\gamma_A'(U_A,\mathbf s_A)=2^{N_A}\,\delta_{\gamma_A,\gamma_A'}\,\delta_{U_\Gamma,V_\gamma}\qquad(\text{式(30) [PRX: (C12)]})
```

である。これを使うと、手順3、4と合わせて

```math
f_{\rho,\sigma}(U)^2\le\frac1{2^{N_A}}\lVert\rho_A-\sigma_A\rVert_2^2\sum_{\gamma_A}\delta_{U_\Gamma,V_\gamma}\,9^{N_\Gamma}\,c_\gamma^2\qquad(\text{式(31) [PRX: (C13)]})
```

**手順6: ユニタリーについて平均する。** 分散は2乗の平均以下で、 $\mathbb E_U[\delta_{U_\Gamma,V_\gamma}]=3^{-N_\Gamma}$ 、 $3^{N_\Gamma}\le3^{N_A}$ を順に使う。

```math
\mathbb V_U[f_{\rho,\sigma}(U)]\le\mathbb E_U[f_{\rho,\sigma}(U)^2]\le\frac{\lVert\rho_A-\sigma_A\rVert_2^2}{2^{N_A}}\sum_{\gamma_A}3^{N_\Gamma}c_\gamma^2\le3^{N_A}\,\lVert O_A^{(1)}\rVert_2^2\,\lVert\rho_A-\sigma_A\rVert_2^2\qquad(\text{式(32) [PRX: (C14)]})
```

(補足)論文の式(32) [PRX: (C14)]の左辺は $\mathbb V_U[f_{\rho,\sigma}(U)^2]$ と印刷されているが、前後の式から $\mathbb V_U[f_{\rho,\sigma}(U)]$ の意味である。

**手順7: ショット雑音の項。** $g_\rho$ の定義(式(17) [PRX: (B10)])の引き算の項を落とすと $g_\rho(U)\le\sum_{\mathbf s}P_\rho(\mathbf s\mid U)w(U,\mathbf s)^2$ である。この右辺を $U$ について平均し、状態について最大を取ったものは、古典シャドウの元論文で「シャドウノルム」の2乗と呼ばれる量で、パウリ測定では Hilbert–Schmidt ノルムで抑えられることが知られている(Huang–Kueng–Preskill の命題3)。

```math
\mathbb E_U[g_\rho(U)]\le\lVert O^{(1)}\rVert_{\rm shadow}^2\le2^{N_A}\,\lVert O_A^{(1)}\rVert_2^2\qquad(\text{式(33) [PRX: (C15)]})
```

**結果。** 手順6と7を式(18) [PRX: (B11)]に入れると

```math
\mathbb V[\hat O]\le\frac{n^2\lVert O_A^{(1)}\rVert_2^2}{N_U}\Bigl(3^{N_A}\lVert\rho_A-\sigma_A\rVert_2^2+\frac{2^{N_A}}{N_M}\Bigr)+\mathcal O\Bigl(\frac1{N_U^2}\Bigr)\qquad(\text{式(6)、式(34) [PRX: (C16)]})
```

が得られる。論文はこれを主要な結果の1つと位置づけている。

### 6.5 式(6)の読み方

1. **prior の良さは、台 $A$ に縮約した状態の Hilbert–Schmidt 距離で入る。** 標準のシャドウ( $\sigma\to0$ )では $\lVert\rho_A-\sigma_A\rVert_2^2$ が $\lVert\rho_A\rVert_2^2=\mathrm{Tr}(\rho_A^2)$ ( $\rho_A$ の純度)に置き換わる。
2. **必要な設定の数が大きく減る条件。** 論文は2つの条件を挙げている。1つめは $\lVert\rho_A-\sigma_A\rVert_2^2\ll\lVert\rho_A\rVert_2^2$ 、2つめは $N_M\gg(2/3)^{N_A}\lVert\rho_A\rVert_2^{-2}$ である。(補足)1つめは「第1項が CRM で大きく減る」条件である。2つめは「標準のシャドウで第1項が第2項より大きい」条件 $3^{N_A}\lVert\rho_A\rVert_2^2\gg2^{N_A}/N_M$ を書き直したもので、これが成り立たないと、もともとショット雑音が誤差の主因で、第1項を減らしても全体はあまり減らない。
3. **$m$ に依存しない。** 式(6)は元々の多コピー推定量( $m=N_U$ )にも当てはまる。
4. **台を小さく取れる場合がある。** $\lVert O_A^{(1)}\rVert_2^2$ は状態に依存する。状態と観測量によっては $O^{(1)}$ の台が $A$ より小さくなり、そのときは式(6)の $A$ 、 $N_A$ を $O^{(1)}$ の台とその大きさに置き換えた、より強い上界が成り立つ。
5. **(補足)式(6)は一般的な代わりに緩い。** $n=1$ でパウリ文字列( $\lVert\gamma_A\rVert_2^2=2^{N_A}$ )に使うと、第2項は $4^{N_A}/(N_UN_M)$ になり、式(4)の $3^{N_A}/(N_UN_M)$ より大きい。式(6)は Cauchy–Schwarz の不等式などを何度か使った上界で、等式ではない。パウリ文字列では5章の等式(式(25) [PRX: (C7)])を使うべきである。

### 6.6 シフト演算子の $O_A^{(1)}$ (付録 E [arXiv: C])

モーメント $p_n=\mathrm{Tr}(\rho_A^n)$ は、 $n$ 個のコピーを巡回的にずらすシフト演算子

```math
\tau_A^{(n)}\lvert\mathbf s_A^{(1)}\rangle\lvert\mathbf s_A^{(2)}\rangle\cdots\lvert\mathbf s_A^{(n)}\rangle=\lvert\mathbf s_A^{(n)}\rangle\lvert\mathbf s_A^{(1)}\rangle\cdots\lvert\mathbf s_A^{(n-1)}\rangle
```

を使って $p_n=\mathrm{Tr}(\tau_A^{(n)}\rho_A^{\otimes n})$ と書ける( $n=2$ が1.1節のスワップ演算子)。このとき式(6)に入る演算子は

```math
O_A^{(1)}=\rho_A^{\,n-1}\qquad(\text{式(39) [PRX: (E5)]})
```

である。

(補足) $n=2$ で確かめる。スワップ演算子は最初から対称なので $O_{\rm sym}=\mathbb S$ で、 $O^{(1)}=\mathrm{Tr}_2[\mathbb S(\mathbb 1\otimes\rho)]$ である。行列要素を計算すると

```math
\langle a\rvert O^{(1)}\lvert b\rangle=\sum_c\langle a\rvert\langle c\rvert\,\mathbb S\,(\mathbb 1\otimes\rho)\lvert b\rangle\lvert c\rangle=\sum_c\langle c\rvert\langle a\rvert(\mathbb 1\otimes\rho)\lvert b\rangle\lvert c\rangle=\sum_c\delta_{cb}\,\langle a\rvert\rho\lvert c\rangle=\langle a\rvert\rho\lvert b\rangle
```

で、 $O^{(1)}=\rho$ となる。

一般の $n$ では、論文は次のように示す。置換 $\pi$ で入れ替えたシフト演算子 $W_\pi^\dagger\tau_nW_\pi$ も、コピーを一巡する巡回置換である。1番目のコピー以外を $\rho$ で埋めて部分トレースを取ると、一巡の途中で $\rho$ の行列要素が $n-1$ 個、鎖のようにつながり、完全性 $\sum_{\mathbf s}\lvert\mathbf s\rangle\langle\mathbf s\rvert=\mathbb 1$ で和を取ると $\rho^{n-1}$ になる(式(37)、(38) [PRX: (E3)、(E4)])。結果は $\pi$ によらないので、 $\pi$ についての平均も $\rho^{n-1}$ である。

(補足)したがって $p_n$ の推定では $\lVert O_A^{(1)}\rVert_2^2=\mathrm{Tr}(\rho_A^{2n-2})=p_{2n-2}$ で、例えば純度( $n=2$ )の上界は

```math
\mathbb V[\hat p_2]\le\frac{4\,p_2}{N_U}\Bigl(3^{N_A}\lVert\rho_A-\sigma_A\rVert_2^2+\frac{2^{N_A}}{N_M}\Bigr)+\mathcal O\Bigl(\frac1{N_U^2}\Bigr)
```

となる。

---

## 7. 応用1 — 実験データによるエンタングルメントの検出(PRX 版の III.A 節。arXiv 版にはない)

論文の中で唯一、実機のデータを使った例である。

(補足)この章の言葉。

- **イオントラップ**は、電場で捕まえたイオン(電気を帯びた原子)の列を量子ビットとして使う装置である。
- **クエンチ**は、系のハミルトニアンを急に切り替え、その後の時間発展を見る実験である。時間とともに、離れた部分どうしが量子的に絡み合っていく。
- **エンタングルメントの検出**とは、2つの部分系の状態が「それぞれ独立な状態を確率で混ぜたもの」では書けないことを、測定データから示すことである。
- **デコヒーレンス**は、装置が外界と触れることで、純粋状態が混合状態になっていく過程である。
- **ジャックナイフ法**は、データの一部を順に抜いて推定をやり直し、その散らばりから誤差を見積もる統計の方法である。

### 7.1 データ

- 10量子ビットのイオントラップ量子シミュレータで行われたクエンチ実験(Brydges ら, *Science* 2019)のランダム測定データ。
- $N_U=500$ 個のランダムなユニタリー、ユニタリーごとに $N_M=150$ 回の測定。
- 元の実験では、時間発展の間に2つの部分系の間にエンタングルメントが生まれることを、このデータから検出している。

### 7.2 prior

実験のユニタリーな時間発展を数値シミュレーションして得た純粋状態 $\sigma=\lvert\psi\rangle\langle\psi\rvert$ を使う。実験の状態はデコヒーレンスの影響で混合状態になっているので、このシミュレーションは近似にすぎない。論文の主張は、**それでも CRM シャドウを使えば、標準のシャドウよりずっと少ない測定でエンタングルメントを検出できる**ということである。

### 7.3 判定に使う量

$p_3$ -PPT 条件と呼ばれる判定法を使う。(補足)部分転置は、2つの部分系のうち片方 $A$ の添字だけについて行列を転置する操作である。 $A$ と $B$ がエンタングルしていない状態では、部分転置した行列も固有値がすべて 0 以上になる(PPT は positive partial transpose の略)。負の固有値があればエンタングルしている。負の固有値があるかどうかを、モーメント $p_2$ と $p_3$ だけから判定できる十分条件が、次の不等式である。部分系 $A$ と $B$ の状態 $\rho_{AB}$ を $A$ について部分転置した行列 $\rho_{AB}^{T_A}$ のモーメント

```math
p_n=\mathrm{Tr}\bigl[(\rho_{AB}^{T_A})^n\bigr]
```

について、 $p_2^2/p_3>1$ なら $A$ と $B$ はエンタングルしている。 $p_n$ は $n$ コピーの観測量として書けるので、6章の推定量で推定できる。

### 7.4 結果(PRX 版 Fig. 1(b))

- 対象は、クエンチの時刻 $t=1$ ms での、部分系 $A=[1,2]$ と $B=[3,4]$ の間のエンタングルメント。
- 実験のデータと古典シミュレーションの両方から、同じ $N_U$ 個のユニタリーだけを取り出した縮小データで $p_2^2/p_3$ を推定し、 $N_U$ を変えて比べる。誤差棒はジャックナイフ法で求めた1標準偏差 $\varepsilon$ である。
- 「推定値から $\varepsilon$ を引いても 1 を超える」ことを検出の基準にすると、**CRM シャドウでは $N_U=53$ 、標準のシャドウでは $N_U=299$** で基準を満たす。
- 必要な測定の数が**約6分の1**になり、論文はこれがデータ取得の時間を数時間節約することに当たると述べている。

この例は、prior が実験の状態の正確なモデルでなくても(ノイズのない純粋状態のシミュレーションでも)役に立つこと、そして過去に取ったデータに後から適用できることを示している。

---

## 8. 応用2 — フォン・ノイマン・エントロピーの推定(III.B 節、付録 D [arXiv: 例1、付録 E])

### 8.1 エントロピーを多項式で近似する

部分系 $A$ のフォン・ノイマン・エントロピー $S=-\mathrm{Tr}(\rho_A\log\rho_A)$ は、エンタングルメントの尺度で、量子多体系の相や相転移の区別に使える。しかし $\rho_A$ の多項式ではないので、多コピー観測量としては直接測れない。そこで $\rho_A$ の固有値 $\lambda$ で $S=-\sum_\lambda\lambda\log\lambda$ と書き、関数 $f(x)=-x\log x$ を区間 $(0,1)$ で多項式

```math
f_{n_{\max}}(x)=\sum_{n=1}^{n_{\max}}a_n\,x^n
```

に最小2乗法で近似する。すると

```math
S_{n_{\max}}=\mathrm{Tr}\bigl[f_{n_{\max}}(\rho_A)\bigr]=\sum_{n=1}^{n_{\max}}a_n\,p_n,\qquad p_n=\mathrm{Tr}(\rho_A^n)\qquad(\text{式(7)})
```

は、モーメント $p_n$ の線形結合になる。各 $p_n$ はシフト演算子の期待値(6.6節)なので、バッチの数を $m=n_{\max}$ としたバッチシャドウの推定量(式(5))で推定できる。論文は、この当てはめが、同じ次数の他の多項式補間より正確な近似を与えると述べている。

### 8.2 係数の求め方(付録 D [arXiv: E])

最小にする量は

```math
I_{n_{\max}}=\int_0^1\bigl[f(x)-f_{n_{\max}}(x)\bigr]^2dx\qquad(\text{式(41) [PRX: (D1)]})
```

である。2乗を展開して積分する。

- (補足) $\int_0^1x^n\cdot x^{n'}\thinspace dx=\frac1{1+n+n'}$ 。
- (補足) $\int_0^1x^k\log x\thinspace dx=-\frac1{(k+1)^2}$ (部分積分)を $k=n+1$ で使うと、 $\int_0^1(-x\log x)\thinspace x^n\thinspace dx=\frac1{(n+2)^2}$ 。

したがって

```math
I_{n_{\max}}=-\sum_{n=1}^{n_{\max}}\frac{2a_n}{(2+n)^2}+\sum_{n,n'=1}^{n_{\max}}\frac{a_na_{n'}}{1+n+n'}+\text{定数}\qquad(\text{式(42) [PRX: (D2)]})
```

で、 $a_n$ で微分して 0 と置くと

```math
\sum_{n'=1}^{n_{\max}}\frac{a_{n'}}{1+n+n'}=\frac1{(2+n)^2}\qquad(\text{式(43) [PRX: (D3)]})
```

という連立1次方程式になる。係数の行列 $A_{n,n'}=1/(1+n+n')$ は Cauchy 行列と呼ばれる形( $1/(x_n+y_{n'})$ 、ここでは $x_n=n+1$ 、 $y_{n'}=n'$ )で、逆行列が閉じた式で書けるので、 $a_n$ も閉じた式で書ける(式(44)、(45) [PRX: (D4)、(D5)])。

$n_{\max}=3$ では、論文に書かれているとおり

```math
f_3(x)=\frac{137}{60}x-4x^2+\frac74x^3
```

である。(補足)この文書で式(43)を解いて確かめた。 $n_{\max}=5$ では $a=(\frac{667}{210},-\frac{32}3,18,-16,\frac{11}2)$ になる。

**近似の誤差の上界。** 区間 $[0,1]$ での $\lvert f(x)-f_{n_{\max}}(x)\rvert$ の最大値を $\alpha_{n_{\max}}$ とする。論文の値は $\alpha_3\approx0.046$ 、 $\alpha_4\approx0.028$ 、 $\alpha_5\approx0.019$ である。 $f(0)=f_{n_{\max}}(0)=0$ なので、0 の固有値は誤差に寄与せず、

```math
\lvert S-S_{n_{\max}}\rvert\le\sum_{\lambda\ne0}\bigl\lvert f(\lambda)-f_{n_{\max}}(\lambda)\bigr\rvert\le\alpha_{n_{\max}}\,\mathrm{rank}(\rho_A)\qquad(\text{式(46) [PRX: (D6)]})
```

となる。最小2乗誤差 $I_{n_{\max}}$ は $n_{\max}$ とともに単調に小さくなる(PRX 版 Fig. 5)。

### 8.3 数値実験の設定

(補足)表の言葉。**イジング鎖**は、隣り合うスピンをそろえようとする項 $Z_iZ_{i+1}$ と、スピンを横に倒そうとする項 $X_i$ を持つ1次元の模型である。2つの項の強さが等しい点が**臨界点**で、そこでは遠く離れたスピンどうしにも相関が残り、エンタングルメントが系の大きさとともに増え続ける。**MPS**(行列積状態)は状態を小さな行列の積で表す方法で、行列の大きさ(**結合次元** $\chi$ )が大きいほど正確である。 $\chi=1$ は、量子ビットどうしが完全に独立な状態(**積状態**)に当たる。**DMRG** は、MPS の中でエネルギー最小の状態を探す計算法である。**MPO**(行列積演算子)は、同じ考え方で密度行列などの演算子を表したものである。

| 項目 | 内容 |
|---|---|
| 状態 | 臨界イジング鎖 $H=-\sum_{i=1}^N(Z_iZ_{i+1}+X_i)$ (開放端、 $Z_{N+1}=0$ )の基底状態 $\lvert G\rangle$ 。DMRG で求め、結合次元 $\chi_G\approx40$ の MPS で表す |
| 部分系 | 鎖の半分。 $N_A=N/2=8$ ( $N=16$ ) |
| prior | $\lvert G\rangle$ を結合次元 $\chi=1,2,3$ に切り詰めた MPS $\lvert\psi_\chi\rangle$ 。 $\sigma$ はその縮約状態で、行列積演算子(MPO)になる |
| 測定 | $N_M=1000$ に固定し、 $N_U$ を変える |
| 評価 | 相対誤差の平均。20回(Fig. 2(b))と50回(Fig. 2(c))の模擬についての平均 |

臨界点なので、半分の部分系のエントロピーは $S=\frac c6\log N_A+\text{定数}$ ( $c=1/2$ はイジング普遍性クラスの中心電荷)のように対数的に増える。

### 8.4 結果

- **多項式近似の質(Fig. 2(a))**: $n_{\max}=3$ でも、 $N_A$ についての対数的な増え方は再現される。 $n_{\max}=5,7$ では $S$ と定量的に合う。
- **統計誤差(Fig. 2(b))**: $\hat S_3$ の統計的な相対誤差 $\mathbb E[\lvert\hat S_3-S_3\rvert]/S_3$ を比べる。
  - $\chi=1$ の prior は積状態で、不正確すぎて標準のシャドウからの改善がない。
  - $\chi=2,3$ では prior が十分に正確で、統計誤差が大きく減る。
- **全体の誤差(Fig. 2(c))**: $\chi=2$ で、真の $S$ に対する相対誤差 $\mathbb E[\lvert\hat S_{n_{\max}}-S\rvert]/S$ を $n_{\max}=3,5,7$ で比べる。これは統計誤差と、多項式近似による系統誤差 $\lvert S_{n_{\max}}-S\rvert$ の両方を含む。
  - 測定が少なく統計誤差が支配的なうちは、 $n_{\max}$ が大きいほど誤差が大きい。論文はこれを式(6)の因子 $n^2$ によるものとしている。
  - 測定が多くなると誤差は系統誤差で頭打ちになり、 $n_{\max}$ が大きい方が有利になる。

**統計誤差を減らせると、より高い次数の近似を使えるようになり、系統誤差も減らせる。** これが CRM のこの例での意義である。

---

## 9. 応用3 — 忠実度の推定と prior の反復改善(III.C 節 [arXiv: 例2])

### 9.1 何を推定するか

用意した状態 $\rho$ と、理論上の純粋状態 $\lvert\psi\rangle$ の間の忠実度

```math
\mathcal F_\psi=\langle\psi\rvert\rho\lvert\psi\rangle=\mathrm{Tr}(O\rho),\qquad O=\lvert\psi\rangle\langle\psi\rvert
```

を推定する。単一コピー( $n=1$ )の観測量だが、 $O$ は全系に作用する。

### 9.2 動機

1. **状態の用意の検証。** 忠実度は、量子デバイスが目標の状態を作れたかどうかの証明になる。大域的な Clifford 測定(量子ビットごとに基底を選ぶ代わりに、全量子ビットをまとめてかき混ぜる回路を掛けてから測る方法。回路が深く、実験では難しい)による古典シャドウなら効率的に推定できるが、**局所的なランダム測定では、必要な測定の数が系の大きさについて指数的に増えうる**。そこを CRM で減らしたい。
2. **良い prior を見つける手段。** 式(6)を、全系に台を持つ観測量( $N_A=N$ )と純粋な prior $\sigma=\lvert\phi\rangle\langle\phi\rvert$ に使うと、**$\mathcal F_\phi\ge1/2$ なら CRM の方が標準のシャドウより分散(の上界)が小さい**。

(補足)2の理由。 $\rho$ も純粋とすると、 $\mathrm{Tr}(\rho^2)=\mathrm{Tr}(\sigma^2)=1$ 、 $\mathrm{Tr}(\rho\sigma)=\mathcal F_\phi$ なので

```math
\lVert\rho-\sigma\rVert_2^2=\mathrm{Tr}(\rho^2)-2\,\mathrm{Tr}(\rho\sigma)+\mathrm{Tr}(\sigma^2)=2(1-\mathcal F_\phi)
```

である。標準のシャドウで対応する量は $\lVert\rho\rVert_2^2=1$ なので、式(6)の第1項が小さくなる条件 $2(1-\mathcal F_\phi)\le1$ は $\mathcal F_\phi\ge1/2$ である。 $\rho$ が混合状態のときは $\lVert\rho-\sigma\rVert_2^2=\mathrm{Tr}(\rho^2)-2\mathcal F_\phi+1$ を $\mathrm{Tr}(\rho^2)$ と比べることになり、同じ条件 $\mathcal F_\phi\ge1/2$ が出る。

### 9.3 prior を反復で選ぶ手順

1. prior $\sigma=\lvert\phi\rangle\langle\phi\rvert$ から始め、 $\mathcal F_\phi$ を CRM シャドウか標準のシャドウで推定する。どちらを使うかは、後処理の段階で経験的な分散を比べて決めればよい。
2. $\mathcal F_\phi$ が閾値 $\mathsf F$ ( $\mathsf F\ge1/2$ )を下回ったら、新しい prior を作る(より多くの古典計算を使ってよい)。1 に戻る。
3. 十分に高い忠実度の prior が見つかったら、それを使って任意の多コピー観測量(他の状態との忠実度を含む)を CRM で推定する。性能は、測った $\mathcal F_\phi$ の値を入れた式(6)で保証される。

**この手順全体が、1つのランダム測定のデータの上で行える。** prior は後処理でしか使わないからである(4.3節)。重点サンプリング(基底を一様にではなく、prior から見て情報の多そうな基底を高い確率で選ぶ方法)による忠実度の推定では、測る設定の選び方が prior に依存するので、prior を変えるたびに測り直しになる。

### 9.4 数値実験

(補足)表の言葉。**ランダム量子回路**は、ランダムに選んだゲート(量子ビットへの操作)を層状に重ねた回路で、深さ $d$ は層の数である。**Haar ランダム**は、ユニタリー全体から偏りなく一様に選ぶことを指す。**脱分極雑音**は、確率 $p$ で量子ビットの状態が完全にでたらめな状態に置き換わる雑音である。これがあると出力は混合状態になる。

| 項目 | 内容 |
|---|---|
| 状態 | 1次元のランダム量子回路の出力。 $N=30$ 量子ビット、深さ $d=6$ 。1量子ビットと隣接2量子ビットの Haar ランダムなゲートの層を交互に重ねる |
| ノイズ | 各ゲートに強さ $p$ の局所的な脱分極雑音。 $p=0$ (純粋状態)と、 $p=10^{-3}$ 、 $10^{-4}$ (混合状態) |
| prior | **ノイズのない回路**の厳密な出力を、結合次元 $\chi$ に切り詰めた MPS $\lvert\phi\rangle$ |
| 測定 | $N_U=15$ 、 $N_M=10^5$ |
| 誤差棒 | ランダムなユニタリーについての平均の標準誤差 |

この prior では、CRM の推定値は MPS の計算で $\chi$ の多項式の時間で求まる。 $\chi$ を上げると忠実度 $\mathcal F_\phi$ は上がるが、計算のコストも増える。

**結果(Fig. 3)**:

- $\chi$ を上げると、推定された忠実度が上がり、CRM の誤差棒が小さくなる。
- $\chi$ が小さいとき(忠実度が低いとき)は、CRM は改善を与えず、誤差棒は標準のシャドウより**大きい**。prior が悪いと損をするという、2章の一般論のとおりである。
- ノイズのある状態( $p>0$ )でも同じ傾向が見られる。このとき純粋な prior $\lvert\phi\rangle$ は、 $\chi$ をいくら上げても混合状態の $\rho$ とは一致しない。

---

## 10. 別の実験から prior を作る(付録 A [arXiv: 付録 D])

### 10.1 想定する状況

$\sigma$ の理論的な記述がなく、代わりに、 $\rho$ に近い状態 $\sigma$ を実現する**別の実験**(companion experiment)がある場合を考える。要求するのは、 $\sigma$ が $\rho$ の近似になっていることだけである。

論文は前提を次のように整理している。もし $\sigma=\rho$ が保証されているなら、両方の実験のデータを全部まとめて標準のシャドウを作る方が効率的だろう。想定しているのは、別の実験が $\rho$ を正確には再現していないが、**ずっと多くのユニタリー( $N_U'\gg N_U$ )で測れる**場合である。装置が速い、あるいは過去の複数の実験のデータを集められる、といった状況がこれに当たる。

### 10.2 作り方

1. 元の実験から、 $N_U$ 個の古典シャドウ $\hat\rho^{(r)}$ を作る( $r=1,\dots,N_U$ )。
2. 別の実験から、**1と同じユニタリー $U^{(r)}$ を使って** $N_U$ 個の古典シャドウ $\hat\sigma^{(r)}$ を作る。
3. 別の実験から、**独立に選んだ** $N_U'$ 個のユニタリーで、追加の古典シャドウ $\hat\sigma^{(r)}$ を作る( $r=N_U+1,\dots,N_U+N_U'$ )。測る順番は問わない。

3つのデータからそれぞれ $m$ 個のバッチシャドウ $\hat\rho^{[t]}$ 、 $\hat\sigma^{[t]}$ 、 $(\hat\sigma')^{[t]}$ を作り、

```math
\hat\rho_\sigma^{[t]}=\hat\rho^{[t]}-\hat\sigma^{[t]}+(\hat\sigma')^{[t]}\qquad(\text{式(40) [PRX: (A1)]})
```

と組み合わせる。

式(2)との対応は次のとおりである。

| 式(2) | 式(40) [PRX: (A1)] | 役割 |
|---|---|---|
| $\hat\rho^{(r)}$ | $\hat\rho^{[t]}$ | 実験のシャドウ |
| $\sigma^{(r)}$ (厳密な確率で計算) | $\hat\sigma^{[t]}$ (**同じユニタリー**で測った別の実験のシャドウ) | 共通の乱数を持つ引き算の相手 |
| $\sigma$ (厳密) | $(\hat\sigma')^{[t]}$ (**独立なユニタリー**で測った別の実験のシャドウ) | 引いた分の平均を足し戻す |

**不偏性**: $\mathbb E[\hat\sigma^{[t]}]=\mathbb E[(\hat\sigma')^{[t]}]=\sigma$ なので、 $\mathbb E[\hat\rho_\sigma^{[t]}]=\rho-\sigma+\sigma=\rho$ である。足し戻す項を引き算の項とは独立なデータから作るのは、同じデータを使うと2つが打ち消し合って意味がなくなるからである。

$N_U'$ が有限であることは追加の統計誤差を生む。 $N_U'\gg1$ の極限では $(\hat\sigma')^{[t]}\to\sigma$ となって式(2)の CRM シャドウと同じになり、式(6)の上界が当てはまる。(補足) $\hat\sigma^{[t]}$ も頻度から作るので、別の実験のショット雑音も加わる。

### 10.3 数値例(PRX 版 Fig. 4)

- 8章と同じ臨界イジング鎖で、 $N=12$ 、 $N_A=6$ 、 $N_M=1000$ 、バッチの数 $m=n_{\max}=3$ 。
- 別の実験は、少しずれたハミルトニアン $H'=H+\sum_i\epsilon_iZ_i$ ( $\epsilon_i$ は $[0,0.02]$ から独立に選ぶ)の基底状態 $\lvert G'\rangle$ を実現するとする。
- $N_U'=4000$ では、誤差は大きく減るが、 $N_U'$ が有限であることによる**頭打ち**が見える。
- $N_U'$ を10倍に増やすと頭打ちの高さが下がり、調べた範囲のすべての $N_UN_M$ で標準のシャドウより良い推定になる。

---

## 11. 結論と展望(IV 節)

論文のまとめ: CRM シャドウは、ランダム測定の後処理に $\rho$ の近似的な知識を取り入れることで、線形・多コピーの観測量の推定を大きく改善する、すぐに使える道具である。

**応用として挙げられているもの**:

- 変分量子アルゴリズムでの勾配の推定
- 物質の量子相の探査
- 量子過程(チャネル)の性質の推定への拡張(PRX 版で追加)

**今後の課題として挙げられているもの**:

- **重点サンプリングや適応的な手法との組み合わせ**(純度 $p_2$ のランダム測定で開発されたもの)
- 補助系を使った測定の方法との組み合わせ
- 量子実験と古典シミュレーションの相関を使って測定誘起のエンタングルメント相転移を調べる最近の研究との組み合わせ、より一般に、古典シミュレーションをモデルに依存しない量子的な性質の学習にどう役立てるか(PRX 版で追加)

論文の Julia のコードは ITensor と PastaQ を使っており、[RandomMeas.jl](https://github.com/bvermersch/RandomMeas.jl) で公開されている。この研究のリポジトリはその派生である。

---

## 12. 読むときの注意(補足)

この章は論文の内容ではなく、論文を使うときに区別しておくべき点をこの文書でまとめたものである。

### 12.1 等式と上界

| 式 | 種類 | 使いどころ |
|---|---|---|
| 式(25) [PRX: (C7)] の左辺 | **等式**(単一コピーのパウリ文字列) | 分散を正確に予測したいとき。標準のシャドウとの比も正確に出る |
| 式(4) | 上界。 $3^{N_A}-1\le3^{N_A}$ と $1-\mathrm{Tr}(\rho O)^2\le1$ を使ったもの | 見通しのよい目安 |
| 式(18) [PRX: (B11)] | 主要項は等式、 $\mathcal O(1/N_U^2)$ は具体的に抑えられていない | 多コピー観測量の一般的な構造 |
| 式(6)、式(34) [PRX: (C16)] | 上界。Cauchy–Schwarz、 $3^{N_\Gamma}\le3^{N_A}$ 、シャドウノルムの上界を使ったもの | 性能の保証。実際の分散はこれより小さいのが普通 |

「忠実度 $1/2$ 以上なら得をする」(9.2節)は式(6)の上界どうしの比較で、分散そのものの比較ではない。

### 12.2 prior の良さを測る量は観測量の種類で違う

| 観測量 | 分散に入る prior の誤差 |
|---|---|
| 単一コピーのパウリ文字列 | その観測量の期待値の外れ $\mathrm{Tr}[O(\rho-\sigma)]$ だけ |
| 多コピー観測量(一般) | 台に縮約した状態の Hilbert–Schmidt 距離 $\lVert\rho_A-\sigma_A\rVert_2$ (上界として) |
| 全系に台を持つ観測量と純粋な prior | 忠実度 $\mathcal F_\phi$ ( $\lVert\rho-\sigma\rVert_2^2=2(1-\mathcal F_\phi)$ を通じて) |

局所的なパウリ観測量では、大域的な忠実度がどれだけ低くても、その観測量の外れが小さければ得をする。この研究の「局所観測量の利得は忠実度ではなく外れで決まる」という結果は、論文の式(25) [PRX: (C7)]のこの性質を系統的に調べたものである([README_crm_paper.md](README_crm_paper.md) 7章)。

### 12.3 CRM が減らすのは設定のくじの分散だけ

5.2節の表のとおり、ショット雑音の項は prior では減らない。したがって、CRM の効果は $N_M$ が大きいほど大きく、1設定1ショット( $N_M=1$ )では限られる。論文の例がどれも $N_M$ を大きく取っている( $150$ 、 $1000$ 、 $10^5$ )のはこのためで、論文自身が「較正に時間がかかり $N_U$ が限られる実験で有用」と位置づけている(5.4節)。

### 12.4 prior が悪いと損をする

推定は常に不偏だが、分散は増えうる。論文の数値例にも、 $\chi=1$ の積状態(8.4節)や、忠実度の低い prior(9.4節)で改善がない、あるいは悪くなる場合が示されている。ただし、標準のシャドウと CRM シャドウは同じデータから作れるので、**経験的な分散を比べて良い方を使えばよい**(9.3節の手順1)。

### 12.5 印刷上の注意

- 式(8) [PRX: (B1)]の前の因子は $N_U/m$ と印刷されているが、平均の因子 $m/N_U$ と読む(6.2節)。
- 式(32) [PRX: (C14)]の左辺は $\mathbb V_U[f_{\rho,\sigma}(U)^2]$ と印刷されているが、 $\mathbb V_U[f_{\rho,\sigma}(U)]$ の意味である(6.4節)。
- arXiv 版と PRX 版では付録の順番が違う(0.3節の表)。arXiv 版の付録 A, B, C, D, E が、PRX 版の付録 B, C, E, A, D に当たる。

### 12.6 論文が扱っていないこと

- どんな物理系でどんな prior が良いか(論文の prior は、厳密な状態を切り詰めた MPS、ノイズのないシミュレーション、別の実験の3種類)。
- 観測量による違い、対称性を破った prior、読み出し誤差への対処。

これらがこの研究の対象で、対応は [README_crm_paper.md](README_crm_paper.md) 9章にまとめてある。

---

## 13. この文書の限界

- 図の内容は、論文の本文と図の説明文に書かれていることを使った。図そのものから読んだのは、PRX 版 Fig. 4 の凡例にある $N_U'N_M$ の2つの値( $4.0\times10^6$ と $4.0\times10^7$ 、10.3節の「10倍」の根拠)だけである。
- 式(33) [PRX: (C15)]のシャドウノルムの上界は、古典シャドウの元論文(Huang–Kueng–Preskill)の結果の引用で、この文書では導いていない。
- 式(10) [PRX: (B3)]の分散の展開は Rath ら(2023)の結果の引用で、この文書では主要項の意味の説明と、展開の係数の確認だけを行った。
- 「(補足)」と書いた導出と数値例は、この文書で補ったものである。5.4節の数値例と8.2節の係数は、この文書で計算した。
