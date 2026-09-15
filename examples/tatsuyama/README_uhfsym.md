# UHF-sym とは何か — スピンの向きを平均した UHF

[README.md](README.md) に戻る / 使っている文書: [README_eigenstate.md](README_eigenstate.md)(3.5〜4)、[README_observables.md](README_observables.md)、[README_experiments.md](README_experiments.md)(§4)

---

## この文書の要点

1. **UHF-sym は、UHF の状態をあらゆるスピンの向きに回して平均した prior** である。UHF が勝手に選んだ「スピンの向き」の情報だけを消し、それ以外はそのまま残す。
2. **1つの波動関数ではなく混合状態**(確率的な混ぜ合わせ)だが、CRM では prior の期待値 $\langle P\rangle_\sigma$ を古典計算できれば十分なので、そのまま使える。**追加の測定はいらない**。
3. **スピン回転で変わらない量(密度、二重占有、オンサイト ZZ、 $\mathbf S_i\cdot\mathbf S_j$ など)は UHF とぴったり同じ値**になり、**向きを持つ量( $S^z_i$ )は 0** になり、**2つのスピンの相関は $x,y,z$ の3方向に均等に配り直される**。
4. 計算は、観測量の「回転で変わらない部分」を取り出すだけで、UHF と同じ手間で済む。
5. **直せるのは対称性の破れによる誤差だけ**である。UHF の相関は古典的なネール状態のものなので、量子的な一重項の相関の不足は直せない。

---

## 1. なぜ必要か: UHF はスピンの向きを勝手に選ぶ

半充填の Hubbard 模型の UHF 解は反強磁性(ネール)状態で、各サイトのスピンが上下交互に向いている。このとき UHF は**スピンの向く軸を1本選んでいる**(この研究では $z$ 軸)。ところが、Hubbard 模型はスピンをどの向きに回しても変わらないので、 $x$ 軸や斜めの軸に向いた UHF もまったく同じエネルギーを持つ。 **$z$ 軸を選んだのは近似の都合であって、物理ではない。**

一方、有限の系の真の基底状態はスピン一重項(全スピン 0)で、**特定の向きを持たない**。1次元 $L=64$ 、 $U=12$ で比べると:

| 量 | 真の状態 | UHF |
|---|---|---|
| サイトのスピン $\langle S^z_i\rangle$ | 0 | **0.486** |
| 隣のスピンの相関 $\langle S^z_iS^z_j\rangle$ | $-0.1275$ | **$-0.2398$** |
| 隣のスピンの相関 $\langle S^x_iS^x_j\rangle$ | $-0.1275$ | **$-0.0034$** |
| 隣の上向き電子の相関 $\langle Z_{i\uparrow}Z_{j\uparrow}\rangle$ | $-0.526$ | **$-0.973$** |

真の状態では $z$ 方向と $x$ 方向の相関が等しいのに、UHF は相関をほとんど $z$ 方向に集めている。**この「向きの偏り」が、スピンの量で UHF prior が損をする主な原因だった**([README_eigenstate.md](README_eigenstate.md) 3.3、[README_observables.md](README_observables.md) 3)。

UHF がこう振る舞うのには理由がある。1枚の Slater 行列式では $d_i=n_i^2/4-\lvert\langle\mathbf S_i\rangle\rvert^2$ が厳密に成り立ち、二重占有を減らす(局所モーメントを作る)にはスピンの向きを決めるしかない([README_observables.md](README_observables.md) 3.1)。**UHF は電荷の量を正しくするために、スピンの向きという本来ない情報を付け足している。** UHF-sym はその付け足した情報だけを取り除く。

---

## 2. 定義: あらゆる向きに回して平均する

系全体のスピンを、軸 $\hat n$ のまわりに角度 $\alpha$ だけ回す演算子を $R=\exp(-i\alpha\,\hat n\cdot\mathbf S_{\rm tot})$ とする。UHF-sym は、UHF の状態 $\sigma_{\rm UHF}$ を**あらゆる回転で回して、均等な重みで平均した状態**である:

```math
\bar\sigma=\int dR\;R\,\sigma_{\rm UHF}\,R^\dagger
```

( $\int dR$ はすべての回転についての一様な平均。) 直感的には、**ネール状態の向きを表す矢印を、全方位に均等に向けて重ね合わせたもの**である。コンパスの針をあらゆる向きに向けて平均すると「平均の向き」は 0 になるが、針の長さは残る。それと同じで、UHF-sym では**スピンの向きの情報は消え、局所モーメントの大きさや隣どうしが逆向きであることは残る**。

$\bar\sigma$ は1つの波動関数ではなく、回した UHF たちを確率的に混ぜた**混合状態**である。量子コンピュータでこの状態を作るのは面倒だが、**CRM では prior の状態を実際に用意する必要はない**。必要なのは各パウリ文字列の期待値 $\langle P\rangle_{\bar\sigma}$ という数だけで、それは古典計算で求まる。したがって UHF-sym を使っても**測定データは同じで、後処理が変わるだけ**である。

---

## 3. 何が変わり、何が変わらないか

期待値の定義から

```math
\langle O\rangle_{\bar\sigma}=\int dR\;\mathrm{Tr}\bigl[O\,R\sigma_{\rm UHF}R^\dagger\bigr]=\int dR\;\bigl\langle R^\dagger O R\bigr\rangle_{\rm UHF}
```

である。つまり**状態を回す代わりに観測量を回して、UHF で期待値を取って平均すればよい**。

### 3.1 回転で変わらない量は、UHF とぴったり同じ値になる

$R^\dagger OR=O$ なら、平均しても $\langle O\rangle_{\bar\sigma}=\langle O\rangle_{\rm UHF}$ である。スピンの向きによらない量はすべてこれに当たる:

- サイトの電子数 $n_i$ 、二重占有 $d_i=n_{i\uparrow}n_{i\downarrow}$ 、オンサイト ZZ $Z_{i\uparrow}Z_{i\downarrow}$
- 電荷の偶奇 $(-1)^{n_i}$ とその積
- 2つのスピンの内積 $\mathbf S_i\cdot\mathbf S_j$
- 両スピンを足したホッピング $\sum_\sigma(c^\dagger_{i\sigma}c_{j\sigma}+{\rm h.c.})$

**UHF-sym は UHF の「良いところ」(局所モーメントによる電荷の量の正しさ)をまったく損なわない。**

### 3.2 向きを持つ量は 0 になる

サイトのスピンは UHF で $\langle\mathbf S_i\rangle=m\,\hat e_z$ ( $m$ は磁化)である。回すと向きが $\hat n$ に変わり、あらゆる $\hat n$ で平均すると

```math
\langle\mathbf S_i\rangle_{\bar\sigma}=m\int\frac{d\hat n}{4\pi}\;\hat n=0
```

になる(単位ベクトルを球面上で平均すると 0)。実際に §2j の計算では、UHF の $\langle S^z_i\rangle=+0.3845$ が平均後に $-6.2\times10^{-17}$ (計算機の精度でゼロ)になっている。

### 3.3 2つのスピンの相関は、3方向に均等に配り直される

2つのスピンの相関 $\langle S^a_iS^b_j\rangle$ ( $a,b=x,y,z$ )を回して平均するには、単位ベクトルの成分の積の球面平均

```math
\int\frac{d\hat n}{4\pi}\;n_an_b=\frac{\delta_{ab}}{3}\qquad\bigl(\text{例えば }n_z^2\text{ の平均は }1/3\bigr)
```

を使う。回転は3方向を入れ替えるだけで、3方向の和 $\mathbf S_i\cdot\mathbf S_j$ は変えないので、どんな状態から出発しても

```math
\langle S^a_iS^b_j\rangle_{\bar\sigma}=\frac{\delta_{ab}}{3}\,\langle\mathbf S_i\cdot\mathbf S_j\rangle_{\rm UHF}
```

になる。**相関の総量は UHF のまま、3方向に均等に配り直される。** 1次元 $L=64$ 、 $U=12$ の隣の結合では:

| | $\langle S^z_iS^z_j\rangle$ | $\langle S^x_iS^x_j\rangle$ | $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ |
|---|---|---|---|
| UHF | $-0.2398$ | $-0.0034$ | $-0.2465$ |
| **UHF-sym** | **$-0.0822$** | **$-0.0822$** | $-0.2465$ (同じ) |
| 真の状態 | $-0.1275$ | $-0.1275$ | $-0.3826$ |

UHF-sym は真の状態と同じく等方的になる。ただし総量 $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ は UHF のまま( $-0.2465$ 、古典的なネール状態の値 $-1/4$ に近い)で、真の値 $-0.3826$ には届かない(5.2)。

### 3.4 この研究で使った観測量の変換表

JW 変換で $Z_{i\uparrow}=1-2n_{i\uparrow}$ 、 $n_{i\uparrow}=n_i/2+S^z_i$ なので、 $Z_{i\uparrow}=(1-n_i)-2S^z_i$ と「電荷の部分」と「スピンの部分」に分けられる。これに 3.1〜3.3 を当てはめると:

| 観測量 | UHF-sym での値 | 1次元 $L=64$ 、 $U=12$ の例(UHF → UHF-sym) |
|---|---|---|
| $n_i$ 、二重占有、 $Z_{i\uparrow}Z_{i\downarrow}$ 、 $(-1)^{n_i}$ | UHF と同じ | $Z_{i\uparrow}Z_{i\downarrow}$ : $-0.9456$ → $-0.9456$ |
| $S^z_i$ | 0 | $0.486$ → $0$ |
| $Z_{i\uparrow}=(1-n_i)-2S^z_i$ | $1-n_i$ (半充填で 0) | $-0.972$ → $0$ |
| $S^z_iS^z_j$ 、 $S^x_iS^x_j$ | $\tfrac13\langle\mathbf S_i\cdot\mathbf S_j\rangle$ | $-0.240$ 、 $-0.003$ → $-0.082$ |
| $Z_{i\uparrow}Z_{j\uparrow}$ | $\langle c_ic_j\rangle+\tfrac43\langle\mathbf S_i\cdot\mathbf S_j\rangle$ ( $c_i=1-n_i$ ) | $-0.973$ → $-0.342$ |
| $c^\dagger_{i\uparrow}c_{j\uparrow}+{\rm h.c.}$ | 上向きと下向きの平均 | 反強磁性の隣の結合では同じ値 |

$Z_{i\uparrow}Z_{j\uparrow}$ の行は、 $Z_{i\uparrow}Z_{j\uparrow}=c_ic_j-2c_iS^z_j-2S^z_ic_j+4S^z_iS^z_j$ と展開すると、スピンを1つだけ含む交差項は 3.2 で 0 に、 $S^z_iS^z_j$ は 3.3 で $\tfrac13\mathbf S_i\cdot\mathbf S_j$ に、電荷だけの $c_ic_j$ は 3.1 でそのまま、となることから出る。 $\lvert y\rvert$ が UHF の 0.973 から 0.342 に下がるのは、 $z$ 方向に集中していた相関が3方向に配り直されて、 $z$ 成分が約3分の1になるためである。

---

## 4. どう計算するか

### 4.1 解析的に: 回転で変わらない部分を取り出す

3.4 の表のとおり、観測量を「回転で変わらない部分」と「向きを持つ部分」に分け、前者だけを UHF で計算すればよい。**必要なのは UHF の期待値だけ**なので、手間は UHF と同じである。この研究の [crm_eigen_uhfsym.py](crm_eigen_uhfsym.py)( $Z_{i\uparrow}Z_{j\uparrow}$ )と [crm_new_observables_analysis.py](crm_new_observables_analysis.py)(二重占有・局所エネルギーの各項)はこの方法を使っている。

### 4.2 数値的に: 回した UHF で計算して、向きについて平均する

回した UHF も、上向きと下向きの軌道を混ぜただけの(一般化された)Slater 行列式である。したがって、回した状態の $2n\times2n$ の相関行列に対して Wick の定理で任意のパウリ文字列の期待値が計算でき、それを向きについて数値的に平均すればよい。既存の実装([crm_1d_hf.jl](crm_1d_hf.jl) の `symavg`、[crm_2d_symuhf.jl](crm_2d_symuhf.jl))は、 $\cos\theta$ について一様な 64 方向で平均している。

**注意(§2i で直した点)**: $z$ 軸から傾ける角 $\theta$ だけで平均すると、 $z$ 軸まわりの回転で変わらない量( $Z$ 型の量やホッピング)には厳密だが、横向きの $S^x$ と $S^y$ を区別する量には足りない。そのような量には方位角 $\varphi$ についての平均も必要である。

2つの方法は一致する: $L=8$ 、 $U=4$ の $\langle S^zS^z\rangle$ で、解析的には $(\langle S^zS^z\rangle+2\langle S^xS^x\rangle)/3=-0.07333$ 、数値的な球面平均では $-0.07331$ (差は数値積分の格子の誤差)である。

---

## 5. CRM の利得にどう効くか

### 5.1 直せる誤差: 対称性の破れ

UHF の外れは「対称性の破れの誤差」と「相関の誤差」の和に分けられる([README_observables.md](README_observables.md) 3.2)。**UHF-sym はこのうち破れの誤差だけを厳密に消す。**

| 観測量 | UHF | **UHF-sym** | 詳細 |
|---|---|---|---|
| オンサイト ZZ(1次元 $L=64$ 、 $U=12$ ) | 472 | 472(変わらない) | 回転で不変 |
| 単一サイト $Z_\uparrow$ ( $U=12$ ) | 0.016 | **1.00** | [README_eigenstate.md](README_eigenstate.md) 2.6 |
| 隣接 $Z_\uparrow Z_\uparrow$ (1次元 $L=64$ 、 $U=12$ ) | 1.38 | **6.78** | 同 3.5 |
| 隣接 $Z_\uparrow Z_\uparrow$ で損をする系(80 系) | 17 | **0** | 同 3.5 |
| 二重占有を4項の和で測る( $U=12$ 、80 系の中央値) | 1.81 | **123** | [README_observables.md](README_observables.md) 5.2 |
| 局所エネルギー( $U=12$ 、5結合の中央値) | 1.65 | **50.9** | 同 5.4 |

### 5.2 直せない誤差: 古典的な相関

- **量子的な一重項の相関は足せない。** 3.3 のとおり $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ は UHF の値のまま( $-0.25$ 付近)で、真の値( $-0.38$ )には届かない。
- **遠くのスピン相関では損をする。** UHF は長距離秩序を持つので、平均しても離れたサイトの相関が $\lvert y\rvert\approx1/3$ のまま残る。真の相関は距離とともに減衰するので、1次元では $r=4$ 付近から損に転じる([README_eigenstate.md](README_eigenstate.md) 4.3)。
- **対称性で値が 0 になる量では、利得は 1 までしか戻らない。** 天井そのものが 1 だからである。
- **相関が強い結合では、かえって悪くなることがある。** 平均で $z$ 成分が約3分の1に下がるため、真の相関がそれより十分強い結合では控えめすぎになる(1次元の強い結合で 9.45 → 3.38)。

### 5.3 平均と射影の違い(UHF-sym は真の一重項ではない)

UHF-sym は回した UHF を**確率的に混ぜた**だけで、干渉がないので、回転で変わらない量はすべて UHF とまったく同じ値のままになる(3.1)。**新しい相関は1つも生まれない。** 一方、回した UHF を**量子的に重ね合わせて**全スピン 0 の成分だけを取り出す方法(スピン射影 HF)では、重ね合わせの干渉によって $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ のような回転で変わらない量も変わり、量子的な一重項の相関に近づく可能性がある。**この研究では射影は計算していない**。5.2 の「直せない誤差」を減らす次の候補である。

---

## 6. いつ UHF-sym を使うべきか

**prior の対称性は、測る状態の対称性に合わせる**のが原則である([README_experiments.md](README_experiments.md) §4)。

| 測る状態 | 向いている prior | 理由 |
|---|---|---|
| 対称性を保つ状態(有限系の厳密な基底状態、ドープして秩序が溶けた状態) | **UHF-sym** | 真の状態に向きがないので、向きの情報は誤差にしかならない |
| 状態自体が対称性を破っている(収束しきっていない2次元の DMRG、対称性を破る外場がある実験) | 共線 UHF のことがある | 状態にも向きがあるので、UHF の向きが合う場合がある |

prior は古典的な後処理なので、**同じ測定データに対して UHF と UHF-sym の両方を試し、観測量ごとに選んでよい**。目安として、観測量が回転で変わらない量だけでできていれば両者は同じ値になる。向きを持つ成分を含むなら、状態が対称なときは UHF-sym を選ぶ。
