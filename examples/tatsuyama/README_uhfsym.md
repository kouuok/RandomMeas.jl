# UHF-sym とは何か — スピンの向きを平均した UHF

[README.md](README.md) に戻る / 使っている文書: [README_eigenstate.md](README_eigenstate.md)(3.5〜4)、[README_observables.md](README_observables.md)、[README_experiments.md](README_experiments.md)(§4) / 続き: [README_phf_mps.md](README_phf_mps.md)、[README_varsym.md](README_varsym.md)、[README_varsym_distance.md](README_varsym_distance.md)

---

## この文書の要点

1. **UHF-sym は、UHF の状態をあらゆるスピンの向きに回して平均した prior** である。UHF が勝手に選んだ「スピンの向き」の情報だけを消し、それ以外はそのまま残す。
2. **1つの波動関数ではなく混合状態**(確率的な混ぜ合わせ)だが、CRM では prior の期待値 $\langle P\rangle_\sigma$ を古典計算できれば十分なので、そのまま使える。**追加の測定はいらない**。
3. **スピン回転で変わらない量**(密度、二重占有、オンサイト ZZ、 $\mathbf S_i\cdot\mathbf S_j$ など)**は UHF とぴったり同じ値**になる。**向きを持つ量**( $S^z_i$ )**は 0** になる。**2つのスピンの相関は $x,y,z$ の3方向に均等に配り直される**。
4. 計算は、観測量の「回転で変わらない部分」を取り出すだけで、UHF と同じ手間で済む。
5. **直せるのは対称性の破れによる誤差だけ**である。UHF の相関は古典的なネール状態のものなので、量子的な一重項の相関の不足は直せない。遠い距離の相関に残る「偽の長距離秩序」も直せない。
6. 直せない誤差を直す候補として、**スピン射影 HF と低い結合次元の MPS** を計算して比べた(5.3、詳しくは [README_phf_mps.md](README_phf_mps.md))。その後、**変分 MPS を UHF-sym と同じように回転平均した prior** が、電荷とスピンの両方で損をしないことが分かった([README_varsym.md](README_varsym.md))。

## 前提として知っておくこと

| 用語・記号 | 意味 |
|---|---|
| UHF | 非制限 Hartree–Fock。上向きと下向きの電子を別々の Slater 行列式で表す平均場。半充填の Hubbard 模型では反強磁性(ネール)解になる |
| $n_{i\sigma}$ 、 $n_i$ | サイト $i$ 、スピン $\sigma$ の電子数と、サイトの電子数 $n_{i\uparrow}+n_{i\downarrow}$ |
| $\mathbf S_i$ | サイト $i$ のスピン。 $S^z_i=\tfrac12(n_{i\uparrow}-n_{i\downarrow})$ |
| $d_i$ | 二重占有 $\langle n_{i\uparrow}n_{i\downarrow}\rangle$ |
| $Z_{i\uparrow}$ | 上向きの軌道の量子ビットのパウリ $Z$ 。Jordan–Wigner 変換で $Z_{i\uparrow}=1-2n_{i\uparrow}$ |
| 利得 $G$ | 標準の古典シャドウの分散 ÷ CRM の分散。損得の境目は $\lvert\Delta\rvert=\lvert x\rvert$ ( $x$ は真の値、 $\Delta$ は prior の外れ) |

---

## 1. なぜ必要か: UHF はスピンの向きを勝手に選ぶ

半充填の Hubbard 模型の UHF 解は反強磁性(ネール)状態で、各サイトのスピンが上下交互に向いている。このとき UHF は**スピンの向く軸を1本選んでいる**(この研究では $z$ 軸)。ところが、Hubbard 模型のハミルトニアンは全スピンをどの向きに回しても変わらないので、 $x$ 軸や斜めの軸に向いた UHF もまったく同じエネルギーを持つ。**$z$ 軸を選んだのは近似の都合であって、物理ではない。**

一方、有限の系の真の基底状態はスピン一重項(全スピン 0)で、**特定の向きを持たない**。一重項は回転で変わらないので、どのサイトでも $\langle\mathbf S_i\rangle=0$ で、相関は $x,y,z$ の3方向で等しい。1次元 $L=64$ 、 $U=12$ で比べると:

| 量 | 真の状態 | UHF |
|---|---|---|
| サイトのスピン $\langle S^z_i\rangle$ | 0 | **0.486** |
| 隣のスピンの相関 $\langle S^z_iS^z_j\rangle$ | $-0.1275$ | **$-0.2398$** |
| 隣のスピンの相関 $\langle S^x_iS^x_j\rangle$ | $-0.1275$ | **$-0.0034$** |
| 隣の上向き電子の相関 $\langle Z_{i\uparrow}Z_{j\uparrow}\rangle$ | $-0.526$ | **$-0.973$** |

真の状態では $z$ 方向と $x$ 方向の相関が等しいのに、UHF は相関をほとんど $z$ 方向に集めている。**この「向きの偏り」が、スピンの量で UHF prior が損をする主な原因だった**([README_eigenstate.md](README_eigenstate.md) 3.3、[README_observables.md](README_observables.md) 3)。

**UHF がこう振る舞う理由。** UHF では上向きと下向きの電子が独立なので、二重占有は $d_i=\langle n_{i\uparrow}\rangle\langle n_{i\downarrow}\rangle$ である。 $\langle n_{i\uparrow}\rangle=\tfrac12\langle n_i\rangle+\langle S^z_i\rangle$ 、 $\langle n_{i\downarrow}\rangle=\tfrac12\langle n_i\rangle-\langle S^z_i\rangle$ を入れると

```math
d_i=\Bigl(\tfrac12\langle n_i\rangle\Bigr)^2-\langle S^z_i\rangle^2
```

となる(スピンが $z$ 方向だけを向く共線の UHF では $\langle S^z_i\rangle^2=\lvert\langle\mathbf S_i\rangle\rvert^2$ )。半充填で $\langle n_i\rangle=1$ なら $d_i=\tfrac14-\lvert\langle\mathbf S_i\rangle\rvert^2$ で、**二重占有を減らす(相互作用 $U$ のエネルギーを下げる)には、磁化 $\langle\mathbf S_i\rangle$ を持つしかない**。真の状態は電子どうしの量子的な相関で二重占有を減らせるが、UHF にはその手段がないので、スピンの向きを決めて代わりにする([README_observables.md](README_observables.md) 3.1)。**UHF は電荷の量を正しくするために、スピンの向きという本来ない情報を付け足している。** UHF-sym はその付け足した情報だけを取り除く。

---

## 2. 定義: あらゆる向きに回して平均する

系全体のスピンを、軸 $\hat n$ のまわりに角度 $\alpha$ だけ回す演算子を $R=\exp(-i\alpha\,\hat n\cdot\mathbf S_{\rm tot})$ とする( $\mathbf S_{\rm tot}=\sum_i\mathbf S_i$ )。UHF-sym は、UHF の状態 $\sigma_{\rm UHF}$ を**あらゆる回転で回して、均等な重みで平均した状態**である:

```math
\bar\sigma=\int dR\;R\,\sigma_{\rm UHF}\,R^\dagger
```

( $\int dR$ はすべての回転についての一様な平均で、 $\int dR\,1=1$ と規格化する。)

直感的には、**ネール状態の向きを表す矢印を、全方位に均等に向けて混ぜ合わせたもの**である。コンパスの針をあらゆる向きに向けて平均すると「平均の向き」は 0 になるが、針の長さは残る。それと同じで、UHF-sym では**スピンの向きの情報は消え、局所モーメントの大きさや隣どうしが逆向きであることは残る**。

$\bar\sigma$ は1つの波動関数ではなく、回した UHF たちを確率的に混ぜた**混合状態**である。量子コンピュータでこの状態を作るのは面倒だが、**CRM では prior の状態を実際に用意する必要はない**。必要なのは各パウリ文字列の期待値 $\langle P\rangle_{\bar\sigma}$ という数だけで、それは古典計算で求まる。したがって UHF-sym を使っても**測定データは同じで、後処理が変わるだけ**である。

---

## 3. 何が変わり、何が変わらないか

期待値の定義とトレースの巡回性 $\mathrm{Tr}[ABC]=\mathrm{Tr}[CAB]$ から

```math
\langle O\rangle_{\bar\sigma}=\int dR\;\mathrm{Tr}\bigl[O\,R\sigma_{\rm UHF}R^\dagger\bigr]=\int dR\;\mathrm{Tr}\bigl[(R^\dagger OR)\,\sigma_{\rm UHF}\bigr]=\int dR\;\bigl\langle R^\dagger O R\bigr\rangle_{\rm UHF}
```

である。つまり**状態を回す代わりに観測量を回して、UHF で期待値を取って平均すればよい**。

### 3.1 回転で変わらない量は、UHF とぴったり同じ値になる

$R^\dagger OR=O$ なら、積分の中身は回転によらず $\langle O\rangle_{\rm UHF}$ なので、平均しても $\langle O\rangle_{\bar\sigma}=\langle O\rangle_{\rm UHF}$ である。スピンの向きによらない量はすべてこれに当たる:

- サイトの電子数 $n_i$ 、二重占有 $n_{i\uparrow}n_{i\downarrow}$ (同じサイトの上下2つの電子はスピンが打ち消し合った一重項なので、回転で変わらない)、オンサイト ZZ $Z_{i\uparrow}Z_{i\downarrow}=1-2n_i+4n_{i\uparrow}n_{i\downarrow}$
- 電荷の偶奇 $(-1)^{n_i}$ とその積
- 2つのスピンの内積 $\mathbf S_i\cdot\mathbf S_j$ (ベクトルの内積は回転で変わらない)
- 両スピンを足したホッピング $\sum_\sigma(c^\dagger_{i\sigma}c_{j\sigma}+{\rm h.c.})$

**UHF-sym は UHF の「良いところ」(局所モーメントによる電荷の量の正しさ)をまったく損なわない。**

### 3.2 向きを持つ量は 0 になる

サイトのスピンは UHF で $\langle\mathbf S_i\rangle=m\,\hat e_z$ ( $m$ は磁化)である。回転 $R$ で回すと向きが単位ベクトル $\hat n$ に変わる。 $R$ をすべての回転について平均することは、 $\hat n$ を球面上で一様に平均することに当たるので

```math
\langle\mathbf S_i\rangle_{\bar\sigma}=m\int\frac{d\hat n}{4\pi}\;\hat n=0
```

になる(単位ベクトルを球面上で平均すると、向きが打ち消し合って 0)。実際に §2j の計算では、UHF の $\langle S^z_i\rangle=+0.3845$ が平均後に $-6.2\times10^{-17}$ (計算機の精度でゼロ)になっている。

### 3.3 2つのスピンの相関は、3方向に均等に配り直される

2つのスピンの相関 $\langle S^a_iS^b_j\rangle$ ( $a,b=x,y,z$ )を回して平均するには、単位ベクトルの成分の積の球面平均

```math
\int\frac{d\hat n}{4\pi}\;n_an_b=\frac{\delta_{ab}}{3}
```

を使う。 $a\ne b$ では $n_a$ の符号が反転した向きと打ち消し合って 0 になり、 $a=b$ では3方向が対等で $n_x^2+n_y^2+n_z^2=1$ なので、それぞれ $1/3$ になる。回転は3方向を混ぜるだけで、3方向の和 $\mathbf S_i\cdot\mathbf S_j$ は変えないので、どんな状態から出発しても

```math
\langle S^a_iS^b_j\rangle_{\bar\sigma}=\frac{\delta_{ab}}{3}\,\langle\mathbf S_i\cdot\mathbf S_j\rangle_{\rm UHF}
```

になる。**相関の総量は UHF のまま、3方向に均等に配り直される。** 1次元 $L=64$ 、 $U=12$ の隣の結合では:

| | $\langle S^z_iS^z_j\rangle$ | $\langle S^x_iS^x_j\rangle$ | $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ |
|---|---|---|---|
| UHF | $-0.2398$ | $-0.0034$ | $-0.2465$ |
| **UHF-sym** | **$-0.0822$** | **$-0.0822$** | $-0.2465$ (同じ) |
| 真の状態 | $-0.1275$ | $-0.1275$ | $-0.3826$ |

UHF の $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ は $-0.2398+2\times(-0.0034)=-0.2465$ で、UHF-sym の各方向は $-0.2465/3=-0.0822$ である。UHF-sym は真の状態と同じく等方的になる。ただし総量 $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ は UHF のまま( $-0.2465$ 、古典的なネール状態の値 $-1/4$ に近い)で、真の値 $-0.3826$ には届かない(5.2)。

### 3.4 この研究で使った観測量の変換表

$n_{i\uparrow}=\tfrac12n_i+S^z_i$ ( $n_i+2S^z_i=2n_{i\uparrow}$ から)なので、

```math
Z_{i\uparrow}=1-2n_{i\uparrow}=(1-n_i)-2S^z_i
```

と「電荷の部分」と「スピンの部分」に分けられる。以下 $c_i\equiv1-n_i$ と書く。これに 3.1〜3.3 を当てはめると:

| 観測量 | UHF-sym での値 | 1次元 $L=64$ 、 $U=12$ の例(UHF → UHF-sym) |
|---|---|---|
| $n_i$ 、二重占有、 $Z_{i\uparrow}Z_{i\downarrow}$ 、 $(-1)^{n_i}$ | UHF と同じ | $Z_{i\uparrow}Z_{i\downarrow}$ : $-0.9456$ → $-0.9456$ |
| $S^z_i$ | 0 | $0.486$ → $0$ |
| $Z_{i\uparrow}=c_i-2S^z_i$ | $\langle c_i\rangle=1-\langle n_i\rangle$ (半充填で 0) | $-0.972$ → $0$ |
| $S^z_iS^z_j$ 、 $S^x_iS^x_j$ | $\tfrac13\langle\mathbf S_i\cdot\mathbf S_j\rangle$ | $-0.240$ 、 $-0.003$ → $-0.082$ |
| $Z_{i\uparrow}Z_{j\uparrow}$ | $\langle c_ic_j\rangle+\tfrac43\langle\mathbf S_i\cdot\mathbf S_j\rangle$ | $-0.973$ → $-0.342$ |
| $c^\dagger_{i\uparrow}c_{j\uparrow}+{\rm h.c.}$ | 上向きと下向きの平均 | 反強磁性の隣の結合では同じ値 |

$Z_{i\uparrow}Z_{j\uparrow}$ の行は次のように出る。展開すると

```math
Z_{i\uparrow}Z_{j\uparrow}=c_ic_j-2c_iS^z_j-2S^z_ic_j+4S^z_iS^z_j
```

で、スピンを1つだけ含む交差項は 3.2 と同じ理由で 0 に、 $S^z_iS^z_j$ は 3.3 で $\tfrac13\mathbf S_i\cdot\mathbf S_j$ に、電荷だけの $c_ic_j$ は 3.1 でそのまま残る。数値では $\langle c_ic_j\rangle\approx-0.013$ 、 $\tfrac43\times(-0.2465)=-0.329$ で、和が $-0.342$ である。 $\lvert y\rvert$ が UHF の 0.973 から 0.342 に下がるのは、 $z$ 方向に集中していた相関が3方向に配り直されて、 $z$ 成分が約3分の1になるためである。

最後の行(ホッピング)は、回転で上向きと下向きが混ざるため、上向きだけのホッピングの値が上下の平均になる。半充填の反強磁性の UHF では、隣の結合で上向きと下向きのホッピングが等しいので値は変わらない。

---

## 4. どう計算するか

### 4.1 解析的に: 回転で変わらない部分を取り出す

3.4 の表のとおり、観測量を「回転で変わらない部分」と「向きを持つ部分」に分け、前者だけを UHF で計算すればよい。**必要なのは UHF の期待値だけ**なので、手間は UHF と同じである。この研究の [crm_eigen_uhfsym.py](crm_eigen_uhfsym.py)( $Z_{i\uparrow}Z_{j\uparrow}$ )と [crm_new_observables_analysis.py](crm_new_observables_analysis.py)(二重占有・局所エネルギーの各項)はこの方法を使っている。

### 4.2 数値的に: 回した UHF で計算して、向きについて平均する

回した UHF も、上向きと下向きの軌道を混ぜただけの(一般化された)Slater 行列式である。したがって、回した状態の $2n\times2n$ の相関行列に対して Wick の定理で任意のパウリ文字列の期待値が計算でき、それを向きについて数値的に平均すればよい。既存の実装([crm_1d_hf.jl](crm_1d_hf.jl) の `symavg`、[crm_2d_symuhf.jl](crm_2d_symuhf.jl))は、 $\cos\theta$ について一様な 64 方向で平均している。

**極角 $\theta$ だけで平均してよい場合。** 共線の UHF は $z$ 軸まわりの回転で変わらない( $z$ 方向の磁化しか持たない)。そのため、回転の効果は「 $z$ 軸をどこへ向けるか」(極角 $\theta$ と方位角 $\varphi$ )だけで決まる。さらに観測量も $z$ 軸まわりの回転で変わらない量( $Z$ だけでできた量やホッピング)なら、 $\varphi$ にもよらないので、 $\theta$ についての平均だけでよい。球面上の一様な平均では $\cos\theta$ が $[-1,1]$ で一様に分布するので、 $\cos\theta$ について一様に点を取る。

**注意(§2i で直した点)**: 横向きの $S^x$ と $S^y$ を区別する量は $z$ 軸まわりの回転で変わるので、 $\theta$ だけの平均では足りず、方位角 $\varphi$ についての平均も必要である。

**2つの方法は一致する**: $L=8$ 、 $U=4$ の $\langle S^zS^z\rangle$ で、解析的には $(\langle S^zS^z\rangle+2\langle S^xS^x\rangle)/3=-0.07333$ 、数値的な球面平均では $-0.07331$ (差は数値積分の格子の誤差)である。

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

**表の読み方**:

- **単一サイト $Z_\uparrow$ が 1.00 で止まる理由。** 真の値は $x=0$ である。利得の式 $G=\dfrac{(3^{\lvert A\rvert}-1)x^2+v_s}{(3^{\lvert A\rvert}-1)\Delta^2+v_s}$ に $x=0$ を入れると、分子は $v_s$ だけになり、 $G\le1$ で、 $\Delta=0$ のときだけ 1 になる。UHF-sym の値は $1-\langle n_i\rangle=0$ なので $\Delta=0$ 、 $G=1$ である。CRM で減らせる分散(基底のくじの揺らぎ $\propto x^2$ )がもともとないので、これ以上は良くならない。
- **二重占有が大きく改善する理由。** 二重占有 $\tfrac14(1-Z_\uparrow-Z_\downarrow+Z_\uparrow Z_\downarrow)$ のうち、UHF の $\langle Z_\uparrow\rangle$ と $\langle Z_\downarrow\rangle$ は磁化の分だけ逆向きに外れる。合計では打ち消し合うが、CRM の分散には項ごとの外れの積が入るので打ち消されない。UHF-sym では両方とも $1-\langle n_i\rangle$ になり、項ごとの外れが消える。

### 5.2 直せない誤差: 古典的な相関

- **量子的な一重項の相関は足せない。** 3.3 のとおり $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ は UHF の値のまま( $-0.25$ 付近)で、真の値( $-0.38$ )には届かない。反平行の積状態 $\lvert\uparrow\downarrow\rangle$ の $\mathbf S_i\cdot\mathbf S_j$ は $-1/4$ 、一重項は $-3/4$ で、真の状態は一重項の成分を多く含む。UHF は積状態の側にとどまる。
- **遠くのスピン相関では損をする。** UHF の相関には、距離によらない $\langle\mathbf S_i\rangle\cdot\langle\mathbf S_j\rangle=\pm m^2$ (長距離秩序)が含まれる。 $\mathbf S_i\cdot\mathbf S_j$ は回転で変わらないので、平均しても 3.4 の表の $\tfrac43\langle\mathbf S_i\cdot\mathbf S_j\rangle$ に $\pm\tfrac43m^2$ (強結合で $m\approx\tfrac12$ なので約 $\tfrac13$ )が残り、離れたサイトでも $\lvert y\rvert\approx1/3$ のままになる。真の相関は距離とともに減衰するので、1次元では $r=4$ 付近から $\lvert\Delta\rvert>\lvert x\rvert$ となって損に転じる([README_eigenstate.md](README_eigenstate.md) 4.3、[README_varsym_distance.md](README_varsym_distance.md))。
- **対称性で値が 0 になる量では、利得は 1 までしか戻らない。** 5.1 の表の読み方のとおり、天井そのものが 1 だからである。
- **相関が強い結合では、かえって悪くなることがある。** 平均で $z$ 成分が約3分の1に下がるため、真の相関がそれより十分強い結合(例えば開放端の端の、ほぼ一重項になっている結合)では控えめすぎになる(1次元の強い結合で 9.45 → 3.38)。

### 5.3 平均と射影の違い、低い結合次元の MPS との比較 → **[README_phf_mps.md](README_phf_mps.md)**

UHF-sym は回した UHF を**確率的に混ぜた**もので、異なる向きの状態の間の干渉がない。そのため、回転で変わらない量はすべて UHF とまったく同じ値のままになる(3.1)。**新しい相関は1つも生まれない。** 5.2 の「直せない誤差」を直せる prior として、次の2つを計算し、同じ表で比べた。

- **スピン射影 HF(PHF)**: 回した UHF を**量子的に重ね合わせて**、全スピン 0 の成分だけを取り出したもの。重ね合わせの干渉で、 $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ のような回転で変わらない量も変わる。
- **低い結合次元の MPS**: 厳密な状態を結合次元 $\chi$ に切り詰めた**切断 MPS**と、結合次元を $\chi$ に制限した DMRG の基底状態である**変分 MPS**。

計算方法、導出、結果の詳しい説明は **[README_phf_mps.md](README_phf_mps.md)** にまとめた。主な数値だけを再掲する(利得は $U=12$ の 20 系の中央値。1次元 $L=64$ の値は中央の結合):

| | UHF | UHF-sym | PHF | 切断 $\chi=4$ | 切断 $\chi=8$ | 変分 $\chi=4$ | 変分 $\chi=16$ |
|---|---|---|---|---|---|---|---|
| 隣接 $Z_\uparrow Z_\uparrow$ の利得 | 1.22 | 6.95 | 14.2 | 29.2 | 29.4 | 15.3 | 28.8 |
| オンサイト ZZ の利得 | 465 | 465 | 513 | 151 | 546 | 553 | 555 |
| 単一サイト $Z_\uparrow$ の利得 | 0.02 | 1.00 | 1.00 | 1.00 | 1.00 | 0.04 | 0.46 |
| 二重占有(4項の和)の利得 | 1.81 | 123 | 127 | 77.8 | 128 | 4.75 | 53.2 |
| $L=64$ の隣の $\langle\mathbf S\cdot\mathbf S\rangle$ (真 $-0.383$ ) | $-0.247$ | $-0.247$ | $-0.262$ | $-0.395$ | $-0.377$ | $-0.407$ | $-0.366$ |
| $L=64$ の相関エネルギーの回収率 | 0% | 0% | 8% | $-162$ % | 91% | 88% | 100% |

(相関エネルギーの回収率は、UHF のエネルギーから厳密なエネルギーまでの差のうち、その prior が取り戻した割合である。負の値は UHF よりエネルギーが高いことを表す。)

- **PHF の効果は $1/L$ で消える。** 射影は UHF の2つの副格子の大きなスピンを一重項に組み直すだけなので、副格子が違う対の相関が距離によらず $2m/L$ だけ強まり、エネルギーは系全体で $2mJ$ しか下がらない( $U=12$ 、 $L=64$ で隣の結合の $L\,\Delta\langle\mathbf S\cdot\mathbf S\rangle$ が予測 $-0.972$ 、実測 $-0.963$ )。射影してから軌道を最適化(VAP)しても変わらない。
- **切断 MPS( $\chi=4$ )はスピン相関が正確だが、電荷のゆらぎを捨てる**(中央の二重占有が 0)。厳密な状態がないと作れない。
- **変分 MPS は電荷のゆらぎとエネルギーが正確だが、低い $\chi$ では UHF と同じくスピンの向きを選ぶ**( $\chi=4$ で中央の磁化 0.30)。そのため、向きを持つ量やそれを含む和では損をする。

**実用上の結論**: 数十〜数百量子ビットの系では、PHF は CRM の prior として UHF-sym からの改善が小さく、計算の手間に見合わない。5.2 の「古典的な相関の不足」は結合ごとの量子的な相関を持つ MPS で直せる。実際に作れる変分 MPS では対称性の破れを別に取り除く必要があり、それは UHF-sym と同じ回転平均でできる。回転平均した変分 MPS( $\chi\ge8$ )は、中央のサイトとその隣の観測量で損が1つもなかった([README_varsym.md](README_varsym.md))。

## 6. いつ UHF-sym を使うべきか

**prior の対称性は、測る状態の対称性に合わせる**のが原則である([README_experiments.md](README_experiments.md) §4)。

| 測る状態 | 向いている prior | 理由 |
|---|---|---|
| 対称性を保つ状態(有限系の厳密な基底状態、ドープして秩序が溶けた状態) | **UHF-sym** | 真の状態に向きがないので、向きの情報は誤差にしかならない |
| 状態自体が対称性を破っている(収束しきっていない2次元の DMRG、対称性を破る外場がある実験) | 共線 UHF のことがある | 状態にも向きがあるので、UHF の向きが合う場合がある |

prior は古典的な後処理なので、**同じ測定データに対して UHF と UHF-sym の両方を試し、観測量ごとに選んでよい**。目安は次のとおりである。

- 観測量が回転で変わらない量だけでできていれば、両者は同じ値になるので、どちらでもよい。
- 向きを持つ成分を含み、測る状態が対称なら、UHF-sym を選ぶ。
- 遠い距離のスピン相関には、どちらも向かない(5.2)。
