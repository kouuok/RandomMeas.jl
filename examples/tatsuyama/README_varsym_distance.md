# 回転平均した prior は遠いスピン相関でも効くか

[README.md](README.md) に戻る / 前の段階: [README_varsym.md](README_varsym.md)(中央のサイトとその隣だけを調べた)

スクリプト:
- [crm_varsym_pairs.jl](crm_varsym_pairs.jl)(clara で実行。保存済みの変分 MPS から、全サイト対について必要な量を求める)
- [crm_varsym_distance.py](crm_varsym_distance.py)(clara で実行。距離ごとの利得を数える)

出力: [crm_varsym_distance.out](crm_varsym_distance.out)

---

## この文書で分かること

- スピンの向きを1つ選んでしまう prior(UHF や低い結合次元の変分 MPS)は、距離によらない「偽の長距離秩序」を持つ。
- その偽の秩序は、スピン回転で平均しても消えない。そのため、回転平均した UHF(UHF-sym)は隣の対では得をするが、距離 $r\ge4$ の対ではすべて損をする( $U\ge4$ )。
- 偽の秩序を差し引いた「連結相関の prior」は、どの状態の期待値でもないが、CRM の推定を偏らせない。損はしないが、遠い対では得もしない。
- 全距離で得をするのは、対称性を保つ量子的な prior(切断 MPS、結合次元の十分大きい変分 MPS の回転平均)である。

## 前提として知っておくこと

| 用語・記号 | 意味 |
|---|---|
| $n_{i\sigma}$ | サイト $i$ 、スピン $\sigma\in\{\uparrow,\downarrow\}$ の電子数(0 か 1) |
| $n_i$ | サイト $i$ の電子数 $n_{i\uparrow}+n_{i\downarrow}$ |
| $\mathbf S_i$ | サイト $i$ のスピン。 $S^z_i=\tfrac12(n_{i\uparrow}-n_{i\downarrow})$ |
| $Z_{i\uparrow}$ | サイト $i$ の上向きスピンの軌道を表す量子ビットのパウリ $Z$ 。Jordan–Wigner 変換で $Z_{i\uparrow}=1-2n_{i\uparrow}$ |
| 利得 $G$ | 標準の古典シャドウの分散 ÷ CRM の分散。 $G>1$ なら CRM が得、 $G<1$ なら損 |
| 損得の境目 | 単一のパウリ文字列では $\lvert\Delta\rvert=\lvert x\rvert$ ( $x$ は真の値、 $\Delta=x-y$ は prior の外れ。元論文 II.C 節)。 $\lvert\Delta\rvert>\lvert x\rvert$ なら損 |
| UHF | 非制限 Hartree–Fock。各サイトのスピンを副格子ごとに上下どちらかに傾けた平均場 |
| UHF-sym | UHF を全スピン回転で平均した混合状態([README_uhfsym.md](README_uhfsym.md)) |
| 変分 MPS | 結合次元を $\chi$ に制限した DMRG の基底状態。低い $\chi$ ではスピンの向きを選ぶことがある([README_phf_mps.md](README_phf_mps.md)) |

---

## 1. 観測量と prior の値

### 1.1 観測量

距離 $r=\lvert j-i\rvert$ (周期境界では短い方の距離)の $Z_{i\uparrow}Z_{j\uparrow}$ を測る。台は2量子ビットで、1つの設定あたりのショット数は $n_m=100$ とする。

- 1次元 $L=64$ 開放端: 端の影響を避けるため、 $17\le i<j\le48$ の対だけを使う。
- 1次元 $L=16$ 周期境界: すべての対を使う。

距離ごとに、利得の中央値と、1% より大きく損をする( $G<0.99$ の)対の割合を数えた。真の値は厳密な基底状態から取った(`crm_eigcorr_*.tsv`)。

### 1.2 $Z_{i\uparrow}$ を電荷とスピンに分ける

$n_{i\uparrow}$ は $n_i$ と $S^z_i$ で書ける:

```math
n_i+2S^z_i=(n_{i\uparrow}+n_{i\downarrow})+(n_{i\uparrow}-n_{i\downarrow})=2n_{i\uparrow}
\quad\Longrightarrow\quad
n_{i\uparrow}=\tfrac12\bigl(n_i+2S^z_i\bigr).
```

これを $Z_{i\uparrow}=1-2n_{i\uparrow}$ に入れると

```math
Z_{i\uparrow}=1-n_i-2S^z_i
```

である。2つのサイトの積は

```math
Z_{i\uparrow}Z_{j\uparrow}=(1-n_i)(1-n_j)\;-\;2(1-n_i)S^z_j\;-\;2S^z_i(1-n_j)\;+\;4S^z_iS^z_j
```

となる。第1項は電荷だけ、第2・3項はスピンの1次、第4項はスピンの2次である。

### 1.3 回転平均した値

全スピンを回転 $R$ で回して平均した混合状態 $\bar\sigma=\int dR\,R\sigma R^\dagger$ での期待値は、観測量の側を回して平均したものに等しい。回転したとき、電荷 $n_i$ は変わらず、スピン $\mathbf S_i$ はベクトルとして回る。

- **スピンの1次の項。** $S^z_j$ は、回転した後には単位ベクトル $\hat{\mathbf u}$ 方向の成分 $\hat{\mathbf u}\cdot\mathbf S_j$ になる。 $\hat{\mathbf u}$ をすべての向きで平均すると、 $\overline{\hat u_a}=0$ なので第2・3項の平均は 0 になる。
- **スピンの2次の項。** $S^z_iS^z_j$ は $(\hat{\mathbf u}\cdot\mathbf S_i)(\hat{\mathbf u}\cdot\mathbf S_j)=\sum_{a,b}\hat u_a\hat u_b\,S^a_iS^b_j$ になる。方向の平均は $\overline{\hat u_a\hat u_b}=\tfrac13\delta_{ab}$ (等方的で、 $\sum_a\hat u_a^2=1$ だから各成分が $1/3$ )なので、平均は $\tfrac13\mathbf S_i\cdot\mathbf S_j$ である。

したがって、回転平均した prior の値は

```math
\langle Z_{i\uparrow}Z_{j\uparrow}\rangle_{\bar\sigma}=\bigl\langle(1-n_i)(1-n_j)\bigr\rangle_\sigma+\tfrac43\langle\mathbf S_i\cdot\mathbf S_j\rangle_\sigma
```

で、元の状態 $\sigma$ の $\langle(1-n_i)(1-n_j)\rangle$ と $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ だけから計算できる。

### 1.4 比べる prior

| prior | $\langle Z_{i\uparrow}Z_{j\uparrow}\rangle$ として使う値 |
|---|---|
| UHF、変分 MPS | 状態の期待値そのもの |
| -sym(回転平均) | $\langle(1-n_i)(1-n_j)\rangle+\tfrac43\langle\mathbf S_i\cdot\mathbf S_j\rangle$ |
| -conn(連結相関) | $\langle(1-n_i)(1-n_j)\rangle+\tfrac43\bigl(\langle\mathbf S_i\cdot\mathbf S_j\rangle-\langle\mathbf S_i\rangle\cdot\langle\mathbf S_j\rangle\bigr)$ |
| 切断 MPS | 状態の期待値そのもの(厳密な状態を切り詰めて作るので、厳密な状態が要る) |

変分 MPS は $S^z$ の総和を保存量として扱うので、 $\langle S^x_i\rangle=\langle S^y_i\rangle=0$ である。したがって $\langle\mathbf S_i\rangle\cdot\langle\mathbf S_j\rangle=\langle S^z_i\rangle\langle S^z_j\rangle$ になる。UHF でも磁化は $z$ 方向だけを向くので同じである。

**照合**: [crm_varsym_pairs.jl](crm_varsym_pairs.jl) が出した変分 MPS の $\langle Z_{i\uparrow}Z_{j\uparrow}\rangle$ は、以前の計算(`crm_mpslow_*.tsv`)の全対 34176 個と完全に一致した。

---

## 2. 結果

### 2.1 1次元 $L=64$ 開放端、 $U=12$

各セルは「利得の中央値 / 1% より大きく損をする対の割合」で、損をする対がないときは割合を省いた。 $r=5$ 〜8 などはその範囲の距離の対をまとめたものである。

| prior | $r=1$ | $r=2$ | $r=3$ | $r=4$ | $r=5$ 〜8 | $r=9$ 〜16 | $r=17$ 〜31 |
|---|---|---|---|---|---|---|---|
| UHF | 3.3 | 0.12 / 100% | 0.10 / 100% | 0.04 / 100% | 0.03 / 100% | 0.02 / 100% | 0.01 / 100% |
| UHF-sym | 4.5 | 3.3 | 2.7 | **0.61 / 100%** | 0.32 / 100% | 0.17 / 100% | 0.13 / 100% |
| UHF-conn | 1.09 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| 切断 $\chi=8$ | 60.1 | 5.2 | 4.5 | 2.2 | 1.6 | 1.16 | 1.03 |
| 変分 $\chi=4$ | 60.0 | 3.1 | 1.8 / 48% | 0.46 / 100% | 0.25 / 100% | 0.13 / 100% | 0.10 / 100% |
| 変分 $\chi=4$ -sym | 31.2 | 5.7 | 3.9 | 2.5 | 1.7 | 0.86 / 78% | 0.63 / 100% |
| 変分 $\chi=8$ -sym | 43.4 | 5.7 | 4.9 | 2.5 | 1.7 | 1.21 | 1.06 / 2% |
| 変分 $\chi=16$ (= -sym) | 56.4 | 5.7 | 5.2 | 2.5 | 1.8 | 1.21 | 1.06 |

**表の読み方**:

- UHF は $r\ge2$ のすべての対で大きく損をする。UHF-sym は $r\le3$ では得をするが、 $r\ge4$ ではすべての対で損をする。
- UHF-conn は損をしないが、 $r\ge2$ では利得がちょうど 1 で、得もしない。
- 切断 MPS と $\chi=16$ の変分 MPS は、最も遠い対まで損をしない。
- 遠い対ほど利得そのものは 1 に近づく。真の相関が小さくなり、CRM で減らせる分散( $(3^{\lvert A\rvert}-1)x^2$ )が小さくなるためで、prior の良し悪しとは別の理由である。
- $U=12$ の $\chi=16$ では、変分 MPS はもともと対称な解( $\langle S^z_i\rangle=0$ )に落ちているので、回転平均しても値が変わらない。

### 2.2 $r\ge2$ の全対での損の割合

各セルは「1次元 $L=64$ 開放端 / $L=16$ 周期境界」での、 $r\ge2$ の対のうち 1% より大きく損をする対の割合である。

| prior | $U=2$ | $U=4$ | $U=8$ | $U=12$ |
|---|---|---|---|---|
| UHF | 97% / 85% | 100% / 100% | 100% / 100% | 100% / 100% |
| UHF-sym | 62% / 0% | 87% / 69% | 87% / 69% | 87% / 69% |
| UHF-conn | 6% / 31% | 6% / 15% | 0% / 0% | 0% / 0% |
| 切断 $\chi=8$ | 0% / 0% | 0% / 0% | 0% / 0% | 0% / 0% |
| 変分 $\chi=4$ | 90% / 57% | 90% / 56% | 90% / 56% | 90% / 58% |
| 変分 $\chi=4$ -sym | 70% / 0% | 54% / 0% | 54% / 0% | 52% / 0% |
| 変分 $\chi=8$ -sym | 13% / 0% | 3% / 0% | 1% / 0% | 1% / 0% |
| 変分 $\chi=16$ -sym | 0% / 0% | 0% / 0% | 0% / 0% | 0% / 0% |

**表の読み方**: 回転平均(-sym)で損が減るのは、変分 MPS では $\chi$ が大きいほど顕著で、 $\chi=16$ では損がなくなる。UHF-sym は $U\ge4$ で開放端の 87%、周期境界の 69% の対で損をし、回転平均では直らない。

---

## 3. 解釈

### 3.1 回転平均で残る偽の秩序

UHF の各サイトは磁化 $\langle\mathbf S_i\rangle=\pm m\,\hat{\mathbf z}$ を持つ(副格子ごとに符号が交互)。どんな状態でも、相関は

```math
\langle\mathbf S_i\cdot\mathbf S_j\rangle
=\underbrace{\langle\mathbf S_i\rangle\cdot\langle\mathbf S_j\rangle}_{\text{UHF では }\pm m^2}
+\underbrace{\bigl\langle(\mathbf S_i-\langle\mathbf S_i\rangle)\cdot(\mathbf S_j-\langle\mathbf S_j\rangle)\bigr\rangle}_{\text{連結部分}}
```

と分けられる(右辺を展開すると左辺に戻る)。UHF では第1項が $\pm m^2$ で、**距離 $r$ がいくら大きくても変わらない**。これが「どこまで行っても揃ったまま」の古典的な Néel 秩序である。

一方、有限の系の真の基底状態は全スピン 0 の一重項で、 $\langle\mathbf S_i\rangle=0$ である。1次元ではスピン相関は距離とともに(べき的に)減衰する。

回転平均は $\langle\mathbf S_i\rangle$ の**向き**を平均するので、 $\langle Z_{i\uparrow}\rangle$ のずれ(1.2 節の1次の項)は消える。しかし $\mathbf S_i\cdot\mathbf S_j$ は回転で変わらない量なので、上の第1項 $\pm m^2$ は回転平均した値にもそのまま残る。1.3 節の式から、UHF-sym の値は真の値に比べて遠い対でおよそ $\tfrac43m^2$ (強結合で $m\to\tfrac12$ なので約 $0.33$ )だけ外れ続ける。

遠い対では真の値 $x$ が小さくなるので、ある距離で $\lvert\Delta\rvert>\lvert x\rvert$ (損得の境目)を越える。2.1 節の表では、それが $r=4$ だった。

[README_eigenstate.md](README_eigenstate.md) の「長距離秩序を持つ prior は遠くの相関で必ず損をする」は、回転平均しても変わらない。**対称性の破れの「向き」は平均で消えるが、秩序の「大きさ」は消えない。**

### 3.2 連結相関の prior

偽の秩序を取り除くには、上の第1項を引いた連結部分だけを使えばよい。これが 1.4 節の「-conn」である。

**この値は、どの状態の期待値でもない。** どんな状態でも $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ と $\langle\mathbf S_i\rangle\cdot\langle\mathbf S_j\rangle$ はその状態から別々に決まり、勝手に一方だけを取り除くことはできないからである。

**それでも CRM の推定は偏らない。** 1つの設定 $u$ での CRM の推定値は、実験の推定値 $X(u)$ から $3^{\lvert A\rvert}\,\mathbf 1[\text{一致}]\,y$ を引き、 $y$ を足し戻したものである。設定について平均すると

```math
\mathbb E\bigl[X-3^{\lvert A\rvert}\mathbf 1[\text{一致}]\,y+y\bigr]
=x-3^{\lvert A\rvert}\cdot3^{-\lvert A\rvert}\,y+y=x
```

となる(基底が一致する確率は $3^{-\lvert A\rvert}$ )。この計算に $y$ が何の期待値であるかは一度も使っていない。したがって、 $y$ はどんな数でもよい。

結果は「損はしないが得もしない」だった。理由は次のとおりである。

1. 平均場(Slater 行列式)の連結部分は、Wick の定理で電子の飛び移りの振幅の2乗(交換項)だけになり、絶縁体では距離とともに速く減衰する。
2. したがって遠い対では、-conn の値はほぼ $\langle(1-n_i)(1-n_j)\rangle\approx0$ (半充填で $n_i\approx1$ )になる。
3. prior の値が 0 なら $\Delta=x$ で、1 節の利得の式から $G=1$ になる(分子と分母が等しくなる)。

真の状態のゆっくり減衰する相関を予言するには、量子的な相関を持った prior が要る。

### 3.3 全距離で効く prior

- **切断 MPS**:厳密な一重項を、 $S^z$ の量子数ごとに上下対称に切り詰めるので、 $\langle S^z_i\rangle$ はほぼ 0 のままで、偽の秩序を持たない([README_phf_mps.md](README_phf_mps.md) の Schmidt 分解の節)。 $\chi\ge8$ で全距離で損がない。
- **変分 MPS の回転平均**: $\chi$ が大きいほど変分 MPS 自体の磁化 $\lvert\langle S^z_i\rangle\rvert$ が小さく、偽の秩序 $\langle S^z_i\rangle\langle S^z_j\rangle$ も小さい。 $\chi=16$ の強結合では対称な解に落ちるので、全距離で損がない。周期境界では $\chi=4$ の回転平均でも損がなかったが、開放端の $\chi=4$ は遠い対で損をする(開放端の変分 MPS の方が磁化が大きいためと思われるが、確かめていない)。

---

## 4. 実務上の意味

1. 遠い相関(相関関数の距離依存、構造因子など)を CRM で測るときは、**対称性を自発的に破った prior を使ってはいけない**。回転平均でも直らない。
2. 平均場しか使えない場合は、連結相関の prior を使えば損はしない。ただし得もしないので、近い対( $r\le3$ )だけ UHF-sym、遠い対は連結相関(または prior なし)と使い分けるのがよい。CRM の prior は観測量ごとに選べるので、これは追加の測定なしでできる。
3. 対称性を保つ量子的な prior(切断 MPS、 $\chi$ の十分大きい変分 MPS)なら全距離で得をする。

---

## 5. 限界

- 1次元の 2 つの系( $L=64$ 開放端、 $L=16$ 周期境界)だけを調べた。2次元では真の状態自体が長距離秩序に近づくので、結論が変わる可能性がある。
- 観測量は $Z_\uparrow Z_\uparrow$ だけである。
