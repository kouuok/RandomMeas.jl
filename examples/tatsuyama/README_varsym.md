# 変分 MPS をスピンの向きについて平均した prior

[README.md](README.md) に戻る / 前の段階: [README_phf_mps.md](README_phf_mps.md)、[README_uhfsym.md](README_uhfsym.md) / 続き: [README_varsym_distance.md](README_varsym_distance.md)

スクリプト:
- [crm_mps_lowchi_sym.jl](crm_mps_lowchi_sym.jl)(clara で実行。保存済みの変分 MPS 320 個から、回転平均に要る量を求める)
- [crm_varsym_compare.py](crm_varsym_compare.py)(clara で実行。利得の表を作る)

出力: [crm_mpslow_sym.tsv](crm_mpslow_sym.tsv)、[crm_varsym_compare.out](crm_varsym_compare.out)

---

## この文書で分かること

[README_phf_mps.md](README_phf_mps.md) の要点4で、変分 MPS(結合次元を $\chi$ に制限した DMRG の基底状態)には次の性質があることが分かっていた。

- 電荷の量(オンサイト ZZ など)では最良の prior である。
- ただし UHF と同じくスピンの向きを1つ選んでしまうので、スピンの量では損をする。

その対策として「UHF-sym と同じ回転平均が使えるはず」と書いたまま、試していなかった。この文書でそれを試した。結果は次のとおりである。

- **回転平均した変分 MPS( $\chi\ge8$ )は、厳密な参照状態を持つ 80 系 × 4 観測量のすべてで損をしない。** 平均しない変分 MPS( $\chi=8$ )は、単一サイト $Z_\uparrow$ で 80 系中 75 系、二重占有(4項の和)で 9 系、隣接 $Z_\uparrow Z_\uparrow$ で 5 系が損をしていた。
- **電荷の量とスピンの量の両方で、UHF-sym より良いか同じである。** $U=12$ の中央値( $n_m=100$ )では、オンサイト ZZ が 555(天井、UHF-sym は 465)、隣接 $Z_\uparrow Z_\uparrow$ が 19.8( $\chi=8$ )と 26.7( $\chi=16$ )(UHF-sym は 6.95、天井は 31.6)、二重占有が 129(UHF-sym は 123)である。
- 隣接スピン相関で UHF-sym に勝つのは、変分 MPS の $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ が UHF のそれより正確だからである。回転平均はスピンの向きの偏りを取り除くが、 $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ は変えない。
- 必要なのは変分 MPS を1回作ることと、期待値をいくつか計算することだけで、厳密な状態は要らない。

## 前提として知っておくこと

| 用語・記号 | 意味 |
|---|---|
| $n_{i\sigma}$ 、 $n_i$ | サイト $i$ 、スピン $\sigma$ の電子数(0 か 1)と、サイトの電子数 $n_{i\uparrow}+n_{i\downarrow}$ |
| $\mathbf S_i$ | サイト $i$ のスピン。 $S^z_i=\tfrac12(n_{i\uparrow}-n_{i\downarrow})$ |
| $Z_{i\uparrow}$ | 上向きスピンの軌道の量子ビットのパウリ $Z$ 。 $Z_{i\uparrow}=1-2n_{i\uparrow}$ |
| オンサイト ZZ | $Z_{i\uparrow}Z_{i\downarrow}$ 。半充填では二重占有率 $d$ と $\langle Z_{i\uparrow}Z_{i\downarrow}\rangle=4d-1$ で対応する |
| 利得 $G$ | 標準の古典シャドウの分散 ÷ CRM の分散( $n_m=100$ )。単一のパウリ文字列では $G=\dfrac{(3^{\lvert A\rvert}-1)x^2+v_s}{(3^{\lvert A\rvert}-1)\Delta^2+v_s}$ 、 $v_s=3^{\lvert A\rvert}(1-x^2)/n_m$ |
| 天井 | prior が完全( $\Delta=0$ )なときの利得 |

---

## 1. 回転平均を観測量の側で計算する

### 1.1 回転平均とは

全スピンを回転 $R$ で一様に回し、すべての回転について平均した混合状態

```math
\bar\sigma=\int dR\;R\,\sigma\,R^\dagger
```

を prior にする。CRM に必要なのは prior での期待値だけなので、混合状態でも構わない。期待値は $\mathrm{Tr}[O\bar\sigma]=\int dR\,\mathrm{Tr}[(R^\dagger OR)\sigma]$ と書けるので、**観測量の側を回して平均したもの**を、元の状態 $\sigma$ で測ればよい。

### 1.2 何が変わり、何が変わらないか

- 電子数 $n_i$ と、その積 $n_in_j$ はスピン回転で変わらない。
- 二重占有 $n_{i\uparrow}n_{i\downarrow}$ も変わらない。同じサイトに上下2つの電子がいる状態は、スピンが打ち消し合った一重項だからである。
- スピン $\mathbf S_i$ はベクトルとして回る。回転した後の $S^z_i$ は、ある単位ベクトル $\hat{\mathbf u}$ の方向の成分 $\hat{\mathbf u}\cdot\mathbf S_i$ になる。すべての回転について平均することは、 $\hat{\mathbf u}$ をすべての向きについて一様に平均することに当たる。

$\hat{\mathbf u}$ の向きの平均は

```math
\overline{\hat u_a}=0,\qquad\overline{\hat u_a\hat u_b}=\tfrac13\,\delta_{ab}
```

である。1つ目は向きが対称に散らばるから、2つ目は等方的で $\sum_a\hat u_a^2=1$ なので各成分の2乗の平均が $1/3$ だからである。したがって

```math
\overline{\langle S^z_i\rangle}=\sum_a\overline{\hat u_a}\,\langle S^a_i\rangle=0,\qquad
\overline{\langle S^z_iS^z_j\rangle}=\sum_{a,b}\overline{\hat u_a\hat u_b}\,\langle S^a_iS^b_j\rangle=\tfrac13\langle\mathbf S_i\cdot\mathbf S_j\rangle
```

となる。 $\mathbf S_i\cdot\mathbf S_j$ そのものは回転で変わらない。

### 1.3 観測量の値

$n_{i\uparrow}=\tfrac12(n_i+2S^z_i)$ ( $n_i+2S^z_i=2n_{i\uparrow}$ から)なので、

```math
Z_{i\uparrow}=1-2n_{i\uparrow}=1-n_i-2S^z_i
```

である。これに 1.2 節を当てはめる。

- **単一サイト $Z_{i\uparrow}$**: $S^z_i$ の平均が 0 なので $\langle Z_{i\uparrow}\rangle_{\bar\sigma}=1-\langle n_i\rangle$ 。
- **隣接 $Z_{i\uparrow}Z_{j\uparrow}$**: 展開すると $(1-n_i)(1-n_j)-2(1-n_i)S^z_j-2S^z_i(1-n_j)+4S^z_iS^z_j$ で、スピンの1次の項は平均で消え、2次の項は $\tfrac13\mathbf S_i\cdot\mathbf S_j$ になる:

```math
\langle Z_{i\uparrow}Z_{j\uparrow}\rangle_{\bar\sigma}=\bigl\langle(1-n_i)(1-n_j)\bigr\rangle+\tfrac43\langle\mathbf S_i\cdot\mathbf S_j\rangle
```

- **オンサイト ZZ**: $Z_{i\uparrow}Z_{i\downarrow}=(1-2n_{i\uparrow})(1-2n_{i\downarrow})=1-2n_i+4n_{i\uparrow}n_{i\downarrow}$ は電子数と二重占有だけで書けるので、回転で変わらない。

いずれも、平均する前の変分 MPS の $\langle n_i\rangle$ 、 $\langle n_in_j\rangle$ 、 $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ 、 $\langle n_{i\uparrow}n_{i\downarrow}\rangle$ から求まる。UHF-sym([README_uhfsym.md](README_uhfsym.md))と同じ式である。

### 1.4 計算

[crm_mps_lowchi_sym.jl](crm_mps_lowchi_sym.jl) は、clara に保存済みの変分 MPS(20 通りの格子 × 4 つの $U$ × $\chi=2,4,8,16$ の 320 個)について、中央のサイトとその隣のこれらの量を MPO で求める。 $\mathbf S_i\cdot\mathbf S_j=S^z_iS^z_j+\tfrac12(S^+_iS^-_j+S^-_iS^+_j)$ の $S^\pm_i$ ( $S^+_i=c^\dagger_{i\uparrow}c_{i\downarrow}$ など)は同じサイトの中で電子を移す演算子なので、Jordan–Wigner の符号は付かない。

**照合**: 同じ出力に入れた、平均する前のオンサイト ZZ・隣接 $Z_\uparrow Z_\uparrow$ ・ $S^z$ ・ $\mathbf S\cdot\mathbf S$ が、以前の計算(`crm_mpslow_*.tsv`)と最大 $10^{-12}$ で一致した。

---

## 2. 結果

データは厳密な参照状態を持つ 80 系(統合表 `crm_edfid_all.tsv` の 20 通りの格子 × 4 つの $U$ 。各 $U$ に 20 系)の、中央のサイトとその隣である。値は利得の中央値( $n_m=100$ )、括弧は 1% より大きく損をする( $G<0.99$ の)系の数で、0 のときは省いた。二重占有は4項の和 $\tfrac14(1-Z_\uparrow-Z_\downarrow+Z_\uparrow Z_\downarrow)$ として測るときの厳密な分散比である。

### 2.1 $U=12$

| 観測量 | UHF | UHF-sym | 変分 $\chi=4$ | 変分 $\chi=4$ -sym | 変分 $\chi=8$ | 変分 $\chi=8$ -sym | 変分 $\chi=16$ | 変分 $\chi=16$ -sym |
|---|---|---|---|---|---|---|---|---|
| オンサイト ZZ | 465 | 465 | 553 | 553 | 555 | 555 | 555 | 555 |
| 隣接 $Z_\uparrow Z_\uparrow$ | 1.22 (7) | 6.95 | 15.3 (1) | 18.0 (1) | 16.8 (1) | **19.8** | 28.8 (1) | **26.7** |
| 単一サイト $Z_\uparrow$ | 0.02 (20) | 1.00 | 0.04 (20) | 1.00 | 0.11 (18) | 1.00 | 0.46 (11) | 1.00 |
| 二重占有 | 1.81 | 123 | 4.75 | 129 | 13.3 | **129** | 53.2 | **129** |

**表の読み方**:

- 天井はオンサイト ZZ が 555、隣接 $Z_\uparrow Z_\uparrow$ が 31.6 である。
- **オンサイト ZZ は回転平均で変わらない**(1.3 節)。変分 MPS は $\chi=4$ からほぼ天井に届いている。
- **単一サイト $Z_\uparrow$ は、どの prior でも 1 を超えない。** 真の値は $x=0$ である(有限系の基底状態は一重項で $\langle S^z_i\rangle=0$ 、半充填で $\langle n_i\rangle=1$ )。利得の式に $x=0$ を入れると $G=v_s/((3^{\lvert A\rvert}-1)\Delta^2+v_s)\le1$ で、 $\Delta=0$ のときだけ 1 になる。回転平均した prior の値は $1-\langle n_i\rangle\approx0$ なので $G\approx1$ (損をしない)になり、平均しない prior は $\pm2\langle S^z_i\rangle$ だけ外れて損をする。
- **隣接 $Z_\uparrow Z_\uparrow$ では、回転平均で損が消える。** $\chi=16$ では平均しない方の中央値(28.8)が平均した方(26.7)より少し大きいが、平均しない方は 1 系で損をしている。 $U=12$ の $\chi=16$ の変分 MPS は、1次元開放端では対称な解に落ちる(平均しても値が変わらない)が、周期境界と2次元では対称性を破ったまま残る([README_phf_mps.md](README_phf_mps.md))。破れた系で、平均しない方の中央値が少し大きくなる理由は調べていない。
- **二重占有は、回転平均で大きく改善する**(4.75〜53.2 → 129)。2.4 節で理由を述べる。

### 2.2 $U=4$

| 観測量 | UHF | UHF-sym | 変分 $\chi=8$ | 変分 $\chi=8$ -sym | 変分 $\chi=16$ | 変分 $\chi=16$ -sym |
|---|---|---|---|---|---|---|
| オンサイト ZZ | 50.3 | 50.3 | 50.1 | 50.1 | 50.6 | 50.6 |
| 隣接 $Z_\uparrow Z_\uparrow$ | 1.89 (5) | 11.2 | 9.71 (1) | 11.7 | 20.8 | 18.9 |
| 単一サイト $Z_\uparrow$ | 0.02 (20) | 1.00 (1) | 0.19 (19) | 1.00 (1) | 0.12 (15) | 1.00 (1) |
| 二重占有 | 1.20 (3) | 28.0 | 8.67 (3) | 28.0 | 5.34 (2) | 28.1 |

**表の読み方**: $U=4$ では UHF-sym と変分 $\chi=8$ -sym がほぼ並ぶ。 $\chi=16$ -sym になると、隣接 $Z_\uparrow Z_\uparrow$ で UHF-sym を上回る(18.9 対 11.2)。

### 2.3 損の数(全 $U$ 、80 系)

| prior | オンサイト ZZ | 隣接 $Z_\uparrow Z_\uparrow$ | 単一サイト $Z_\uparrow$ | 二重占有 |
|---|---|---|---|---|
| UHF | 0 | 17 | 77 | 13 |
| UHF-sym | 0 | 0 | 2 | 0 |
| 変分 $\chi=4$ | 1 | 3 | 80 | 6 |
| 変分 $\chi=4$ -sym | 1 | 1 | 7 | 1 |
| 変分 $\chi=8$ | 0 | 5 | 75 | 9 |
| **変分 $\chi=8$ -sym** | **0** | **0** | **2** | **0** |
| 変分 $\chi=16$ | 0 | 2 | 53 | 4 |
| **変分 $\chi=16$ -sym** | **0** | **0** | **2** | **0** |

**表の読み方**:

- 単一サイト $Z_\uparrow$ の 2 系は、UHF-sym でも損をしている。2.1 節のとおり、この観測量は $G\le1$ しかありえないので、prior の $1-\langle n_i\rangle$ がわずかにずれるだけで(例えば二部格子でない格子で $\langle n_i\rangle\ne1$ になると)損になる。
- $\chi=4$ -sym にはまだ損が残る。 $\chi=4$ の変分 MPS は局所解に落ちやすく、回転で変わらないオンサイト ZZ でも 1 系( $U=2$ )が損をしている([README_phf_mps.md](README_phf_mps.md) の変分 MPS の作り方の節)。回転平均では直せない種類の誤りである。

### 2.4 なぜ二重占有が大きく改善するのか

二重占有の4項のうち、 $Z_\uparrow$ と $Z_\downarrow$ の2項について、平均しない変分 MPS は $1-\langle n_i\rangle\mp2\langle S^z_i\rangle$ と、磁化の分だけ逆向きに外れる。合計では打ち消し合うが、CRM の分散には項ごとの外れの積が入る([README_experiments.md](README_experiments.md) §3 の和の観測量の厳密な式)ので、打ち消されない。回転平均で $\langle S^z_i\rangle$ の寄与が消えると、項ごとの外れも小さくなる。

### 2.5 なぜ隣接 $Z_\uparrow Z_\uparrow$ で UHF-sym に勝つのか

1.3 節のとおり、回転平均した prior の隣接 $Z_\uparrow Z_\uparrow$ は、 $\langle(1-n_i)(1-n_j)\rangle$ (電荷の量)と $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ (スピンの量)で決まる。2つの prior の差は主に後者の正確さから来る。出力の [C] 節では、相対誤差 $\lvert\text{prior}/\text{真}-1\rvert$ の中央値は次のとおりだった。

| $U$ | UHF(= UHF-sym) | 変分 $\chi=4$ | 変分 $\chi=8$ | 変分 $\chi=16$ |
|---|---|---|---|---|
| 2 | 0.245 | 0.190 | 0.187 | 0.073 |
| 4 | 0.224 | 0.212 | 0.131 | 0.088 |
| 8 | 0.327 | 0.211 | 0.123 | 0.040 |
| 12 | 0.349 | 0.209 | 0.107 | 0.049 |

- **UHF の $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ は、古典的な Néel(反平行の積状態)の値である。** 反平行の積状態 $\lvert\uparrow\downarrow\rangle$ の $\mathbf S_i\cdot\mathbf S_j$ は $-\tfrac14$ で、隣どうしが量子的な一重項を組んだ値 $-\tfrac34$ に比べて相関が弱い。強結合ほど真の状態は一重項に近づくので、誤差が大きくなる( $U=12$ で 35%)。これが UHF-sym の隣接 $Z_\uparrow Z_\uparrow$ が天井に届かない理由である([README_prior_survey.md](README_prior_survey.md) 結果2)。
- **変分 MPS は結合ごとに量子的な相関を持てる**ので、 $\chi$ を上げるほど $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ が正確になる。平均する前は、その正確さがスピンの向きの偏りに隠れていた。回転平均で偏りだけを取り除くと、その正確さがそのまま利得に出る。

---

## 3. 実務上の意味

1. **電荷とスピンの両方を測るなら、変分 MPS( $\chi\ge8$ )を回転平均したものが、この研究で試した prior のうち厳密な状態を使わないものでは最も良い。** 中央のサイトとその隣では損が1つもなく、隣接スピン相関で UHF-sym の 2〜3 倍得をする( $U\ge8$ )。
2. 回転平均のコストは、変分 MPS の期待値をいくつか計算するだけである。観測量ごとに次のものがあればよい。
   - $\langle n_i\rangle$ 、 $\langle n_in_j\rangle$ 、 $\langle n_{i\uparrow}n_{i\downarrow}\rangle$ (電荷の量)
   - $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ (スピンの量)
3. 回転平均は、真の状態がスピン回転の対称性を持つことを前提にしている(半充填の有限系の基底状態は全スピン 0 の一重項)。2次元の DMRG の参照状態のように、測る状態の側が対称性を破っているときは、項ごとの $\langle Z_q\rangle$ が合わなくなり、和の観測量で損をすることがある([README_experiments.md](README_experiments.md) §4 の二重占有の退行と同じ理由)。
4. **遠い距離のスピン相関では、回転平均でも直らない誤りが残る**( $\chi$ が小さいとき)。スピンの向きを選んだ状態は距離によらない $\langle\mathbf S_i\rangle\cdot\langle\mathbf S_j\rangle$ を持ち、これは回転で変わらないからである([README_varsym_distance.md](README_varsym_distance.md))。

---

## 4. 限界

- 中央のサイトとその隣の観測量だけを調べた。距離の遠い対は [README_varsym_distance.md](README_varsym_distance.md) で調べた。ブロックの電荷の偶奇 $(-1)^{N_B}$ は調べていないが、電子数だけで決まる量なので回転平均で変わらない。
- 系は厳密な参照を持つ 80 系(16 サイト以下の2次元と、1次元 $L\le64$ )だけである。
- 変分 MPS が局所解に落ちる問題([README_phf_mps.md](README_phf_mps.md))はそのまま残る。 $\chi=4$ ではその影響が見える。
