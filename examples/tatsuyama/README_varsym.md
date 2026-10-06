# 変分 MPS をスピンの向きについて平均した prior

[README.md](README.md) に戻る / 前の段階: [README_phf_mps.md](README_phf_mps.md)、[README_uhfsym.md](README_uhfsym.md)

スクリプト:
- [crm_mps_lowchi_sym.jl](crm_mps_lowchi_sym.jl)(clara で実行。保存済みの変分 MPS 320 個から回転平均に要る量を求める)
- [crm_varsym_compare.py](crm_varsym_compare.py)(clara で実行。利得の表)

出力: [crm_mpslow_sym.tsv](crm_mpslow_sym.tsv)、[crm_varsym_compare.out](crm_varsym_compare.out)

---

## 0. 要点

[README_phf_mps.md](README_phf_mps.md) の要点4で、変分 MPS(結合次元を $\chi$ に制限した DMRG の基底状態)には次の性質があった。

- 電荷の量では最良の prior である。
- ただし UHF と同じくスピンの向きを選んでしまう。

その対策として「UHF-sym と同じ回転平均が使えるはず」と書いて、試していなかった。この文書でそれを試した。

- **回転平均した変分 MPS( $\chi=8$ 以上)は、厳密な 80 系 × 4 観測量のすべてで損をしない。**
  - 平均しない変分 MPS( $\chi=8$ )は、単一サイト $Z_\uparrow$ で 80 系中 75 系、二重占有(4項の和)で 9 系、隣接 $Z_\uparrow Z_\uparrow$ で 5 系が損をしていた。
- **電荷の量とスピンの量の両方で、UHF-sym より良い。** $U=12$ の中央値( $n_m=100$ )では次のとおり。
  - オンサイト ZZ:555(天井)。UHF-sym は 465。
  - 隣接 $Z_\uparrow Z_\uparrow$ :19.8( $\chi=8$ )、26.7( $\chi=16$ )。UHF-sym は 6.95、天井は 31.6。
  - 二重占有:129。UHF-sym は 123。
- 隣接 $Z_\uparrow Z_\uparrow$ で UHF-sym に勝つ理由は、変分 MPS の $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ が UHF より正確だからである(相対誤差 $U=12$ で UHF 0.35 に対し変分 $\chi=8$ で 0.11、 $\chi=16$ で 0.05)。回転平均はスピンの向きの偏りを消すが、 $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ は変えない。したがって平均した後の質は $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ の正確さで決まる。
- 必要なのは変分 MPS の DMRG 1回と、期待値の計算だけである。厳密な状態は要らない(切断 MPS との違い)。
  - 切断 MPS $\chi=4$ の隣接 $Z_\uparrow Z_\uparrow$ は 29.2 で、 $\chi=16$-sym の 26.7 とほぼ同じである。ただし切断 MPS は厳密な状態がないと作れない。

---

## 1. 回転平均を観測量の側で計算する

全スピンを一様に回して平均した混合状態

```math
\bar\sigma=\int d\Omega\,R(\Omega)\,\sigma\,R(\Omega)^\dagger
```

での期待値は、観測量の側を回して平均したものに等しい。サイト $i$ のスピン $\mathbf S_i$ はベクトルとして回るので、方向について平均すると

```math
\overline{\langle S^z_i\rangle}=0,\qquad
\overline{\langle S^z_iS^z_j\rangle}=\tfrac13\langle\mathbf S_i\cdot\mathbf S_j\rangle
```

となり、 $n_i$ 、 $n_i n_j$ 、 $n_{i\uparrow}n_{i\downarrow}$ 、 $\mathbf S_i\cdot\mathbf S_j$ は変わらない。 $Z_{i\uparrow}=1-2n_{i\uparrow}=1-n_i-2S^z_i$ を使うと

```math
\langle Z_{i\uparrow}\rangle_{\bar\sigma}=1-\langle n_i\rangle,\qquad
\langle Z_{i\uparrow}Z_{j\uparrow}\rangle_{\bar\sigma}=\bigl\langle(1-n_i)(1-n_j)\bigr\rangle+\tfrac43\langle\mathbf S_i\cdot\mathbf S_j\rangle
```

で、オンサイト ZZ $=1-2n_i+4n_{i\uparrow}n_{i\downarrow}$ は回転で変わらない。いずれも元の(平均する前の)変分 MPS の $\langle n_i\rangle$ 、 $\langle n_in_j\rangle$ 、 $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ から求まる。UHF-sym([README_uhfsym.md](README_uhfsym.md))と同じ式である。

[crm_mps_lowchi_sym.jl](crm_mps_lowchi_sym.jl) は、clara に保存済みの変分 MPS(80 系 × 4 つの $U$ × $\chi=2,4,8,16$ の 320 個)について、中央のサイトとその隣のこれらの量を MPO で求める。 $S^\pm_i=c^\dagger_{i\uparrow}c_{i\downarrow}$ などは同じサイトの演算子なので、JW の符号は付かない。

**照合**: 同じ出力に入れたオンサイト ZZ・隣接 $Z_\uparrow Z_\uparrow$ ・ $S^z$ ・ $\mathbf S\cdot\mathbf S$ (平均する前の値)が、以前の計算(`crm_mpslow_*.tsv`)と最大 $10^{-12}$ で一致した。

---

## 2. 結果

厳密な 80 系(統合表 `crm_edfid_all.tsv` で厳密な参照状態を持つ系、各 $U$ に 20 系)の中央のサイトと隣である。値は利得の中央値、括弧は 1% より大きく損をする系の数( $n_m=100$ )。二重占有は4項の和として測るときの厳密な分散比である。

### 2.1 $U=12$

| 観測量 | UHF | UHF-sym | 変分 $\chi=4$ | 変分 $\chi=4$ -sym | 変分 $\chi=8$ | 変分 $\chi=8$ -sym | 変分 $\chi=16$ | 変分 $\chi=16$ -sym |
|---|---|---|---|---|---|---|---|---|
| オンサイト ZZ | 465 | 465 | 553 | 553 | 555 | 555 | 555 | 555 |
| 隣接 $Z_\uparrow Z_\uparrow$ | 1.22 (7) | 6.95 | 15.3 (1) | 18.0 (1) | 16.8 (1) | **19.8** | 28.8 (1) | **26.7** |
| 単一サイト $Z_\uparrow$ | 0.02 (20) | 1.00 | 0.04 (20) | 1.00 | 0.11 (18) | 1.00 | 0.46 (11) | 1.00 |
| 二重占有 | 1.81 | 123 | 4.75 | 129 | 13.3 | **129** | 53.2 | **129** |

天井は、オンサイト ZZ が 555、隣接 $Z_\uparrow Z_\uparrow$ が 31.6 である。単一サイト $Z_\uparrow$ は真の値が 0 なので、どの prior でも 1 を超えない。

### 2.2 $U=4$

| 観測量 | UHF | UHF-sym | 変分 $\chi=8$ | 変分 $\chi=8$ -sym | 変分 $\chi=16$ | 変分 $\chi=16$ -sym |
|---|---|---|---|---|---|---|
| オンサイト ZZ | 50.3 | 50.3 | 50.1 | 50.1 | 50.6 | 50.6 |
| 隣接 $Z_\uparrow Z_\uparrow$ | 1.89 (5) | 11.2 | 9.71 (1) | 11.7 | 20.8 | 18.9 |
| 単一サイト $Z_\uparrow$ | 0.02 (20) | 1.00 (1) | 0.19 (19) | 1.00 (1) | 0.12 (15) | 1.00 (1) |
| 二重占有 | 1.20 (3) | 28.0 | 8.67 (3) | 28.0 | 5.34 (2) | 28.1 |

$U=4$ では UHF-sym と変分 $\chi=8$ -sym がほぼ並ぶ。 $\chi=16$ -sym で隣接 $Z_\uparrow Z_\uparrow$ が UHF-sym を上回る。

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

単一サイト $Z_\uparrow$ は、UHF-sym でも 2 系が損をしている。真の値がほぼ 0 で $G\le1$ しかありえない観測量なので、prior の $1-\langle n_i\rangle$ がわずかにずれるだけで損になる。

$\chi=4$ -sym にはまだ損が残る。 $\chi=4$ の変分 MPS は局所解に落ちやすく、オンサイト ZZ でも1系( $U=2$ )が損をしている([README_phf_mps.md](README_phf_mps.md) の変分 MPS の作り方の節)。

### 2.4 なぜ隣接 $Z_\uparrow Z_\uparrow$ で UHF-sym に勝つのか

回転平均した prior の隣接 $Z_\uparrow Z_\uparrow$ は、 $\langle(1-n_i)(1-n_j)\rangle$ と $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ で決まる。前者は電荷の量で、2つの prior の差は主に後者の $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ の正確さから来る。出力の [C] 節では、相対誤差 $\lvert\text{prior}/\text{真}-1\rvert$ の中央値は次のとおりだった。

| $U$ | UHF(= UHF-sym) | 変分 $\chi=4$ | 変分 $\chi=8$ | 変分 $\chi=16$ |
|---|---|---|---|---|
| 2 | 0.245 | 0.190 | 0.187 | 0.073 |
| 4 | 0.224 | 0.212 | 0.131 | 0.088 |
| 8 | 0.327 | 0.211 | 0.123 | 0.040 |
| 12 | 0.349 | 0.209 | 0.107 | 0.049 |

- UHF の $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ は古典的な Néel(反平行の積状態)の値で、隣どうしの一重項の量子的な相関を持たない。強結合ほど誤差が大きい( $U=12$ で 35%)。これが UHF-sym の隣接 $Z_\uparrow Z_\uparrow$ が天井に届かない理由である([README_prior_survey.md](README_prior_survey.md) 結果2)。
- 変分 MPS は結合ごとに量子的な相関を持つので、 $\chi$ を上げるほど $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ が正確になる。平均する前は、その正確さがスピンの向きの偏りに隠れていた。回転平均で偏りだけを取り除くと、その正確さがそのまま利得に出る。

---

## 3. 実務上の意味

1. **電荷とスピンの両方を測るなら、変分 MPS( $\chi\ge8$ )を回転平均したものが、この研究で試した prior のうち厳密な状態を使わないものでは最も良い。** 損は1つもなく、UHF-sym より隣接スピン相関で 2〜3 倍得をする( $U\ge8$ )。
2. 回転平均のコストは、変分 MPS の期待値をいくつか計算するだけである。観測量ごとに次のものがあればよい。
   - $\langle n_i\rangle$ と $\langle n_in_j\rangle$ (電荷の量)
   - $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ (スピンの量)
3. 回転平均は、真の状態がスピン回転の対称性を持つことを前提にしている(半充填の有限系の基底状態は全スピン 0 の一重項)。ただし、2次元の DMRG の参照状態のように実験の状態の側が対称性を破っているときは、項ごとの $\langle Z_q\rangle$ が合わなくなり、和の観測量で損をすることがある([README_experiments.md](README_experiments.md) §4 の二重占有の退行と同じ理由)。

---

## 4. 限界

- 中央のサイトと隣の観測量だけを調べた。距離の遠い対は [README_varsym_distance.md](README_varsym_distance.md) で調べた(回転平均では偽の長距離秩序が残るので、 $\chi$ が小さいと遠い対で損をする)。ブロックの電荷の偶奇は調べていない。ブロックの偶奇は電荷の量なので、回転平均で変わらない。
- 系は厳密な参照を持つ 80 系(16 サイト以下の2次元と、1次元 $L\le64$ )だけである。
- 変分 MPS が局所解に落ちる問題([README_phf_mps.md](README_phf_mps.md))はそのまま残る。 $\chi=4$ ではその影響が見える。
