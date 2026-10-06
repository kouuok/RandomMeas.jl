# 回転平均した prior は遠いスピン相関でも効くか

[README.md](README.md) に戻る / 前の段階: [README_varsym.md](README_varsym.md)(中央のサイトと隣だけを調べた)

スクリプト:
- [crm_varsym_pairs.jl](crm_varsym_pairs.jl)(clara で実行。保存済みの変分 MPS から全サイト対の量を求める)
- [crm_varsym_distance.py](crm_varsym_distance.py)(clara で実行。距離ごとの利得)

出力: [crm_varsym_distance.out](crm_varsym_distance.out)

---

## 0. 要点

- **スピンの向きを選んだ prior は、偽の長距離秩序を持つ。** $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ に $\langle\mathbf S_i\rangle\cdot\langle\mathbf S_j\rangle=\pm m^2$ が距離によらず含まれる。 $\mathbf S_i\cdot\mathbf S_j$ は回転で変わらないので、**回転平均してもこの部分は残る**。
- その結果、UHF-sym は隣( $r=1$ )では得をするが、 $r\ge4$ の対ではすべて損をする( $U\ge4$ 、1次元 $L=64$ 開放端と $L=16$ 周期境界の両方)。
- **偽の秩序を引いた「連結相関の prior」は損をしないが、遠い対では得もしない**(利得 1.00)。CRM はどんな数値を prior に使っても偏らないので、物理的な状態の期待値でない prior も使える。ただし平均場の連結相関は数サイトで消えてしまい、真の状態のゆっくり減衰する相関を予言できない。
- **全距離で得をするのは、対称性を保つ prior である。** 切断 MPS( $\chi\ge8$ )と、回転平均した変分 MPS( $\chi=16$ 、および周期境界の $\chi\ge4$ )は、 $r=1$ から最も遠い対まで損がほとんどない。

---

## 1. 観測量と prior

観測量は距離 $r$ の $Z_{i\uparrow}Z_{j\uparrow}$ (台は 2 量子ビット、 $n_m=100$ )。1次元 $L=64$ 開放端では端の影響を避けて $17\le i<j\le48$ の対、 $L=16$ 周期境界では全対を使い、距離ごとに利得の中央値と損をする対の割合を数えた。真の値は厳密な基底状態(`crm_eigcorr_*.tsv`)である。

$Z_{i\uparrow}=1-n_i-2S^z_i$ から、比べる prior の値は次のとおりである。

| prior | $\langle Z_{i\uparrow}Z_{j\uparrow}\rangle$ の値 |
|---|---|
| UHF、変分 MPS | 状態の期待値そのもの |
| -sym(回転平均) | $\langle(1-n_i)(1-n_j)\rangle+\tfrac43\langle\mathbf S_i\cdot\mathbf S_j\rangle$ |
| -conn(連結相関) | $\langle(1-n_i)(1-n_j)\rangle+\tfrac43\bigl(\langle\mathbf S_i\cdot\mathbf S_j\rangle-\langle\mathbf S_i\rangle\cdot\langle\mathbf S_j\rangle\bigr)$ |
| 切断 MPS | 状態の期待値そのもの(厳密な状態が要る) |

変分 MPS は $S^z$ の総和を保つので $\langle S^x_i\rangle=\langle S^y_i\rangle=0$ で、 $\langle\mathbf S_i\rangle\cdot\langle\mathbf S_j\rangle=\langle S^z_i\rangle\langle S^z_j\rangle$ である。

**照合**: [crm_varsym_pairs.jl](crm_varsym_pairs.jl) が出した変分 MPS の $Z_\uparrow Z_\uparrow$ は、以前の計算(`crm_mpslow_*.tsv`)の全対 34176 個と完全に一致した。

---

## 2. 結果

### 2.1 1次元 $L=64$ 開放端、 $U=12$

利得の中央値 / 1% より大きく損をする対の割合。

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

$U=12$ の $\chi=16$ では、変分 MPS はもともと対称な解に落ちていて、回転平均しても値は変わらない。

### 2.2 $r\ge2$ の全対での損の割合

1次元 $L=64$ 開放端 / $L=16$ 周期境界。

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

---

## 3. 解釈

### 3.1 回転平均で残る偽の秩序

UHF は副格子ごとに磁化 $\pm m$ を持つ。その $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ は

```math
\langle\mathbf S_i\cdot\mathbf S_j\rangle_{\rm UHF}=\underbrace{\langle\mathbf S_i\rangle\cdot\langle\mathbf S_j\rangle}_{=\pm m^2\ (\text{距離によらない})}+\bigl(\text{連結部分}\bigr)
```

で、第1項が「どこまで行っても揃ったまま」の古典的な Néel 秩序である。有限の系の真の基底状態は一重項で $\langle\mathbf S_i\rangle=0$ 、1次元では相関が距離とともに減衰する。回転平均は $\langle\mathbf S_i\rangle$ の向きを平均するが、 $\mathbf S_i\cdot\mathbf S_j$ は回転で変わらないので第1項は消えない。したがって UHF-sym は遠い対で $\frac43m^2$ だけ外れ続け、真の値が小さくなる $r\ge4$ で損得の境目を越える。

[README_eigenstate.md](README_eigenstate.md) の「長距離秩序を持つ prior は遠くの相関で必ず損をする」は、回転平均しても変わらない。対称性の破れの「向き」は平均で消えるが、「秩序の大きさ」は消えない。

### 3.2 連結相関の prior

第1項を引いた値は、もはやどの状態の期待値でもない(どの状態でも $\langle\mathbf S_i\cdot\mathbf S_j\rangle$ と $\langle\mathbf S_i\rangle\cdot\langle\mathbf S_j\rangle$ は別々に決まる)。それでも CRM の推定は偏らない。CRM に必要なのは各設定での prior の条件付き期待値 $3^{\lvert A\rvert}\mathbf 1[\text{一致}]\,y$ と足し戻す値 $y$ だけで、 $y$ が何であっても期待値が打ち消し合うからである。

結果は「損はしないが得もしない」だった。平均場の連結相関は短い距離で急に減衰し、prior の値はほぼ $\langle(1-n_i)(1-n_j)\rangle\approx0$ になる。外れ $\Delta$ がほぼ真の値そのものになり、利得は 1 に近づく。真の状態の相関はべき的にゆっくり減衰するので、それを予言するには量子的な相関を持った prior が要る。

### 3.3 全距離で効く prior

- **切断 MPS**:厳密な一重項を切り詰めるので、対称性を保ち、偽の秩序を持たない。 $\chi\ge8$ で全距離で損がない。
- **変分 MPS の回転平均**: $\chi$ が大きいほど変分 MPS 自体の磁化が小さく、偽の秩序も小さい。 $\chi=16$ では強結合で対称な解に落ちるので、全距離で損がない。周期境界では $\chi=4$ の回転平均でも損がなかったが、開放端の $\chi=4$ は遠い対で損をする(開放端の変分 MPS の方が磁化が大きいためと思われるが、確かめていない)。

---

## 4. 実務上の意味

1. 遠い相関(相関関数の距離依存、構造因子など)を CRM で測るときは、**対称性を自発的に破った prior を使ってはいけない**。回転平均でも直らない。
2. 平均場しか使えない場合は、連結相関の prior を使えば損はしない。ただし得もしないので、近い対( $r\le3$ )だけ UHF-sym、遠い対は連結相関(または prior なし)と使い分けるのがよい。CRM の prior は観測量ごとに選べるので、これは追加の測定なしでできる。
3. 対称性を保つ量子的な prior(切断 MPS、 $\chi$ の十分大きい変分 MPS)なら全距離で得をする。

---

## 5. 限界

- 1次元の 2 つの系( $L=64$ 開放端、 $L=16$ 周期境界)だけである。2次元では真の状態自体が長距離秩序に近づくので、結論が変わる可能性がある。
- 観測量は $Z_\uparrow Z_\uparrow$ だけである。
