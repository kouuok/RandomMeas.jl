"""
読み出し誤差に対する CRM の利得 (README_readout_noise.md の数値の出典)

誤差モデル: 測定した各量子ビットの結果が、独立に確率 p で反転する(対称なビット反転)。
台 |A| のパウリ文字列 P の1ショットの値 (±1) の平均は f x (f = (1-2p)^|A|) になる。

比べる推定量(どれも同じ量 f x を推定するので、偏りは同じ。最後に f で割る読み出し補正も共通):
  標準        : 3^|A| m                        (m: 一致した基底でのショット平均)
  CRM (そのまま): 3^|A| (m - y) + y'           prior の値 y を誤差なしで引く。y' は足し戻す値
  CRM (誤差込み): 3^|A| (m - f y) + f y         prior にも同じ誤差モデルを掛ける

分散はどれも元論文の式(25) [PRX: (C7)] の形で、x -> f x、Δ -> f x - y (そのまま) または f (x - y) (誤差込み) と
置き換えたものになる。[0] でビット反転を直接まねたモンテカルロと照合する。

データ: 1次元 L=64 開放端の DMRG 基底状態(厳密な参照) と prior の値
  crm_eigcorr_W1L64_cylinder_U*.tsv  (ZupZdn = オンサイト ZZ、ZupZup = 隣接 Z↑Z↑)
  crm_eigparity_W1L64_cylinder_U*.tsv (P = 電荷の偶奇 (-1)^{N_B}、ブロック ℓ サイト = 台 2ℓ)
"""
import csv
import numpy as np

rng = np.random.default_rng(20261006)
PS = (0.0, 0.001, 0.002, 0.005, 0.01, 0.02)


def gain(x, y, nA, nm):
    K = 3.0**nA - 1
    vs = 3.0**nA * (1 - x * x) / nm
    return (K * x * x + vs) / (K * (x - y) ** 2 + vs)


def g_unaware(x, y, nA, nm, p):
    f = (1 - 2 * p) ** nA
    return gain(f * x, y, nA, nm)


def g_aware(x, y, nA, nm, p):
    f = (1 - 2 * p) ** nA
    return gain(f * x, f * y, nA, nm)


# ------------------------------------------------------------
# [0] ビット反転を直接まねたモンテカルロでの確認
# ------------------------------------------------------------
def mc_check(x, y, nA, nm, p, NU=4_000_000):
    # 基底が一致しないと推定量は 0 (標準) か定数なので、一致した基底のショットだけをまねる
    nmatch = rng.binomial(NU, 3.0**-nA)
    # 真の固有値の積 (±1、平均 x) を各ショットで引き、台の各量子ビットを確率 p で反転させる
    s = np.where(rng.random((nmatch, nm)) < (1 + x) / 2, 1.0, -1.0)
    flips = rng.binomial(nA, p, size=(nmatch, nm))
    s *= np.where(flips % 2 == 1, -1.0, 1.0)
    m = s.mean(axis=1)
    f = (1 - 2 * p) ** nA
    z = np.zeros(NU - nmatch)
    std = np.concatenate([3.0**nA * m, z])
    un = np.concatenate([3.0**nA * (m - y), z]) + y
    aw = np.concatenate([3.0**nA * (m - f * y), z]) + f * y
    return std.mean(), un.mean(), aw.mean(), std.var() / un.var(), std.var() / aw.var()


print("[0] 利得の式(x -> f x)とビット反転のモンテカルロの照合 (n_u = 4000000 基底)")
print("    x       y     |A| n_m   p     | 平均(標準, そのまま, 誤差込み) 真の f x | G そのまま MC / 式 | G 誤差込み MC / 式")
for (x, y, nA, nm, p) in [(-0.93, -0.945, 2, 100, 0.01), (0.6, 0.55, 4, 10, 0.02),
                          (-0.9, -0.92, 4, 20, 0.02), (0.93, 0.94, 6, 100, 0.005), (0.5, 0.45, 3, 1, 0.05)]:
    ms, mu, ma, gu, ga = mc_check(x, y, nA, nm, p)
    f = (1 - 2 * p) ** nA
    print(f"  {x:+.3f} {y:+.3f}  {nA:2d} {nm:4d} {p:.3f} | {ms:+.4f} {mu:+.4f} {ma:+.4f}  {f*x:+.4f} |"
          f" {gu:8.2f} / {g_unaware(x,y,nA,nm,p):8.2f} | {ga:8.2f} / {g_aware(x,y,nA,nm,p):8.2f}")


# ------------------------------------------------------------
# データの読み込み
# ------------------------------------------------------------
def load(U):
    d = {}
    for r in csv.DictReader(open(f"crm_eigcorr_W1L64_cylinder_U{U}.tsv"), delimiter="\t"):
        i, j = int(r["i"]), int(r["j"])
        if r["kind"] == "ZupZdn" and i == 32 and j == 32:
            d[("ZZ onsite", r["prior"])] = (float(r["true"]), float(r["prior_val"]), 2)
        if r["kind"] == "ZupZup" and i == 32 and j == 33:
            d[("隣接 Z↑Z↑", r["prior"])] = (float(r["true"]), float(r["prior_val"]), 2)
    for r in csv.DictReader(open(f"crm_eigparity_W1L64_cylinder_U{U}.tsv"), delimiter="\t"):
        if r["kind"] == "P" and r["start"] == "17" and r["l"] in ("4", "8", "16"):
            l = int(r["l"])
            d[(f"P_{l} (台{2*l})", r["prior"])] = (float(r["true"]), float(r["prior_val"]), 2 * l)
    return d


OBS = ["ZZ onsite", "隣接 Z↑Z↑", "P_4 (台8)", "P_8 (台16)", "P_16 (台32)"]
PRIORS = ["UHF", "chi8", "chi32"]

for U in ("12.0", "4.0"):
    d = load(U)
    for nm in (100, 1):
        print(f"\n[1] L=64 開放端, U={U}, n_m={nm}: 利得 G (そのまま / 誤差込み)")
        print("    観測量           prior  x        y/x   | " + "  ".join(f"p={p:<5}" + " " * 9 for p in PS))
        for name in OBS:
            for pr in PRIORS:
                if (name, pr) not in d:
                    continue
                x, y, nA = d[(name, pr)]
                cells = "  ".join(f"{g_unaware(x,y,nA,nm,p):6.1f}/{g_aware(x,y,nA,nm,p):6.1f}  " for p in PS)
                print(f"    {name:14s} {pr:6s} {x:+.4f} {y/x:5.2f} | {cells}")

# ------------------------------------------------------------
# [2] 損得の境目: そのままの prior で G<1 になる p
#    G<1 ⇔ |f x - y| > |f x| ⇔ (y≈x のとき) f < y/(2x) ⇔ (1-2p)^|A| < y/(2x)
# ------------------------------------------------------------
print("\n[2] そのままの prior で損(G<1)に転じる p (n_m=100 での数値解 / 式 (1-2p)^|A| = y/(2x) の解)")
d = load("12.0")
for name in OBS:
    x, y, nA = d[(name, "UHF")]
    ps = np.linspace(0, 0.2, 200001)
    g = np.array([g_unaware(x, y, nA, 100, p) for p in ps[::100]])
    idx = np.argmax(g < 1)
    p_num = ps[::100][idx] if g[idx] < 1 else float("nan")
    r = y / (2 * x)
    p_app = (1 - r ** (1 / nA)) / 2 if 0 < r < 1 else float("nan")
    print(f"    {name:14s} UHF: |A|={nA:2d}  y/x={y/x:.3f}  p* 数値 {p_num:.4f}  式 {p_app:.4f}")

# ------------------------------------------------------------
# [3] 誤差込みの prior の天井 (Δ=0): 1 + (3^|A|-1) f^2 x^2 n_m / (3^|A| (1 - f^2 x^2))
# ------------------------------------------------------------
print("\n[3] 誤差込みの prior の天井 (Δ=0, n_m=100, U=12 の x)")
for name in OBS:
    x, _, nA = d[(name, "UHF")]
    print(f"    {name:14s} |A|={nA:2d} x={x:+.4f}: " +
          "  ".join(f"p={p}: {gain((1-2*p)**nA*x, (1-2*p)**nA*x, nA, 100):7.1f}" for p in PS))

# ------------------------------------------------------------
# [4] 誤差率の見積もりがずれた場合: prior に掛ける f を p の推定値 p^ で作る
#    推定量 3^|A| (m - f^ y) + f^ y は f x を偏りなく推定し、外れは Δ' = f x - f^ y
# ------------------------------------------------------------
print("\n[4] 誤差率の見積もり p^ がずれたときの利得 (U=12, UHF, n_m=100, 真の p=0.01)")
p = 0.01
for name in OBS:
    x, y, nA = d[(name, "UHF")]
    f = (1 - 2 * p) ** nA
    cells = []
    for ph in (0.0, 0.005, 0.008, 0.009, 0.01, 0.011, 0.012, 0.015):
        fh = (1 - 2 * ph) ** nA
        cells.append(f"p^={ph}: {gain(f*x, fh*y, nA, 100):6.1f}")
    print(f"    {name:14s} |A|={nA:2d}: " + "  ".join(cells))
