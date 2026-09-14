"""HF prior の CRM 利得が観測量によって変わる理由を、既存データで分解して確かめる。

  [1] G = G_max / (1 + R),  R = ε²(G_max - 1): 天井で決まるか、prior の誤差で決まるか
  [2] Slater 行列式の恒等式 d = n²/4 - |<S>|²(二重占有を減らすには局所スピンの期待値が要る)
  [3] スピン相関の誤差の向き: 真の状態は等方的、UHF は z 軸方向に偏る
  [4] Δ = Δ_相関 + Δ_破れ(Δ_破れ = y_sym - y は群平均で消える部分)の符号と大きさ
  [5] 隣接 Z↑Z↑ の y/x = (y/y_sym)·(y_sym/x): 成分の集中 × 古典的な相関と量子的な相関の比
  [6] 天井の U 依存: 電荷は U とともに伸び、スピンは頭打ちになる
"""
import csv, collections, math, statistics as st
from crm_eigen_uhfsym import lat_edges, solve_uhf, wick, gain

V = collections.defaultdict(dict)
for r in csv.DictReader(open("crm_edfid_all.tsv"), delimiter="\t"):
    if r["exact"] == "yes":
        V[(int(r["W"]), r["geometry"], int(r["LX"]), float(r["U"]), r["prior"])][r["observable"]] = (float(r["true"]), float(r["prior_val"]))
systems = sorted({k[:4] for k in V})
bip = lambda k: not (k[0] == 4 and k[1] == "torus")          # 4×3 トーラスだけが非二部格子
Gmax = lambda x, nA=2: gain(x, x, nA)

# ---- [1] -------------------------------------------------------------------------
print("[1] G = G_max/(1+R), R = ε²(G_max-1)。R≪1 なら天井で決まり、R≫1 なら prior の誤差で決まる(UHF, 1次元 L=64)")
print(f"  {'U':>3s} | {'オンサイト ZZ: G_max':>20s} {'ε':>6s} {'R':>7s} {'G':>6s} | {'隣接 Z↑Z↑: G_max':>18s} {'ε':>6s} {'R':>6s} {'G':>5s} | 単一サイト Z↑: G_max  G")
for U in (2.0, 4.0, 8.0, 12.0):
    d = V[(1, "cylinder", 64, U, "UHF")]; cells = []
    for obs in ("ZZ onsite", "ZZ up-up nb"):
        x, y = d[obs]; gm = Gmax(x); eps = abs(x-y)/abs(x); R = eps**2*(gm-1)
        assert abs(gm/(1+R) - gain(x, y)) < 1e-9*gain(x, y)
        cells.append(f"{gm:20.1f} {eps:6.3f} {R:7.3f} {gain(x,y):6.2f}" if obs == "ZZ onsite" else f"{gm:18.1f} {eps:6.3f} {R:6.1f} {gain(x,y):5.2f}")
    xz = 1 - d["n"][0] - 2*d["Sz"][0]; yz = 1 - d["n"][1] - 2*d["Sz"][1]
    print(f"  {U:3.0f} | " + " | ".join(cells) + f" | {Gmax(xz,1):8.2f}  {gain(xz,yz,1):.4f}")

# ---- [2] -------------------------------------------------------------------------
print("\n[2] 二重占有 d と局所スピン: Slater 行列式(共線 UHF)では d = n²/4 - <S^z>² が厳密")
mx = 0; gaps = []
for k in systems:
    for p in ("UHF",):
        d = V[k + (p,)]
        n, sz, zz = d["n"][1], d["Sz"][1], d["ZZ onsite"][1]
        dd = (zz - 1 + 2*n)/4
        if bip(k): mx = max(mx, abs(dd - (n*n/4 - sz*sz)))
    n, sz, zz = d["n"][0], d["Sz"][0], d["ZZ onsite"][0]
    gaps.append((k, (zz - 1 + 2*n)/4, n*n/4 - sz*sz))
print(f"  UHF(二部格子 76 系): |d - (n²/4 - <S^z>²)| の最大 {mx:.1e}")
for k, dt, sl in [g for g in gaps if g[0][:3] == (1, "cylinder", 64)]:
    print(f"  真の状態 L=64 U={k[3]:4.0f}: d = {dt:.4f}, <S>=0 なので n²/4 - |<S>|² = {sl:.4f}  → 恒等式からのずれ {sl-dt:.4f}")

# ---- [3] -------------------------------------------------------------------------
print("\n[3] 隣接スピン相関の等方性(統合表の SzSz nb, SxSx nb)")
iso = 0; sign_ok = 0; nsys = 0; ex = None
for k in systems:
    d = V[k + ("UHF",)]
    (tz, uz), (tx, ux) = d["SzSz nb"], d["SxSx nb"]
    iso = max(iso, abs(tz - tx)); nsys += 1
    if (uz - tz) * (ux - tx) < 0: sign_ok += 1
    if k == (1, "cylinder", 8, 4.0): ex = (tz, tx, uz, ux)
print(f"  真の状態: |<SzSz> - <SxSx>| の最大 {iso:.1e}(80 系)")
print(f"  UHF の誤差の符号が z 成分と x 成分で逆: {sign_ok}/{nsys} 系")
tz, tx, uz, ux = ex
print(f"  例 L=8 U=4: 真 {tz:+.5f} / {tx:+.5f}   UHF {uz:+.5f} / {ux:+.5f}   群平均 {(uz+2*ux)/3:+.5f}")

# ---- [4], [5] 中央の結合 --------------------------------------------------------------
print("\n[4][5] 隣接 Z↑Z↑(中央の結合): Δ = Δ_相関(x - y_sym) + Δ_破れ(y_sym - y)、 y/x = (y/y_sym)·(y_sym/x)")
print(f"  {'系':26s} | {'x':>7s} {'y_sym':>7s} {'y':>7s} | {'Δ_相関':>7s} {'Δ_破れ':>7s} {'符号':4s} | {'y/y_sym':>7s} {'y_sym/x':>7s} {'y/x':>5s}")
opp = 0; tot = 0; rows45 = []
for k in systems:
    W, geo, LX, U = k; n = W*LX; e = lat_edges(LX, W, geo == "torus")
    c0 = (max(1, -(-LX//2))-1)*W + max(1, -(-W//2)); nb = next(b for (a, b) in e if a == c0)
    Gu, Gd = solve_uhf(e, n, LX, W, U, n//2, n//2); w = wick(Gu, Gd, c0-1, nb-1)
    x, y = V[k + ("UHF",)]["ZZ up-up nb"]; ys = w["zz_uu_sym"]
    dc, ds = x - ys, ys - y; tot += 1; opp += dc*ds < 0
    rows45.append((k, x, ys, y, dc, ds))
    if U == 12.0 and bip(k) and (W, geo) in ((1, "torus"), (2, "torus"), (4, "cylinder")) or k == (1, "cylinder", 64, 12.0):
        print(f"  {str(k):26s} | {x:+7.3f} {ys:+7.3f} {y:+7.3f} | {dc:+7.3f} {ds:+7.3f} {'逆' if dc*ds<0 else '同':4s} | {y/ys:7.2f} {ys/x:7.2f} {y/x:5.2f}")
print(f"  中央の結合で Δ_相関 と Δ_破れ の符号が逆(打ち消し合う)系: {opp}/{tot}")
for U in (2.0, 4.0, 8.0, 12.0):
    r_ = [a for a in rows45 if a[0][3] == U and bip(a[0])]
    print(f"  U={U:4.0f}: y/y_sym の中央値 {st.median(a[3]/a[2] for a in r_):.2f}   y_sym/x の中央値 {st.median(a[2]/a[1] for a in r_):.2f}")

# ---- [4] 距離ごと -------------------------------------------------------------------------
print("\n[4] 距離 r ごとの符号(1次元の全系、UHF と UHF-sym の対)")
P = collections.defaultdict(dict)
for r in csv.DictReader(open("crm_eigcorr_pairs.tsv"), delimiter="\t"):
    if r["W"] == "1" and r["prior"] in ("UHF", "UHF-sym"):
        P[(r["geometry"], r["LX"], r["U"], r["i"], r["j"])][r["prior"]] = (int(r["r"]), float(r["x"]), float(r["y"]))
by = collections.defaultdict(lambda: [0, 0, [], []])
for v in P.values():
    (r_, x, yu), (_, _, ys) = v["UHF"], v["UHF-sym"]
    key = r_ if r_ <= 4 else "5+"
    b = by[key]; b[1] += 1; b[0] += (x - ys)*(ys - yu) < 0
    b[2].append(abs(x - ys)); b[3].append(abs(ys - yu))
for key in (1, 2, 3, 4, "5+"):
    o, t, dc, ds = by[key]
    print(f"  r={str(key):>2s}: 符号が逆 {o:4d}/{t:4d} ({100*o/t:3.0f}%)   |Δ_相関| 中央値 {st.median(dc):.3f}   |Δ_破れ| 中央値 {st.median(ds):.3f}")

# ---- [6] --------------------------------------------------------------------------------
print("\n[6] 天井の U 依存(1次元 L=64、真の状態だけで決まる)")
prev = None
for U in (2.0, 4.0, 8.0, 12.0):
    d = V[(1, "cylinder", 64, U, "UHF")]
    xo, xn = d["ZZ onsite"][0], d["ZZ up-up nb"][0]; dd = (xo + 1)/4
    s = f"  U={U:4.0f}: 二重占有 d={dd:.4f}  オンサイト G_max={Gmax(xo):6.1f}  |  隣接 |x|={abs(xn):.3f}  隣接 G_max={Gmax(xn):5.1f}"
    if prev: s += f"   d の指数 d∝U^-a: a={math.log(prev[1]/dd)/math.log(U/prev[0]):.2f}"
    print(s); prev = (U, dd)
xr = {U: V[(1, "torus", 16, U, "UHF")]["ZZ up-up nb"][0] for U in (2.0, 4.0, 8.0, 12.0)}
print("  周期 L=16 の隣接 |x|:", {int(U): round(abs(v), 3) for U, v in xr.items()}, "→ U→∞ のハイゼンベルク環 L=16 は 0.595")
