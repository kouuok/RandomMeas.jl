"""スピン射影 HF(PHF、UHF を射影したもの)を prior にしたときの CRM の利得を、UHF・UHF-sym・RHF・MPS と比べる。

  [A] 統合表の厳密な 80 系: オンサイト ZZ、隣接 Z↑Z↑、単一サイト Z↑、二重占有(4項の和)、<S_i·S_j>
  [B] 射影の効果の系の大きさ依存(1次元開放端)
  [C] 距離 r の Z↑Z↑(crm_eigcorr_pairs.tsv の対)
  [D] ブロックの電荷の偶奇 P_ℓ(crm_eigparity_*.tsv)
  [E] 変分エネルギー
"""
import csv, collections, glob, statistics as st
import numpy as np
from crm_phf import PHF, solve_uhf_orbitals
from crm_eigen_uhfsym import lat_edges, solve_uhf, rhf_G, wick, gain
from crm_pauli_sum_variance import variances

NM = 100; NB = 64
V = collections.defaultdict(dict)
for r in csv.DictReader(open("crm_edfid_all.tsv"), delimiter="\t"):
    if r["exact"] == "yes":
        V[(int(r["W"]), r["geometry"], int(r["LX"]), float(r["U"]), r["prior"])][r["observable"]] = (float(r["true"]), float(r["prior_val"]))
systems = sorted({k[:4] for k in V})

CACHE = {}
def models(W, geo, LX, U):
    key = (W, geo, LX, U)
    if key not in CACHE:
        n = W*LX; e = lat_edges(LX, W, geo == "torus")
        Pu, Pd = solve_uhf_orbitals(e, n, LX, W, U, n//2, n//2)
        Gu, Gd = Pu @ Pu.T, Pd @ Pd.T
        CACHE[key] = (e, Gu, Gd, rhf_G(e, n, n//2), PHF(Pu, Pd, nbeta=NB))
    return CACHE[key]

def center(W, LX, e):
    c0 = (max(1, -(-LX//2))-1)*W + max(1, -(-W//2)); nb = next(b for (a, b) in e if a == c0)
    return c0 - 1, nb - 1                       # 0 始まり

# ---------------- [A] ----------------
res = collections.defaultdict(list); ss_rows = []
for k in systems:
    W, geo, LX, U = k; e, Gu, Gd, Gr, ph = models(*k)
    i, j = center(W, LX, e)
    d = V[k + ("UHF",)]
    x_on, x_nb = d["ZZ onsite"][0], d["ZZ up-up nb"][0]
    x_z = 1 - d["n"][0] - 2*d["Sz"][0]
    wu = wick(Gu, Gd, i, j); wr = wick(Gr, Gr, i, j)
    pri = {
        "UHF":     dict(on=wu["zz_on"], nb=wu["zz_uu"], z=1-2*Gu[i, i], zd=1-2*Gd[i, i]),
        "UHF-sym": dict(on=wu["zz_on"], nb=wu["zz_uu_sym"], z=1-Gu[i, i]-Gd[i, i], zd=1-Gu[i, i]-Gd[i, i]),
        "PHF":     dict(on=ph.ZupZdn_onsite(i), nb=ph.ZupZup(i, j), z=ph.Zup(i), zd=ph.Zup(i)),
        "RHF":     dict(on=wr["zz_on"], nb=wr["zz_uu"], z=0.0, zd=0.0)}
    ex = {"II": 1.0, "ZI": x_z, "IZ": 1 - d["n"][0] + 2*d["Sz"][0], "ZZ": x_on}
    for lab, p in pri.items():
        res[("オンサイト ZZ", U, lab)].append(gain(x_on, p["on"], 2))
        res[("隣接 Z↑Z↑", U, lab)].append(gain(x_nb, p["nb"], 2))
        res[("単一サイト Z↑", U, lab)].append(gain(x_z, p["z"], 1))
        ey = {"II": 1.0, "ZI": p["z"], "IZ": p["zd"], "ZZ": p["on"]}
        vs, vc = variances([(-0.25, "ZI"), (-0.25, "IZ"), (0.25, "ZZ")], ex, ey, NM)
        res[("二重占有 4項", U, lab)].append(vs/vc)
    ss_rows.append((k, 3*d["SzSz nb"][0], wu["ss"], ph.SS(i, j)))

print("[A] 統合表の厳密な 80 系(中央値 [最小]、損をする系の数)")
for obs in ("オンサイト ZZ", "隣接 Z↑Z↑", "単一サイト Z↑", "二重占有 4項"):
    print(f"  {obs}")
    for U in (2.0, 4.0, 8.0, 12.0):
        cells = []
        for lab in ("UHF", "UHF-sym", "PHF", "RHF"):
            v = res[(obs, U, lab)]
            cells.append(f"{lab} {st.median(v):7.2f}[{min(v):6.2f}] 損{sum(g < 1-1e-9 for g in v):2d}")
        print(f"    U={U:4.0f}: " + "  ".join(cells))

print("\n[A'] 隣接の <S_i·S_j>: 真 / UHF(= UHF-sym)/ PHF(1次元開放端)")
for (k, st_, su, sp) in ss_rows:
    if k[0] == 1 and k[1] == "cylinder" and k[2] in (8, 16, 64):
        print(f"    L={k[2]:2d} U={k[3]:4.0f}: 真 {st_:+.4f}  UHF {su:+.4f}  PHF {sp:+.4f}   PHF が真に近づいた割合 {(sp-su)/(st_-su):+.2f}")

# ---------------- [B] 射影の効果の系の大きさ依存 ----------------
print("\n[B] 射影による隣接 <S·S> の変化 ΔSS = PHF - UHF と、真の値との差(1次元)")
print(f"  {'系':18s} {'U':>3s} | {'ΔSS':>8s} {'L·ΔSS':>7s} | {'真 - UHF':>9s}")
for U in (8.0, 12.0):
    for (geo, LX) in [("cylinder", L) for L in (8, 10, 12, 14, 16, 24, 32, 48, 64)] + [("torus", L) for L in (8, 12, 16)]:
        k = (1, geo, LX, U)
        if k + ("UHF",) not in V: continue
        e, Gu, Gd, Gr, ph = models(*k); i, j = center(1, LX, e)
        su = wick(Gu, Gd, i, j)["ss"]; sp = ph.SS(i, j); stv = 3*V[k + ("UHF",)]["SzSz nb"][0]
        print(f"  {geo+' L='+str(LX):18s} {U:3.0f} | {sp-su:+8.4f} {LX*(sp-su):+7.3f} | {stv-su:+9.4f}")

# ---------------- [C] 距離 r の Z↑Z↑ ----------------
print("\n[C] 距離 r の Z↑Z↑(中央値の G。損 = G<1 の対の割合)")
pairs = collections.defaultdict(dict)
for r in csv.DictReader(open("crm_eigcorr_pairs.tsv"), delimiter="\t"):
    if r["W"] == "1" and r["prior"] in ("UHF", "UHF-sym", "RHF", "chi8"):
        pairs[(r["geometry"], int(r["LX"]), float(r["U"]), int(r["r"]), int(r["i"]), int(r["j"]))][r["prior"]] = (float(r["x"]), float(r["y"]))
for (geo, LX) in (("cylinder", 64), ("torus", 16)):
    for U in (8.0, 12.0):
        e, Gu, Gd, Gr, ph = models(1, geo, LX, U)
        print(f"  {geo} L={LX} U={U:.0f}")
        rs = sorted({k[3] for k in pairs if k[:3] == (geo, LX, U)})
        for r_ in [r for r in rs if r <= 6 or r in (8, 16, 31)]:
            ks = [k for k in pairs if k[:4] == (geo, LX, U, r_)]
            cells = []
            for lab in ("UHF", "UHF-sym", "PHF", "RHF", "chi8"):
                gs = []
                for k in ks:
                    x = pairs[k]["UHF"][0]
                    y = ph.ZupZup(k[4]-1, k[5]-1) if lab == "PHF" else pairs[k][lab][1]
                    gs.append(gain(x, y, 2))
                cells.append(f"{lab} {st.median(gs):7.2f} 損{sum(g<1 for g in gs)/len(gs):4.0%}")
            print(f"    r={r_:2d}: " + "  ".join(cells))

# ---------------- [D] ブロックの電荷の偶奇 ----------------
print("\n[D] ブロックの電荷の偶奇 P_ℓ(1次元 L=64、始点 17)")
par = collections.defaultdict(dict)
for fn in glob.glob("crm_eigparity_W1L64_cylinder_U*.tsv"):
    for r in csv.DictReader(open(fn), delimiter="\t"):
        if r["kind"] == "P" and r["start"] == "17":
            par[(float(r["U"]), int(r["l"]))][r["prior"]] = (float(r["true"]), float(r["prior_val"]))
for U in (8.0, 12.0):
    e, Gu, Gd, Gr, ph = models(1, "cylinder", 64, U)
    row = []
    for l in (1, 2, 4, 8, 16, 32):
        x, yu = par[(U, l)]["UHF"]; yp = ph.P_block(list(range(16, 16+l)))
        row.append(f"ℓ={l}: UHF {gain(x, yu, 2*l):6.1f} PHF {gain(x, yp, 2*l):6.1f}")
    print(f"  U={U:.0f}: " + "  ".join(row))

# ---------------- [E] 変分エネルギー ----------------
print("\n[E] 1サイトあたりのエネルギー(1次元)")
def energy(ph_or_G, e, n, U, kind):
    if kind == "PHF":
        kin = sum(-2*ph_or_G.hop(a-1, b-1) for a, b in e)
        pot = sum(U*(ph_or_G.ZupZdn_onsite(i) - 1 + 2*(1 - ph_or_G.one_minus_n(i)))/4 for i in range(n))
    else:
        Gu, Gd = ph_or_G
        kin = sum(-(Gu[a-1, b-1] + Gu[b-1, a-1] + Gd[a-1, b-1] + Gd[b-1, a-1]) for a, b in e)
        pot = sum(U*Gu[i, i]*Gd[i, i] for i in range(n))
    return (kin + pot)/n
Eref = {}
for r in csv.DictReader(open("crm_edfid_all.tsv"), delimiter="\t"):
    Eref[(int(r["W"]), r["geometry"], int(r["LX"]), float(r["U"]))] = float(r["Eref"])
for (geo, LX) in (("cylinder", 8), ("cylinder", 16), ("cylinder", 64), ("torus", 16)):
    for U in (4.0, 8.0, 12.0):
        k = (1, geo, LX, U); e, Gu, Gd, Gr, ph = models(*k); n = LX
        eu = energy((Gu, Gd), e, n, U, "UHF"); ep = energy(ph, e, n, U, "PHF"); ex = Eref[k]/n
        print(f"  {geo} L={LX:2d} U={U:4.0f}: 厳密 {ex:+.5f}  UHF {eu:+.5f}  PHF {ep:+.5f}   相関エネルギーのうち射影で得た割合 {(ep-eu)/(ex-eu):.3f}   L×(E_PHF-E_UHF) {n*(ep-eu):+.3f}")
