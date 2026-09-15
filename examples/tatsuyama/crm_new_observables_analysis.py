"""物理的に自然な観測量で、HF を使った CRM の利得を厳密に評価する(和の分散は crm_pauli_sum_variance の式)。

  [A] 二重占有 D = n↑n↓ = (1 - Z↑ - Z↓ + Z↑Z↓)/4 を4項の和として測る場合(統合表の厳密な80系)
"""
import csv, collections, statistics as st
from crm_pauli_sum_variance import variances
from crm_eigen_uhfsym import gain

NM = 100
V = collections.defaultdict(dict)
for r in csv.DictReader(open("crm_edfid_all.tsv"), delimiter="\t"):
    if r["exact"] == "yes":
        V[(int(r["W"]), r["geometry"], int(r["LX"]), float(r["U"]), r["prior"])][r["observable"]] = (float(r["true"]), float(r["prior_val"]))
systems = sorted({k[:4] for k in V})
bip = lambda k: not (k[0] == 4 and k[1] == "torus")

def site_paulis(d, which):          # which: 0 = 真の状態, 1 = prior
    n, sz, zz = d["n"][which], d["Sz"][which], d["ZZ onsite"][which]
    return {"II": 1.0, "ZI": 1 - n - 2*sz, "IZ": 1 - n + 2*sz, "ZZ": zz}

TERMS_D = [(-0.25, "ZI"), (-0.25, "IZ"), (0.25, "ZZ")]
def gain_sum(terms, ex, ey):
    vs, vc = variances(terms, ex, ey, NM); return vs/vc

res = collections.defaultdict(list)
for k in systems:
    ex = site_paulis(V[k + ("UHF",)], 0)
    priors = {}
    pu = site_paulis(V[k + ("UHF",)], 1)
    priors["UHF"] = pu
    priors["UHF-sym"] = {"II": 1.0, "ZI": (pu["ZI"] + pu["IZ"])/2, "IZ": (pu["ZI"] + pu["IZ"])/2, "ZZ": pu["ZZ"]}   # 階数0の部分だけ残す
    priors["RHF"] = {"II": 1.0, "ZI": 0.0, "IZ": 0.0, "ZZ": 0.0}                                                  # 一様密度 1、d = 1/4
    for c in ("chi2", "chi4", "chi8"):
        if k + (c,) in V: priors[c] = site_paulis(V[k + (c,)], 1)
    for lab, ey in priors.items():
        res[("D 4項", k[3], lab)].append(gain_sum(TERMS_D, ex, ey))
        res[("ZZ 1項", k[3], lab)].append(gain(ex["ZZ"], ey["ZZ"], 2))

print("[A] 二重占有: 4項の和として測る場合と、Z↑Z↓ だけで測る場合(80系の中央値 [最小])")
print(f"  {'U':>3s} | " + " | ".join(f"{lab:>8s}: 4項 / ZZ だけ" for lab in ("UHF", "UHF-sym", "RHF", "chi4", "chi8")))
for U in (2.0, 4.0, 8.0, 12.0):
    cells = []
    for lab in ("UHF", "UHF-sym", "RHF", "chi4", "chi8"):
        a, b = res[("D 4項", U, lab)], res[("ZZ 1項", U, lab)]
        cells.append(f"{st.median(a):7.2f}[{min(a):5.2f}] / {st.median(b):6.1f}")
    print(f"  {U:3.0f} | " + " | ".join(cells))

# 1例の内訳: どの項が CRM の分散を支配しているか
k = (1, "cylinder", 64, 12.0)
ex = site_paulis(V[k + ("UHF",)], 0); ey = site_paulis(V[k + ("UHF",)], 1)
print(f"\n  内訳 {k}: 真の値 Z↑={ex['ZI']:+.4f} Z↓={ex['IZ']:+.4f} Z↑Z↓={ex['ZZ']:+.4f} / UHF Z↑={ey['ZI']:+.4f} Z↓={ey['IZ']:+.4f} Z↑Z↓={ey['ZZ']:+.4f}")
for sub, lab in (([(-0.25, "ZI"), (-0.25, "IZ")], "Z↑ と Z↓ の項だけ"), ([(0.25, "ZZ")], "Z↑Z↓ の項だけ"), (TERMS_D, "4項すべて")):
    vs, vc = variances(sub, ex, ey, NM)
    print(f"    {lab:18s}: 標準 {vs:.5f}  CRM(UHF) {vc:.5f}")

# =====================================================================================
#  [B][C][D] clara の crm_new_observables.jl の出力(中央の結合と距離 r)
# =====================================================================================
import glob
import numpy as np
from crm_eigen_uhfsym import lat_edges, solve_uhf, rhf_G

T = collections.defaultdict(dict)      # (LX, geo, U) -> {(state, kind, i, j, name): value}
HV = {}
for fn in sorted(glob.glob("crm_newobs_W1L*_*_U*.tsv")):
    for r in csv.DictReader(open(fn), delimiter="\t"):
        k = (int(r["LX"]), r["geometry"], float(r["U"]))
        T[k][(r["state"], r["kind"], int(r["i"]), int(r["j"]), r["name"])] = float(r["value"]); HV[k] = float(r["Hvar"])
if T:
    print("\n系:", [(k, f"{HV[k]:.0e}") for k in sorted(T)])

MASK = {"Za": "ZIII", "Zb": "IZII", "Zc": "IIZI", "Zd": "IIIZ"}
def zlabel(name):                        # "Zabd" -> "ZZIZ"
    return "".join("Z" if q in name[1:] else "I" for q in "abcd")

def bond_dicts(t, state, i, full):
    """4量子ビット (a=i↑, b=i↓, c=j↑, d=j↓) のパウリ文字列 → 期待値。full なら積に要る量もすべて入れる。"""
    g = lambda nm: t[(state, "bond", i, i+1, nm)]
    e = {"IIII": 1.0, "XZXI": g("hopU"), "YZYI": g("hopU"), "IXZX": g("hopD"), "IYZY": g("hopD")}
    for nm in ("Za", "Zb", "Zc", "Zd", "Zab", "Zcd"):
        e[zlabel(nm)] = g(nm)
    if full:
        for m in range(1, 16):
            nm = "Z" + "".join("abcd"[q] for q in range(4) if (m >> q) & 1); e[zlabel(nm)] = g(nm)
        e["XIXI"] = e["YIYI"] = g("hopU_Zb"); e["XZXZ"] = e["YZYZ"] = g("hopU_Zd")
        e["IXIX"] = e["IYIY"] = g("hopD_Zc"); e["ZXZX"] = e["ZYZY"] = g("hopD_Za")
    return e

def uhfsym_of(e):
    s = dict(e); hop = (e["XZXI"] + e["IXZX"])/2
    for p in ("XZXI", "YZYI", "IXZX", "IYZY"): s[p] = hop
    zi = (e["ZIII"] + e["IZII"])/2; zj = (e["IIZI"] + e["IIIZ"])/2
    s["ZIII"] = s["IZII"] = zi; s["IIZI"] = s["IIIZ"] = zj
    return s

KIN = [(-0.5, "XZXI"), (-0.5, "YZYI"), (-0.5, "IXZX"), (-0.5, "IYZY")]
def energy_terms(U):
    return KIN + [(-U/8, "ZIII"), (-U/8, "IZII"), (U/8, "ZZII"), (-U/8, "IIZI"), (-U/8, "IIIZ"), (U/8, "IIZZ")]

if T:
    outB = collections.defaultdict(list)
    for (LX, geo, U), t in sorted(T.items()):
        n = LX; e_ = lat_edges(LX, 1, geo == "torus")
        Gr = rhf_G(e_, n, n//2)
        bonds = sorted({kk[2] for kk in t if kk[1] == "bond"})
        for i in bonds:
            ex = bond_dicts(t, "true", i, True)
            priors = {"UHF": bond_dicts(t, "UHF", i, False)}
            priors["UHF-sym"] = uhfsym_of(priors["UHF"])
            hr = 2*Gr[i-1, i]
            priors["RHF"] = {"IIII": 1.0, "XZXI": hr, "YZYI": hr, "IXZX": hr, "IYZY": hr, "ZIII": 0.0, "IZII": 0.0, "IIZI": 0.0, "IIIZ": 0.0, "ZZII": 0.0, "IIZZ": 0.0}
            for c in ("chi4", "chi8", "chi32"):
                if ("chi4" if c == "chi4" else c, "bond", i, i+1, "hopU") in t: priors[c] = bond_dicts(t, c, i, False)
            for lab, ey in priors.items():
                vs, vc = variances(KIN, ex, ey, NM); outB[("運動 4項", U, lab)].append(vs/vc)
                outB[("hop↑ 1項", U, lab)].append(gain(ex["XZXI"], ey["XZXI"], 3))
                vs, vc = variances(energy_terms(U), ex, ey, NM); outB[("結合エネルギー", U, lab)].append(vs/vc)
                vs, vc = variances(energy_terms(U)[4:], ex, ey, NM); outB[("相互作用 6項", U, lab)].append(vs/vc)
            outB[("|hop↑| 真", U, "")].append(abs(ex["XZXI"]))
    for obs in ("hop↑ 1項", "運動 4項", "相互作用 6項", "結合エネルギー"):
        print(f"\n[B/C] {obs}(中央の結合、全系の中央値 [最小])")
        for U in (2.0, 4.0, 8.0, 12.0):
            if not outB.get((obs, U, "UHF")): continue
            cells = []
            for lab in ("UHF", "UHF-sym", "RHF", "chi4", "chi8", "chi32"):
                v = outB.get((obs, U, lab))
                cells.append(f"{lab} {st.median(v):7.2f}[{min(v):6.2f}]" if v else "")
            extra = f"  |hop↑|={st.median(outB[('|hop↑| 真', U, '')]):.3f}" if obs == "hop↑ 1項" else ""
            print(f"  U={U:4.0f}: " + "  ".join(c for c in cells if c) + extra)

    # ---- [D] 距離 r: 電荷の偶奇の2点相関と、ホッピング ----
    print("\n[D] 距離 r(1次元 L=64 開放端と周期 L=16)")
    for key in [k for k in sorted(T) if k[0] in (64, 16)]:
        LX, geo, U = key; t = T[key]; n = LX; e_ = lat_edges(LX, 1, geo == "torus")
        Gu, Gd = solve_uhf(e_, n, LX, 1, U, n//2, n//2); Gr = rhf_G(e_, n, n//2)
        dist = sorted({(kk[2], kk[3]) for kk in t if kk[1] == "dist"})
        print(f"  {geo} L={LX} U={U:.0f}")
        print(f"    {'r':>2s} | {'FF 真':>7s} {'UHF':>7s} {'G_UHF':>7s} {'G_RHF':>6s} {'G_χ4':>6s} | {'hop 真':>7s} {'G_UHF':>6s} {'G_sym':>6s} {'G_RHF':>6s} {'G_χ8':>6s}")
        for (i, j) in dist:
            r_ = j - i
            if r_ > 6 and r_ not in (8, 12, 16, 31): continue
            B = [i-1, j-1]; I2 = np.eye(2)
            ff = lambda st_: t[(st_, "dist", i, j, "FF")]
            x = ff("true"); yu = ff("UHF")
            yr = np.linalg.det(I2 - 2*Gr[np.ix_(B, B)])**2
            g4 = gain(x, ff("chi4"), 4) if ("chi4", "dist", i, j, "FF") in t else float("nan")
            hx = t[("true", "dist", i, j, "hopU")]; hu = t[("UHF", "dist", i, j, "hopU")]
            hsym = (2*Gu[i-1, j-1] + 2*Gd[i-1, j-1])/2; hr = 2*Gr[i-1, j-1]
            nA = 2*r_ + 1
            g8 = gain(hx, t[("chi8", "dist", i, j, "hopU")], nA) if ("chi8", "dist", i, j, "hopU") in t else float("nan")
            print(f"    {r_:2d} | {x:+7.3f} {yu:+7.3f} {gain(x, yu, 4):7.2f} {gain(x, yr, 4):6.2f} {g4:6.2f} | {hx:+7.3f} {gain(hx, hu, nA):6.2f} {gain(hx, hsym, nA):6.2f} {gain(hx, hr, nA):6.2f} {g8:6.2f}")
