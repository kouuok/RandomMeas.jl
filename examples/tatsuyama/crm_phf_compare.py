"""スピン射影 HF(PHF、UHF を射影したもの)を prior にしたときの CRM の利得を、UHF・UHF-sym・RHF・低い結合次元の MPS と比べる。

  MPS は2種類:
    切断 MPS: 厳密な状態を χ に切り詰めたもの(crm_edfid_all.tsv / crm_eigcorr_pairs.tsv / crm_eigparity_*.tsv の chi*)。
              厳密な状態がないと作れないので、「χ で表せるほぼ最良の prior」の目安。
    変分 MPS: 結合次元を χ に制限した DMRG の基底状態(crm_mps_lowchi.jl、clara の crm_mpslow_*.tsv)。
              厳密な状態を知らなくても作れる。

  [0] MPS のデータの照合
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
CH = (2, 4, 8, 16)
HF = ("UHF", "UHF-sym", "PHF", "RHF")
MPS = tuple(f"{m}χ{c}" for m in ("切断", "変分") for c in CH)
V = collections.defaultdict(dict); Eref = {}
for r in csv.DictReader(open("crm_edfid_all.tsv"), delimiter="\t"):
    Eref[(int(r["W"]), r["geometry"], int(r["LX"]), float(r["U"]))] = float(r["Eref"])
    if r["exact"] == "yes":
        V[(int(r["W"]), r["geometry"], int(r["LX"]), float(r["U"]), r["prior"])][r["observable"]] = (float(r["true"]), float(r["prior_val"]))
systems = sorted({k[:4] for k in V})

# crm_mps_lowchi.jl の出力: (W, geo, LX, U, 方法, χ) -> {(量, i, j): 値}。方法は exact / trunc / var
MV = collections.defaultdict(dict)
for fn in sorted(glob.glob("crm_mpslow_W*L*_*_U*.tsv")):
    for r in csv.DictReader(open(fn), delimiter="\t"):
        MV[(int(r["W"]), r["geometry"], int(r["LX"]), float(r["U"]), r["method"], int(r["chi"]))][(r["quantity"], int(r["i"]), int(r["j"]))] = float(r["value"])
CENTRE = ("ZZ onsite", "ZZ up-up nb", "SzSz nb", "SxSx nb", "DoubleOcc", "n", "Sz")
def centre_vals(d): return {q: v for (q, i, j), v in d.items() if q in CENTRE}

def mps_values(k):
    """系 k の MPS prior の中央の値 {ラベル: {観測量: 値}}。"""
    out = {}
    for c in CH:
        if k + (f"chi{c}",) in V:
            out[f"切断χ{c}"] = {o: v[1] for o, v in V[k + (f"chi{c}",)].items()}
        if k + ("var", c) in MV:
            out[f"変分χ{c}"] = centre_vals(MV[k + ("var", c)])
    return out

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

# ---------------- [0] MPS のデータの照合 ----------------
print("[0] MPS のデータの照合")
w_tr = 0.0; n_tr = 0; w_e = 0.0; w_z = 0.0; w_p = 0.0; below = 0; nvar = 0
spread = collections.defaultdict(list); agree = collections.defaultdict(list)
for key, d in MV.items():
    W, geo, LX, U, m, c = key; k = (W, geo, LX, U); n = W*LX
    if m == "trunc" and k + (f"chi{c}",) in V:              # 統合表の切断 MPS と同じ値になるはず
        for o, v in centre_vals(d).items():
            w_tr = max(w_tr, abs(v - V[k + (f"chi{c}",)][o][1])); n_tr += 1
    if m == "exact":
        w_e = max(w_e, abs(d[("E", 0, 0)] - Eref[k]))
    if m != "var": continue
    nvar += 1
    i, j = center(W, LX, lat_edges(LX, W, geo == "torus")); cv = centre_vals(d)
    if ("Zup", i+1, i+1) in d:                               # 相関行列から組んだ値と MPO の値
        w_z = max(w_z, abs(d[("ZupZdn", i+1, i+1)] - cv["ZZ onsite"]),
                  abs(d[("ZupZup", min(i, j)+1, max(i, j)+1)] - cv["ZZ up-up nb"]),
                  abs(d[("Zup", i+1, i+1)] - (1 - cv["n"] - 2*cv["Sz"])))
    for (q, s0, l), v in d.items():                          # ℓ=1 の電荷の偶奇はオンサイト ZZ
        if q == "P" and l == 1 and ("ZupZdn", s0, s0) in d:
            w_p = max(w_p, abs(v - d[("ZupZdn", s0, s0)]))
    Es = [v for (q, a, b), v in d.items() if q == "E_start"]
    spread[c].append((max(Es) - min(Es))/n); agree[c].append(sum(e_ < min(Es) + 1e-6*n for e_ in Es))
    below += d[("E", 0, 0)] < Eref[k] - 1e-8
print(f"  切断 MPS の中央の値 vs 統合表: {n_tr} 値、最大差 {w_tr:.1e}   厳密なエネルギー vs 統合表: 最大差 {w_e:.1e}")
print(f"  変分 MPS {nvar} 個: 相関行列 vs MPO の最大差 {w_z:.1e}、P(ℓ=1) vs オンサイト ZZ {w_p:.1e}、厳密なエネルギーを下回ったもの {below}")
for c in CH:
    if spread[c]:
        print(f"    χ={c:2d}: 出発点ごとのエネルギーの幅(1サイトあたり)中央値 {st.median(spread[c]):.1e} 最大 {max(spread[c]):.1e}   "
              f"最良から 1e-6 以内に入った出発点の数 中央値 {st.median(agree[c])} 最小 {min(agree[c])}")

# ---------------- [A] ----------------
res = collections.defaultdict(list); ss_rows = []; grp = collections.defaultdict(list); nb_gain = {}
for k in systems:
    W, geo, LX, U = k; e, Gu, Gd, Gr, ph = models(*k)
    i, j = center(W, LX, e)
    d = V[k + ("UHF",)]
    x_on, x_nb = d["ZZ onsite"][0], d["ZZ up-up nb"][0]
    x_z = 1 - d["n"][0] - 2*d["Sz"][0]
    wu = wick(Gu, Gd, i, j); wr = wick(Gr, Gr, i, j)
    pri = {
        "UHF":     dict(on=wu["zz_on"], nb=wu["zz_uu"], z=1-2*Gu[i, i], zd=1-2*Gd[i, i], ss=wu["ss"]),
        "UHF-sym": dict(on=wu["zz_on"], nb=wu["zz_uu_sym"], z=1-Gu[i, i]-Gd[i, i], zd=1-Gu[i, i]-Gd[i, i], ss=wu["ss"]),
        "PHF":     dict(on=ph.ZupZdn_onsite(i), nb=ph.ZupZup(i, j), z=ph.Zup(i), zd=ph.Zup(i), ss=ph.SS(i, j)),
        "RHF":     dict(on=wr["zz_on"], nb=wr["zz_uu"], z=0.0, zd=0.0, ss=wr["ss"])}
    for lab, p in mps_values(k).items():
        pri[lab] = dict(on=p["ZZ onsite"], nb=p["ZZ up-up nb"], z=1-p["n"]-2*p["Sz"], zd=1-p["n"]+2*p["Sz"],
                        ss=p["SzSz nb"] + 2*p["SxSx nb"])
    ex = {"II": 1.0, "ZI": x_z, "IZ": 1 - d["n"][0] + 2*d["Sz"][0], "ZZ": x_on}
    for lab, p in pri.items():
        res[("オンサイト ZZ", U, lab)].append(gain(x_on, p["on"], 2))
        res[("隣接 Z↑Z↑", U, lab)].append(gain(x_nb, p["nb"], 2)); nb_gain[(k, lab)] = res[("隣接 Z↑Z↑", U, lab)][-1]
        res[("単一サイト Z↑", U, lab)].append(gain(x_z, p["z"], 1))
        ey = {"II": 1.0, "ZI": p["z"], "IZ": p["zd"], "ZZ": p["on"]}
        vs, vc = variances([(-0.25, "ZI"), (-0.25, "IZ"), (0.25, "ZZ")], ex, ey, NM)
        res[("二重占有 4項", U, lab)].append(vs/vc)
        if U >= 4:
            grp[("16 サイト以下" if W*LX <= 16 else "24 サイト以上", lab)].append(gain(x_nb, p["nb"], 2))
    ss_rows.append((k, 3*d["SzSz nb"][0], {lab: p["ss"] for lab, p in pri.items()}))

def cell(v):
    s = f"{st.median(v):7.2f}[{min(v):6.2f}] 損{sum(g < 0.99 for g in v):2d}"
    return s if len(v) == 20 else s + f" (n={len(v)})"
print("\n[A] 統合表の厳密な 80 系(中央値 [最小]、1% より大きく損をする系の数)")
for obs in ("オンサイト ZZ", "隣接 Z↑Z↑", "単一サイト Z↑", "二重占有 4項"):
    print(f"  {obs}")
    for U in (2.0, 4.0, 8.0, 12.0):
        print(f"    U={U:4.0f}: " + "  ".join(f"{lab} {cell(res[(obs, U, lab)])}" for lab in HF))
        for m in ("切断", "変分"):
            labs = [f"{m}χ{c}" for c in CH if res.get((obs, U, f"{m}χ{c}"))]
            if labs:
                print("            " + "  ".join(f"{lab} {cell(res[(obs, U, lab)])}" for lab in labs))

print("\n[A2] 隣接 Z↑Z↑ の利得の中央値を系の大きさで分ける(U=4〜12)")
for g in ("16 サイト以下", "24 サイト以上"):
    print(f"  {g}: " + "  ".join(f"{lab} {st.median(grp[(g, lab)]):.2f}" for lab in HF + MPS if grp.get((g, lab)))
          + f"   (系の数 {len(grp[(g, 'UHF')])})")
    rat = [p/s for p, s in zip(grp[(g, "PHF")], grp[(g, "UHF-sym")])]
    print(f"    系ごとの PHF / UHF-sym: 中央値 {st.median(rat):.3f}  範囲 {min(rat):.2f}〜{max(rat):.2f}  1 未満の系 {sum(r < 1 for r in rat)}")

print("\n[A3] MPS がスピン回転の対称性を破っているか: 中央サイトの |<S^z>|(中央値 [最小, 最大])と、隣接 Z↑Z↑ で損をする系")
for U in (4.0, 12.0):
    for m, mm in (("切断", "chi"), ("変分", "var")):
        cells = []
        for c in CH:
            if mm == "chi":
                sz = [abs(V[k + (f"chi{c}",)]["Sz"][1]) for k in systems if k[3] == U and k + (f"chi{c}",) in V]
            else:
                sz = [abs(centre_vals(MV[k + ("var", c)])["Sz"]) for k in systems if k[3] == U and k + ("var", c) in MV]
            if sz:
                cells.append(f"χ{c} {st.median(sz):.3f} [{min(sz):.3f}, {max(sz):.3f}]")
        print(f"  U={U:4.0f} {m}: " + "  ".join(cells))
print("  変分 MPS の中央 |<S^z>| を系ごとに(U=12、χ=4 / 8 / 16)")
for k in systems:
    if k[3] == 12.0 and all(k + ("var", c) in MV for c in (4, 8, 16)):
        print(f"    {str(k[:3]):24s} " + " / ".join(f"{abs(centre_vals(MV[k + ('var', c)])['Sz']):.3f}" for c in (4, 8, 16)))
for lab in MPS:
    lose = [k[:3] for (k, l), g in nb_gain.items() if l == lab and k[3] == 12.0 and g < 1]
    if lose:
        print(f"  隣接 Z↑Z↑ で損をする系(U=12、{lab}): {lose}")
print("  中央サイトの二重占有 d(真 / UHF / 切断χ4 / 切断χ8 / 変分χ4 / 変分χ8)")
for k in [(1, "cylinder", 16, 4.0), (1, "cylinder", 16, 12.0), (1, "cylinder", 64, 4.0), (1, "cylinder", 64, 12.0)]:
    t, yu = V[k + ("UHF",)]["DoubleOcc"]
    print(f"    {str(k):28s} 真 {t:.4f}  UHF {yu:.4f}  " + "  ".join(
        f"{lab} {mps_values(k)[lab]['DoubleOcc']:.4f}" for lab in ("切断χ4", "切断χ8", "変分χ4", "変分χ8")))

print("\n[A'] 隣接の <S_i·S_j>(1次元開放端)。括弧は UHF から真の値までの差のうち埋まった割合")
for (k, sst, ssd) in ss_rows:
    if k[0] == 1 and k[1] == "cylinder" and k[2] in (8, 16, 64):
        su = ssd["UHF"]
        print(f"    L={k[2]:2d} U={k[3]:4.0f}: 真 {sst:+.4f}  UHF {su:+.4f}  " +
              "  ".join(f"{lab} {ssd[lab]:+.4f}({(ssd[lab]-su)/(sst-su):+.2f})" for lab in ("PHF",) + MPS if lab in ssd))

# ---------------- [B] 射影の効果の系の大きさ依存 ----------------
print("\n[B] 射影による隣接 <S·S> の変化 ΔSS = PHF - UHF と、真の値との差(1次元)。回転子の描像の予測は L·ΔSS = -2m")
print(f"  {'系':18s} {'U':>3s} | {'ΔSS':>8s} {'L·ΔSS':>7s} {'-2m':>7s} | {'真 - UHF':>9s}")
for U in (4.0, 8.0, 12.0):
    for (geo, LX) in [("cylinder", L) for L in (8, 10, 12, 14, 16, 24, 32, 48, 64)] + [("torus", L) for L in (8, 12, 16)]:
        k = (1, geo, LX, U)
        if k + ("UHF",) not in V: continue
        e, Gu, Gd, Gr, ph = models(*k); i, j = center(1, LX, e)
        su = wick(Gu, Gd, i, j)["ss"]; sp = ph.SS(i, j); stv = 3*V[k + ("UHF",)]["SzSz nb"][0]
        m = abs(Gu[i, i] - Gd[i, i])/2
        print(f"  {geo+' L='+str(LX):18s} {U:3.0f} | {sp-su:+8.4f} {LX*(sp-su):+7.3f} {-2*m:+7.3f} | {stv-su:+9.4f}")
print("  距離 r の対(1次元開放端 L=64、中央のサイトから): L·ΔSS。回転子の描像では奇数 r(副格子が違う)で -2m、偶数 r で 0")
for U in (4.0, 8.0, 12.0):
    e, Gu, Gd, Gr, ph = models(1, "cylinder", 64, U); c = 31
    print(f"    U={U:4.0f}: " + "  ".join(f"r={r}: {64*(ph.SS(c, c+r) - wick(Gu, Gd, c, c+r)['ss']):+.3f}" for r in range(1, 7)))

# ---------------- [C] 距離 r の Z↑Z↑ ----------------
print("\n[C] 距離 r の Z↑Z↑(中央値の G。損 = G<1 の対の割合)")
pairs = collections.defaultdict(dict)                 # 切断 MPS はこのファイルにある χ=4, 8 だけ
for r in csv.DictReader(open("crm_eigcorr_pairs.tsv"), delimiter="\t"):
    if r["W"] == "1":
        pairs[(r["geometry"], int(r["LX"]), float(r["U"]), int(r["r"]), int(r["i"]), int(r["j"]))][r["prior"]] = (float(r["x"]), float(r["y"]))
for (geo, LX) in (("cylinder", 64), ("torus", 16)):
    for U in (8.0, 12.0):
        e, Gu, Gd, Gr, ph = models(1, geo, LX, U)
        print(f"  {geo} L={LX} U={U:.0f}")
        rs = sorted({k[3] for k in pairs if k[:3] == (geo, LX, U)})
        for r_ in [r for r in rs if r <= 6 or r in (8, 16, 31)]:
            ks = [k for k in pairs if k[:4] == (geo, LX, U, r_)]
            cells = []
            for lab in HF + MPS:
                gs = []
                for k in ks:
                    x = pairs[k]["UHF"][0]; ii, jj = k[4], k[5]
                    if lab == "PHF":           y = ph.ZupZup(ii-1, jj-1)
                    elif lab.startswith("切断"): y = pairs[k].get("chi" + lab[3:], (None, None))[1]
                    elif lab.startswith("変分"): y = MV.get((1, geo, LX, U, "var", int(lab[3:])), {}).get(("ZupZup", ii, jj))
                    else:                      y = pairs[k][lab][1]
                    if y is None: break
                    gs.append(gain(x, y, 2))
                if len(gs) == len(ks):
                    cells.append(f"{lab} {st.median(gs):6.2f} 損{sum(g<1 for g in gs)/len(gs):4.0%}")
            print(f"    r={r_:2d}: " + "  ".join(cells[:4]) + "\n          " + "  ".join(cells[4:]))

# ---------------- [D] ブロックの電荷の偶奇 ----------------
print("\n[D] ブロックの電荷の偶奇 P_ℓ(1次元 L=64、始点 17)")
par = collections.defaultdict(dict)
for fn in glob.glob("crm_eigparity_W1L64_cylinder_U*.tsv"):
    for r in csv.DictReader(open(fn), delimiter="\t"):
        if r["kind"] == "P" and r["start"] == "17":
            par[(float(r["U"]), int(r["l"]))][r["prior"]] = (float(r["true"]), float(r["prior_val"]))
for U in (8.0, 12.0):
    e, Gu, Gd, Gr, ph = models(1, "cylinder", 64, U)
    for l in (1, 2, 4, 8, 16, 32):
        x, yu = par[(U, l)]["UHF"]; yp = ph.P_block(list(range(16, 16+l)))
        cells = [f"UHF {gain(x, yu, 2*l):6.1f}", f"PHF {gain(x, yp, 2*l):6.1f}"]
        for c in CH:
            if f"chi{c}" in par[(U, l)]:
                cells.append(f"切断χ{c} {gain(x, par[(U, l)][f'chi{c}'][1], 2*l):6.1f}")
        for c in CH:
            y = MV.get((1, "cylinder", 64, U, "var", c), {}).get(("P", 17, l))
            if y is not None:
                cells.append(f"変分χ{c} {gain(x, y, 2*l):6.1f}")
        print(f"  U={U:.0f} ℓ={l:2d}: " + "  ".join(cells))

# ---------------- [E] 変分エネルギー ----------------
print("\n[E] 1サイトあたりのエネルギー(1次元)。回収率 = UHF から厳密な値までの相関エネルギーのうち取り戻した割合")
def energy(ph_or_G, e, n, U, kind):
    if kind == "PHF":
        kin = sum(-2*ph_or_G.hop(a-1, b-1) for a, b in e)
        pot = sum(U*(ph_or_G.ZupZdn_onsite(i) - 1 + 2*(1 - ph_or_G.one_minus_n(i)))/4 for i in range(n))
    else:
        Gu, Gd = ph_or_G
        kin = sum(-(Gu[a-1, b-1] + Gu[b-1, a-1] + Gd[a-1, b-1] + Gd[b-1, a-1]) for a, b in e)
        pot = sum(U*Gu[i, i]*Gd[i, i] for i in range(n))
    return (kin + pot)/n
for (geo, LX) in (("cylinder", 8), ("cylinder", 16), ("cylinder", 64), ("torus", 16)):
    for U in (4.0, 8.0, 12.0):
        k = (1, geo, LX, U); e, Gu, Gd, Gr, ph = models(*k); n = LX
        eu = energy((Gu, Gd), e, n, U, "UHF"); ep = energy(ph, e, n, U, "PHF"); ex = Eref[k]/n
        i, _ = center(1, LX, e); m = abs(Gu[i, i] - Gd[i, i])/2
        print(f"  {geo} L={LX:2d} U={U:4.0f}: 厳密 {ex:+.5f}  UHF {eu:+.5f}  PHF {ep:+.5f}   相関エネルギーのうち射影で得た割合 {(ep-eu)/(ex-eu):.3f}   "
              f"L×(E_PHF-E_UHF) {n*(ep-eu):+.3f}  回転子の予測 -2mJ {-2*m*4/U:+.3f}")
        cells = []
        for m, mm in (("切断", "trunc"), ("変分", "var")):
            for c in CH:
                d = MV.get((1, geo, LX, U, mm, c))
                if d:
                    cells.append(f"{m}χ{c} {(d[('E', 0, 0)]/n - eu)/(ex - eu):+.3f}")
        if cells:
            print("      回収率: " + "  ".join(cells))
