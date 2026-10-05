"""変分 MPS をスピンの向きについて平均した prior(変分 MPS-sym)の利得を、UHF・UHF-sym・変分 MPS と比べる。

  変分 MPS(結合次元 χ に制限した DMRG、crm_mps_lowchi.jl)は電荷の量では良い prior だが、スピンの向きを選ぶ
  (README_phf_mps.md の要点4)。UHF-sym と同じ回転平均 σ̄ = ∫dΩ R σ R† を掛けると、
      <Z_{i↑}>       -> 1 - <n_i>
      <Z_{i↑}Z_{j↑}> -> <(1-n_i)(1-n_j)> + (4/3)<S_i·S_j>
      <Z_{i↑}Z_{i↓}>, <n_i n_j>, <S_i·S_j> は不変
  回転平均に要る量は clara で crm_mps_lowchi_sym.jl が保存済みの状態から求めた(crm_mpslow_sym.tsv)。

  真の値は統合表 crm_edfid_all.tsv の厳密な 80 系(中央のサイトと隣)。利得は n_m=100。
  二重占有は 4 項の和として測るときの厳密な分散比(crm_pauli_sum_variance.py)。
"""
import csv, collections, glob, statistics as st
from crm_eigen_uhfsym import lat_edges, solve_uhf, wick, gain
from crm_pauli_sum_variance import variances

NM = 100
CH = (4, 8, 16)

# 真の値(厳密な 80 系)
V = collections.defaultdict(dict)
for r in csv.DictReader(open("crm_edfid_all.tsv"), delimiter="\t"):
    if r["exact"] == "yes":
        V[(int(r["W"]), r["geometry"], int(r["LX"]), float(r["U"]), r["prior"])][r["observable"]] = (float(r["true"]), float(r["prior_val"]))
systems = sorted({k[:4] for k in V})

# 変分 MPS の中央の値(crm_mps_lowchi.jl の出力)
MV = collections.defaultdict(dict)
for fn in sorted(glob.glob("crm_mpslow_W*L*_*_U*.tsv")):
    for r in csv.DictReader(open(fn), delimiter="\t"):
        if r["method"] == "var":
            MV[(int(r["W"]), r["geometry"], int(r["LX"]), float(r["U"]), int(r["chi"]))][r["quantity"]] = float(r["value"])

# 回転平均に要る量(crm_mps_lowchi_sym.jl の出力)
S = {}
for r in csv.DictReader(open("crm_mpslow_sym.tsv"), delimiter="\t"):
    S[(int(r["W"]), r["geometry"], int(r["LX"]), float(r["U"]), int(r["chi"]))] = {k: float(v) for k, v in r.items() if k not in ("geometry",)}

def center(W, LX, e):
    c0 = (max(1, -(-LX//2))-1)*W + max(1, -(-W//2)); nb = next(b for (a, b) in e if a == c0)
    return c0 - 1, nb - 1

# ---------------- [0] 照合 ----------------
w = 0.0; nchk = 0
for k, s in S.items():
    d = MV.get(k)
    if d is None: continue
    w = max(w, abs(s["ZZ_onsite"] - d["ZZ onsite"]), abs(s["ZupZup_nb"] - d["ZZ up-up nb"]), abs(s["Sz_i"] - d["Sz"]),
            abs(s["SS"] - (d["SzSz nb"] + 2*d["SxSx nb"])))
    nchk += 1
print(f"[0] crm_mpslow_sym.tsv と crm_mpslow_*.tsv の照合: {nchk} 状態、オンサイト ZZ・隣接 Z↑Z↑・S^z・S·S の最大差 {w:.1e}")

# ---------------- [A] 利得 ----------------
res = collections.defaultdict(list)
for k in systems:
    W, geo, LX, U = k
    e = lat_edges(LX, W, geo == "torus"); n = W*LX
    i, j = center(W, LX, e)
    d = V[k + ("UHF",)]
    x_on, x_nb = d["ZZ onsite"][0], d["ZZ up-up nb"][0]
    x_z = 1 - d["n"][0] - 2*d["Sz"][0]
    ex = {"II": 1.0, "ZI": x_z, "IZ": 1 - d["n"][0] + 2*d["Sz"][0], "ZZ": x_on}
    Gu, Gd = solve_uhf(e, n, LX, W, U, n//2, n//2)
    wu = wick(Gu, Gd, i, j)
    pri = {"UHF": dict(on=wu["zz_on"], nb=wu["zz_uu"], z=1-2*Gu[i, i], zd=1-2*Gd[i, i]),
           "UHF-sym": dict(on=wu["zz_on"], nb=wu["zz_uu_sym"], z=1-Gu[i, i]-Gd[i, i], zd=1-Gu[i, i]-Gd[i, i])}
    for c in CH:
        mv = MV.get(k + (c,)); sv = S.get(k + (c,))
        if mv is None or sv is None: continue
        pri[f"変分χ{c}"] = dict(on=mv["ZZ onsite"], nb=mv["ZZ up-up nb"],
                               z=1-mv["n"]-2*mv["Sz"], zd=1-mv["n"]+2*mv["Sz"])
        cc = 1 - sv["n_i"] - sv["n_j"] + sv["nn"]
        pri[f"変分χ{c}-sym"] = dict(on=sv["ZZ_onsite"], nb=cc + (4/3)*sv["SS"], z=1-sv["n_i"], zd=1-sv["n_i"])
    for lab, p in pri.items():
        res[("オンサイト ZZ", U, lab)].append(gain(x_on, p["on"], 2))
        res[("隣接 Z↑Z↑", U, lab)].append(gain(x_nb, p["nb"], 2))
        res[("単一サイト Z↑", U, lab)].append(gain(x_z, p["z"], 1))
        ey = {"II": 1.0, "ZI": p["z"], "IZ": p["zd"], "ZZ": p["on"]}
        vs, vc = variances([(-0.25, "ZI"), (-0.25, "IZ"), (0.25, "ZZ")], ex, ey, NM)
        res[("二重占有 4項", U, lab)].append(vs/vc)

def cell(v):
    if not v: return "—"
    return f"{st.median(v):7.2f}[{min(v):6.2f}] 損{sum(g < 0.99 for g in v):2d}/{len(v)}"
LABS = ["UHF", "UHF-sym"] + [f"変分χ{c}{s}" for c in CH for s in ("", "-sym")]
print("\n[A] 厳密な 80 系(中央値 [最小]、1% より大きく損をする系の数 / 系の数)。n_m=100")
for obs in ("オンサイト ZZ", "隣接 Z↑Z↑", "単一サイト Z↑", "二重占有 4項"):
    print(f"  {obs}")
    for U in (2.0, 4.0, 8.0, 12.0):
        print(f"    U={U:4.0f}: " + "  ".join(f"{lab} {cell(res[(obs, U, lab)])}" for lab in LABS[:4]))
        print("            " + "  ".join(f"{lab} {cell(res[(obs, U, lab)])}" for lab in LABS[4:]))

print("\n[B] 全 U をまとめた損の数(1% より大きく損をする系 / 系の数)")
for lab in LABS:
    cells = []
    for obs in ("オンサイト ZZ", "隣接 Z↑Z↑", "単一サイト Z↑", "二重占有 4項"):
        v = [g for U in (2.0, 4.0, 8.0, 12.0) for g in res[(obs, U, lab)]]
        cells.append(f"{obs} {sum(g < 0.99 for g in v)}/{len(v)}")
    print(f"  {lab:12s} " + "  ".join(cells))

# ---------------- [C] 隣の <S_i·S_j> の相対誤差 ----------------
# 回転平均は <S_i·S_j> を変えないので、隣接 Z↑Z↑ の利得を決めるのは prior の <S_i·S_j> の正確さである。
print("\n[C] 隣の <S_i·S_j> の相対誤差 |prior/真 - 1|(中央値、厳密な 80 系のうち各 U の 20 系)")
T = {k: 3*V[k + ("UHF",)]["SzSz nb"][0] for k in systems}      # 真の状態はスピン回転不変なので S·S = 3<SzSz>
for U in (2.0, 4.0, 8.0, 12.0):
    cells = []
    eu = []
    for k in systems:
        if k[3] != U: continue
        W, geo, LX, _ = k; e = lat_edges(LX, W, geo == "torus"); n = W*LX; i, j = center(W, LX, e)
        Gu, Gd = solve_uhf(e, n, LX, W, U, n//2, n//2)
        eu.append(abs(wick(Gu, Gd, i, j)["ss"]/T[k] - 1))
    cells.append(f"UHF(=UHF-sym) {st.median(eu):.3f}")
    for c in CH:
        ev = [abs(S[k + (c,)]["SS"]/T[k] - 1) for k in systems if k[3] == U and k + (c,) in S]
        cells.append(f"変分χ{c} {st.median(ev):.3f}")
    print(f"  U={U:4.0f}: " + "  ".join(cells))
