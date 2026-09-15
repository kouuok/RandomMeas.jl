"""ブロックの偶奇の文字列(crm_eigen_parity.jl の出力)で、利得の大きい観測量の特徴を確かめる。

  P_ℓ = (-1)^{N_B}   ブロック内部の二重占有・空サイトの対は偶奇を変えず、境界をまたぐ対だけが変える
                     → ブロックを長くしても固有状態からの漏れ η は増えないはず。スピン回転で不変。
  Q_ℓ = (-1)^{N↑_B}  スピンのもつれで揺らぎ、UHF はネール状態で N↑_B をほぼ固定する → 過信するはず。

Slater 行列式では <Π_{i∈B}(1-2n_iσ)> = det(I - 2G^σ_B) なので、UHF と RHF の値は軌道から直接出る。
"""
import csv, glob, collections, statistics as st
import numpy as np
from crm_eigen_uhfsym import lat_edges, solve_uhf, rhf_G, gain

R = collections.defaultdict(dict)   # (LX, geo, U, kind, start, l) -> {prior: (x, y)}
meta = {}
for fn in sorted(glob.glob("crm_eigparity_W1L*_*_U*.tsv")):
    for r in csv.DictReader(open(fn), delimiter="\t"):
        k = (int(r["LX"]), r["geometry"], float(r["U"]), r["kind"], int(r["start"]), int(r["l"]))
        R[k][r["prior"]] = (float(r["true"]), float(r["prior_val"])); meta[k[:3]] = float(r["Hvar"])
systems = sorted(meta)
print("系:", [(s, f"{meta[s]:.0e}") for s in systems])

# ---- 検証: UHF の MPS の値と Wick(行列式)の値 ------------------------------------------------
worst = 0.0; WICK = {}
for (LX, geo, U) in systems:
    n = LX; e = lat_edges(LX, 1, geo == "torus")
    Gu, Gd = solve_uhf(e, n, LX, 1, U, n//2, n//2); Gr = rhf_G(e, n, n//2)
    for k, d in R.items():
        if k[:3] != (LX, geo, U): continue
        kind, s0, l = k[3:]; B = list(range(s0-1, s0-1+l)); I = np.eye(l)
        du, dd = np.linalg.det(I - 2*Gu[np.ix_(B, B)]), np.linalg.det(I - 2*Gd[np.ix_(B, B)])
        dr = np.linalg.det(I - 2*Gr[np.ix_(B, B)])
        wu = du*dd if kind == "P" else du
        wr = dr*dr if kind == "P" else dr
        WICK[k] = (wu, wr)
        worst = max(worst, abs(wu - d["UHF"][1]))
print(f"UHF: MPS と Wick(行列式)の最大差 {worst:.1e}")

def row(k, p):
    x = R[k]["UHF"][0]
    y = WICK[k][1] if p == "RHF" else R[k][p][1]
    nA = 2*k[5] if k[3] == "P" else k[5]
    return x, y, gain(x, y, nA)

for kind, title in (("P", "電荷の偶奇 P_ℓ = (-1)^{N_B}(台 2ℓ 量子ビット)"), ("Q", "上向きの偶奇 Q_ℓ = (-1)^{N↑_B}(台 ℓ 量子ビット)")):
    print(f"\n==== {title} ====")
    for (LX, geo, U) in systems:
        if not (LX in (64, 16)): continue
        starts = sorted({k[4] for k in R if k[:4] == (LX, geo, U, kind)})
        s0 = starts[0]
        print(f"  {geo} L={LX} U={U:.0f}  始点 {s0}")
        print(f"    {'ℓ':>2s} {'|x|':>6s} {'η_ρ':>6s} {'G_max':>8s} | {'UHF: y/x':>8s} {'δ':>6s} {'G':>8s} | {'RHF: G':>7s} | {'χ4: G':>7s} {'χ8: G':>7s}")
        for l in sorted({k[5] for k in R if k[:5] == (LX, geo, U, kind, s0)}):
            k = (LX, geo, U, kind, s0, l)
            if l > 4 and l % 4 and l != max({kk[5] for kk in R if kk[:5] == (LX, geo, U, kind, s0)}): continue
            x, y, g = row(k, "UHF"); nA = 2*l if kind == "P" else l
            er = (1-abs(x))/2; es = (1-abs(y))/2; dl = (es-er)/er if er > 1e-12 else float("nan")
            _, _, gr = row(k, "RHF")
            g4 = row(k, "chi4")[2] if "chi4" in R[k] else float("nan"); g8 = row(k, "chi8")[2] if "chi8" in R[k] else float("nan")
            print(f"    {l:2d} {abs(x):6.3f} {er:6.3f} {gain(x,x,nA):8.1f} | {y/x if abs(x)>1e-3 else float("nan"):8.2f} {dl:+6.2f} {g:8.2f} | {gr:7.2f} | {g4:7.2f} {g8:7.2f}")
