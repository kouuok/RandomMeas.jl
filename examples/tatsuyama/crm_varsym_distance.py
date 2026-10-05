"""回転平均した prior は遠くのスピン相関でも効くか、と「連結相関の prior」の効果。

  観測量: 距離 r の Z↑Z↑ (= <(1-n_i)(1-n_j)> + 4<S^z_i S^z_j> + (S^z の 1 次の項))。台は 2 量子ビット、n_m=100。
  prior:
    UHF, 切断 MPS χ=4/8/16(crm_eigcorr_*.tsv の値)
    UHF-sym       回転平均: <(1-n_i)(1-n_j)> + (4/3)<S_i·S_j>
    UHF-conn      連結相関: <(1-n_i)(1-n_j)> + (4/3)(<S_i·S_j> - <S_i>·<S_j>)
    変分χ, 変分χ-sym, 変分χ-conn(crm_varsym_pairs.tsv、clara で保存済みの状態から)
  「-conn」は物理的な状態の期待値ではない。CRM はどんな数値を prior に使っても偏らないので、それでもよい。

  系: 1 次元 L=64 開放端(端の影響を避けて 17 ≤ i < j ≤ 48 の対)と L=16 周期境界(全対)。
"""
import csv, collections, statistics as st
from crm_eigen_uhfsym import lat_edges, solve_uhf, wick, gain

CH = (4, 8, 16)
SYS = [("cylinder", 64), ("torus", 16)]
US = (2.0, 4.0, 8.0, 12.0)
RBINS = [(1, 1), (2, 2), (3, 3), (4, 4), (5, 8), (9, 16), (17, 31)]

def dist(geo, L, i, j):
    r = abs(j - i)
    return min(r, L - r) if geo == "torus" else r

def in_bulk(geo, L, i, j):
    return True if geo == "torus" else (17 <= min(i, j) and max(i, j) <= 48)

# 真の値と、UHF・切断 MPS の値
T = collections.defaultdict(dict)
for geo, L in SYS:
    for U in US:
        for r in csv.DictReader(open(f"crm_eigcorr_W1L{L}_{geo}_U{U}.tsv"), delimiter="\t"):
            if r["kind"] != "ZupZup": continue
            i, j = int(r["i"]), int(r["j"])
            if not in_bulk(geo, L, i, j): continue
            d = T[(geo, L, U, i, j)]
            d["true"] = float(r["true"]); d[r["prior"]] = float(r["prior_val"])

# 変分 MPS の対の量
VP = {}
for r in csv.DictReader(open("crm_varsym_pairs.tsv"), delimiter="\t"):
    VP[(r["geometry"], int(r["LX"]), float(r["U"]), int(r["chi"]), int(r["i"]), int(r["j"]))] = {k: float(v) for k, v in r.items() if k != "geometry"}

# [0] 照合: crm_varsym_pairs.tsv の Z↑Z↑ を、以前の計算 crm_mpslow_*.tsv の全対の値と比べる
import glob
MV = {}
for fn in glob.glob("crm_mpslow_W1L64_cylinder_U*.tsv") + glob.glob("crm_mpslow_W1L16_torus_U*.tsv"):
    for r in csv.DictReader(open(fn), delimiter="\t"):
        if r["method"] == "var" and r["quantity"] == "ZupZup":
            MV[(r["geometry"], int(r["LX"]), float(r["U"]), int(r["chi"]), int(r["i"]), int(r["j"]))] = float(r["value"])
w0 = max((abs(v["ZupZup"] - MV[k]) for k, v in VP.items() if k in MV), default=float("nan"))
n0 = sum(k in MV for k in VP)
wn = max(abs(v["n_i"] - 1) for v in VP.values())
print(f"[0] 変分 MPS の対 {len(VP)} 個。以前の Z↑Z↑ との照合 {n0} 対で最大差 {w0:.1e}。max|n_i - 1| = {wn:.1e}")

res = collections.defaultdict(list)
for geo, L in SYS:
    e = lat_edges(L, 1, geo == "torus")
    for U in US:
        Gu, Gd = solve_uhf(e, L, L, 1, U, L//2, L//2)
        for (g_, L_, U_, i, j), d in T.items():
            if (g_, L_, U_) != (geo, L, U) or "true" not in d: continue
            x = d["true"]; r = dist(geo, L, i, j)
            w = wick(Gu, Gd, i-1, j-1)
            szi = 0.5*(Gu[i-1, i-1] - Gd[i-1, i-1]); szj = 0.5*(Gu[j-1, j-1] - Gd[j-1, j-1])
            pri = {"UHF": d.get("UHF"), "UHF-sym": w["zz_uu_sym"],
                   "UHF-conn": w["cc"] + (4/3)*(w["ss"] - szi*szj)}
            for c in CH:
                pri[f"切断χ{c}"] = d.get(f"chi{c}")
                v = VP.get((geo, L, U, c, i, j))
                if v is None: continue
                cc = 1 - v["n_i"] - v["n_j"] + v["nn"]
                pri[f"変分χ{c}"] = v["ZupZup"]
                pri[f"変分χ{c}-sym"] = cc + (4/3)*v["SS"]
                pri[f"変分χ{c}-conn"] = cc + (4/3)*(v["SS"] - v["Sz_i"]*v["Sz_j"])
            for lab, y in pri.items():
                if y is None: continue
                res[(geo, L, U, r, lab)].append(gain(x, y, 2))

LABS = ["UHF", "UHF-sym", "UHF-conn"] + [f"切断χ{c}" for c in CH] + \
       [f"変分χ{c}{s}" for c in CH for s in ("", "-sym", "-conn")]

def cell(v):
    if not v: return "   —   "
    return f"{st.median(v):6.2f}/{sum(g < 0.99 for g in v)/len(v):4.0%}"

for geo, L in SYS:
    for U in US:
        print(f"\n[A] {geo} L={L} U={U:.0f}: 距離 r の Z↑Z↑ の利得の中央値 / 1% より大きく損をする対の割合(n_m=100)")
        rb = [b for b in RBINS if b[0] <= (L//2 if geo == "torus" else 31)]
        print("  " + f"{'prior':14s}" + "".join(f"{('r='+str(a) if a == b else f'r={a}-{b}'):>14s}" for a, b in rb))
        for lab in LABS:
            cells = []
            for a, b in rb:
                v = [g for r in range(a, b+1) for g in res.get((geo, L, U, r, lab), [])]
                cells.append(f"{cell(v):>14s}")
            if any(c.strip() != "—" for c in cells):
                print(f"  {lab:14s}" + "".join(cells))

# [B] 全距離をまとめた損の割合(U ごと、r ≥ 2)
print("\n[B] r ≥ 2 の全対での損の割合(L=64 開放端 / L=16 周期境界)")
for U in US:
    cells = []
    for lab in LABS:
        parts = []
        for geo, L in SYS:
            v = [g for (g_, L_, U_, r, l), gs in res.items() if (g_, L_, U_, l) == (geo, L, U, lab) and r >= 2 for g in gs]
            parts.append(f"{sum(g < 0.99 for g in v)/len(v):4.0%}" if v else "  — ")
        cells.append(f"{lab} {'/'.join(parts)}")
    print(f"  U={U:4.0f}: " + "  ".join(cells))
