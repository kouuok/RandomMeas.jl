"""距離 r の Z↑Z↑ 相関で「固有状態への近さ」と利得の関係を連続的に確かめる。

入力: crm_eigcorr_W*L*_*_U*.tsv(crm_eigen_corr.jl、clara)
  1. 検証: 中央の値が統合表と一致するか / UHF の全サイト対が Python の Wick(移植した UHF)と一致するか
  2. UHF は長距離秩序を持つので |y| はどの距離でも大きいまま。真の |x| は減衰する。
     1.4 の境目 |x| = |y|/2 を距離のどこで横切るか。
  3. SU(2) 群平均した UHF(y_sym = <c c> + 4/3 <S·S>)は |y|≈1/3 で、|x|<1/6 から損をすると予想した。
"""
import csv, glob, collections, math, statistics as st, re
import numpy as np
from crm_eigen_uhfsym import lat_edges, solve_uhf, wick, gain

def load():
    D = collections.defaultdict(dict)       # (W,geo,LX,U) -> {(prior,kind,i,j): (x,y)}
    meta = {}
    for fn in sorted(glob.glob("crm_eigcorr_W*L*_*_U*.tsv")):
        for r in csv.DictReader(open(fn), delimiter="\t"):
            k = (int(r["W"]), r["geometry"], int(r["LX"]), float(r["U"]))
            D[k][(r["prior"], r["kind"], int(r["i"]), int(r["j"]))] = (float(r["true"]), float(r["prior_val"]))
            meta[k] = (float(r["Hvar"]), int(r["chi_ref"]))
    return D, meta

if __name__ == "__main__":
    D, meta = load()
    VARTOL = 1e-8          # 厳密な ρ とみなすエネルギー分散の上限(統合表と同じ基準)
    print("読み込んだ系:", len(D))
    for k in sorted(D):
        bad = meta[k][0] > VARTOL
        print(f"  {k}  分散 {meta[k][0]:.1e}  χ {meta[k][1]}" + ("   ← 厳密でないので除外" if bad else ""))
        if bad: del D[k]
    print("解析に使う系:", len(D))

    # ---- 1a. 統合表との照合(中央のサイトと結合) ---------------------------------
    M = {}
    for r in csv.DictReader(open("crm_edfid_all.tsv"), delimiter="\t"):
        if r["observable"] in ("ZZ onsite", "ZZ up-up nb"):
            M[(int(r["W"]), r["geometry"], int(r["LX"]), float(r["U"]), r["prior"], r["observable"])] = (float(r["true"]), float(r["prior_val"]))
    worst = 0.0; nchk = 0
    for k, d in D.items():
        W, geo, LX, U = k
        e = lat_edges(LX, W, geo == "torus")
        c0 = (max(1, -(-LX//2))-1)*W + max(1, -(-W//2)); nb = next(b for (a, b) in e if a == c0)
        for (p, kind, i, j), (x, y) in d.items():
            for obs, (ii, jj, kk) in (("ZZ onsite", (c0, c0, "ZupZdn")), ("ZZ up-up nb", (min(c0,nb), max(c0,nb), "ZupZup"))):
                if (kind, i, j) == (kk, ii, jj) and (k + (p, obs)) in M:
                    mx, my = M[k + (p, obs)]; worst = max(worst, abs(x-mx), abs(y-my)); nchk += 1
    print(f"\n[1a] 統合表との照合: {nchk} 値、最大差 {worst:.1e}")

    # ---- 1b. UHF の全サイト対を Wick と照合、ついでに y_sym を作る ----------------------
    SYM = {}; worst_all = collections.defaultdict(float)
    for k, d in D.items():
        W, geo, LX, U = k; n = W*LX
        Gu, Gd = solve_uhf(lat_edges(LX, W, geo == "torus"), n, LX, W, U, n//2, n//2)
        for (p, kind, i, j), (x, y) in d.items():
            if p != "UHF": continue
            if kind == "Zup":
                wz = 1 - 2*Gu[i-1, i-1]; worst_all[kind] = max(worst_all[kind], abs(wz - y))
            elif kind == "ZupZup":
                w = wick(Gu, Gd, i-1, j-1); worst_all[kind] = max(worst_all[kind], abs(w["zz_uu"] - y))
                SYM[(k, i, j)] = w["zz_uu_sym"]
            elif kind == "ZupZdn":
                wz = (1-2*Gu[i-1, i-1])*(1-2*Gd[j-1, j-1]) - (0 if i == j else 0)   # ↑↓ は無相関(同一サイトも異サイトも)
                worst_all[kind] = max(worst_all[kind], abs(wz - y))
    print("[1b] UHF 全サイト対: MPS と Wick の最大差", {kk: f"{v:.1e}" for kk, v in worst_all.items()})

    # ---- 2. 距離 r ごとの利得 -------------------------------------------------------
    def pairs_by_r(k, d):
        W, geo, LX, U = k
        out = collections.defaultdict(list)
        for (p, kind, i, j) in d:
            if p != "UHF" or kind != "ZupZup": continue
            xi, yi = divmod(i-1, W); xj, yj = divmod(j-1, W)
            if W == 1:
                lo, hi = LX//4, 3*LX//4                         # 開放端の影響を避けて中央の半分だけ使う
                if not (lo <= xi < hi and lo <= xj < hi): continue
                r = abs(xi-xj)
            else:
                dy = abs(yi-yj); dy = min(dy, W-dy) if W > 2 else dy
                r = abs(xi-xj) + dy
            out[r].append((i, j))
        return out

    PRI = ["UHF", "UHF-sym", "chi4", "chi8", "chi16"]
    allpts = []                                                 # (系, r, prior, x, y, G)
    for k, d in sorted(D.items()):
        byr = pairs_by_r(k, d)
        print(f"\n[2] {k}   (各 r でサイト対の中央値。損 = G<1 の対の割合)")
        print(f"  {'r':>2s} {'対':>3s} {'|x|':>6s} | " + " | ".join(f"{p:>7s} y/x      G  損" for p in PRI))
        for r in sorted(byr):
            cells = []; xs = []
            for p in PRI:
                ratios, gs = [], []
                for (i, j) in byr[r]:
                    x = d[("UHF", "ZupZup", i, j)][0]
                    if p == "UHF-sym": y = SYM[(k, i, j)]
                    elif (p, "ZupZup", i, j) in d: y = d[(p, "ZupZup", i, j)][1]
                    else: continue
                    g = gain(x, y); ratios.append(y/x); gs.append(g); allpts.append((k, r, p, x, y, g, i, j))
                    if p == "UHF": xs.append(abs(x))
                cells.append(f"{st.median(ratios):5.2f} {st.median(gs):6.2f} {sum(g < 1-1e-9 for g in gs)/len(gs):4.0%}" if gs else " "*17)
            print(f"  {r:2d} {len(byr[r]):3d} {st.median(xs):6.3f} | " + " | ".join(cells))

    # ---- 3. 予想の検定 ------------------------------------------------------------------
    bad = sum((g < 1-1e-9) != (not (0 < y/x < 2)) for (_, _, _, x, y, g, _, _) in allpts if abs(g-1) > 1e-9)
    print(f"\n[3a] 全 {len(allpts)} 点で G<1 ⟺ not(0<y/x<2) の破れ: {bad}")
    for p in ("UHF", "UHF-sym"):
        pts = [(abs(x), abs(y), g) for (_, _, q, x, y, g, _, _) in allpts if q == p and x*y > 0]
        lose = [a for a in pts if a[2] < 1]
        pred = sum((g < 1) == (ax < ay/2) for ax, ay, g in pts)
        opp = sum(1 for (_, _, q, x, y, g, _, _) in allpts if q == p and x*y < 0)
        print(f"[3b] {p:7s}: 同符号 {len(pts)} 対、損 {len(lose)}。'|x|<|y|/2 なら損' の的中 {pred}/{len(pts)}"
              f"  損の対の |y| 中央値 {st.median([a[1] for a in lose]) if lose else float('nan'):.3f}  逆符号 {opp} 対")

    # ---- 4. 図のためのサイト対ごとの表 ----------------------------------------------------
    with open("crm_eigcorr_pairs.tsv", "w") as f:
        f.write("W\tgeometry\tLX\tU\tr\ti\tj\tprior\tx\ty\tG\n")
        for (k, r, p, x, y, g, i, j) in allpts:
            f.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{k[3]:.1f}\t{r}\t{i}\t{j}\t{p}\t{x:.10f}\t{y:.10f}\t{g:.8e}\n")
    print("書き出し: crm_eigcorr_pairs.tsv", len(allpts), "行")
