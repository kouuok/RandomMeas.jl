"""UHF prior を SU(2) で群平均すると、隣接 Z↑Z↑ の「固有状態への寄りすぎ」は解消するか。

Z_{i↑}Z_{j↑} = c_i c_j - 2 c_i S^z_j - 2 S^z_i c_j + 4 S^z_i S^z_j   (c_i = 1 - n_i)
スピン回転の Haar 平均はスカラー c_i c_j を残し、1階テンソルの交差項を 0 にし、
S^z_i S^z_j を (1/3) S_i·S_j に置き換える。したがって群平均した prior の値は厳密に
    y_sym = <c_i c_j>_σ + (4/3) <S_i·S_j>_σ
で、求積はいらない。UHF は Slater 行列式なので右辺は Wick の定理で軌道から直接出る。

crm_2d_ed_fid_table.jl の lat_edges / solve_uhf を移植し、まず統合表の UHF の値
(ZZ onsite, ZZ up-up nb, SzSz nb, SxSx nb, n, Sz)を再現することを確かめてから使う。
"""
import csv, collections, math
import numpy as np

def lat_edges(LX, W, pbc):
    sidx = lambda x, y: (x-1)*W + y
    e = []
    for x in range(1, LX+1):            # Julia の `for x in 1:Lx, y in 1:Wd` は x が外側
        for y in range(1, W+1):
            s = sidx(x, y)
            if x < LX:                e.append(tuple(sorted((s, sidx(x+1, y)))))
            elif pbc and LX > 2:      e.append(tuple(sorted((s, sidx(1, y)))))
            if W > 2:                 e.append(tuple(sorted((s, sidx(x, 1 if y == W else y+1)))))
            elif W == 2 and y == 1:   e.append(tuple(sorted((s, sidx(x, 2)))))
    out = []
    for t in e:
        if t not in out: out.append(t)
    return out

def solve_uhf(edges, n, LX, W, U, Nup, Ndn, t=1.0, iters=8000, tol=1e-13):
    sidx = lambda x, y: (x-1)*W + y
    T = np.zeros((n, n))
    for a, b in edges: T[a-1, b-1] = T[b-1, a-1] = -t
    fa = (Nup+Ndn)/(2*n); nup = np.zeros(n); ndn = np.zeros(n)
    for x in range(1, LX+1):
        for y in range(1, W+1):
            s = sidx(x, y)-1; g = 1 if (x+y) % 2 == 1 else -1
            nup[s] = fa + 0.4*g; ndn[s] = fa - 0.4*g
    nup = np.clip(nup, 0.02, 0.98); ndn = np.clip(ndn, 0.02, 0.98)
    nup *= Nup/nup.sum(); ndn *= Ndn/ndn.sum()
    for _ in range(iters):
        _, Vu = np.linalg.eigh(T + np.diag(U*ndn)); _, Vd = np.linalg.eigh(T + np.diag(U*nup))
        nu = (Vu[:, :Nup]**2).sum(1); nd = (Vd[:, :Ndn]**2).sum(1)
        if max(abs(nu-nup).max(), abs(nd-ndn).max()) < tol:
            nup, ndn = nu, nd; break
        nup = 0.5*nu + 0.5*nup; ndn = 0.5*nd + 0.5*ndn
    return Vu[:, :Nup] @ Vu[:, :Nup].T, Vd[:, :Ndn] @ Vd[:, :Ndn].T     # 相関行列 <c†_i c_j>

def wick(Gu, Gd, i, j):
    """UHF(↑↓ を混ぜない Slater 行列式)での期待値。i≠j は隣接対、添字は 0 始まり。"""
    nu, nd = np.diag(Gu), np.diag(Gd)
    z = lambda k, G: 1 - 2*G[k, k]
    zz_on = z(i, Gu)*z(i, Gd)                                       # 同じサイトの ↑↓ は無相関
    zz_uu = z(i, Gu)*z(j, Gu) - 4*Gu[i, j]**2
    szsz = 0.25*((nu[i]-nd[i])*(nu[j]-nd[j]) - Gu[i, j]**2 - Gd[i, j]**2)
    sxsx = -0.5*Gu[i, j]*Gd[i, j]
    ninj = (nu[i]+nd[i])*(nu[j]+nd[j]) - Gu[i, j]**2 - Gd[i, j]**2
    cc = 1 - (nu[i]+nd[i]) - (nu[j]+nd[j]) + ninj
    return dict(zz_on=zz_on, zz_uu=zz_uu, szsz=szsz, sxsx=sxsx, n=nu[i]+nd[i], sz=0.5*(nu[i]-nd[i]),
                cc=cc, ss=szsz + 2*sxsx, zz_uu_sym=cc + (4/3)*(szsz + 2*sxsx))

def gain(x, y, nA=2, nm=100):
    K = 3.0**nA - 1; vs = 3.0**nA*(1 - x*x)/nm
    return (K*x*x + vs)/(K*(x-y)**2 + vs)

if __name__ == "__main__":
    V = collections.defaultdict(dict)
    for r in csv.DictReader(open("crm_edfid_all.tsv"), delimiter="\t"):
        if r["exact"] == "yes" and r["prior"] == "UHF":
            V[(int(r["W"]), r["geometry"], int(r["LX"]), float(r["U"]))][r["observable"]] = (float(r["true"]), float(r["prior_val"]))
    key = {"ZZ onsite": "zz_on", "ZZ up-up nb": "zz_uu", "SzSz nb": "szsz", "SxSx nb": "sxsx", "n": "n", "Sz": "sz"}
    worst = collections.defaultdict(float); rows = []
    for (W, geo, LX, U), d in sorted(V.items()):
        n = W*LX; edges = lat_edges(LX, W, geo == "torus")
        c0 = (max(1, -(-LX//2)) - 1)*W + max(1, -(-W//2))
        nb = next(b for (a, b) in edges if a == c0)
        Gu, Gd = solve_uhf(edges, n, LX, W, U, n//2, n//2)
        w = wick(Gu, Gd, c0-1, nb-1)
        for obs, k in key.items():
            diff = abs(w[k] - d[obs][1]); worst[obs] = max(worst[obs], diff)
            if diff > 1e-6: print(f"  不一致 {(W,geo,LX,U)} {obs}: Wick {w[k]:+.8f}  表 {d[obs][1]:+.8f}")
        rows.append(((W, geo, LX, U), d, w))
    print("統合表の UHF 値の再現(最大差):", {k: f"{v:.1e}" for k, v in worst.items()})

    # ---- §2i の求積値との照合(L=8 OBC, U=4: UHF-sym の SzSz = -0.07331) -----------
    for (k, d, w) in rows:
        if k == (1, "cylinder", 8, 4.0):
            print(f"§2i 照合 L=8 U=4: (SzSz+2SxSx)/3 = {w['ss']/3:+.5f}  (§2i の球面求積 -0.07331、真値 -0.06980)")

    # ---- 群平均で隣接 Z↑Z↑ の位置と利得がどう変わるか ------------------------------
    cls = lambda x, y: "逆符号" if x*y < 0 else ("過信" if abs(y) > abs(x) else "控えめ")
    print(f"\n{'系':26s} | {'x':>7s} | {'UHF y':>7s} {'y/x':>5s} {'G':>6s} | {'sym y':>7s} {'y/x':>5s} {'G':>6s}")
    cnt = collections.Counter(); lose = collections.Counter(); ratio = []
    for (k, d, w) in rows:
        x = d["ZZ up-up nb"][0]; yu = d["ZZ up-up nb"][1]; ys = w["zz_uu_sym"]
        gu, gs = gain(x, yu), gain(x, ys)
        cnt[(cls(x, yu), cls(x, ys))] += 1; lose["UHF"] += gu < 1 - 1e-9; lose["sym"] += gs < 1 - 1e-9
        ratio.append(gs/gu)
        if k[3] == 12.0 or k == (1, "cylinder", 8, 4.0):
            print(f"{str(k):26s} | {x:+7.4f} | {yu:+7.4f} {yu/x:5.2f} {gu:6.2f} | {ys:+7.4f} {ys/x:5.2f} {gs:6.2f}")
    print("\n分類の遷移 (UHF → UHF-sym):", dict(cnt))
    print(f"損をする系: UHF {lose['UHF']}/80  →  UHF-sym {lose['sym']}/80")
    ratio.sort()
    print(f"G_sym/G_UHF: 中央値 {ratio[len(ratio)//2]:.2f}  最小 {ratio[0]:.2f}  最大 {ratio[-1]:.2f}")
