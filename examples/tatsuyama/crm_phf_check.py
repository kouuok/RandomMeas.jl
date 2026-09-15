"""crm_phf.PHF を、フォック空間での厳密な S=0 射影と照合する(小さい系)。"""
import numpy as np
from crm_phf import PHF, solve_uhf_orbitals
from crm_eigen_uhfsym import lat_edges

def fock_phf(Pu, Pd, n):
    Nu, Nd = Pu.shape[1], Pd.shape[1]
    ups = [m for m in range(1 << n) if bin(m).count("1") == Nu]; dns = [m for m in range(1 << n) if bin(m).count("1") == Nd]
    basis = [(u, d) for u in ups for d in dns]; idx = {b: k for k, b in enumerate(basis)}
    occ = lambda m: [i for i in range(n) if (m >> i) & 1]
    psi = np.array([np.linalg.det(Pu[occ(u)]) * np.linalg.det(Pd[occ(d)]) for (u, d) in basis])
    ups1 = [m for m in range(1 << n) if bin(m).count("1") == Nu+1]; dns1 = [m for m in range(1 << n) if bin(m).count("1") == Nd-1]
    basis1 = [(u, d) for u in ups1 for d in dns1]; idx1 = {b: k for k, b in enumerate(basis1)}
    Sp = np.zeros((len(basis1), len(basis)))
    for k, (u, d) in enumerate(basis):
        for i in range(n):
            if (d >> i) & 1 and not (u >> i) & 1:
                sgn = (-1)**(bin(u).count("1") + bin(d & ((1 << i)-1)).count("1") + bin(u & ((1 << i)-1)).count("1"))
                Sp[idx1[(u | (1 << i), d & ~(1 << i))], k] += sgn
    ev, V = np.linalg.eigh(Sp.T @ Sp); P0 = V[:, ev < 1e-8]
    phi = P0 @ (P0.T @ psi); phi /= np.linalg.norm(phi)
    def z(modes):
        vals = np.array([np.prod([1 - 2*(((u >> p) & 1) if p < n else ((d >> (p-n)) & 1)) for p in modes]) for (u, d) in basis])
        return float(np.sum(vals*phi**2))
    def hop_up(i, j):
        out = np.zeros(len(basis))
        for k, (u, d) in enumerate(basis):
            if (u >> j) & 1 and not (u >> i) & 1:
                lo, hi = sorted((i, j)); between = bin(u & (((1 << hi)-1) & ~((1 << (lo+1))-1))).count("1")
                out[idx[((u & ~(1 << j)) | (1 << i), d)]] += (-1)**between*phi[k]
        return float(2*phi @ out)
    return z, hop_up

for (n, U) in ((6, 4.0), (8, 8.0), (8, 2.0)):
    e = lat_edges(n, 1, False); Pu, Pd = solve_uhf_orbitals(e, n, n, 1, U, n//2, n//2)
    z, h = fock_phf(Pu, Pd, n); i, j, k = n//2 - 1, n//2, n//2 + 1
    print(f"L={n} U={U}")
    for nb in (48, 96, 192):
        ph = PHF(Pu, Pd, nbeta=nb)
        rows = [("Z_i↑", ph.Zup(i), z([i])), ("Z_i↑Z_i↓", ph.ZupZdn_onsite(i), z([i, n+i])),
                ("Z_i↑Z_j↑(隣)", ph.ZupZup(i, j), z([i, j])), ("Z_i↑Z_k↑(距離2)", ph.ZupZup(i-1, k), z([i-1, k])),
                ("F_iF_j", ph.FF(i, j), z([i, n+i, j, n+j])), ("P_3", ph.P_block([i-1, i, j]), z([i-1, n+i-1, i, n+i, j, n+j])),
                ("hop↑(隣)", ph.hop(i, j), h(i, j))]
        print(f"  nβ={nb:3d}(使った β {ph.nbeta_used}): 最大差 {max(abs(a-b) for _, a, b in rows):.1e}")
    for lab, a, b in rows: print(f"    {lab:16s} 行列式 {a:+.10f}  厳密 {b:+.10f}")
