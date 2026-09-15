"""スピン射影 HF(PHF): UHF の行列式を全スピン 0 に射影した状態の期待値。

UHF(共線、N↑ = N↓)は S^z_tot = 0 の固有状態である。スピン回転で変わらない(スカラーの)観測量 O なら
    <O>_PHF = <Φ|O P|Φ> / <Φ|P|Φ>,   P|Φ> = (1/2) ∫_0^π sin β  e^{-iβS_y}|Φ> dβ
で厳密に求まる(左の <Φ| と O が S_z を保つので、z 軸まわりの回転の積分は自明になる)。
スカラーでない観測量は、全スピン 0 の状態では「回転で変わらない部分」しか期待値に効かないので、それに置き換える:
    Z_{i↑} → 1 - n_i ,   Z_{i↑}Z_{j↑} → <c_i c_j> + (4/3)<S_i·S_j> ,   c†_{i↑}c_{j↑}+h.c. → 上向きと下向きの平均
遷移行列要素は一般化 Wick の定理: g(p,q) = <Φ|c†_p c_q|Φ_β>/<Φ|Φ_β> = [Φ_β (Φ^T Φ_β)^{-1} Φ^T]_{qp}
重なりが <Φ|Φ> の 1e-14 未満の β は捨てる(局所的な量への寄与はその程度に小さい)。
スピン軌道の並びは (1↑..n↑, 1↓..n↓)。
"""
import numpy as np

def solve_uhf_orbitals(edges, n, LX, W, U, Nup, Ndn, t=1.0, iters=8000, tol=1e-13):
    """crm_eigen_uhfsym.solve_uhf と同じ反復で、占有軌道そのものを返す。"""
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
            break
        nup = 0.5*nu + 0.5*nup; ndn = 0.5*nd + 0.5*ndn
    return Vu[:, :Nup], Vd[:, :Ndn]

class PHF:
    def __init__(self, Phi_up, Phi_dn, nbeta=96, rel=1e-14):
        n = Phi_up.shape[0]; self.n = n
        Nu, Nd = Phi_up.shape[1], Phi_dn.shape[1]
        Phi = np.zeros((2*n, Nu + Nd)); Phi[:n, :Nu] = Phi_up; Phi[n:, Nu:] = Phi_dn
        x, w = np.polynomial.legendre.leggauss(nbeta)            # x = cos β(重み sin β dβ)
        self.G = []; ws = []; lds = []; sgs = []
        for xb, wb in zip(x, w):
            b = np.arccos(xb); c, s = np.cos(b/2), np.sin(b/2)
            Pb = np.vstack([c*Phi[:n] - s*Phi[n:], s*Phi[:n] + c*Phi[n:]])     # e^{-iβS_y}|Φ>
            M = Phi.T @ Pb
            sg, ld = np.linalg.slogdet(M)
            if sg == 0: continue
            self.G.append((Pb @ np.linalg.solve(M, Phi.T)).T)    # g[p, q] = <c†_p c_q>_tr
            ws.append(wb/2); lds.append(ld); sgs.append(sg)
        lds = np.array(lds); keep = lds - lds.max() > np.log(rel)
        self.G = [g for g, k in zip(self.G, keep) if k]
        self.wt = np.array(ws)[keep]*np.array(sgs)[keep]*np.exp(lds[keep] - lds.max())
        self.wt /= self.wt.sum()                                  # <Φ|P|Φ> で規格化した重み
        self.nbeta_used = int(keep.sum())

    def up(self, i): return i
    def dn(self, i): return self.n + i

    def _avg(self, f):
        return float(np.sum([w*f(g) for w, g in zip(self.wt, self.G)]))

    # ---- スカラーの量 ----
    def zstring(self, modes):
        """Π_{p∈modes}(1-2n_p)。スカラーな組(サイトごとに ↑↓ がそろう電荷の偶奇など)に使う。"""
        m = list(modes)
        return self._avg(lambda g: np.linalg.det(np.eye(len(m)) - 2*g[np.ix_(m, m)]))

    def one_minus_n(self, i):
        return self._avg(lambda g: 1 - g[self.up(i), self.up(i)] - g[self.dn(i), self.dn(i)])

    def cc(self, i, j):
        """<(1-n_i)(1-n_j)>(i≠j)。"""
        u, d = self.up, self.dn
        def f(g):
            ops = {}
            Z = lambda S: np.linalg.det(np.eye(len(S)) - 2*g[np.ix_(S, S)])
            return 0.25*(Z([u(i), u(j)]) + Z([u(i), d(j)]) + Z([d(i), u(j)]) + Z([d(i), d(j)]))
        return self._avg(f)

    def SS(self, i, j):
        """<S_i·S_j>(i≠j) = S^zS^z + (S+S- + S-S+)/2。"""
        u, d = self.up, self.dn
        def f(g):
            Z = lambda S: np.linalg.det(np.eye(len(S)) - 2*g[np.ix_(S, S)])
            szsz = (Z([d(i), d(j)]) - Z([d(i), u(j)]) - Z([u(i), d(j)]) + Z([u(i), u(j)]))/16
            flip = lambda a, b, c, dd: g[a, b]*g[c, dd] - g[a, dd]*g[c, b]   # <c†_a c_b c†_c c_d>
            spsm = flip(u(i), d(i), d(j), u(j))      # S+_i S-_j = c†_{i↑}c_{i↓} c†_{j↓}c_{j↑}
            smsp = flip(d(i), u(i), u(j), d(j))      # S-_i S+_j
            return szsz + 0.5*(spsm + smsp)
        return self._avg(f)

    def hop(self, i, j):
        """上向きと下向きの平均 (1/2)Σ_σ <c†_{iσ}c_{jσ} + h.c.>。"""
        u, d = self.up, self.dn
        return self._avg(lambda g: 0.5*(g[u(i), u(j)] + g[u(j), u(i)] + g[d(i), d(j)] + g[d(j), d(i)]))

    # ---- 観測量(全スピン 0 の状態での期待値) ----
    def Zup(self, i):          return self.one_minus_n(i)
    def ZupZdn_onsite(self, i): return self.zstring([self.up(i), self.dn(i)])
    def ZupZup(self, i, j):    return self.cc(i, j) + (4/3)*self.SS(i, j)
    def FF(self, i, j):        return self.zstring([self.up(i), self.dn(i), self.up(j), self.dn(j)])
    def P_block(self, sites):  return self.zstring([p for s in sites for p in (self.up(s), self.dn(s))])

# ---------------------------------------------------------------------------------------------
#  射影してから変分する PHF(VAP): 射影したエネルギーを、共線(N↑ = N↓)の軌道について最小化する
# ---------------------------------------------------------------------------------------------
def phf_energy(phf, edges, U, t=1.0):
    n = phf.n
    ai = np.array([a-1 for a, b in edges]); bj = np.array([b-1 for a, b in edges]); ii = np.arange(n)
    def f(g):
        kin = -t*np.sum(g[ai, bj] + g[bj, ai] + g[n+ai, n+bj] + g[n+bj, n+ai])
        pot = U*np.sum(g[ii, ii]*g[n+ii, n+ii] - g[ii, n+ii]*g[n+ii, ii])
        return kin + pot
    return phf._avg(f)

def vap(edges, n, LX, W, U, nbeta=48, maxiter=400, t=1.0):
    """UHF の軌道から出発し、Thouless の回転 (占有 + 仮想·K) を QR で直交化して射影エネルギーを最小化する。"""
    from scipy.optimize import minimize
    Nu = Nd = n//2
    T = np.zeros((n, n))
    for a, b in edges: T[a-1, b-1] = T[b-1, a-1] = -t
    Pu0, Pd0 = solve_uhf_orbitals(edges, n, LX, W, U, Nu, Nd)
    def full_basis(P):
        Q, _ = np.linalg.qr(np.hstack([P, np.random.default_rng(0).normal(size=(n, n - P.shape[1]))]))
        Q[:, :P.shape[1]] = P; Q, _ = np.linalg.qr(Q); return Q[:, :P.shape[1]]*np.sign(np.sum(Q[:, :P.shape[1]]*P, 0)), Q[:, P.shape[1]:]
    Ou, Vu = full_basis(Pu0); Od, Vd = full_basis(Pd0)
    nk = (n - Nu)*Nu
    def orbs(x):
        Ku = x[:nk].reshape(n - Nu, Nu); Kd = x[nk:].reshape(n - Nd, Nd)
        Pu, _ = np.linalg.qr(Ou + Vu @ Ku); Pd, _ = np.linalg.qr(Od + Vd @ Kd); return Pu, Pd
    def E(x):
        Pu, Pd = orbs(x); return phf_energy(PHF(Pu, Pd, nbeta=nbeta), edges, U, t)/n
    x0 = np.zeros(2*nk)
    r = minimize(E, x0, method="L-BFGS-B", options=dict(maxiter=maxiter, gtol=1e-8))
    Pu, Pd = orbs(r.x)
    return Pu, Pd, r.fun, E(x0), r
