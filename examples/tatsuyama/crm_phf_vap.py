"""射影してから変分する PHF(VAP)を、UHF を射影しただけの PHF(PAV)と比べる(1次元開放端)。
初期点は UHF と、磁化を弱めた UHF(U/2 で解いた軌道)の2通り。"""
import time, csv
import numpy as np
from scipy.optimize import minimize
from crm_phf import PHF, phf_energy, solve_uhf_orbitals
from crm_eigen_uhfsym import lat_edges

Eref = {}
for r in csv.DictReader(open("crm_edfid_all.tsv"), delimiter="\t"):
    if r["W"] == "1" and r["geometry"] == "cylinder": Eref[(int(r["LX"]), float(r["U"]))] = float(r["Eref"])

def run_vap(n, U, start_U):
    e = lat_edges(n, 1, False); Nu = n//2
    Pu0, Pd0 = solve_uhf_orbitals(e, n, n, 1, start_U, Nu, Nu)
    def E(x):
        Pu, _ = np.linalg.qr(x[:n*Nu].reshape(n, Nu)); Pd, _ = np.linalg.qr(x[n*Nu:].reshape(n, Nu))
        return phf_energy(PHF(Pu, Pd, nbeta=48), e, U)/n
    r = minimize(E, np.concatenate([Pu0.ravel(), Pd0.ravel()]), method="L-BFGS-B", options=dict(maxiter=3000, gtol=1e-9, ftol=1e-15))
    Pu, _ = np.linalg.qr(r.x[:n*Nu].reshape(n, Nu)); Pd, _ = np.linalg.qr(r.x[n*Nu:].reshape(n, Nu))
    return r.fun, Pu, Pd, r

if __name__ == "__main__":
    out = {}
    for n in (8, 16):
        for U in (4.0, 8.0, 12.0):
            e = lat_edges(n, 1, False)
            Pu, Pd = solve_uhf_orbitals(e, n, n, 1, U, n//2, n//2)
            Gu, Gd = Pu @ Pu.T, Pd @ Pd.T
            Euhf = (sum(-(Gu[a-1,b-1]+Gu[b-1,a-1]+Gd[a-1,b-1]+Gd[b-1,a-1]) for a, b in e) + U*sum(Gu[i,i]*Gd[i,i] for i in range(n)))/n
            Epav = phf_energy(PHF(Pu, Pd, nbeta=48), e, U)/n
            t0 = time.time(); best = None
            for sU in (U, U/2):
                Ev, Pv, Qv, r = run_vap(n, U, sU)
                if best is None or Ev < best[0]: best = (Ev, Pv, Qv, sU, r.nit)
            Ev, Pv, Qv, sU, nit = best; Ex = Eref[(n, U)]/n
            m_uhf = abs(np.diag(Gu)[n//2] - np.diag(Gd)[n//2])/2; m_vap = abs(np.diag(Pv @ Pv.T)[n//2] - np.diag(Qv @ Qv.T)[n//2])/2
            print(f"L={n:2d} U={U:4.0f}: 厳密 {Ex:+.5f}  UHF {Euhf:+.5f}  PAV {Epav:+.5f}  VAP {Ev:+.5f}(初期 U={sU:g}, {nit} 反復)"
                  f"  相関エネルギーの回収率 PAV {(Epav-Euhf)/(Ex-Euhf):.3f} VAP {(Ev-Euhf)/(Ex-Euhf):.3f}  中央の磁化 UHF {m_uhf:.3f} VAP {m_vap:.3f}  ({time.time()-t0:.0f} 秒)")
            out[(n, U)] = (Pv, Qv)
    np.save("crm_phf_vap_orbitals.npy", {f"{k[0]}_{k[1]}": v for k, v in out.items()}, allow_pickle=True)
