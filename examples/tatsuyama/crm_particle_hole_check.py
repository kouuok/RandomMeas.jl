"""半充填の二部格子で <n_i> = 1 がサイトごとに厳密に成り立つことと、それが崩れる条件を厳密対角化で確かめる。

1次元 Hubbard 鎖(開放端)。粒子数 (N↑, N↓) を固定したセクターで基底状態を求め、
サイトごとの電子数 <n_i>、二重占有 <d_i>、結合の運動エネルギー <K_ij> を出す。
  (a) 半充填・二部格子                 → <n_i> = 1 が機械精度で成り立つはず
  (b) ドープ(半充填から外す)          → Friedel 振動でサイトごとに変わるはず
  (c) 次近接ホッピング t'(二部格子でない) → 半充填でも変わるはず
  (d) 端に外場ポテンシャル             → 変わるはず
  (e) 奇数長のリング(二部格子でないが並進対称) → 基底状態が縮退していると、1つの固有ベクトルは並進対称性を破りうる。
      縮退した多様体で平均すれば並進対称性だけで <n_i> = 1 に戻る
さらに粒子正孔変換が対称性であることを、その帰結で確かめる: Π† H Π = H - U(N̂ - L) なので、
(N↑, N↓) のセクターと (L-N↑, L-N↓) のセクターのスペクトルは U(L - N↑ - N↓) だけずれて一致するはず。
二部格子でない(t' がある)とこの関係は崩れる。
"""
import itertools
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla


def basis(L, Nu, Nd):
    ups = [sum(1 << i for i in c) for c in itertools.combinations(range(L), Nu)]
    dns = [sum(1 << i for i in c) for c in itertools.combinations(range(L), Nd)]
    states = [(u, d) for u in ups for d in dns]
    return states, {s: k for k, s in enumerate(states)}


def hop_sign(b, i, j):
    """c†_i c_j を同じスピンのビット列 b に作用させたときの符号(i と j の間の占有数の偶奇)。"""
    lo, hi = min(i, j), max(i, j)
    mask = ((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1)
    return -1 if bin(b & mask).count("1") % 2 else 1


def hamiltonian(L, Nu, Nd, bonds, U, eps=None):
    """bonds: [(i, j, t_ij)]。H = -Σ t_ij (c†_i c_j + h.c.) + U Σ n↑n↓ + Σ eps_i n_i"""
    states, idx = basis(L, Nu, Nd)
    rows, cols, vals = [], [], []
    for k, (u, d) in enumerate(states):
        diag = U * bin(u & d).count("1")
        if eps is not None:
            diag += sum(eps[i] * (((u >> i) & 1) + ((d >> i) & 1)) for i in range(L))
        rows.append(k); cols.append(k); vals.append(diag)
        for (i, j, t) in bonds:
            for a, b in ((i, j), (j, i)):                     # c†_a c_b
                for spin in (0, 1):
                    bits = u if spin == 0 else d
                    if (bits >> b) & 1 and not (bits >> a) & 1:
                        nb = bits ^ (1 << b) ^ (1 << a)
                        s = hop_sign(bits, a, b)
                        new = (nb, d) if spin == 0 else (u, nb)
                        rows.append(idx[new]); cols.append(k); vals.append(-t * s)
    n = len(states)
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n)), states, idx


def ground(H):
    w, v = sla.eigsh(H, k=2, which="SA")
    o = np.argsort(w)
    return w[o], v[:, o[0]]


def observables(psi, states, L, bonds, idx):
    p = np.abs(psi) ** 2
    n = np.zeros(L); dbl = np.zeros(L)
    for c, (u, d) in zip(p, states):
        for i in range(L):
            n[i] += c * (((u >> i) & 1) + ((d >> i) & 1))
            dbl[i] += c * (((u >> i) & 1) & ((d >> i) & 1))
    K = {}
    for (i, j, t) in bonds:                                   # <Σσ c†_i c_j + h.c.>
        val = 0.0
        for k, (u, d) in enumerate(states):
            if psi[k] == 0: continue
            for spin in (0, 1):
                bits = u if spin == 0 else d
                for a, b in ((i, j), (j, i)):
                    if (bits >> b) & 1 and not (bits >> a) & 1:
                        nb = bits ^ (1 << b) ^ (1 << a)
                        new = (nb, d) if spin == 0 else (u, nb)
                        val += psi[idx[new]] * hop_sign(bits, a, b) * psi[k]
        K[(i, j)] = val
    return n, dbl, K


def chain(L, t=1.0, tp=0.0, pbc=False):
    b = [(i, i + 1, t) for i in range(L - 1)]
    if pbc: b.append((0, L - 1, t))
    if tp: b += [(i, i + 2, tp) for i in range(L - 2)]
    return b


def report(title, L, Nu, Nd, bonds, U, eps=None):
    H, states, idx = hamiltonian(L, Nu, Nd, bonds, U, eps)
    w, psi = ground(H)
    n, dbl, K = observables(psi, states, L, bonds, idx)
    nn = [K[(i, i + 1)] for i in range(L - 1) if (i, i + 1) in K]
    print(f"\n{title}   (次の固有値との差 {w[1]-w[0]:.3f})")
    print("  <n_i> : " + " ".join(f"{x:.4f}" for x in n) + f"   max|n_i - 平均| = {np.max(np.abs(n - n.mean())):.1e}")
    print("  <d_i> : " + " ".join(f"{x:.4f}" for x in dbl))
    print("  隣の結合 <K> : " + " ".join(f"{x:.4f}" for x in nn))
    return H, states, idx, psi, n


if __name__ == "__main__":
    L, U = 8, 4.0
    print(f"1次元 Hubbard 鎖 L={L}、U={U}(t=1)、開放端")
    H, states, idx, psi, n = report("(a) 半充填・二部格子 N↑=N↓=4", L, 4, 4, chain(L), U)
    report("(b) ドープ N↑=N↓=3(半充填から電子を2つ抜く)", L, 3, 3, chain(L), U)
    report("(c) 半充填 + 次近接ホッピング t'=0.3(二部格子でない)", L, 4, 4, chain(L, tp=0.3), U)
    eps = np.zeros(L); eps[0] = 0.5
    report("(d) 半充填 + 左端のサイトだけに外場 0.5", L, 4, 4, chain(L), U, eps)
    report("(e) 半充填・奇数長のリング L=7、N↑=4 N↓=3(二部格子でないが並進対称)", 7, 4, 3, chain(7, pbc=True), U)
    H7, st7, idx7 = hamiltonian(7, 4, 3, chain(7, pbc=True), U)
    w7, v7 = sla.eigsh(H7, k=8, which="SA"); o = np.argsort(w7); w7, v7 = w7[o], v7[:, o]
    g = int(np.sum(w7 < w7[0] + 1e-8))
    navg = np.mean([observables(v7[:, k], st7, 7, [], idx7)[0] for k in range(g)], axis=0)
    print(f"  縮退した {g} 個の基底状態で平均した <n_i> : " + " ".join(f"{x:.4f}" for x in navg)
          + f"   max|n_i - 1| = {np.max(np.abs(navg - 1)):.1e}")

    # 粒子正孔変換の帰結: (N↑,N↓) と (L-N↑,L-N↓) のスペクトルが U(L-N↑-N↓) ずれて一致する
    print("\n粒子正孔変換の確認: セクター (3,3) と (5,5) の低い固有値(後者から 2U を引く)")
    for label, tp in (("二部格子(t'=0)", 0.0), ("t'=0.3", 0.3)):
        e33 = np.sort(sla.eigsh(hamiltonian(L, 3, 3, chain(L, tp=tp), U)[0], k=4, which="SA")[0])
        e55 = np.sort(sla.eigsh(hamiltonian(L, 5, 5, chain(L, tp=tp), U)[0], k=4, which="SA")[0]) - 2*U
        print(f"  {label}: (3,3) " + " ".join(f"{x:+.8f}" for x in e33))
        print(f"  {'':{len(label)}s}  (5,5)-2U " + " ".join(f"{x:+.8f}" for x in e55) + f"   最大差 {np.max(np.abs(e33-e55)):.1e}")
