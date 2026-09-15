"""パウリ文字列の和 O = Σ_k c_k P_k を局所ランダムパウリ測定(1基底あたり n_m ショット)で推定するときの、
標準シャドウと CRM の1ユニットあたりの分散の厳密な式と、その乱数による検証。

  両立する対(重なる量子ビットで同じパウリ)だけが同じ基底で同時に当たる(確率 3^{-|A_k ∪ A_l|})。
  Var_std = Σ_{k,l 両立} c_k c_l 3^{|A_k∩A_l|} [ <P_k P_l>/n_m + (1-1/n_m) x_k x_l ] - (Σ c_k x_k)²
  Var_CRM = Σ_{k,l 両立} c_k c_l 3^{|A_k∩A_l|} [ (<P_k P_l> - x_k x_l)/n_m + Δ_k Δ_l ] - (Σ c_k Δ_k)²
  (x_k = <P_k>_ρ, Δ_k = x_k - <P_k>_σ。k=l なら単一パウリの式に戻る)
"""
import itertools
import numpy as np

I2 = np.eye(2); X = np.array([[0, 1], [1, 0]]); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1, -1])
PM = {"I": I2, "X": X, "Y": Y, "Z": Z}

def compatible(p, q):
    return all(a == "I" or b == "I" or a == b for a, b in zip(p, q))

def overlap(p, q):
    return sum(1 for a, b in zip(p, q) if a != "I" and b != "I")

def product(p, q):
    """両立する2つのパウリ文字列の積(重なりは同じ文字なので恒等になる)。"""
    return "".join(b if a == "I" else (a if b == "I" else "I") for a, b in zip(p, q))

def variances(terms, ex, ey, nm):
    """terms: [(c_k, P_k)]、ex/ey: パウリ文字列 → ρ/σ での期待値(積 P_kP_l も ex に含める)。"""
    sx = sum(c*ex[p] for c, p in terms); sd = sum(c*(ex[p] - ey[p]) for c, p in terms)
    vs = vc = 0.0
    for (ck, pk), (cl, pl) in itertools.product(terms, terms):
        if not compatible(pk, pl): continue
        w = ck*cl*3.0**overlap(pk, pl); pp = ex[product(pk, pl)]
        vs += w*(pp/nm + (1 - 1/nm)*ex[pk]*ex[pl])
        vc += w*((pp - ex[pk]*ex[pl])/nm + (ex[pk] - ey[pk])*(ex[pl] - ey[pl]))
    return vs - sx**2, vc - sd**2

if __name__ == "__main__":
    rng = np.random.default_rng(7)
    n = 3
    def rand_rho():
        A = rng.normal(size=(2**n, 2**n)) + 1j*rng.normal(size=(2**n, 2**n)); R = A @ A.conj().T; return R/np.trace(R)
    op = lambda p: np.array(np.kron(np.kron(PM[p[0]], PM[p[1]]), PM[p[2]]))
    ev = lambda R, p: float(np.real(np.trace(R @ op(p))))
    rho, sig = rand_rho(), rand_rho()
    # 重なり・両立・非両立がすべて含まれる和
    terms = [(0.7, "ZZI"), (-0.4, "ZIZ"), (0.3, "XZI"), (0.5, "IIZ"), (-0.2, "YYX"), (0.6, "ZII")]
    allp = {"".join(t) for t in itertools.product("IXYZ", repeat=n)}
    ex = {p: ev(rho, p) for p in allp}; ey = {p: ev(sig, p) for p in allp}
    nm = 5
    vs_th, vc_th = variances(terms, ex, ey, nm)
    # 乱数で測定を模擬: 各量子ビットの基底を一様に選び、その基底での出力分布から n_m ショット
    Hd = np.array([[1, 1], [1, -1]])/np.sqrt(2); Sdg = np.diag([1, -1j])
    rot = {"Z": I2, "X": Hd, "Y": Hd @ Sdg}
    NU = 400_000; est_s = np.empty(NU); est_c = np.empty(NU)
    bases = rng.integers(0, 3, size=(NU, n)); letters = "XYZ"
    cache = {}
    for u in range(NU):
        b = "".join(letters[i] for i in bases[u])
        if b not in cache:
            Ub = np.kron(np.kron(rot[b[0]], rot[b[1]]), rot[b[2]])
            pr = np.clip(np.real(np.diag(Ub @ rho @ Ub.conj().T)), 0, None); cache[b] = pr/pr.sum()
        shots = rng.choice(2**n, size=nm, p=cache[b])
        bits = (shots[:, None] >> np.arange(n)[::-1]) & 1          # 量子ビット0が最上位
        s_acc = c_acc = 0.0
        for c, p in terms:
            supp = [q for q in range(n) if p[q] != "I"]
            if all(b[q] == p[q] for q in supp):
                m = np.mean(np.prod(1 - 2*bits[:, supp], axis=1))
                s_acc += c*3**len(supp)*m; c_acc += c*3**len(supp)*(m - ey[p])
            c_acc += 0  # 足し戻す定数 Σ c_k y_k は分散に効かない
        est_s[u] = s_acc; est_c[u] = c_acc
    print(f"平均: 標準 {est_s.mean():+.4f}(真 {sum(c*ex[p] for c,p in terms):+.4f})  CRM+Σcy {est_c.mean()+sum(c*ey[p] for c,p in terms):+.4f}")
    print(f"標準シャドウの分散: 乱数 {est_s.var():.4f}  式 {vs_th:.4f}  比 {est_s.var()/vs_th:.4f}")
    print(f"CRM の分散:         乱数 {est_c.var():.4f}  式 {vc_th:.4f}  比 {est_c.var()/vc_th:.4f}")
    print(f"(乱数の分散の統計誤差は約 {np.sqrt(2/NU)*100:.1f}%)")
