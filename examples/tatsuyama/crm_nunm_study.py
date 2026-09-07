"""n_u と n_m が利得に与える影響を、実際にシャドウをサンプリングして測る。

利得法則では n_u は約分されて消える(標準・CRM 双方の分散が 1/n_u に比例)。
本当にそうか、そして n_m はどう効くかを、1D L=12 OBC U=8 の onsite ZZ で確かめる。
台=サイト内2量子ビットなので、ρ_A の周辺分布だけで厳密にシミュレートできる。
"""
import numpy as np, csv, math
rng = np.random.default_rng(4271)

R = [x for x in csv.DictReader(open('crm_edfid_all.tsv'), delimiter='\t')]
s = [x for x in R if x['W']=='1' and x['LX']=='12' and x['geometry']=='cylinder'
     and x['U']=='8.0' and x['observable']=='ZZ onsite']
ZZ_rho = float(s[0]['true'])
d_rho = (ZZ_rho + 1)/4
rho = np.array([d_rho, .5-d_rho, .5-d_rho, d_rho])      # (空, ↑, ↓, 二重)
ZZ = np.array([1., -1., -1., 1.])
PRI = {x['prior']: ZZ_rho - float(x['Delta']) for x in s}   # <ZZ>_σ
NA, GRID = 2, [(nu, nm) for nu in (25, 100, 400, 1600) for nm in (1, 10, 100, 1000)]

def law(Ps, nm):
    a = 3.0**NA - 1; vs = 3.0**NA*(1-ZZ_rho**2)/nm
    return (a*ZZ_rho**2 + vs)/(a*(ZZ_rho-Ps)**2 + vs)

def one_run(nu, nm, Ps):
    """1回分の推定値を、標準シャドウと CRM シャドウで返す。"""
    # 各設定で 2 量子ビットの基底を独立に引く。ZZ に寄与するのは両方 Z のときだけ
    match = rng.random(nu) < 1/9
    est_s = np.zeros(nu); est_c = np.zeros(nu)
    n_hit = int(match.sum())
    if n_hit:
        idx = rng.choice(4, size=(n_hit, nm), p=rho)
        m = (ZZ[idx]).mean(axis=1)          # 各設定のショット平均
        est_s[match] = 9.0*m
        est_c[match] = 9.0*(m - Ps)
    est_c += Ps
    return est_s.mean(), est_c.mean()

REP = 4000
print(f"1D L=12 OBC U=8, onsite ZZ:  <ZZ>_ρ = {ZZ_rho:.4f}")
print(f"\n{'prior':<8}{'n_u':>6}{'n_m':>6}{'総測定':>9}{'G(法則)':>10}{'G(実測)':>10}{'比':>7}{'誤差':>8}")
for p in ("chi4", "UHF", "chi32"):
    Ps = PRI[p]
    for nu, nm in GRID:
        if nu*nm > 40000: continue
        out = np.array([one_run(nu, nm, Ps) for _ in range(REP)])
        vs_, vc_ = out[:,0].var(ddof=1), out[:,1].var(ddof=1)
        gemp = vs_/vc_; gl = law(Ps, nm)
        se = gemp*math.sqrt(2/(REP-1))*math.sqrt(2)   # 分散比の粗い標準誤差
        print(f"{p:<8}{nu:>6}{nm:>6}{nu*nm:>9}{gl:>10.2f}{gemp:>10.2f}{gemp/gl:>7.2f}{se/gemp:>7.1%}")
    print()
