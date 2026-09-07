"""総予算 B = n_u * n_m を固定し、配分だけ変えたときの実誤差を測る。

利得 G は n_m とともに伸びるが、G は「標準シャドウに対する比」でしかない。
実際に効くのは CRM 推定量の絶対誤差なので、そちらを直接見る。
1D L=12 OBC U=8 の onsite ZZ(台=サイト内2量子ビット)。
"""
import numpy as np, csv, math
rng = np.random.default_rng(90211)
R = [x for x in csv.DictReader(open('crm_edfid_all.tsv'), delimiter='\t')]
s = [x for x in R if x['W']=='1' and x['LX']=='12' and x['geometry']=='cylinder'
     and x['U']=='8.0' and x['observable']=='ZZ onsite']
ZZ_rho = float(s[0]['true']); d_rho = (ZZ_rho + 1)/4
rho = np.array([d_rho, .5-d_rho, .5-d_rho, d_rho])
ZZ  = np.array([1., -1., -1., 1.])
PRI = {x['prior']: ZZ_rho - float(x['Delta']) for x in s}
NA = 2

def law(Ps, nm):
    a = 3.0**NA - 1; vs = 3.0**NA*(1-ZZ_rho**2)/nm
    return (a*ZZ_rho**2 + vs)/(a*(ZZ_rho-Ps)**2 + vs)

def run(nu, nm, Ps, rep):
    out_s = np.empty(rep); out_c = np.empty(rep)
    for r in range(rep):
        match = rng.random(nu) < 1/9
        n_hit = int(match.sum())
        es = np.zeros(nu); ec = np.zeros(nu)
        if n_hit:
            idx = rng.choice(4, size=(n_hit, nm), p=rho)
            m = ZZ[idx].mean(axis=1)
            es[match] = 9.0*m; ec[match] = 9.0*(m - Ps)
        out_s[r] = es.mean(); out_c[r] = ec.mean() + Ps
    return out_s.std(ddof=1), out_c.std(ddof=1)

B = 10000
print(f"総予算 B = {B} を固定。<ZZ>_ρ = {ZZ_rho:.4f}\n")
print(f"{'prior':<7}{'n_m':>7}{'n_u':>7}{'G(法則)':>10}{'G(実測)':>10}"
      f"{'CRM の誤差':>13}{'標準の誤差':>13}{'CRM 誤差の比':>13}")
for p in ("chi4", "UHF"):
    Ps = PRI[p]; base = None
    for nm in (1, 10, 100, 1000, 10000):
        nu = B // nm
        rep = 3000 if nu*nm <= 20000 else 1200
        sd_s, sd_c = run(nu, nm, Ps, rep)
        if base is None: base = sd_c
        print(f"{p:<7}{nm:>7}{nu:>7}{law(Ps,nm):>10.1f}{(sd_s/sd_c)**2:>10.1f}"
              f"{sd_c:>13.5f}{sd_s:>13.5f}{sd_c/base:>12.2f}x")
    print()
