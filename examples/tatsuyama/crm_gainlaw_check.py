"""利得法則を、厳密な ρ を持つ全ての系でサンプリング検証する。

単一 Pauli 列 P(固有値 ±1)では、基底が一致した設定でのショットは
「平均 <P> の ±1 ベルヌーイ」に従う。つまり <P> さえ分かれば
シャドウ測定を厳密に模擬でき、ED も RDM も要らない。

§1 は8量子ビット1系での検証だった。ここでは 16〜128量子ビット、
U=2〜12、prior 4種、観測量2種で、実測利得と法則値を突き合わせる。
"""
import numpy as np, csv, math, sys
rng = np.random.default_rng(31337)
f = float
NA, NM_LIST, NU = 2, (1, 10, 100), 400
REP = 3000

R = [x for x in csv.DictReader(open('crm_edfid_all.tsv'), delimiter='\t')
     if x['exact'] == 'yes' and x['single_pauli'] == 'true' and x['dim'] == '1D'
     and x['geometry'] == 'cylinder']

def law(P, d, nm):
    a = 3.0**NA - 1; vs = 3.0**NA*(1-P*P)/nm
    return (a*P*P + vs)/(a*d*d + vs)

def simulate(P, Ps, nm, rep):
    """標準シャドウと CRM シャドウの推定値を rep 回。"""
    # 各設定で基底が P と一致する確率は 3^-|A|
    hit = rng.random((rep, NU)) < 3.0**(-NA)
    n = int(hit.sum())
    # 一致した設定のショット平均(±1, 平均 P の nm 回平均)
    m = np.zeros((rep, NU))
    if n:
        s = np.where(rng.random((n, nm)) < (1+P)/2, 1.0, -1.0)
        m[hit] = s.mean(axis=1)
    amp = 3.0**NA
    est_s = (amp*m*hit).mean(axis=1)
    est_c = (amp*(m - Ps)*hit).mean(axis=1) + Ps
    return est_s, est_c

out = []
for L in ("16","32","64"):
    for U in ("2.0","4.0","8.0","12.0"):
        for o in ("ZZ onsite","ZZ up-up nb"):
            for p in ("chi2","chi8","chi32","UHF"):
                x = [y for y in R if y['LX']==L and y['U']==U and y['observable']==o and y['prior']==p]
                if not x: continue
                x = x[0]; P = f(x['true']); Ps = P - f(x['Delta'])
                for nm in NM_LIST:
                    es, ec = simulate(P, Ps, nm, REP)
                    g_emp = es.var(ddof=1)/ec.var(ddof=1)
                    g_law = law(P, P-Ps, nm)
                    out.append((int(L), f(U), o, p, nm, g_law, g_emp, g_emp/g_law))

print(f"検証点数 {len(out)}  (n_u={NU}, 反復={REP})")
print(f"\n{'L':>4}{'qubit':>6}{'U':>4}{'観測量':<13}{'prior':<7}{'n_m':>5}{'G(法則)':>11}{'G(実測)':>11}{'比':>7}")
for L,U,o,p,nm,gl,ge,r in out:
    if L == 64 and nm == 100:
        print(f"{L:>4}{2*L:>6}{U:>4.0f}{o:<13}{p:<7}{nm:>5}{gl:>11.2f}{ge:>11.2f}{r:>7.3f}")
rs = np.array([r for *_,r in out])
print(f"\n全 {len(rs)} 点の 実測/法則:  中央値 {np.median(rs):.4f}   範囲 [{rs.min():.3f}, {rs.max():.3f}]")
print(f"  |比-1| の 95 パーセンタイル: {np.percentile(np.abs(rs-1),95):.3f}")
print(f"  (分散比の統計誤差は反復 {REP} 回で約 {math.sqrt(4/(REP-1)):.3f})")
big = [(L,U,o,p,nm,gl,ge,r) for L,U,o,p,nm,gl,ge,r in out if abs(r-1) > 0.15]
print(f"\n比が 15% 以上ずれた点: {len(big)}")
for t in big[:8]: print(f"   L={t[0]} U={t[1]:.0f} {t[2]} {t[3]} n_m={t[4]}  法則 {t[5]:.2f} 実測 {t[6]:.2f} 比 {t[7]:.3f}")
