"""clara の結果を取り込み、統合表を作り直して系サイズ依存を集計する。

crm_2d_edfid_W{W}L{LX}_{geo}_U{U}.tsv を全部読み、恒等成分の補正
(c0, single_pauli, P_fluc, eps, G, G_max)を掛けてから 1 つにまとめる。
"""
import csv, glob, math, sys, statistics as st
f = float
SPEC = {"ZZ onsite":(0.0,True), "ZZ up-up nb":(0.0,True),
        "SzSz nb":(0.0,False), "SxSx nb":(0.0,False),
        "DoubleOcc":(0.25,False), "n":(1.0,False), "Sz":(0.0,False)}
NM = 100

def gain(P, d, nA):
    a = 3.0**nA - 1; vs = 3.0**nA*(1-P*P)/NM
    return (a*P*P + vs)/(a*d*d + vs) if vs > 0 else float('nan')

rows = []
for fn in sorted(glob.glob("crm_2d_edfid_W*_U*.tsv")):
    for x in csv.DictReader(open(fn), delimiter='\t'):
        c0, sp = SPEC[x['observable']]
        P = f(x['true']); nA = int(x['nA']); d = f(x['Delta']); Pf = P - c0
        x['dim'] = "1D" if x['W'] == '1' else "2D"
        x['exact'] = 'yes' if abs(f(x['Hvar'])) < 1e-8 else 'no'
        x['c0'] = f"{c0:.4f}"; x['single_pauli'] = str(sp).lower(); x['P_fluc'] = f"{Pf:.10f}"
        x['eps'] = f"{abs(d)/abs(Pf):.6e}" if abs(Pf) > 1e-12 else "NaN"
        if sp:
            x['G'] = f"{gain(Pf,d,nA):.6f}"; x['G_max'] = f"{gain(Pf,0.0,nA):.6f}"
        else:
            x['G'] = "NaN"; x['G_max'] = "NaN"
        rows.append(x)

COLS = ['dim','W','LX','nsites','geometry','bipartite','U','Eref','Hvar','exact',
        'prior','global_fid','observable','nA','c0','single_pauli','P_fluc',
        'F_supp','D_supp','true','prior_val','Delta','eps','G','G_max']
with open("crm_edfid_all.tsv","w") as io:
    io.write("\t".join(COLS)+"\n")
    for x in rows: io.write("\t".join(str(x.get(c,"")) for c in COLS)+"\n")
nex = sum(1 for x in rows if x['exact']=='yes')
print(f"統合: {len(rows)}行 / 厳密 {nex}行 / 系 {len({(x['W'],x['LX'],x['geometry'],x['U']) for x in rows})}")

# --- 系サイズ依存 ---
print("\n【1D OBC の系サイズ依存 — UHF, ZZ onsite】")
for U in ("2.0","4.0","8.0","12.0"):
    s = [x for x in rows if x['dim']=='1D' and x['geometry']=='cylinder' and x['U']==U
         and x['prior']=='UHF' and x['observable']=='ZZ onsite' and x['exact']=='yes']
    if len(s) < 2: continue
    s.sort(key=lambda y: int(y['LX']))
    gf = [f(x['global_fid']) for x in s]; fs = [f(x['F_supp']) for x in s]; g = [f(x['G']) for x in s]
    Ls = [x['LX'] for x in s]
    print(f"  U={U.replace('.0',''):<3} L={','.join(Ls)}")
    print(f"       大域F {gf[0]:.5f} → {gf[-1]:.3e}  ({gf[0]/gf[-1]:>8.1f}倍の崩壊)")
    print(f"       F_台  {fs[0]:.5f} → {fs[-1]:.5f}  (変動 {max(fs)-min(fs):.1e})")
    print(f"       G     {g[0]:.1f} → {g[-1]:.1f}  ({max(g)/min(g):.3f}倍)")

# --- 平均場 vs MPS ---
print("\n【平均場は MPS を上回るか — 1D OBC, ZZ onsite】")
print(f"{'U':>4}{'L':>5}{'chi32 の F':>12}{'G(chi32)':>10}{'G(UHF)':>10}{'UHF 勝ち':>9}")
for U in ("2.0","4.0","8.0","12.0"):
    for L in sorted({x['LX'] for x in rows if x['dim']=='1D'}, key=int):
        s = [x for x in rows if x['dim']=='1D' and x['geometry']=='cylinder' and x['U']==U
             and x['LX']==L and x['observable']=='ZZ onsite' and x['exact']=='yes']
        d = {x['prior']: x for x in s}
        if 'UHF' not in d or 'chi32' not in d: continue
        gu = f(d['UHF']['G']); gc = f(d['chi32']['G'])
        print(f"{U.replace('.0',''):>4}{L:>5}{f(d['chi32']['global_fid']):>12.6f}"
              f"{gc:>10.1f}{gu:>10.1f}{'○' if gu>gc else '×':>7}")
