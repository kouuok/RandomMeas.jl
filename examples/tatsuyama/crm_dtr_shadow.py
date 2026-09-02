"""台のRDMをシャドウ(標準/CRM)から再構成し、トレース距離が測れるかを調べる。

rho_A, sigma_A は実データ(L=16, U=8 半充填, site 2, prior chi_p=4)。
各設定でランダムな局所Pauli基底を引き、NM ショット測る。16個のPauli期待値を
標準シャドウとCRMシャドウの両方で推定し、RDMを再構成してトレース距離を出す。
README 準備4f 用。
"""
import numpy as np
rng=np.random.default_rng(20260902)
I2=np.eye(2); X=np.array([[0,1],[1,0]],complex); Y=np.array([[0,-1j],[1j,0]]); Z=np.diag([1,-1]).astype(complex)
P1={'I':I2,'X':X,'Y':Y,'Z':Z}
EV={b:np.linalg.eigh(P1[b]) for b in 'XYZ'}
d_rho,d_sig=0.09924590,0.09023233
rho=np.diag([d_rho,.5-d_rho,.5-d_rho,d_rho]).astype(complex)
sig=np.diag([d_sig,.5-d_sig,.5-d_sig,d_sig]).astype(complex)
S=[a+b for a in 'IXYZ' for b in 'IXYZ']
op=lambda s: np.kron(P1[s[0]],P1[s[1]])
sigP={s:np.real(np.trace(op(s)@sig)) for s in S}
Dtrue=.5*np.sum(np.abs(np.linalg.eigvalsh(rho-sig))); DZZ=np.real(np.trace(op('ZZ')@(rho-sig)))
BP=[(a,b) for a in 'XYZ' for b in 'XYZ']
PROB={};SG={}
for a,b in BP:
    U=np.kron(EV[a][1].conj().T,EV[b][1].conj().T)
    p=np.real(np.diag(U@rho@U.conj().T)); PROB[(a,b)]=np.clip(p,0,None)/p.sum()
    SG[(a,b)]=(EV[a][0],EV[b][0])

def run(NU,NM,crm):
    tot={s:0.0 for s in S}
    cnt=rng.multinomial(NU,[1/9]*9)
    for (a,b),n in zip(BP,cnt):
        if n==0: continue
        idx=rng.choice(4,size=(n,NM),p=PROB[(a,b)])
        s0=SG[(a,b)][0][idx//2]; s1=SG[(a,b)][1][idx%2]
        m0,m1,m01=s0.mean(1),s1.mean(1),(s0*s1).mean(1)
        for s,f,mv in ((a+'I',3,m0),('I'+b,3,m1),(a+b,9,m01)):
            tot[s]+= np.sum(f*(mv-sigP[s])) if crm else np.sum(f*mv)
    hat={}
    for s in S:
        if s=='II': hat[s]=1.0
        else: hat[s]=tot[s]/NU + (sigP[s] if crm else 0.0)
    R=sum(hat[s]*op(s) for s in S)/4; R=(R+R.conj().T)/2
    return .5*np.sum(np.abs(np.linalg.eigvalsh(R-sig))), hat['ZZ']-sigP['ZZ']

print(f"真値  D_A = {Dtrue:.4e}   Δ_ZZ = {DZZ:.4e}   |Δ|/2D = {abs(DZZ)/(2*Dtrue):.4f}")
print(f"総測定数 = N_U×N_M\n")
import sys
print(f"{"N_U":>8}{'N_M':>6}{'方式':>6}{'D̂ 平均':>11}{'D̂ 偏り':>9}{'D̂ 相対誤差':>11}{'Δ̂ 相対誤差':>12}{'Δ の利得':>9}")
REP=60
for NU,NM in ((100,100),(1000,100),(10000,100)):
    res={}
    for crm in (False,True):
        o=np.array([run(NU,NM,crm) for _ in range(REP)])
        D,Dl=o[:,0],o[:,1]; res[crm]=(D,Dl)
        g="" if not crm else f"{(res[False][1].std()/Dl.std())**2:>9.1f}"
        print(f"{NU:>8}{NM:>6}{'CRM' if crm else '標準':>7}{D.mean():>11.3e}{D.mean()/Dtrue:>8.1f}倍"
              f"{D.std()/Dtrue:>10.1%}{Dl.std()/abs(DZZ):>11.1%}{g}")
    print()
