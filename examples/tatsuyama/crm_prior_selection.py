"""段階2: 有限の測定データだけで prior を選べるか。

同じ1回分のランダム測定データ(元論文の反復手順と同じ設定)を後処理し、
7つの候補 prior それぞれについて Δ̂, D̂_A, F̂_A を推定して利得を予測する。
真の利得と突き合わせ、順位相関・選択の正しさ・予測誤差を測る。

対象: L=16, U=8 半充填, site 2, ZZ onsite (単一Pauli列、台=2量子ビット)
README 準備4g 用。
"""
import numpy as np
rng=np.random.default_rng(7217)
I2=np.eye(2); X=np.array([[0,1],[1,0]],complex); Y=np.array([[0,-1j],[1j,0]]); Z=np.diag([1,-1]).astype(complex)
P1={'I':I2,'X':X,'Y':Y,'Z':Z}
EV={b:np.linalg.eigh(P1[b]) for b in 'XYZ'}
S=[a+b for a in 'IXYZ' for b in 'IXYZ']
op=lambda s: np.kron(P1[s[0]],P1[s[1]])
NM=100; NA=2

d_rho=0.09924590
PRIORS={"chi2":0.0,"chi4":0.09023233,"chi8":0.098702,"chi16":0.099212,
        "chi32":0.099244,"UHF":0.104691,"UHFsym":0.104691}
# UHF は S^z 対称性を破る: 台の RDM が d だけでは決まらない
UHF_M=0.762310      # データの D_A(UHF)=0.3866 に合わせた局所磁化(F_A も 0.7006 対 0.7004 で一致)
MAG={k:(UHF_M if k=="UHF" else 0.0) for k in PRIORS}
def state(d,m):
    # diag(d, 1/2-d+m/2, 1/2-d-m/2, d) を規格化(m は S^z 方向のずれ)
    v=np.array([d,.5-d+m/2,.5-d-m/2,d]); v=np.clip(v,1e-12,None); return np.diag(v/v.sum()).astype(complex)
rho=state(d_rho,0.0)
SIG={k:state(PRIORS[k],MAG[k]) for k in PRIORS}
PP={k:{s:float(np.real(np.trace(op(s)@SIG[k]))) for s in S} for k in PRIORS}
ZZ_rho=float(np.real(np.trace(op('ZZ')@rho)))
Dtrue={k:.5*float(np.sum(np.abs(np.linalg.eigvalsh(rho-SIG[k])))) for k in PRIORS}
Dl_true={k:ZZ_rho-PP[k]['ZZ'] for k in PRIORS}
def gain(P,dl): 
    vs=3**NA*(1-P*P)/NM; a=3**NA-1
    return (a*P*P+vs)/(a*min(abs(dl),2.0)**2+vs)
Gtrue={k:gain(ZZ_rho,Dl_true[k]) for k in PRIORS}

BP=[(a,b) for a in 'XYZ' for b in 'XYZ']
PROB={};SG={}
for a,b in BP:
    U=np.kron(EV[a][1].conj().T,EV[b][1].conj().T)
    p=np.real(np.diag(U@rho@U.conj().T)); PROB[(a,b)]=np.clip(p,0,None)/p.sum()
    SG[(a,b)]=(EV[a][0],EV[b][0])

def dataset(NU):
    """1回分のデータ。文字列ごとに (Σ shot平均, 一致設定数, 係数) を返す"""
    acc={s:[0.0,0] for s in S}
    for (a,b),n in zip(BP,rng.multinomial(NU,[1/9]*9)):
        if n==0: continue
        idx=rng.choice(4,size=(n,NM),p=PROB[(a,b)])
        s0=SG[(a,b)][0][idx//2]; s1=SG[(a,b)][1][idx%2]
        for s,mv in ((a+'I',s0.mean(1)),('I'+b,s1.mean(1)),(a+b,(s0*s1).mean(1))):
            acc[s][0]+=float(mv.sum()); acc[s][1]+=n
    return acc

def estimate(acc,NU,k,crm):
    """prior k での全16成分の推定"""
    hat={}
    for s in S:
        if s=='II': hat[s]=1.0; continue
        f=3**sum(1 for c in s if c!='I'); Ss,ns=acc[s]
        hat[s]=f*Ss/NU + (PP[k][s]*(1-f*ns/NU) if crm else 0.0)
    R=sum(hat[s]*op(s) for s in S)/4; R=(R+R.conj().T)/2
    ev=np.linalg.eigvalsh(R-SIG[k])
    D=.5*float(np.sum(np.abs(ev)))
    # 忠実度(再構成した R を物理的な状態に射影してから)
    w,V=np.linalg.eigh(R); w=np.clip(w,0,None); w=w/w.sum(); Rp=(V*w)@V.conj().T
    sq=np.linalg.eigvalsh(Rp); sq=np.clip(sq,0,None)
    A=np.linalg.cholesky(Rp+1e-15*np.eye(4)) if False else None
    ws,Vs=np.linalg.eigh(Rp); ws=np.clip(ws,0,None)
    sr=(Vs*np.sqrt(ws))@Vs.conj().T
    M=sr@SIG[k]@sr; mw=np.clip(np.linalg.eigvalsh(M),0,None)
    F=float(np.sum(np.sqrt(mw))**2)
    return hat['ZZ']-PP[k]['ZZ'], D, F, hat['ZZ']

def spearman(a,b):
    ra=np.argsort(np.argsort(a)).astype(float); rb=np.argsort(np.argsort(b)).astype(float)
    ra-=ra.mean(); rb-=rb.mean()
    return float(ra@rb/np.sqrt((ra@ra)*(rb@rb)))

keys=list(PRIORS); gt=np.array([Gtrue[k] for k in keys])
best=keys[int(np.argmax(gt))]
print("真値")
print(f"{'prior':<9}{'Δ_ZZ':>12}{'D_A':>12}{'F_A':>11}{'G':>9}")
for k in keys:
    Fk=None
    print(f"{k:<9}{Dl_true[k]:>12.3e}{Dtrue[k]:>12.3e}{'':>11}{Gtrue[k]:>9.2f}")
print(f"最良 prior = {best} (G={gt.max():.2f})\n")

REP=200
print(f"{'N_U':>7}{'予測に使う量':>14}{'順位相関':>10}{'最良を当てる率':>14}{'G の相対誤差(中央)':>18}")
for NU in (100,1000,10000):
    coll={"Δ̂ から":[], "D̂_A から":[], "F̂_A から":[]}
    hit={k:0 for k in coll}; sp={k:[] for k in coll}
    for _ in range(REP):
        acc=dataset(NU)
        gd,gD,gF=[],[],[]
        for k in keys:
            dl,D,F,zz=estimate(acc,NU,k,True)
            gd.append(gain(zz,dl)); gD.append(gain(zz,2*D))
            gF.append(gain(zz,2*np.sqrt(max(1-F,0))))
        for nm,g in (("Δ̂ から",gd),("D̂_A から",gD),("F̂_A から",gF)):
            g=np.array(g); sp[nm].append(spearman(g,gt))
            if keys[int(np.argmax(g))]==best or Gtrue[keys[int(np.argmax(g))]]>=0.95*gt.max(): hit[nm]+=1
            coll[nm].append(np.median(np.abs(g-gt)/gt))
    for nm in coll:
        print(f"{NU:>7}{nm:>15}{np.mean(sp[nm]):>10.3f}{100*hit[nm]/REP:>13.0f}%{np.median(coll[nm]):>17.1%}")
    print()
