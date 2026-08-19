import re, math
CHIS=[2,4,6,8,12,16,24,32,48,64,96,128,192,256]
rows=[]; cfg=None
for ln in open('crm_chistar.log'):
    m=re.match(r'U=([\d.]+) δ=([\d.]+) L=\s*(\d+)\s*:', ln)
    if m: cfg=(float(m.group(1)),float(m.group(2)),int(m.group(3))); continue
    if cfg is None: continue
    m=re.match(r'(\S.*?)\s{2,}((?:[\d.eE+-]+|NaN)(?:\s+(?:[\d.eE+-]+|NaN)){13})', ln)
    if not m: continue
    name=m.group(1).strip()
    vals=[float('nan') if v=='NaN' else float(v) for v in m.group(2).split()]
    if len(vals)!=14: continue
    rows.append((cfg,name,vals))

def chistar(eps,target):
    # χ >= maxlinkdim (NaN) は切断が厳密 ⇒ 目標を満たすとみなす
    ok=True; best=float('nan')
    for k in range(len(CHIS)-1,-1,-1):
        e=eps[k]
        if not math.isnan(e) and e>=target: ok=False
        if ok: best=CHIS[k]
    return best

cfgs=[(1.0,0.0,"U=1 half"),(2.0,0.0,"U=2 half"),(4.0,0.0,"U=4 half"),(8.0,0.0,"U=8 half"),
      (4.0,0.125,"U=4 δ=1/8"),(8.0,0.125,"U=8 δ=1/8")]
Ls=[16,32,64,128]
obs=["DoubleOcc","ZZ onsite","SzSz r=1","ZZ up-up r=1","hop up r=1"]
print("χ*(ε<1e-2) の L 依存 【NaN=厳密表現 として修正】")
print(f"{'config':16s} {'observable':14s}" + "".join(f"{'L='+str(l):>9s}" for l in Ls) + f"{'L128/L16':>10s}")
out=open('crm_chistar_fixed.tsv','w'); out.write("U\tdoping\tL\tobservable\ttarget\tchi_star\n")
for U,d,lab in cfgs:
    for o in obs:
        line=f"{lab:16s} {o:14s}"; v=[]
        for L in Ls:
            hit=[r for r in rows if r[0]==(U,d,L) and r[1]==o]
            if not hit: line+=f"{'-':>9s}"; v.append(float('nan')); continue
            cs=chistar(hit[0][2],1e-2); v.append(cs)
            line+=f"{('>256' if math.isnan(cs) else int(cs)):>9}"
            for tg in (1e-1,1e-2,1e-3):
                c2=chistar(hit[0][2],tg)
                out.write(f"{U}\t{d}\t{L}\t{o}\t{tg}\t{'NaN' if math.isnan(c2) else int(c2)}\n")
        r = v[3]/v[0] if (not math.isnan(v[0]) and not math.isnan(v[3]) and v[0]>0) else float('nan')
        line+=f"{r:>10.2f}" if not math.isnan(r) else f"{'-':>10s}"
        print(line)
out.close()
