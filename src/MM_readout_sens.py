# -*- coding: utf-8 -*-
"""Axis (4) readout-choice robustness. Attractors of the (unperturbed) dynamical
model are found ONCE per condition and cached; we then re-score living cells under
different readout definitions and compare the 12-regimen ranking to the baseline
linear 8-output residual by Spearman rho, also tracking whether the best regimen
and the top-3 set are preserved. Apoptotic cells contribute 0 under every readout."""
import random, numpy as np
from MM_network import FAULTS, N_FAULTS
from MM_rule_perturb import nxt, _params, REGIMENS, faults, IDX, iAPOP, spearman

def out8(x):
    ANTI=1 if (x[19] or x[20]) else 0; SURV=1 if (ANTI and not x[21]) else 0
    # MYC, CCND2, E2F, NFKB, ANTIAPOP, VEGF(=STAT3), CDK46, SURVIVAL
    return (x[12],x[13],x[18],x[10],ANTI,x[0],x[17],SURV)
NAMES=['MYC','CCND2','E2F','NFKB','ANTIAPOP','VEGF','CDK46','SURVIVAL']

def cache_attractors(n_init,seed):
    """returns living_cycles[j][f] = list over inits of (list of out8 vectors) or None if apoptotic."""
    random.seed(seed); cache={}
    for j,(nm,vec) in enumerate(REGIMENS):
        for f in faults:
            P,ovi=_params([f],vec); cells=[]
            for _ in range(n_init):
                st=[random.getrandbits(1) for _ in range(22)]; st[iAPOP]=0
                for i,v in ovi.items(): st[i]=v
                st=tuple(st); seen={}; path=[]
                while st not in seen:
                    seen[st]=len(path); path.append(st); st=nxt(st,P,{})
                    if ovi:
                        l=list(st)
                        for i,v in ovi.items(): l[i]=v
                        st=tuple(l)
                cyc=path[seen[st]:]
                if any(x[iAPOP] for x in cyc): cells.append(None)          # dead -> 0
                else: cells.append([out8(x) for x in cyc])
            cache[(j,f)]=cells
    return cache

def rank(cache,score):
    """score(outvec)->float in [0,1]; residual_cell=mean over cycle; per-regimen mean over faults/inits."""
    res=np.zeros(len(REGIMENS))
    for j in range(len(REGIMENS)):
        r=0.0
        for f in faults:
            cells=cache[(j,f)]; s=0.0
            for c in cells:
                s += 0.0 if c is None else sum(score(o) for o in c)/len(c)
            r+=s/len(cells)
        res[j]=r/len(faults)
    return res

if __name__=='__main__':
    import time,csv; t=time.time()
    NINIT=80; SEED=0
    cache=cache_attractors(NINIT,SEED)
    base=rank(cache, lambda o: sum(o)/8.0)
    order=np.argsort(base)                      # ascending residual (best first)
    best0=order[0]; top3_0=set(order[:3].tolist())
    def summarize(name,score):
        r=rank(cache,score); rho=spearman(base,r); od=np.argsort(r)
        best=od[0]; top3=set(od[:3].tolist())
        return name,rho,(best==best0),(top3==top3_0),r
    variants=[]
    for k in range(8):
        variants.append((f'drop {NAMES[k]}', (lambda o,k=k: sum(o[i] for i in range(8) if i!=k)/7.0)))
    variants += [
        ('squared (on/8)^2', lambda o:(sum(o)/8.0)**2),
        ('binary @0.5',      lambda o:1.0 if sum(o)/8.0>=0.5 else 0.0),
        ('proliferation-only(4)', lambda o:(o[0]+o[1]+o[2]+o[6])/4.0),
        ('survival-axis(4)',      lambda o:(o[3]+o[4]+o[5]+o[7])/4.0),
    ]
    print(f"baseline best regimen = {REGIMENS[best0][0]}; top-3 = {[REGIMENS[i][0] for i in order[:3]]}")
    print(f"\n{'readout variant':24s} {'rho':>6s} {'best?':>6s} {'top3?':>6s}")
    rows=[]; rhos=[]
    for nm,sc in variants:
        name,rho,bsame,t3same,r=summarize(nm,sc); rhos.append(rho); rows.append((name,rho,bsame,t3same))
        print(f"{name:24s} {rho:6.3f} {str(bsame):>6s} {str(t3same):>6s}")
    rhos=np.array(rhos)
    loo=rhos[:8]
    print(f"\nleave-one-out (8): min rho={loo.min():.3f}  median={np.median(loo):.3f}  best&top3 preserved all={all(r[2] and r[3] for r in rows[:8])}")
    print(f"all variants: min rho={rhos.min():.3f}  median={np.median(rhos):.3f}")
    with open('MM_readout_sens.csv','w',newline='') as fh:
        w=csv.writer(fh); w.writerow(['readout_variant','rho_vs_baseline','best_preserved','top3_preserved'])
        for name,rho,b,t3 in rows: w.writerow([name,round(rho,3),b,t3])
    print(f"runtime {time.time()-t:.1f}s -> MM_readout_sens.csv")
