# -*- coding: utf-8 -*-
"""Axis (1) update-scheme robustness: synchronous vs random-asynchronous update.
Paired comparison from identical initial states. Reports:
 (a) death/survival OUTCOME agreement (sync vs async), overall %,
 (b) regimen-ranking preservation (Spearman rho on mean residual and on apop prop),
 (c) apoptosis reachability preservation,
 plus the fixed-point coincidence check (a sync fixed point is an async fixed point).
Uses the optimized step() (same as MM_attractor_ci) so any sync/async difference is
purely the update scheme, not the implementation."""
import random, itertools
import numpy as np
from MM_network import FAULTS, N_FAULTS

NODES=['STAT3','PI3K','AKT','MTOR','ERK','BCAT','IGF1R','RANK','NIK','TNF','NFKB',
       'IRF4','MYC','CCND2','MDM2','P53','P21','CDK46','E2F','MCL1','BCLXL','APOP']
IDX={n:i for i,n in enumerate(NODES)}; iAPOP=IDX['APOP']

def step(s,Btz,Dex,Len,antiCD38,PTEN,Len_active,bim_del):
    STAT3,PI3K,AKT,MTOR,ERK,BCAT,IGF1R,RANK,NIK,TNF,NFKB,IRF4,MYC,CCND2,MDM2,P53,P21,CDK46,E2F,MCL1,BCLXL,APOP=s
    nPI3K=1 if (IGF1R or (not PTEN)) else 0
    nNFKB=1 if ((TNF or RANK or AKT or NIK) and (not Btz) and (not Dex)) else 0
    nIRF4=1 if ((not Len_active) and (NFKB or MYC)) else 0
    nMYC =1 if (IRF4 and (ERK or BCAT or STAT3 or NFKB or MTOR)) else 0
    nCCND2=1 if (MYC or ERK or BCAT or NFKB) else 0
    nMDM2=1 if (AKT or ERK) else 0
    nP53 =1 if (not MDM2) else 0
    nP21 =1 if (P53 or Len_active) else 0
    nCDK46=1 if (CCND2 and (not P21)) else 0
    nMCL1=1 if (STAT3 or NFKB or IRF4) else 0
    nBCLXL=1 if (NFKB or STAT3) else 0
    NOXA=Btz; BIM=1 if ((Dex and not (STAT3 or AKT)) and (not bim_del)) else 0
    guarded=(MCL1 and not (NOXA or BIM)) or (BCLXL and not BIM); pro=NOXA or BIM or antiCD38
    nAPOP=1 if (APOP or (pro and not guarded)) else 0
    return (0,nPI3K,PI3K,AKT,IGF1R,0,0,0,RANK,0,nNFKB,nIRF4,nMYC,nCCND2,nMDM2,nP53,nP21,nCDK46,CDK46,nMCL1,nBCLXL,nAPOP)

def outs_sum(x):
    ANTI=1 if (x[19] or x[20]) else 0; SURV=1 if (ANTI and not x[21]) else 0
    return x[12]+x[13]+x[18]+x[10]+ANTI+x[0]+x[17]+SURV   # 8 oncogenic outputs

def _params(fault_ids,drugs):
    ov={}
    for fid in fault_ids:
        node,val,_=FAULTS[fid]; ov[node]=val
    Isa,Dara,Btz,Dex,Len=[int(bool(d)) for d in drugs]; antiCD38=1 if (Isa or Dara) else 0
    PTEN=ov.get('PTEN',1); Len_active=1 if (Len and ov.get('CRBN',1)) else 0; bim_del=(ov.get('BIM',1)==0)
    ovi={IDX[n]:v for n,v in ov.items() if n in IDX}
    return (Btz,Dex,Len,antiCD38,PTEN,Len_active,bim_del),ovi

def sync_from(init,P,ovi):
    st=tuple(init); seen=set(); path=[]
    while st not in seen:
        seen.add(st); path.append(st)
        st=step(st,*P)
        if ovi:
            l=list(st)
            for i,v in ovi.items(): l[i]=v
            st=tuple(l)
    cyc=path[path.index(st):]
    apop=any(x[iAPOP] for x in cyc)
    if apop: return True,0.0,len(cyc)
    resid=sum(outs_sum(x) for x in cyc)/(8*len(cyc))
    return False,resid,len(cyc)

def async_from(init,P,ovi,burn,window):
    s=list(init); free=[i for i in range(22) if i not in ovi]
    for _ in range(burn):
        i=random.choice(free); full=step(tuple(s),*P); s[i]=full[i]
        if s[iAPOP]==1: return True,0.0
    acc=0.0
    for _ in range(window):
        i=random.choice(free); full=step(tuple(s),*P); s[i]=full[i]
        if s[iAPOP]==1: return True,0.0
        acc+=outs_sum(tuple(s))/8.0
    return False,acc/window

def spearman(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    ra=np.argsort(np.argsort(a)); rb=np.argsort(np.argsort(b))
    return float(np.corrcoef(ra,rb)[0,1])

REGIMENS=[('Untreated',[0,0,0,0,0]),('Isatuximab',[1,0,0,0,0]),('Daratumumab',[0,1,0,0,0]),
    ('Bortezomib',[0,0,1,0,0]),('Dexamethasone',[0,0,0,1,0]),('Lenalidomide',[0,0,0,0,1]),
    ('Vd',[0,0,1,1,0]),('Rd',[0,0,0,1,1]),('VRd',[0,0,1,1,1]),('Dara+Rd',[0,1,0,1,1]),
    ('D-VRd',[0,1,1,1,1]),('Isa-VRd',[1,0,1,1,1])]
faults=list(range(1,N_FAULTS+1))

def run(n_init=40,burn=250,window=120,seed=0):
    random.seed(seed)
    agree=0; tot=0; both_live_dabs=[]; both_live_pairs=[]
    sync_res=np.zeros(len(REGIMENS)); async_res=np.zeros(len(REGIMENS))
    sync_ap =np.zeros(len(REGIMENS)); async_ap =np.zeros(len(REGIMENS))
    death_sync_cond=0; death_kept=0
    for j,(nm,vec) in enumerate(REGIMENS):
        sres=ares=sap=aap=0.0
        for f in faults:
            P,ovi=_params([f],vec)
            s_ap_c=a_ap_c=0; s_r=a_r=0.0; s_rn=a_rn=0
            for _ in range(n_init):
                init=[random.getrandbits(1) for _ in range(22)]; init[iAPOP]=0
                for i,v in ovi.items(): init[i]=v
                sa,sr,_=sync_from(init,P,ovi)
                aa,ar=async_from(init,P,ovi,burn,window)
                tot+=1; agree+=(sa==aa)
                if sa: s_ap_c+=1
                else:  s_r+=sr; s_rn+=1
                if aa: a_ap_c+=1
                else:  a_r+=ar; a_rn+=1
                if (not sa) and (not aa): both_live_dabs.append(abs(sr-ar)); both_live_pairs.append((sr,ar))
            # per-condition (regimen,fault) aggregates
            s_apc=s_ap_c/n_init; a_apc=a_ap_c/n_init
            sap+=s_apc; aap+=a_apc
            sres+=(s_r/s_rn if s_rn else 0.0); ares+=(a_r/a_rn if a_rn else 0.0)
            if s_apc>=0.5:
                death_sync_cond+=1; death_kept+=(a_apc>=0.5)
        sync_res[j]=sres/len(faults); async_res[j]=ares/len(faults)
        sync_ap[j]=sap/len(faults);   async_ap[j]=aap/len(faults)
    dabs=np.array(both_live_dabs); pr=np.array(both_live_pairs)
    r_live=float(np.corrcoef(pr[:,0],pr[:,1])[0,1]) if len(pr)>2 else float('nan')
    return dict(
        outcome_agree=agree/tot, n_pairs=tot,
        rho_residual=spearman(sync_res,async_res), rho_apop=spearman(sync_ap,async_ap),
        live_mad=float(dabs.mean()), live_r=r_live, n_live=len(dabs),
        death_kept=death_kept, death_cond=death_sync_cond,
        sync_res=sync_res, async_res=async_res, sync_ap=sync_ap, async_ap=async_ap)

def stability():
    print("\n=== seed stability (n_init=40, burn=250, window=120) ===")
    keys=['outcome_agree','rho_residual','rho_apop','live_r','live_mad']
    acc={k:[] for k in keys}; dk=[]
    for sd in range(5):
        R=run(n_init=40,burn=250,window=120,seed=sd)
        for k in keys: acc[k].append(R[k])
        dk.append(R['death_kept']/R['death_cond'])
    for k in keys:
        v=np.array(acc[k]); print(f"  {k:14s} mean={v.mean():.3f}  range[{v.min():.3f},{v.max():.3f}]")
    dk=np.array(dk); print(f"  death_kept_frac  mean={dk.mean():.3f}  range[{dk.min():.3f},{dk.max():.3f}]")

def fixedpoint_check():
    """A synchronous fixed point is an asynchronous fixed point (definitional);
    confirmed on real fixed points collected from synchronous runs."""
    random.seed(7); fp=set(); cyc=[]
    for j,(nm,vec) in enumerate(REGIMENS):
        for f in faults:
            P,ovi=_params([f],vec)
            for _ in range(120):
                init=[random.getrandbits(1) for _ in range(22)]; init[iAPOP]=0
                for i,v in ovi.items(): init[i]=v
                st=tuple(init); seen={}; path=[]
                while st not in seen:
                    seen[st]=len(path); path.append(st); st=step(st,*P)
                    if ovi:
                        l=list(st)
                        for i,v in ovi.items(): l[i]=v
                        st=tuple(l)
                c=path[seen[st]:]; cyc.append(len(c))
                if len(c)==1: fp.add((c[0],j,f))
    viol=0
    for s,j,f in fp:
        P,ovi=_params([f],REGIMENS[j][1])
        for i in range(22):
            if i in ovi: continue
            if step(s,*P)[i]!=s[i]: viol+=1; break
    cl=np.array(cyc)
    print("\n=== fixed points / cycles (synchronous) ===")
    print(f"  fixed points={100*np.mean(cl==1):.1f}%  2-cycles={100*np.mean(cl==2):.1f}%  >2={100*np.mean(cl>2):.1f}%  max_len={cl.max()}")
    print(f"  distinct fixed points collected={len(fp)}; async-update violations={viol} (0 => all are async fixed points)")

if __name__=='__main__':
    import time, csv
    t=time.time(); R=run(n_init=40,burn=250,window=120,seed=0)
    print(f"runtime {time.time()-t:.1f}s  paired inits={R['n_pairs']}")
    print(f"(a) outcome agreement (sync vs async) = {100*R['outcome_agree']:.1f}%  | live r={R['live_r']:.3f} mad={R['live_mad']:.3f}")
    print(f"(b) ranking Spearman rho: residual={R['rho_residual']:.3f} apop={R['rho_apop']:.3f}")
    print(f"(c) death-dominant kept = {R['death_kept']}/{R['death_cond']}")
    with open('MM_async_robust.csv','w',newline='') as fh:
        w=csv.writer(fh); w.writerow(['regimen','sync_residual','async_residual','sync_apop','async_apop'])
        for j,(nm,_) in enumerate(REGIMENS):
            w.writerow([nm,round(R['sync_res'][j],4),round(R['async_res'][j],4),round(R['sync_ap'][j],4),round(R['async_ap'][j],4)])
    stability(); fixedpoint_check()
