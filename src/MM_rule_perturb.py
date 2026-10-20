# -*- coding: utf-8 -*-
"""Axis (3) rule/structure-perturbation robustness.
For the dynamical model we (a) remove feedback edges one at a time
(MYC->IRF4, IRF4->MYC, IRF4->MCL1, APOP absorbing latch) and (b) flip OR<->AND at
fused convergence nodes, one perturbation at a time. For each we recompute the
12-regimen ranking (mean attractor residual over the 22 single faults) and the
multistability fraction (conditions with >1 attractor), and compare to baseline
by Spearman rho. Honest goal: reveal which edges are load-bearing."""
import random, numpy as np
from MM_network import FAULTS, N_FAULTS

NODES=['STAT3','PI3K','AKT','MTOR','ERK','BCAT','IGF1R','RANK','NIK','TNF','NFKB',
       'IRF4','MYC','CCND2','MDM2','P53','P21','CDK46','E2F','MCL1','BCLXL','APOP']
IDX={n:i for i,n in enumerate(NODES)}; iAPOP=IDX['APOP']
OR=lambda *a:1 if any(a) else 0; AND=lambda *a:1 if all(a) else 0

def nxt(s,P,cfg):
    Btz,Dex,Len,antiCD38,PTEN,Len_active,bim_del=P
    (STAT3,PI3K,AKT,MTOR,ERK,BCAT,IGF1R,RANK,NIK,TNF,NFKB,IRF4,MYC,CCND2,
     MDM2,P53,P21,CDK46,E2F,MCL1,BCLXL,APOP)=s
    fl=cfg.get('flip',set()); rm=cfg.get('rm',set())
    def orand(node,args): return AND(*args) if node in fl else OR(*args)
    nPI3K=orand('PI3K',(IGF1R,(not PTEN)))
    nAKT=PI3K; nMTOR=AKT; nERK=IGF1R; nNIK=RANK
    orNFKB=orand('NFKB',(TNF,RANK,AKT,NIK)); nNFKB=AND(orNFKB,(not Btz),(not Dex))
    # IRF4: base ~Len_active & (NFKB|MYC)
    irf4_src = (NFKB,) if 'MYC->IRF4' in rm else (NFKB,MYC)
    nIRF4 = AND((not Len_active), orand('IRF4',irf4_src)) if len(irf4_src)>1 else AND((not Len_active),NFKB)
    # MYC: base IRF4 & (ERK|BCAT|STAT3|NFKB|MTOR)
    orMYC=orand('MYC',(ERK,BCAT,STAT3,NFKB,MTOR))
    nMYC = orMYC if 'IRF4->MYC' in rm else AND(IRF4,orMYC)
    nCCND2=orand('CCND2',(MYC,ERK,BCAT,NFKB))
    nMDM2=orand('MDM2',(AKT,ERK)); nP53=1 if not MDM2 else 0
    nP21=orand('P21',(P53,Len_active)); nCDK46=(OR(CCND2,(not P21)) if 'CDK46' in fl else AND(CCND2,(not P21))); nE2F=CDK46
    mcl_src=(STAT3,NFKB) if 'IRF4->MCL1' in rm else (STAT3,NFKB,IRF4)
    nMCL1=orand('MCL1',mcl_src); nBCLXL=orand('BCLXL',(NFKB,STAT3))
    NOXA=Btz; BIM=AND(Dex,(not (STAT3 or AKT)),(not bim_del))
    free_MCL1=AND(MCL1,(not (NOXA or BIM))); free_BCLXL=AND(BCLXL,(not BIM))
    guarded=orand('guarded',(free_MCL1,free_BCLXL)); pro=orand('pro',(NOXA,BIM,antiCD38))
    commit=AND(pro,(not guarded))
    nAPOP = commit if 'APOP_latch' in rm else OR(APOP,commit)
    # env-driven nodes forced to 0 (IL6=WNT=IGF1=RANKL=TNF_in=0); faults re-pinned by ovi
    return (0,nPI3K,nAKT,nMTOR,nERK,0,0,0,nNIK,0,nNFKB,nIRF4,nMYC,
            nCCND2,nMDM2,nP53,nP21,nCDK46,nE2F,nMCL1,nBCLXL,nAPOP)

def outs_sum(x):
    ANTI=OR(x[19],x[20]); SURV=AND(ANTI,(not x[21]))
    return x[12]+x[13]+x[18]+x[10]+ANTI+x[0]+x[17]+SURV

def _params(fault_ids,drugs):
    ov={}
    for fid in fault_ids: node,val,_=FAULTS[fid]; ov[node]=val
    Isa,Dara,Btz,Dex,Len=[int(bool(d)) for d in drugs]; antiCD38=OR(Isa,Dara)
    PTEN=ov.get('PTEN',1); Len_active=AND(Len,ov.get('CRBN',1)); bim_del=(ov.get('BIM',1)==0)
    ovi={IDX[n]:v for n,v in ov.items() if n in IDX}
    return (Btz,Dex,Len,antiCD38,PTEN,Len_active,bim_del),ovi

def analyze(fault_ids,drugs,cfg,n_init):
    P,ovi=_params(fault_ids,drugs); resid=0.0; attr=set(); apop=0
    for _ in range(n_init):
        st=[random.getrandbits(1) for _ in range(22)]; st[iAPOP]=0
        for i,v in ovi.items(): st[i]=v
        st=tuple(st); seen={}; path=[]
        while st not in seen:
            seen[st]=len(path); path.append(st); st=nxt(st,P,cfg)
            if ovi:
                l=list(st)
                for i,v in ovi.items(): l[i]=v
                st=tuple(l)
        cyc=path[seen[st]:]; attr.add(frozenset(cyc))
        if any(x[iAPOP] for x in cyc): apop+=1; continue
        resid+=sum(outs_sum(x) for x in cyc)/(8*len(cyc))
    ap=apop/n_init
    return resid/n_init, len(attr), ap

REGIMENS=[('Untreated',[0,0,0,0,0]),('Isatuximab',[1,0,0,0,0]),('Daratumumab',[0,1,0,0,0]),
    ('Bortezomib',[0,0,1,0,0]),('Dexamethasone',[0,0,0,1,0]),('Lenalidomide',[0,0,0,0,1]),
    ('Vd',[0,0,1,1,0]),('Rd',[0,0,0,1,1]),('VRd',[0,0,1,1,1]),('Dara+Rd',[0,1,0,1,1]),
    ('D-VRd',[0,1,1,1,1]),('Isa-VRd',[1,0,1,1,1])]
faults=list(range(1,N_FAULTS+1))

def spearman(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    ra=np.argsort(np.argsort(a)); rb=np.argsort(np.argsort(b))
    return float(np.corrcoef(ra,rb)[0,1])

def screen(cfg,n_init,seed):
    random.seed(seed); res=np.zeros(len(REGIMENS)); fate_ms=0; attr_ms=0; tot=0
    for j,(nm,vec) in enumerate(REGIMENS):
        r=0.0
        for f in faults:
            rr,na,ap=analyze([f],vec,cfg,n_init); r+=rr
            fate_ms+=(0.0<ap<1.0); attr_ms+=(na>1); tot+=1
        res[j]=r/len(faults)
    return res, fate_ms/tot, attr_ms/tot

PERTURB=[
 ('baseline',dict()),
 ('rm MYC->IRF4 (open loop)',{'rm':{'MYC->IRF4'}}),
 ('rm IRF4->MYC (open loop)',{'rm':{'IRF4->MYC'}}),
 ('rm IRF4->MCL1',{'rm':{'IRF4->MCL1'}}),
 ('rm APOP latch',{'rm':{'APOP_latch'}}),
 ('flip PI3K OR->AND',{'flip':{'PI3K'}}),
 ('flip NFKB OR->AND',{'flip':{'NFKB'}}),
 ('flip MYC OR->AND',{'flip':{'MYC'}}),
 ('flip CCND2 OR->AND',{'flip':{'CCND2'}}),
 ('flip MDM2 OR->AND',{'flip':{'MDM2'}}),
 ('flip P21 OR->AND',{'flip':{'P21'}}),
 ('flip CDK46 AND->OR',{'flip':{'CDK46'}}),
 ('flip MCL1 OR->AND',{'flip':{'MCL1'}}),
 ('flip BCLXL OR->AND',{'flip':{'BCLXL'}}),
 ('flip IRF4 OR->AND',{'flip':{'IRF4'}}),
 ('flip guarded OR->AND',{'flip':{'guarded'}}),
 ('flip pro OR->AND',{'flip':{'pro'}}),
]

if __name__=='__main__':
    import time; t=time.time()
    NINIT=80; SEED=0
    base_res,base_fms,base_ams=screen(PERTURB[0][1],NINIT,SEED)
    print(f"baseline fate-multistability = {100*base_fms:.1f}%  (attractor-count>1 = {100*base_ams:.1f}%)  [rt {time.time()-t:.1f}s]")
    print(f"\n{'perturbation':26s} {'rho_rank':>8s} {'fate_ms%':>9s} {'note':>10s}")
    print(f"{'baseline':26s} {1.000:8.3f} {100*base_fms:9.1f}")
    rhos=[]
    for nm,cfg in PERTURB[1:]:
        res,fms,ams=screen(cfg,NINIT,SEED); rho=spearman(base_res,res); rhos.append(rho)
        note='FATE-MS COLLAPSE' if fms< max(0.02, 0.25*base_fms) else ''
        print(f"{nm:26s} {rho:8.3f} {100*fms:9.1f}  {note}")
    rhos=np.array(rhos)
    print(f"\nsingle-perturbation ranking Spearman rho: median={np.median(rhos):.3f}  "
          f"min={rhos.min():.3f}  frac>=0.9={np.mean(rhos>=0.9):.2f}")
    import csv
    with open('MM_rule_perturb.csv','w',newline='') as fh:
        w=csv.writer(fh); w.writerow(['perturbation','rho_rank_vs_baseline','fate_multistab_pct'])
        w.writerow(['baseline',1.000,round(100*base_fms,1)])
        for (nm,cfg),rho in zip(PERTURB[1:],rhos):
            _,fms,_=screen(cfg,NINIT,SEED); w.writerow([nm,round(rho,3),round(100*fms,1)])
    print(f"runtime {time.time()-t:.1f}s  -> MM_rule_perturb.csv")
