# -*- coding: utf-8 -*-
"""Monte-Carlo confidence intervals for the attractor residual per regimen.
Runs R independent replicates (different seeds); reports mean and 95% CI.
Demonstrates that Isatuximab and Daratumumab are statistically indistinguishable
(the model encodes both as one antiCD38 node)."""
import random, itertools, csv
import numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
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

def analyze_fast(fault_ids,drugs,n_init):
    ov={}
    for fid in fault_ids:
        node,val,_=FAULTS[fid]; ov[node]=val
    Isa,Dara,Btz,Dex,Len=[int(bool(d)) for d in drugs]; antiCD38=1 if (Isa or Dara) else 0
    PTEN=ov.get('PTEN',1); Len_active=1 if (Len and ov.get('CRBN',1)) else 0; bim_del=(ov.get('BIM',1)==0)
    ovi={IDX[n]:v for n,v in ov.items() if n in IDX}
    resid=0.0
    for _ in range(n_init):
        st=[random.getrandbits(1) for _ in range(22)]; st[iAPOP]=0
        for i,v in ovi.items(): st[i]=v
        st=tuple(st); seen=set(); path=[]
        while st not in seen:
            seen.add(st); path.append(st)
            st=step(st,Btz,Dex,Len,antiCD38,PTEN,Len_active,bim_del)
            if ovi:
                l=list(st)
                for i,v in ovi.items(): l[i]=v
                st=tuple(l)
        cyc=path[path.index(st):]
        if any(x[iAPOP] for x in cyc): continue
        tot=0
        for x in cyc:
            ANTI=1 if (x[19] or x[20]) else 0; SURV=1 if (ANTI and not x[21]) else 0
            tot+=x[12]+x[13]+x[18]+x[10]+ANTI+x[0]+x[17]+SURV
        resid+=tot/(8*len(cyc))
    return resid/n_init

REGIMENS=[('Untreated',[0,0,0,0,0]),('Isatuximab',[1,0,0,0,0]),('Daratumumab',[0,1,0,0,0]),
    ('Bortezomib',[0,0,1,0,0]),('Dexamethasone',[0,0,0,1,0]),('Lenalidomide',[0,0,0,0,1]),
    ('Vd',[0,0,1,1,0]),('Rd',[0,0,0,1,1]),('VRd',[0,0,1,1,1]),('Dara+Rd',[0,1,0,1,1]),
    ('D-VRd',[0,1,1,1,1]),('Isa-VRd',[1,0,1,1,1])]
faults=list(range(1,N_FAULTS+1))
R=60; NINIT=60
M=np.zeros((R,len(REGIMENS)))
for r in range(R):
    random.seed(1000+r)
    for j,(_,vec) in enumerate(REGIMENS):
        M[r,j]=np.mean([analyze_fast([f],vec,NINIT) for f in faults])
mean=M.mean(0); lo=np.percentile(M,2.5,0); hi=np.percentile(M,97.5,0); sd=M.std(0)

with open('MM_attractor_ci.csv','w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['regimen','mean_residual','ci95_lo','ci95_hi','sd',f'R={R}',f'n_init={NINIT}'])
    for j,(nm,_) in enumerate(REGIMENS): w.writerow([nm,round(mean[j],4),round(lo[j],4),round(hi[j],4),round(sd[j],4)])

# Isa vs Dara
iI,iD=1,2
diff=M[:,iI]-M[:,iD]
print(f"Isatuximab : mean={mean[iI]:.4f}  95% CI [{lo[iI]:.4f}, {hi[iI]:.4f}]")
print(f"Daratumumab: mean={mean[iD]:.4f}  95% CI [{lo[iD]:.4f}, {hi[iD]:.4f}]")
print(f"Isa-Dara difference: mean={diff.mean():.4f}  95% CI [{np.percentile(diff,2.5):.4f}, {np.percentile(diff,97.5):.4f}]  (includes 0: {np.percentile(diff,2.5)<=0<=np.percentile(diff,97.5)})")
print(f"CIs overlap: {not (hi[iI]<lo[iD] or hi[iD]<lo[iI])}")

# figure
order=np.argsort(mean)
fig,ax=plt.subplots(figsize=(9,6))
for k,j in enumerate(order):
    nm=REGIMENS[j][0]; isacd=nm in ('Isatuximab','Daratumumab')
    c=('#C0392B' if isacd else '#1F3A5F')
    ax.errorbar(mean[j],k,xerr=[[max(0,mean[j]-lo[j])],[max(0,hi[j]-mean[j])]],fmt='o',color=c,
                ecolor=c,elinewidth=1.6,capsize=3,ms=6)
    ax.text(hi[j]+0.006,k,nm,va='center',fontsize=9,color=c,fontweight=('bold' if isacd else 'normal'))
ax.set_yticks([]); ax.set_xlabel('attractor residual burden  (mean over 22 single faults; lower = better)')
ax.set_title(f'Per-regimen residual with 95% CI  (R={R} independent Monte-Carlo replicates, n_init={NINIT})\n'
             'Isatuximab and Daratumumab (red) overlap — the model encodes them identically (antiCD38)')
ax.set_xlim(-0.02,max(hi)+0.06); ax.grid(axis='x',alpha=0.25)
ax.axhspan(min(np.where(np.isin(order,[iI,iD]))[0])-0.4,max(np.where(np.isin(order,[iI,iD]))[0])+0.4,color='#C0392B',alpha=0.05)
fig.tight_layout(); fig.savefig('MM_attractor_ci.png',dpi=150,bbox_inches='tight')
print('\nsaved MM_attractor_ci.png / MM_attractor_ci.csv')
