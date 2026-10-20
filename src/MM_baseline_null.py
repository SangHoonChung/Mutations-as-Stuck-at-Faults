# -*- coding: utf-8 -*-
"""§3.6(ii) Baselines and null models -- EXECUTED.
(1) drug-count baseline: does ranking skill reduce to "more drugs = better"?
(2) randomized-rule null: same topology, operators/signs shuffled.
All regimen scores use a SHARED SEED so that same-target agents are bit-identical."""
import itertools, random
import numpy as np
from MM_attractor import analyze
from MM_network import FAULTS, N_FAULTS

faults=list(range(1,N_FAULTS+1))
AG=['antiCD38','Btz','Dex','Len']              # effective agents (Isa==Dara)
PAT=list(itertools.product([0,1],repeat=4))
def vec(p): a,b,d,l=p; return [a,0,b,d,l]      # -> [Isa,Dara,Btz,Dex,Len]
def nm(p):  return '+'.join(x for x,on in zip(AG,p) if on) or 'Untreated'

def model_scores(seed=0, n_init=120):
    """mean attractor residual over the 22 single faults, SHARED seed per condition."""
    out=np.zeros(len(PAT))
    for j,p in enumerate(PAT):
        s=0.0
        for f in faults:
            random.seed(seed)                   # shared seed => identical inputs give identical outputs
            s+=analyze([f],vec(p),n_init=n_init)['residual']
        out[j]=s/len(faults)
    return out

def spearman(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float)
    ra=np.argsort(np.argsort(a)); rb=np.argsort(np.argsort(b))
    return float(np.corrcoef(ra,rb)[0,1])

M=model_scores(seed=0)
K=np.array([sum(p) for p in PAT],float)         # drug count
# baseline predicts: more drugs -> lower residual. Correlate model residual with -count.
rho_all=spearman(M,-K)
print("="*64)
print("(1) DRUG-COUNT BASELINE")
print("="*64)
print(f"Spearman rho(model residual, drug count)  = {spearman(M,K):+.3f}")
print(f"  (baseline 'more drugs = better' agreement: rho(model, -count) = {rho_all:+.3f})")
print()
print("Within-stratum discrimination (drug count is CONSTANT here -> baseline has zero resolving power):")
for k in range(0,5):
    idx=[j for j in range(len(PAT)) if sum(PAT[j])==k]
    if len(idx)<2: 
        if idx: print(f"  {k} drug(s): n=1  ({nm(PAT[idx[0]])} = {M[idx[0]]:.3f})")
        continue
    vals=M[idx]; order=np.argsort(vals)
    span=vals.max()-vals.min()
    best=nm(PAT[idx[order[0]]]); worst=nm(PAT[idx[order[-1]]])
    print(f"  {k} drug(s): n={len(idx)}  model spread = {span:.3f}  "
          f"(best {best} {vals.min():.3f} .. worst {worst} {vals.max():.3f})")
print()
print("Key cases the drug-count baseline CANNOT call (same count, model separates them):")
one=[(nm(p),M[j]) for j,p in enumerate(PAT) if sum(p)==1]
for n_,v in sorted(one,key=lambda t:t[1]): print(f"   1 drug: {n_:10s} {v:.3f}")
i3=[j for j,p in enumerate(PAT) if sum(p)==3]; i4=[j for j,p in enumerate(PAT) if sum(p)==4]
print(f"\n  best 3-drug = {M[i3].min():.3f}   vs   4-drug = {M[i4].min():.3f}  "
      f"-> {'3-drug ties/beats 4-drug (baseline says 4>3)' if M[i3].min()<=M[i4].min()+1e-9 else '4-drug better'}")
np.save('MM_real_scores.npy',M)   # consumed by MM_null_model.py
