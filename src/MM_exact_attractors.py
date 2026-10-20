# -*- coding: utf-8 -*-
"""EXACT attractor enumeration vs Monte-Carlo sampling.
Every attractor lies in the image of the transition map, where the environment
nodes are forced; the reachable subspace is 2^17, so the attractor set can be
enumerated EXHAUSTIVELY rather than estimated from n_init sampled initial states.
Start states are restricted to APOP=0 (alive), matching the paper's sampling."""
import time, random
import numpy as np
from MM_rule_perturb import nxt, _params, IDX, REGIMENS, faults
from MM_attractor import analyze

FREE=['PI3K','AKT','MTOR','ERK','NIK','NFKB','IRF4','MYC','CCND2','MDM2',
      'P53','P21','CDK46','E2F','MCL1','BCLXL']          # APOP fixed 0 at start
FREEI=[IDX[n] for n in FREE]
ENV=[IDX[n] for n in ['STAT3','BCAT','IGF1R','RANK','TNF']]
iAPOP=IDX['APOP']

def exact(fault_ids, drugs):
    """exact set of attractors reachable from ANY alive (APOP=0) initial state"""
    P,ovi=_params(fault_ids,drugs)
    base=[0]*22
    for i,v in ovi.items(): base[i]=v
    for i in ENV:
        if i not in ovi: base[i]=0
    seen={}; attractors=set()
    free=[i for i in FREEI if i not in ovi]
    for bits in range(1<<len(free)):
        s=list(base)
        for k,i in enumerate(free): s[i]=(bits>>k)&1
        if iAPOP not in ovi: s[iAPOP]=0
        st=tuple(s)
        if st in seen: continue
        path=[]; idx={}
        while st not in idx and st not in seen:
            idx[st]=len(path); path.append(st)
            st=nxt(st,P,{})
            if ovi:
                l=list(st)
                for i,v in ovi.items(): l[i]=v
                st=tuple(l)
        if st in idx:
            cyc=frozenset(path[idx[st]:]); attractors.add(cyc)
            for x in path: seen[x]=cyc
        else:
            c=seen[st]
            for x in path: seen[x]=c
    return attractors, len(free)

if __name__=='__main__':
    t=time.time()
    ex_ms=0; sm_ms=0; tot=0; rows=[]
    n_ex=[]; n_sm=[]
    for rn,vec in REGIMENS:
        for f in faults:
            A,nfree=exact([f],vec)
            random.seed(0)
            r=analyze([f],vec,n_init=120)
            ne=len(A); ns=r['n_attractors']
            ex_ms+=(ne>1); sm_ms+=(ns>1); tot+=1
            n_ex.append(ne); n_sm.append(ns)
            rows.append((rn,f,ne,ns))
    n_ex=np.array(n_ex); n_sm=np.array(n_sm)
    print("="*66)
    print("EXACT attractor enumeration  vs  Monte-Carlo sampling (n_init=120)")
    print("="*66)
    print(f"conditions: {tot}   (state subspace enumerated per condition: 2^16 = 65,536 alive states)")
    print(f"multistable (>1 attractor)   EXACT: {100*ex_ms/tot:5.1f}%    SAMPLED: {100*sm_ms/tot:5.1f}%")
    print(f"mean #attractors             EXACT: {n_ex.mean():5.2f}      SAMPLED: {n_sm.mean():5.2f}")
    print(f"max  #attractors             EXACT: {n_ex.max():5d}      SAMPLED: {n_sm.max():5d}")
    print(f"conditions where sampling UNDERCOUNTS: {100*np.mean(n_sm<n_ex):.1f}%")
    print(f"conditions where they agree          : {100*np.mean(n_sm==n_ex):.1f}%")
    print(f"\nruntime {time.time()-t:.0f}s")
    import csv
    with open('MM_exact_attractors.csv','w',newline='') as fh:
        w=csv.writer(fh); w.writerow(['regimen','fault','n_attractors_exact','n_attractors_sampled'])
        w.writerows(rows)
