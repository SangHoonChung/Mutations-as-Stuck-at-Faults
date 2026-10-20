# -*- coding: utf-8 -*-
"""Regenerate Table I with a SHARED SEED per condition, so that regimens which the
model encodes identically (Isa vs Dara -> the same antiCD38 input) receive
bit-identical scores instead of Monte-Carlo noise presented as a ranking."""
import itertools, random
import numpy as np
from MM_attractor import analyze
from MM_network import FAULTS, N_FAULTS

faults=list(range(1,N_FAULTS+1))
TABLE=[('Btz+Dex+Len (VRd)',       [0,0,1,1,1], 3),
       ('Isa+Btz+Dex+Len',         [1,0,1,1,1], 4),
       ('Dara+Btz+Dex+Len',        [0,1,1,1,1], 4),
       ('Isa+Btz+Dex',             [1,0,1,1,0], 3),
       ('Btz+Dex',                 [0,0,1,1,0], 2)]

def mean_resid(drugs, fsets, n_init=120):
    s=0.0
    for fs in fsets:
        random.seed(0)                      # SHARED seed across regimens
        s+=analyze(list(fs),drugs,n_init=n_init)['residual']
    return s/len(fsets)

singles=[(f,) for f in faults]
doubles=list(itertools.combinations(faults,2))
triples=list(itertools.combinations(faults,3))
print(f"fault sets: single={len(singles)} double={len(doubles)} triple={len(triples)}")
print(f"\n{'Regimen':22s} {'n':>2s} {'single':>7s} {'double':>7s} {'triple':>7s}")
res={}
for name,vecd,k in TABLE:
    s1=mean_resid(vecd,singles); s2=mean_resid(vecd,doubles); s3=mean_resid(vecd,triples)
    res[name]=(k,s1,s2,s3)
    print(f"{name:22s} {k:2d} {s1:7.3f} {s2:7.3f} {s3:7.3f}")
a=res['Isa+Btz+Dex+Len']; b=res['Dara+Btz+Dex+Len']
print("\n=== the contradiction the reviewer flagged ===")
print(f"  Isa-quad  : {a[1]:.4f} / {a[2]:.4f} / {a[3]:.4f}")
print(f"  Dara-quad : {b[1]:.4f} / {b[2]:.4f} / {b[3]:.4f}")
print(f"  identical now? {a[1:]==b[1:]}  -> the two rows should be MERGED as 'anti-CD38+VRd'")
import csv
with open('MM_table1_shared_seed.csv','w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['regimen','n_drugs','single','double','triple'])
    for n,(k,s1,s2,s3) in res.items(): w.writerow([n,k,round(s1,4),round(s2,4),round(s3,4)])
