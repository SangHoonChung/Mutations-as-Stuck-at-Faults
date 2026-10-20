# -*- coding: utf-8 -*-
"""Enabled analysis: automatic test-pattern generation (ATPG) for combination
discovery. A drug combination is a test pattern; a fault is 'detected' (the
malignant cell is killed) when the pattern drives the absorbing APOP state from
every sampled initial condition (apop_propensity = 1). Because the drug-input
space is tiny (4 effective agents: antiCD38, Btz, Dex, Len -> 16 patterns), we
run EXHAUSTIVE ATPG: the search is complete, not heuristic. We report the minimal
detecting pattern per fault, the undetectable (intrinsically resistant) faults,
and the minimal single pattern that covers all detectable faults."""
import itertools, numpy as np
from MM_attractor import analyze
from MM_network import FAULTS, FAULT_LABELS, N_FAULTS

AGENTS=['antiCD38','Btz','Dex','Len']          # Isa==Dara -> one antiCD38 agent
def to_drugvec(pat):                            # pat=(a,b,d,l) -> [Isa,Dara,Btz,Dex,Len]
    a,b,d,l=pat; return [a,0,b,d,l]
PATTERNS=list(itertools.product([0,1],repeat=4))
faults=list(range(1,N_FAULTS+1))
NINIT=120; TAU=1.0                              # detected = kill from ALL sampled inits

# apop propensity matrix: faults x patterns
AP=np.zeros((len(faults),len(PATTERNS)))
for i,f in enumerate(faults):
    for j,pat in enumerate(PATTERNS):
        AP[i,j]=analyze([f],to_drugvec(pat),n_init=NINIT)['apop_propensity']

def card(pat): return sum(pat)
def name(pat): return '+'.join(a for a,on in zip(AGENTS,pat) if on) or 'none'

# minimal detecting pattern per fault
print(f"{'fault':34s} {'min kill pattern':22s} {'#drugs':>6s} {'maxAP':>6s}")
undetectable=[]; per_fault_min=[]
for i,f in enumerate(faults):
    det=[j for j in range(len(PATTERNS)) if AP[i,j]>=TAU]
    if not det:
        undetectable.append(f); per_fault_min.append(None)
        print(f"{FAULT_LABELS[f]:34s} {'-- UNDETECTABLE --':22s} {'':>6s} {AP[i].max():6.2f}")
    else:
        jmin=min(det,key=lambda j:card(PATTERNS[j])); per_fault_min.append(jmin)
        print(f"{FAULT_LABELS[f]:34s} {name(PATTERNS[jmin]):22s} {card(PATTERNS[jmin]):6d} {AP[i].max():6.2f}")

detectable=[f for f in faults if f not in undetectable]
print(f"\ndetectable faults: {len(detectable)}/{len(faults)}   undetectable: {len(undetectable)}")
print("undetectable (intrinsic resistance):", [FAULT_LABELS[f] for f in undetectable])

# minimal SINGLE pattern covering all detectable faults
best=None
for j,pat in enumerate(PATTERNS):
    covered=sum(1 for i,f in enumerate(faults) if f in detectable and AP[i,j]>=TAU)
    if covered==len(detectable):
        if best is None or card(pat)<card(PATTERNS[best]): best=j
print(f"\nminimal single pattern covering ALL {len(detectable)} detectable faults:",
      (name(PATTERNS[best])+f"  ({card(PATTERNS[best])} drugs)") if best is not None else "NONE (no single pattern covers all)")

# coverage of each pattern (how many faults each kills) - top patterns
cov=[(name(PATTERNS[j]),card(PATTERNS[j]),int((AP[:,j]>=TAU).sum())) for j in range(len(PATTERNS))]
cov.sort(key=lambda t:(-t[2],t[1]))
print("\ntop patterns by fault coverage (kills = apop=1):")
for nm,k,c in cov[:6]: print(f"  {nm:22s} drugs={k} kills={c}/{len(faults)}")

# distribution of minimal cardinality among detectable
mins=[card(PATTERNS[j]) for j in per_fault_min if j is not None]
print(f"\nminimal-kill cardinality among detectable faults: "
      f"1-drug={mins.count(1)}, 2-drug={mins.count(2)}, 3-drug={mins.count(3)}, 4-drug={mins.count(4)}")

import csv
with open('MM_atpg.csv','w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['fault','min_kill_pattern','n_drugs','max_apop_propensity','detectable'])
    for i,f in enumerate(faults):
        j=per_fault_min[i]
        w.writerow([FAULT_LABELS[f], (name(PATTERNS[j]) if j is not None else 'UNDETECTABLE'),
                    (card(PATTERNS[j]) if j is not None else ''), round(AP[i].max(),3), j is not None])
print("\n-> MM_atpg.csv")
