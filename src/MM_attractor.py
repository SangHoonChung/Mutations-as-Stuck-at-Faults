# -*- coding: utf-8 -*-
"""
MM_attractor.py  -- dynamical (attractor-based) version of the MM Boolean model.

Differences vs the single-pass DAG (MM_network.py):
  * Feedback restored: IRF4<->MYC positive autoregulatory loop (core MM dependency);
    MCL1 is sustained by IRF4 (so survival is coupled to the loop).
  * Apoptosis is an ABSORBING latch: APOP(t+1) = APOP(t) OR (pro AND not guarded)
    -> commitment to death is irreversible (a death fixed point).
  * State is evolved by SYNCHRONOUS updates from sampled initial conditions until it
    reaches an attractor (fixed point or cycle). Multistability => the same condition
    can reach different attractors depending on the initial state.

Readouts per (faults, drugs):
  apop_propensity = fraction of sampled initial states reaching an apoptotic attractor
  residual        = mean over inits of [0 if apoptotic else (malignant outputs ON / 8)]
                    (dead cells contribute 0; living cells contribute their program load)
  n_attractors    = number of distinct attractors found (multistability indicator)

Faults reuse MM_network.FAULTS (same 22) so results are comparable.
"""
import random
from MM_network import FAULTS, FAULT_LABELS, N_FAULTS

random.seed(0)
NODES = ['STAT3','PI3K','AKT','MTOR','ERK','BCAT','IGF1R','RANK','NIK','TNF','NFKB',
         'IRF4','MYC','CCND2','MDM2','P53','P21','CDK46','E2F','MCL1','BCLXL','APOP']
IDX = {n:i for i,n in enumerate(NODES)}
OUTS = ['MYC','CCND2','E2F','NFKB','ANTIAPOP','VEGF','CDK46','SURVIVAL']

def _next(s, drugs, ov, Len_active, bim_del):
    Isa,Dara,Btz,Dex,Len = drugs
    antiCD38 = Isa or Dara
    g = lambda n: s[IDX[n]]
    def F(name, v):
        return ov[name] if name in ov else int(bool(v))
    IL6=IGF1=RANKL=WNT=HYPOXIA=TNF_in=0
    PTEN = ov.get('PTEN',1)
    n = [0]*len(NODES)
    n[IDX['IGF1R']] = F('IGF1R', IGF1)
    n[IDX['RANK']]  = F('RANK',  RANKL)
    n[IDX['TNF']]   = F('TNF',   TNF_in)
    n[IDX['STAT3']] = F('STAT3', IL6)
    n[IDX['PI3K']]  = F('PI3K',  g('IGF1R') or IL6 or (not PTEN))
    n[IDX['AKT']]   = F('AKT',   g('PI3K'))
    n[IDX['MTOR']]  = F('MTOR',  g('AKT'))
    n[IDX['ERK']]   = F('ERK',   g('IGF1R'))
    n[IDX['BCAT']]  = F('BCAT',  WNT)
    n[IDX['NIK']]   = F('NIK',   g('RANK'))
    n[IDX['NFKB']]  = F('NFKB', (g('TNF') or g('RANK') or g('AKT') or g('NIK')) and (not Btz) and (not Dex))
    # ---- IRF4 <-> MYC positive feedback loop ----
    n[IDX['IRF4']] = F('IRF4', (not Len_active) and (g('NFKB') or g('MYC')))
    n[IDX['MYC']]  = F('MYC',  g('IRF4') and (g('ERK') or g('BCAT') or g('STAT3') or g('NFKB') or g('MTOR')))
    n[IDX['CCND2']]= F('CCND2', g('MYC') or g('ERK') or g('BCAT') or g('NFKB'))
    n[IDX['MDM2']] = (g('AKT') or g('ERK'))
    n[IDX['P53']]  = F('P53', not g('MDM2'))
    n[IDX['P21']]  = (g('P53') or Len_active)
    n[IDX['CDK46']]= F('CDK46', g('CCND2') and (not g('P21')))
    n[IDX['E2F']]  = F('E2F', g('CDK46'))
    n[IDX['MCL1']] = F('MCL1', g('STAT3') or g('NFKB') or g('IRF4'))   # IRF4 sustains MCL1
    n[IDX['BCLXL']]= F('BCLXL', g('NFKB') or g('STAT3'))
    # ---- apoptosis (absorbing latch) ----
    NOXA = Btz
    BIM  = (Dex and not (g('STAT3') or g('AKT'))) and (not bim_del)
    aCD38 = antiCD38
    free_MCL1  = g('MCL1') and not (NOXA or BIM)
    free_BCLXL = g('BCLXL') and not BIM
    guarded = free_MCL1 or free_BCLXL
    pro = NOXA or BIM or aCD38
    n[IDX['APOP']] = int(g('APOP') or (pro and not guarded))
    return tuple(int(bool(x)) for x in n)

def _setup(fault_set, drugs):
    ov = {}
    for fid in fault_set:
        if fid in FAULTS:
            node,val,_ = FAULTS[fid]; ov[node]=val
    Isa,Dara,Btz,Dex,Len = [int(bool(d)) for d in drugs]
    Len_active = int(Len and ov.get('CRBN',1))
    bim_del = (ov.get('BIM',1)==0)
    return ov, (Isa,Dara,Btz,Dex,Len), Len_active, bim_del

def find_attractor(init, drugs2, ov, Len_active, bim_del):
    seen={}; s=init; path=[]
    while s not in seen:
        seen[s]=len(path); path.append(s)
        s=_next(s, drugs2, ov, Len_active, bim_del)
    start=seen[s]
    return path[start:]            # list of states in the attractor cycle

def _outputs(state):
    g=lambda n: state[IDX[n]]
    ANTIAPOP = g('MCL1') or g('BCLXL')
    VEGF     = g('STAT3')          # HYPOXIA=0
    SURVIVAL = ANTIAPOP and not g('APOP')
    return [g('MYC'),g('CCND2'),g('E2F'),g('NFKB'),int(ANTIAPOP),int(VEGF),g('CDK46'),int(SURVIVAL)]

def analyze(fault_set, drugs, n_init=120):
    ov, drugs2, Len_active, bim_del = _setup(fault_set, drugs)
    apop=0; resid=0.0; attractors=set()
    for _ in range(n_init):
        init=[random.randint(0,1) for _ in NODES]; init[IDX['APOP']]=0      # start alive
        for fnode,fval in ov.items():
            if fnode in IDX: init[IDX[fnode]]=fval
        att=find_attractor(tuple(init), drugs2, ov, Len_active, bim_del)
        attractors.add(frozenset(att))
        is_apop = any(st[IDX['APOP']]==1 for st in att)
        if is_apop:
            apop+=1
        else:
            load=sum(sum(_outputs(st)) for st in att)/(8*len(att))
            resid+=load
    return {'apop_propensity':apop/n_init, 'residual':resid/n_init, 'n_attractors':len(attractors)}

if __name__=='__main__':
    cases=[
        ('baseline (healthy, untreated)', [0], [0,0,0,0,0]),
        ('NFKB fault, untreated',         [5], [0,0,0,0,0]),
        ('NFKB fault + Len',              [5], [0,0,0,0,1]),
        ('NFKB fault + VRd',              [5], [0,0,1,1,1]),
        ('ERK fault, untreated (loop)',   [3], [0,0,0,0,0]),
        ('ERK fault + Len',               [3], [0,0,0,0,1]),
        ('MCL1-high fault + Btz',         [17],[0,0,1,0,0]),
        ('MCL1-high fault + Dex',         [17],[0,0,0,1,0]),
        ('BCLXL-high fault + VRd',        [18],[0,0,1,1,1]),
    ]
    print(f"{'case':34s} {'apop':>6s} {'residual':>9s} {'#attr':>6s}")
    for lab,f,d in cases:
        r=analyze(f,d)
        print(f"{lab:34s} {r['apop_propensity']:6.2f} {r['residual']:9.3f} {r['n_attractors']:6d}")
