# -*- coding: utf-8 -*-
"""
MM_network_gates.py  (v3)
The MM Boolean network (v3, intrinsic, no immune compartment) using ONLY logic
gates AND/OR/NOT, plus fault MUXes built from those gates. Verified identical to
MM_network.py.
"""
from MM_network import FAULTS

def AND(*a): return int(all(bool(x) for x in a))
def OR(*a):  return int(any(bool(x) for x in a))
def NOT(a):  return int(not bool(a))

NODE_FAULT = {node: (fid, val) for fid, (node, val, _lab) in FAULTS.items()}

def _ov(node, L, active):
    if node in NODE_FAULT:
        fid, v = NODE_FAULT[node]
        sel = 1 if fid in active else 0
        return OR(L, sel) if v == 1 else AND(L, NOT(sel))
    return L

def mm_outputs_gates(fault_set, drugs):
    active = set(fault_set)
    Isa, Dara, Btz, Dex, Len = [int(bool(d)) for d in drugs]
    antiCD38 = OR(Isa, Dara)

    IL6 = 0; IGF1 = 0; RANKL = 0; WNT = 0; HYPOXIA = 0; TNF_in = 0
    CRBN = _ov('CRBN', 1, active)
    PTEN = _ov('PTEN', 1, active)

    TNF   = _ov('TNF',   TNF_in, active)
    IGF1R = _ov('IGF1R', IGF1,   active)
    RANK  = _ov('RANK',  RANKL,  active)

    STAT3 = _ov('STAT3', IL6, active)
    PI3K  = _ov('PI3K',  OR(IGF1R, IL6, NOT(PTEN)), active)
    AKT   = _ov('AKT',   PI3K, active)
    MTOR  = _ov('MTOR',  AKT, active)
    ERK   = _ov('ERK',   IGF1R, active)
    BCAT  = _ov('BCAT',  WNT, active)

    NIK  = _ov('NIK', RANK, active)
    NFKB = _ov('NFKB', AND(OR(TNF, RANK, AKT, NIK), NOT(Btz), NOT(Dex)), active)

    Len_active = AND(Len, CRBN)
    IKZF = NOT(Len_active)
    IRF4 = _ov('IRF4', IKZF, active)
    MYC  = _ov('MYC', AND(OR(ERK, BCAT, STAT3, NFKB, MTOR), IRF4), active)

    CCND2 = _ov('CCND2', OR(MYC, ERK, BCAT, NFKB), active)
    MDM2  = OR(AKT, ERK)
    P53   = _ov('P53', NOT(MDM2), active)
    P21   = OR(P53, Len_active)
    CDK46 = _ov('CDK46', AND(CCND2, NOT(P21)), active)
    E2F   = _ov('E2F', CDK46, active)

    MCL1  = _ov('MCL1',  OR(STAT3, NFKB), active)
    BCLXL = _ov('BCLXL', OR(NFKB, STAT3), active)
    NOXA  = Btz
    BIM   = _ov('BIM', AND(Dex, NOT(STAT3), NOT(AKT)), active)
    aCD38 = antiCD38
    free_MCL1  = AND(MCL1, NOT(NOXA), NOT(BIM))
    free_BCLXL = AND(BCLXL, NOT(BIM))
    guarded = OR(free_MCL1, free_BCLXL)
    pro     = OR(NOXA, BIM, aCD38)
    APOP    = AND(pro, NOT(guarded))

    ANTIAPOP = OR(MCL1, BCLXL)
    VEGF     = OR(STAT3, HYPOXIA)
    SURVIVAL = AND(ANTIAPOP, NOT(APOP))

    return [MYC, CCND2, E2F, NFKB, ANTIAPOP, VEGF, CDK46, SURVIVAL]


if __name__ == '__main__':
    import itertools
    from MM_network import mm_outputs, N_FAULTS
    combos = list(itertools.product([0, 1], repeat=5))
    mism = ck = 0
    for f in range(0, N_FAULTS + 1):
        for c in combos:
            ck += 1
            if mm_outputs([f], list(c)) != mm_outputs_gates([f], list(c)): mism += 1
    for f1, f2 in itertools.combinations(range(1, N_FAULTS + 1), 2):
        for c in combos:
            ck += 1
            if mm_outputs([f1, f2], list(c)) != mm_outputs_gates([f1, f2], list(c)): mism += 1
    print(f'checked {ck} cases, mismatches = {mism}  ->  GATE==ORIGINAL: {mism==0}')
