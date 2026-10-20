# -*- coding: utf-8 -*-
"""
MM_network.py  (v3 - intrinsic, no immune compartment)
Boolean-network core for a Multiple Myeloma (MM) drug-combination screen.

Design choices (per user request):
  * NO immune compartment (no ADCC/CDC). anti-CD38 is modeled ONLY as its
    cell-intrinsic, Fc-independent direct pro-apoptotic effect (subject to the
    same BH3 anti-apoptotic guards). Its dominant clinical mechanism (ADCC/CDC)
    is cell-EXTRINSIC and is explicitly OUT OF SCOPE for this intrinsic model.
  * NO "death zeroes all outputs" collapse rule. Oncogenic outputs are turned
    off only to the extent a drug acts on the corresponding network node, so the
    score reflects pathway suppression (this restores a meaningful spread rather
    than many trivial zeros). Apoptosis is reflected in the SURVIVAL output.
  * Bortezomib vs Dexamethasone are still separated by BH3 logic
    (NOXA->Mcl-1 ; BIM->Bcl-xL/Bcl-2, blunted by IL-6/AKT).

5 drugs (binary): x1 Isatuximab, x2 Daratumumab, x3 Bortezomib, x4 Dexamethasone, x5 Lenalidomide
A "fault" forces one node to its oncogenic value. fault 0 = baseline (healthy, all outputs 0).
measure = ((b+c)^2)/((a+b+c+d)^2) with ideal output all zeros -> (oncogenic_outputs_ON / 8)^2.
Every rule is an explicit, editable modeling choice; not clinically validated.
"""

# id -> (node, forced value, biological label)
FAULTS = {
    1:  ('STAT3', 1, 'STAT3 constitutive (IL-6/JAK)'),
    2:  ('AKT',   1, 'AKT activation'),
    3:  ('ERK',   1, 'RAS/MAPK (KRAS/NRAS/BRAF)'),
    4:  ('BCAT',  1, 'Wnt / beta-catenin'),
    5:  ('NFKB',  1, 'NF-kB constitutive'),
    6:  ('NIK',   1, 'TRAF3 loss / NIK amp (noncanon.)'),
    7:  ('RANK',  1, 'RANKL axis (bone)'),
    8:  ('TNF',   1, 'NF-kB upstream (TNF)'),
    9:  ('MYC',   1, 'MYC rearrangement'),
    10: ('IRF4',  1, 'IRF4 high'),
    11: ('CCND2', 1, 'CCND t(11;14)/t(14;16)/t(4;14)'),
    12: ('CDK46', 1, 'CDK4/6 activation'),
    13: ('E2F',   1, 'RB1 loss'),
    14: ('MTOR',  1, 'mTOR activation'),
    15: ('PI3K',  1, 'PI3K activation'),
    16: ('IGF1R', 1, 'IGF-1 axis'),
    17: ('MCL1',  1, '1q21 gain / MCL1 amp'),
    18: ('BCLXL', 1, 'BCL-XL/BCL-2 high'),
    19: ('P53',   0, 'TP53 loss (del17p)'),
    20: ('PTEN',  0, 'PTEN loss'),
    21: ('BIM',   0, 'BIM (BCL2L11) deletion'),
    22: ('CRBN',  0, 'CRBN loss (Len resistance)'),
}
N_FAULTS = len(FAULTS)
FAULT_LABELS = {fid: lab for fid, (n, v, lab) in FAULTS.items()}
OUTPUT_NAMES = ['MYC', 'CCND2', 'E2F', 'NFKB', 'ANTIAPOP', 'VEGF', 'CDK46', 'SURVIVAL']


def mm_outputs(fault_set, drugs):
    Isa, Dara, Btz, Dex, Len = [int(bool(d)) for d in drugs]
    antiCD38 = int(Isa or Dara)

    forced = {}
    for fid in fault_set:
        if fid in FAULTS:
            node, val, _ = FAULTS[fid]
            forced[node] = val

    def F(name, val):
        return forced.get(name, int(bool(val)))

    IL6 = 0; IGF1 = 0; RANKL = 0; WNT = 0; HYPOXIA = 0; TNF_in = 0
    CRBN = F('CRBN', 1)
    PTEN = F('PTEN', 1)

    TNF   = F('TNF',   TNF_in)
    IGF1R = F('IGF1R', IGF1)
    RANK  = F('RANK',  RANKL)

    STAT3 = F('STAT3', IL6)
    PI3K  = F('PI3K',  IGF1R or IL6 or (not PTEN))
    AKT   = F('AKT',   PI3K)
    MTOR  = F('MTOR',  AKT)
    ERK   = F('ERK',   IGF1R)
    BCAT  = F('BCAT',  WNT)

    NIK = F('NIK', RANK)
    nfkb_drivers = (TNF or RANK or AKT or NIK)
    NFKB = F('NFKB', nfkb_drivers and (not Btz) and (not Dex))   # Btz: proteasome; Dex: transrepression

    Len_active = int(Len and CRBN)
    IKZF = int(not Len_active)
    IRF4 = F('IRF4', IKZF)
    MYC  = F('MYC', (ERK or BCAT or STAT3 or NFKB or MTOR) and IRF4)

    CCND2 = F('CCND2', MYC or ERK or BCAT or NFKB)
    MDM2  = int(AKT or ERK)
    P53   = F('P53', not MDM2)
    P21   = int(P53 or Len_active)
    CDK46 = F('CDK46', CCND2 and (not P21))
    E2F   = F('E2F', CDK46)

    # ---- apoptosis: BH3 logic (Btz vs Dex), intrinsic anti-CD38 ----
    MCL1  = F('MCL1',  STAT3 or NFKB)
    BCLXL = F('BCLXL', NFKB or STAT3)
    NOXA  = int(Btz)                              # neutralizes MCL1 only
    BIM   = F('BIM', Dex and (not (STAT3 or AKT)))  # neutralizes MCL1+BCLXL; blunted by IL-6/AKT
    aCD38_apop = int(antiCD38)                    # direct (Fc-independent) apoptosis only
    free_MCL1  = int(MCL1 and not (NOXA or BIM))
    free_BCLXL = int(BCLXL and not BIM)
    guarded = int(free_MCL1 or free_BCLXL)
    pro     = int(NOXA or BIM or aCD38_apop)
    APOP    = int(pro and not guarded)

    ANTIAPOP = int(MCL1 or BCLXL)
    VEGF     = int(STAT3 or HYPOXIA)
    SURVIVAL = int(ANTIAPOP and not APOP)         # apoptosis reflected here

    return [int(MYC), int(CCND2), int(E2F), int(NFKB),
            int(ANTIAPOP), int(VEGF), int(CDK46), int(SURVIVAL)]


def mm_measure(fault_set, drugs):
    out1 = mm_outputs(fault_set, drugs)
    out2 = [0] * 8
    a = b = c = d = 0
    for i in range(8):
        if   out1[i] == 1 and out2[i] == 1: a += 1
        elif out1[i] == 0 and out2[i] == 1: b += 1
        elif out1[i] == 1 and out2[i] == 0: c += 1
        else:                               d += 1
    return ((b + c) ** 2) / ((a + b + c + d) ** 2)
