# -*- coding: utf-8 -*-
"""
MM_network_partial.py  (v4 - "intermediate" anti-CD38)
Keeps the v3 intrinsic network UNCHANGED, but represents anti-CD38 as PARTIAL,
population-level cytoreduction (ADCC), not as a per-cell binary kill:

    burden       = (oncogenic outputs ON) / 8           # intrinsic, WITHOUT anti-CD38
    if anti-CD38 present:  burden *= (1 - k)             # ADCC removes a fraction k of cells
    measure      = burden ** 2

  k = K_BASE (0.55) single-agent;  K_LEN (0.75) when lenalidomide is present
      (lenalidomide enhances NK/ADCC and up-regulates CD38 -> Dara+Len synergy).

So anti-CD38 is neither useless (v3) nor a magic bullet (v2 immune): it gives a
real but partial benefit, strongest combined with lenalidomide. k is tunable.

This module reuses the SAME 22 faults as v3 (MM_network.FAULTS), so results are
directly comparable to MM_drug_screen_results.csv on an identical fault axis.
"""
from MM_network import (FAULTS, FAULT_LABELS, N_FAULTS, OUTPUT_NAMES,
                        mm_outputs as _v3_outputs)

K_BASE = 0.55     # anti-CD38 ADCC cytoreduction fraction (single agent)
K_LEN  = 0.75     # boosted when lenalidomide is present

def mm_outputs_partial(fault_set, drugs):
    """Intrinsic 8 outputs WITHOUT anti-CD38's direct effect (anti-CD38 acts via ADCC scaling)."""
    Isa, Dara, Btz, Dex, Len = [int(bool(d)) for d in drugs]
    return _v3_outputs(fault_set, [0, 0, Btz, Dex, Len])

def adcc_factor(drugs):
    Isa, Dara, Btz, Dex, Len = [int(bool(d)) for d in drugs]
    if not (Isa or Dara):
        return 1.0
    return 1.0 - (K_LEN if Len else K_BASE)

def mm_measure_partial(fault_set, drugs):
    outs = mm_outputs_partial(fault_set, drugs)
    burden = (sum(outs) / 8.0) * adcc_factor(drugs)
    return burden ** 2
