# -*- coding: utf-8 -*-
"""Verilog <-> Python equivalence check for the exported FSM (hdl/MM_attractor.v).

Draws N random (state, fault set, drug) cases, computes the next state with the
Python reference (MM_attractor._next), and writes them as test vectors. A
self-checking testbench (hdl/MM_attractor_equiv_tb.v) loads each state into the
Verilog register, applies one clock, and compares the register with the expected
next state.

Usage (from src/):
    python3 MM_verilog_equiv.py            # 20,000 cases (paper setting)
    python3 MM_verilog_equiv.py 5000       # custom N
Requires Icarus Verilog (iverilog, vvp) on PATH to run the simulation; without it,
the vectors are still written to hdl/equiv_vectors.txt.
Expected output: "checked 20000 cases, mismatches = 0"
"""
import os, random, shutil, subprocess, sys
from MM_network import FAULTS
from MM_attractor import NODES, IDX, _next, _setup

N = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
SEED = 2026
HDL = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'hdl')
VEC = os.path.join(HDL, 'equiv_vectors.txt')

def bits(vec):            # list[0/1] indexed by node -> string with bit 0 rightmost
    return ''.join(str(int(b)) for b in reversed(vec))

def main():
    rng = random.Random(SEED)
    fids = list(FAULTS)
    lines = []
    for _ in range(N):
        fault_set = rng.sample(fids, rng.choice([0, 1, 1, 2, 3]))
        drugs = [rng.randint(0, 1) for _ in range(5)]          # Isa Dara Btz Dex Len
        ov, drugs2, Len_active, bim_del = _setup(fault_set, drugs)
        state = [rng.randint(0, 1) for _ in NODES]
        sa_en = [0] * len(NODES); sa_val = [0] * len(NODES)
        for node, val in ov.items():
            if node in IDX:                                     # state-node fault
                state[IDX[node]] = val; sa_en[IDX[node]] = 1; sa_val[IDX[node]] = val
        nxt = _next(tuple(state), drugs2, ov, Len_active, bim_del)
        ctrl = list(drugs2) + [int(ov.get('PTEN', 1)), int(ov.get('CRBN', 1)), int(bim_del)]
        # one line = ctrl(8) | state(22) | sa_en(22) | sa_val(22) | expected(22) = 96 bits
        lines.append(''.join(map(str, ctrl)) + bits(state) + bits(sa_en)
                     + bits(sa_val) + bits(nxt))
    with open(VEC, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    print(f'wrote {N} vectors -> {os.path.relpath(VEC)}')

    if not (shutil.which('iverilog') and shutil.which('vvp')):
        print('iverilog/vvp not found: install Icarus Verilog to run the simulation.')
        return
    sim = os.path.join(HDL, 'mm_equiv')
    subprocess.run(['iverilog', f'-DNVEC={N}', '-o', sim,
                    os.path.join(HDL, 'MM_attractor.v'),
                    os.path.join(HDL, 'MM_attractor_equiv_tb.v')], check=True)
    out = subprocess.run(['vvp', sim], cwd=HDL, capture_output=True, text=True).stdout
    for line in out.splitlines():
        if 'mismatch' in line:
            print(line)
    os.remove(sim)

if __name__ == '__main__':
    main()
