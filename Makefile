# Reproduce the paper's computational results (run from repo root).
PY=python3
SRC=src

.PHONY: all equiv screen exact ci table1 baseline sens atpg figures verilog clean
all: equiv screen ci table1 baseline sens atpg

equiv:    ; cd $(SRC) && $(PY) MM_network_gates.py
screen:   ; cd $(SRC) && $(PY) MM_attractor_screen.py
exact:    ; cd $(SRC) && $(PY) MM_exact_attractors.py
ci:       ; cd $(SRC) && $(PY) MM_attractor_ci.py
table1:   ; cd $(SRC) && $(PY) MM_table1_fix.py
baseline: ; cd $(SRC) && $(PY) MM_baseline_null.py && $(PY) MM_null_model.py
sens:     ; cd $(SRC) && $(PY) MM_async_robust.py && $(PY) MM_rule_perturb.py && $(PY) MM_readout_sens.py
atpg:     ; cd $(SRC) && $(PY) MM_atpg.py
figures: screen
	cd $(SRC) && $(PY) make_fig1_minimal.py && $(PY) make_compare_compact.py

verilog:
	cd $(SRC) && $(PY) MM_verilog_equiv.py
	cd hdl && iverilog -o mm_tb MM_attractor.v MM_attractor_tb.v && vvp mm_tb

clean:
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
	rm -f hdl/mm_tb hdl/mm_equiv hdl/*.vvp hdl/equiv_vectors.txt
	rm -f src/*.csv src/*.png src/*.svg src/*.npy src/circuit.dot
