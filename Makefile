.PHONY: src/fine.a cfine/fine.a cppfine/fine.a pyfine/fine.so


all: pyfine/fine.so

slatec/slatec.a:
	cd slatec; make

symb_diffs/heston_psi.f90 symb_diffs/heston_psi_grad.f90: symb_diffs/heston.py symb_diffs/sympy_utils.py
	cd symb_diffs; python3 heston.py

ffine/fine.a: slatec/slatec.a symb_diffs/heston_psi.f90 symb_diffs/heston_psi_grad.f90
	cd ffine; make

cfine/fine.a: ffine/fine.a
	cd cfine; make fine.a

cppfine/fine.a: cfine/fine.a
	cd cppfine; make fine.a

pyfine/fine.so: cppfine/fine.a
	cd pyfine; make fine.so
