this fork
=========

This is Curtis Peterson's personal fork of the MILC code. It 
contains a functioning HMC code that runs with HISQ w/o rooting.
To compile it, copy `Makefile` into `ks_imp_dyn` and edit as 
you need for configuring on whatever machine that you're running
one. Then run `make` in `ks_imp_dyn` with either `su3_hmc1_hisq_pppa`
or `su3_hmc1_hisq_apbc` as your make target. The former will compile
with anti-periodic boundary conditions in the time direction, while
the latter will compile with anti-periodic boundary conditions in 
all four directions. 

milc_qcd
========

MILC collaboration code for lattice QCD calculations

The MILC Code is a body of high performance research software written
in C (with some C++) for doing SU(3) lattice gauge theory on high
performance computers as well as single-processor workstations.  A
wide variety of applications are included.  See NEWS for recent
important changes.
