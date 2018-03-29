spmom
=====

A "method of moments" spray simulation code


--

directories:

/src  [source codes (*.f90)]
/bin  [binary and makefile]
/lib  [LAPACK and LAPACK-95 (not used)]
/cas  [case grid and output]



source:

divides into

grd_*
cfd_*
solv_*
mom_*
spr_*

The first three are required to simulate continuum flow,
the last two are for spray modelling.

The hierarchy of the files is (for the commonly edited files):

main.f90
cfd_solve.f90
cfd_define.f90
cfd_setup.f90

Real precision is defined in

precision.f90


Grids must be formatted in the neutral format of Gambit,
where boundary conditions must be named as

inlet_N (= 10 + N)
outlet_N (= 20 + N)
wall_N (= 30 + N)
symmetry_N (= 40 + N)

where N is the user defined index of the boundary, starting at N=1.

Output files are formatted for Tecplot.


To run the code, change directory to cfd/bin
edit the makefile if necessary
Type <make> (or <make clean> then <make>)
Type <./spmom>


Dependencies on lapack are found in
mom_math.f90
mom_closure.f90
spr_hydro.f90
and are used by the Splines and Max. Entropy solvers
