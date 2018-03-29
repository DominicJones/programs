To build, the following are required:
-------------------------------------
GNU Make (make)
GNU Fortran (gfortran)
GNU C (gcc)


Mesh generation and field plotting:
-----------------------------------
Gmsh v2.5 (gmsh)

Four types of boundary conditions in the mesh file are identified:
physical type = 11 - 19  inlet condition
physical type = 21 - 29  outlet condition
physical type = 31 - 39  wall condition
physical type = 41 - 49  symmetry condition

Movable walls must be indexed as 31
Slip walls must be indexed as 39


To build standard program:
--------------------------
$ cd bin/gpde
$ make clean
$ make [OPTIONS]

OPTIONS:
  BUILD=(<blank>|DEBUG|PROFILE|RELEASE)


demo run:
---------
*** NOTE: Boundary conditions are hard-coded ***

$ cd bin/gpde
$ ./gpde ../cas/duct/duct_m.msh ../mshf/gmsh.fmt ./

$ cd ../cas/duct
$ gmsh duct_m.msh pres.msh vel.msh

