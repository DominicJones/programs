## pcns
A pressure-correction navier-stokes solver
* SIMPLE pressure-velocity coupling scheme
* unstructured meshes
* 2D or 3D geometries


## Requirements
### build tools
* g++ [sudo apt-get install g++]
* gfortran [sudo apt-get install gfortran]
* gmake

### libraries
* LAPACK [sudo apt-get install libblas-dev liblapack-dev]
* openmpi (parallel only) [sudo apt-get install libblas-dev libopenmpi-dev]

### user tools
* [gmsh](http://geuz.org/gmsh/) [sudo apt-get install gmsh]


## Build
### serial
```
  $ make
```

### parallel
```
  $ make mpi=yes
```


## Run
### serial
```
  $ cd bin/
  $ ./flowlab.exe ../examples/cavity_2d/cavity_2d.json
  $ cd ../examples/cavity_2d/
  $ gmsh cavity_2d.geo cavity_2d.msh vel.msh pres.msh
```

### parallel
```
  $ cd bin/
  $ mpirun -n 4 ./flowlab.exe ../examples/cavity_2d/cavity_2d.json part
  $ mpirun -n 4 ./flowlab.exe ../examples/cavity_2d/cavity_2d.json
  $ cd ../examples/cavity_2d/
  $ gmsh cavity_2d.geo cavity_2d.msh_p* vel.msh_p* pres.msh_p*
```
