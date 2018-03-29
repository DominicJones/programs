#include "parms_intf.h"
#ifdef USE_MPI
#include "parms.h"
#include "mpi.h"
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <climits>


#ifdef USE_MPI
parms_Map PAmap;
parms_Mat PAmat;
parms_PC  PApc;
parms_Solver PAsol;
parms_Timer  PAtmr;
double *PAres;
#endif


void parms_init_ (int* mp)
{
#ifdef USE_MPI
  int mpi_size, mpi_rank;
  int i, j, jj0, jjn, jj, m;
  int *part, *dist;
  int offset, dof;
  int n_levels, max_iter;
  double drop_tol;
  char buf[256];

  m = *mp;
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);

  offset = 1;
  dof = 1;

  // create pARMS map object
  if (mpi_size > 1) {
    dist = (int*) malloc ((mpi_size+1) * sizeof(int));
    part = (int*) malloc (m * sizeof(int));

    // Set ownership range for each process
    MPI_Allgather (
      &m, 1, MPI_INT, &dist[1], 1,
      MPI_INT, MPI_COMM_WORLD);

    dist[0] = offset;
    for (i = 0; i < mpi_size; i++) {
      dist[i+1] += dist[i];
    }

    // Set process ownership of the nodes
    for (i = 0; i < m; i++) {
      part[i] = mpi_rank + offset;
    }

    parms_MapCreateFromDist (
      &PAmap, dist, part, MPI_COMM_WORLD,
      offset, dof, NONINTERLACED);

    free (part);
    free (dist);
  }
  else {
    parms_MapCreateFromLocal (&PAmap, m, offset);
  }

  // Create a distributed matrix based on the parms_Map created above.
  parms_MatCreate (&PAmat, PAmap);

  // Create a preconditioner object based on the matrix
  parms_PCCreate (&PApc, PAmat);

//   parms_PCSetType (PApc, PCRAS);
  parms_PCSetType (PApc, PCBJ);
//   parms_PCSetType (PApc, PCSCHUR);

//   parms_PCSetILUType (PApc, PCARMS);
//   parms_PCSetILUType (PApc, PCILU0);
  parms_PCSetILUType (PApc, PCILUK);
//   parms_PCSetILUType (PApc, PCILUT);

  n_levels = 3;
  parms_PCSetNlevels (PApc, n_levels);

//   max_iter = 5; //25;
//   parms_PCSetInnerMaxits (PApc, max_iter);

//   drop_tol = 0.00001;
//   parms_PCSetTol (PApc, &drop_tol);


  // Create a solver object based on the matrix and the preconditioner
  parms_SolverCreate (&PAsol, PAmat, PApc);

  max_iter = 25;
  sprintf (buf, "%d", max_iter);
  parms_SolverSetParam (PAsol, MAXITS, buf);

  drop_tol = 0.01;
  sprintf (buf, "%g", drop_tol);
  parms_SolverSetParam (PAsol, DTOL, buf);


  // initialise residual
  PAres = (double*) malloc (m * sizeof(double));
#endif
}


void parms_solv_ (int* mp, int* im, int* ia, int* ja,
                  double* aij, double* rhs, double* sol,
                  double* res0)
{
#ifdef USE_MPI
  int mpi_size, mpi_rank;
  int i, j, jj0, jjn, jj, m;
  double *res, res1;

  m = *mp;
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);

  // Insert entries into the matrix
  parms_MatSetValues (PAmat, m, im, ia, ja, aij, INSERT);

  // After calling the following statement, no entries can be inserted
  // into the matrix
  parms_MatSetup (PAmat);

  // Setup the preconditioning matrix
  parms_PCSetup (PApc);

  parms_MatVec (PAmat, sol, PAres);
  parms_VecAXPY (PAres, rhs, -1., PAmap);
  parms_VecGetNorm2 (PAres, res0, PAmap);

  // Solve the linear system of equations
  parms_SolverApply (PAsol, rhs, sol);
  parms_SolverGetResidualNorm2 (PAsol, rhs, sol, &res1);
  parms_MatReset (PAmat, SAME_NONZERO_STRUCTURE);
#endif
}
