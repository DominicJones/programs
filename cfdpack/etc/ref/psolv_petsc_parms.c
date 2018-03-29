#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "mshf.h"
#include "mesh.h"
#include "mpi.h"
#include "map.h"
#include "graph.h"

#include "parms.h"
#include "petscksp.h"

void write_gmsh_field (struct mesh_t* mesh, char* filename, int n_phi, double* phi);


int main (int argc, char* argv[]) {
  char fieldfile[256], meshfile[256], mshffile[256], buf[256];
  int i, j, k, l, m, n;
  int ii0, iin, ii, jj0, jjn, jj, kk0, kkn, kk, ll0, lln, ll;
  int n_dim, n_elm, n_fce, n_vrt, n_phi, ty, idx, dim;
  int igbl, iloc, idmn, ibnd, ighst, ipart, ios;
  int *fce_vrt_ptr, *fce_vrt_lst, *fce_elms;

  int mpi_err, mpi_size, mpi_rank, mpi_tag = 100;
  MPI_Comm mpi_comm_world;
  MPI_Status* mpi_stat;
  MPI_Request* mpi_rqst;

  struct mesh_t mesh;
  struct ghst_t ghst, locl;
  struct mshf_t mshf;

  int key, value;
  struct Map *pair;
  struct Map *ismap = NULL;

  double *phi, *phi_locl, *rhs;

  parms_Map PAmap;
  parms_Mat PAmat;
  parms_PC  PApc;
  parms_Solver PAsol;
  parms_Timer  PAtmr;

  Vec PSsol, PSrhs;
  Mat PSmat;
  KSP PSksp;
  PC  PSpc;
  PetscErrorCode pts_err;
  PetscInt PSm;


//   mpi_err = MPI_Init (&argc, &argv);
  pts_err = PetscInitialize (&argc, &argv, (char *)0, (char *)0);
  mpi_comm_world = PETSC_COMM_WORLD;

  mpi_err = MPI_Comm_size (PETSC_COMM_WORLD, &mpi_size);
  mpi_err = MPI_Comm_rank (PETSC_COMM_WORLD, &mpi_rank);

  mpi_stat = malloc ((mpi_size*2) * sizeof(MPI_Status));
  mpi_rqst = malloc ((mpi_size*2) * sizeof(MPI_Request));


  if (argc != 4) {
    fprintf (stderr, "Incorrect arguments\n"
	"Enter <mesh file> <mesh format file> <max. num. dimensions>\n");
    exit (1);
  }

  strcpy (meshfile, argv[1]);
  strcat (meshfile, "_p");
  sprintf (buf, "%d", mpi_rank);
  strcat (meshfile, buf);

  printf ("Reading %s\n", meshfile);
  read_imf_mesh (&mesh, 0, meshfile, 's', 'V', mpi_size, mpi_rank);
  read_imf_mesh (&mesh, 0, meshfile, 's', 'P', mpi_size, mpi_rank);
  read_imf_mesh (&mesh, 0, meshfile, 's', 'E', mpi_size, mpi_rank);
  read_imf_mesh (&mesh, &ghst, meshfile, 's', 'G', mpi_size, mpi_rank);

  printf ("[%d] No. vertices %d\n",mpi_rank,mesh.n_vrt);
  printf ("[%d] No. elements %d\n",mpi_rank,mesh.n_elm);
  printf ("[%d] No. ghosts %d\n",mpi_rank,ghst.num);

  mesh.dmn_part_ptr = calloc ((mpi_size+1), sizeof(int));

//   for (i=0; i<mpi_size; ++i) {
//     printf ("Part %d  Element range [%d %d)\n", mpi_rank, mesh.elm_part_ptr[i], mesh.elm_part_ptr[i+1]);
//   }

  ghst.part_idx = malloc (ghst.num * sizeof(int));
  for (i=0; i<ghst.num; ++i) {
    igbl = ghst.elm_idx[i];
    for (j=0; j<mpi_size; ++j) {
      if (j == mpi_rank) continue;
      if (igbl >= mesh.elm_part_ptr[j] && igbl < mesh.elm_part_ptr[j+1]) {
	ghst.part_idx[i] = j;
	break;
      }
    }
    printf ("[%d] Ghost %d is at [idx %d, part %d]\n",mpi_rank,i,ghst.elm_idx[i],ghst.part_idx[i]);
  }


  strcat (meshfile, "_x");
  printf ("Reading %s\n", meshfile);
  read_imf_mesh (&mesh, 0, meshfile, 's', 'F', mpi_size, mpi_rank);


  strcpy (mshffile, argv[2]);

  printf ("Reading %s\n", mshffile);
  read_mesh_format (&mshf, mshffile);


  n_dim = atoi (argv[3]);


  mesh.n_bfce = 0;
  for (i=0; i<mesh.n_elm; ++i) {
    ty = mesh.elm_types[i];
    idx = mshf.elm_idx[ty];
    dim = mshf.elm_dim[idx];
//     printf ("%d %d %d %d\n",i,ty,idx,dim);
    if (dim == n_dim - 1) ++mesh.n_bfce;
  }

  printf ("[%d] No. boundary elements %d\n",mpi_rank,mesh.n_bfce);


  // create mapping from global element index to local element index
  iloc = 0;
  ibnd = mesh.n_elm - mesh.n_bfce;
  ighst = mesh.n_elm;

  for (i=0; i<mesh.n_elm; ++i) {
    ty = mesh.elm_types[i];
    idx = mshf.elm_idx[ty];
    dim = mshf.elm_dim[idx];
    igbl = i + mesh.elm_part_ptr[mpi_rank];
    key = igbl;

    switch (n_dim - dim) {
    case 0:
      value = iloc;
      map_add (&ismap, key, value);
//       printf ("[%d] ismap elm: %d -> %d\n",mpi_rank,igbl,iloc);
      ++iloc;
      break;
    case 1:
      value = ibnd;
      map_add (&ismap, key, value);
//       printf ("[%d] ismap bnd: %d -> %d\n",mpi_rank,igbl,ibnd);
      ++ibnd;
      break;
    }
  }

  for (i=0; i<ghst.num; ++i) {
    igbl = ghst.elm_idx[i];
    key = igbl;
    value = ighst;
    map_add (&ismap, key, value);
//     printf ("[%d] ismap ghst: %d -> %d\n",mpi_rank,igbl,ighst);
    ++ighst;
  }


  // allocate workspace vectors
  n = mesh.n_elm + ghst.num;
  phi = calloc (n, sizeof(FLOAT));

  n = mesh.n_elm - mesh.n_bfce;
  rhs = calloc (n, sizeof(FLOAT));


  // initialise local and boundary values
  for (i=0; i<mesh.n_elm; ++i) {
    phi[i] = mpi_rank + 0.1;
  }



  // INITIALISE GHOST ELEMENT VALUES




  // how many element values are required from each partition?
  ghst.n_elm = calloc (mpi_size, sizeof(int));

  for (i=0; i<ghst.num; ++i) {
    igbl = ghst.elm_idx[i];
    for (j=0; j<mpi_size; ++j) {
      if (j == mpi_rank) continue;
      if (igbl >= mesh.elm_part_ptr[j] && igbl < mesh.elm_part_ptr[j+1]) {
	++ghst.n_elm[j];
      }
    }
  }

  locl.n_elm = calloc (mpi_size, sizeof(int));

  mpi_err = MPI_Barrier (PETSC_COMM_WORLD);

  for (i=0; i<mpi_size; ++i) {
    if (i == mpi_rank) continue;
    mpi_err = MPI_Irecv(&locl.n_elm[i], 1, MPI_INT, i, MPI_ANY_TAG,
		       PETSC_COMM_WORLD, &mpi_rqst[i]);
  }

  for (i=0; i<mpi_size; ++i) {
    if (i == mpi_rank) continue;
    mpi_err = MPI_Send(&ghst.n_elm[i], 1, MPI_INT, i, mpi_tag,
			PETSC_COMM_WORLD);
  }

  mpi_err = MPI_Barrier (PETSC_COMM_WORLD);




  // which elements are they (what are their indices)?
  locl.num = 0;
  for (i=0; i<mpi_size; ++i) {
    locl.num += locl.n_elm[i];
  }

  locl.elm_idx = calloc (locl.num, sizeof(int));
  phi_locl = calloc (locl.num, sizeof(double));

  mpi_err = MPI_Barrier (PETSC_COMM_WORLD);

  ii = 0;
  for (i=0; i<mpi_size; ++i) {
    if (i == mpi_rank || locl.n_elm[i] == 0) continue;
    mpi_err = MPI_Irecv(&locl.elm_idx[ii], locl.n_elm[i], MPI_INT, i, MPI_ANY_TAG,
		       PETSC_COMM_WORLD, &mpi_rqst[i]);
    ii += locl.n_elm[i];
  }

  ii = 0;
  for (i=0; i<mpi_size; ++i) {
    if (i == mpi_rank || ghst.n_elm[i] == 0) continue;
    mpi_err = MPI_Send(&ghst.elm_idx[ii], ghst.n_elm[i], MPI_INT, i, mpi_tag,
		       PETSC_COMM_WORLD);
    ii += ghst.n_elm[i];
  }

  mpi_err = MPI_Barrier (PETSC_COMM_WORLD);




  // send the requested elmement values
  for (i=0; i<locl.num; ++i) {
    igbl = locl.elm_idx[i];
    pair = map_find (&ismap, igbl);
    iloc = pair->value;
    phi_locl[i] = phi[iloc];
  }

  mpi_err = MPI_Barrier (PETSC_COMM_WORLD);

  ii = mesh.n_elm;
  for (i=0; i<mpi_size; ++i) {
    if (i == mpi_rank || ghst.n_elm[i] == 0) continue;
    mpi_err = MPI_Irecv(&phi[ii], ghst.n_elm[i], MPI_DOUBLE, i, MPI_ANY_TAG,
		       PETSC_COMM_WORLD, &mpi_rqst[mpi_size+i]);
    ii += ghst.n_elm[i];
  }

  ii = 0;
  for (i=0; i<mpi_size; ++i) {
    if (i == mpi_rank || locl.n_elm[i] == 0) continue;
    mpi_err = MPI_Send(&phi_locl[ii], locl.n_elm[i], MPI_DOUBLE, i, mpi_tag,
		       PETSC_COMM_WORLD);
    ii += locl.n_elm[i];
  }

  mpi_err = MPI_Barrier (PETSC_COMM_WORLD);


  // test communications
  for (i=0; i<ghst.num; ++i) {
    igbl = ghst.elm_idx[i];

    jj = -1;
    for (j=0; j<mpi_size; ++j) {
      if (j == mpi_rank || ghst.n_elm[i] == 0) continue;
      if (igbl >= mesh.elm_part_ptr[j] && igbl < mesh.elm_part_ptr[j+1]) {
	jj = j; break;
      }
    }

    ii = -1;
    pair = map_find (&ismap, igbl);
    if (pair) ii = pair->value;

//     printf ("[%d] ghost comm: %d %d %d %d %g\n",mpi_rank,i,igbl,jj,ii,phi[ii]);
  }



//   for (i=0; i<mesh.n_fce; ++i) {
//     jj0 = 2*i;
//     jjn = 2*(i+1);
//     printf ("[%d] tup %d, val[",mpi_rank,i);
//     for (jj=jj0; jj<jjn; ++jj) {
//       j = mesh.fce_elms[jj];
//       printf (" %d",j);
//     }
//     printf ("]\n");
//   }

  int m_grph, *grph_ia, *grph_ja;
  tuple_to_graph (2, mesh.n_fce, mesh.fce_elms, &m_grph, &grph_ia, &grph_ja);


  int *ord;
  ord = malloc (m_grph * sizeof(int));

  // internal faces and ghost faces *then* boundary faces
  idmn = 0; ibnd = m_grph - 1;
  for (i=0; i<m_grph; ++i) {
    if (grph_ja[2*i+1] > 0) {
      ord[i] = idmn; ++idmn;
    }
    else {
      ord[i] = ibnd; --ibnd;
    }
  }
  n = idmn;
//   printf ("[%d] DMN %d, BND %d\n",mpi_rank, idmn, (m_grph-1)-ibnd);
  reorder_graph (m_grph, grph_ia, grph_ja, ord);


  // internal faces *then* ghost faces
  idmn = 0; ibnd = n - 1;
  for (i=0; i<m_grph; ++i) {
    if (grph_ja[2*i] > 0) {
      ord[i] = idmn; ++idmn;
    }
    else {
      ord[i] = ibnd; --ibnd;
    }
  }

  reorder_graph (n, grph_ia, grph_ja, ord);


  // subtract the '+1' offset from local domain indices
  for (i=0; i<grph_ia[m_grph]; ++i) {
    if (grph_ja[i] > 0) --grph_ja[i];
  }


//   for (i=0; i<n; ++i) {
//     jj0 = grph_ia[i];
//     jjn = grph_ia[i+1];
//     printf ("[%d] G row %d, val[",mpi_rank,i);
//     for (jj=jj0; jj<jjn; ++jj) {
//       j = grph_ja[jj];
//       printf (" %d",j);
//     }
//     printf ("]\n");
//   }


  int m_dual, *dual_ia, *dual_ja;
  graph_to_dual (n, grph_ia, grph_ja, &m_dual, &dual_ia, &dual_ja);


//   for (i=0; i<m_dual; ++i) {
//     jj0 = dual_ia[i];
//     jjn = dual_ia[i+1];
//     printf ("[%d] D row %d, val[",mpi_rank,i);
//     for (jj=jj0; jj<jjn; ++jj) {
//       j = dual_ja[jj];
//       printf (" %d",j);
//     }
//     printf ("]\n");
//   }


  int m_mat, *mat_ia, *mat_ja;
  graph_to_matrix (n, grph_ia, grph_ja, &m_mat, &mat_ia, &mat_ja);


//   for (i=0; i<m_mat; ++i) {
//     jj0 = mat_ia[i];
//     jjn = mat_ia[i+1];
//     printf ("[%d] M row %d, val[",mpi_rank,i);
//     for (jj=jj0; jj<jjn; ++jj) {
//       j = mat_ja[jj];
//       printf (" %d",j);
//     }
//     printf ("]\n");
//   }


  
  // construct global domain element index offset list
  mesh.dmn_part_ptr[mpi_rank+1] = mesh.n_elm - mesh.n_bfce;
  mpi_err = MPI_Gather (&mesh.dmn_part_ptr[mpi_rank+1], 1, MPI_INT,
		        &mesh.dmn_part_ptr[mpi_rank+1], 1, MPI_INT,
		        0, PETSC_COMM_WORLD);
  
  if (mpi_rank == 0) {
    for (i=0; i<mpi_size; ++i) {
      mesh.dmn_part_ptr[i+1] += mesh.dmn_part_ptr[i];
    }
  }

  mpi_err = MPI_Bcast (mesh.dmn_part_ptr, mpi_size+1, MPI_INT,
		       0, PETSC_COMM_WORLD);


  for (i=0; i<mpi_size; ++i) {
    printf ("[%d] DMN: [%d %d)\n",mpi_rank,mesh.dmn_part_ptr[i],mesh.dmn_part_ptr[i+1]);
  }


  // convert element matrix indices to global indices
  for (i=0; i<m_mat; ++i) {
    jj0 = mat_ia[i];
    jjn = mat_ia[i+1];
    for (jj=jj0; jj<jjn; ++jj) {
      j = mat_ja[jj];
      if (j < 0) {
	igbl = -j -1;
	pair = map_find (&ismap, igbl);
	if (pair) {
	  iloc = pair->value - mesh.n_elm;
	  ipart = ghst.part_idx[iloc];
	  ios = mesh.elm_part_ptr[ipart] - mesh.dmn_part_ptr[ipart];
	} else {
	  printf ("[E] Cannot find pair from key %d\n", igbl);
	}
// 	printf ("[%d] Ghost: G %d  L %d  P %d  OS %d\n",mpi_rank,igbl,iloc,ipart,ios);
	j = igbl - ios;
      }
      else {
	iloc = j;
	igbl = iloc + mesh.dmn_part_ptr[mpi_rank];
	j = igbl;
      }
      mat_ja[jj] = j;
    }
  }



  int *imat;
  imat = malloc (m_mat * sizeof(int));

  for (i=0; i<m_mat; ++i) {
    jj0 = mat_ia[i];
    jjn = mat_ia[i+1];
    ios = mesh.dmn_part_ptr[mpi_rank];
    imat[i] = i + ios;
    printf ("[%d] N row %d, val[",mpi_rank,i+ios);
    for (jj=jj0; jj<jjn; ++jj) {
      j = mat_ja[jj];
      printf (" %d",j);
    }
    printf ("]\n");
  }


  // initialise matrix values
  double *mat_a, val;
  mat_a = calloc (mat_ia[m_mat], sizeof(double));

  for (ii=0; ii<m_mat; ++ii) {
    i = imat[ii];
    jj0 = mat_ia[ii];
    jjn = mat_ia[ii+1];
    for (jj=jj0; jj<jjn; ++jj) {
      j = mat_ja[jj];
      val = -1.; if (i == j) val = 4.;
      mat_a[jj] = val;
    }
  }



  //------------------------------------------------------------//
  // PETSc METHOD
  //------------------------------------------------------------//


//   /* create a matrix object */
//   MatCreateMPIAIJ (PETSC_COMM_WORLD, m_mat, m_mat,
// 		   PETSC_DETERMINE, PETSC_DETERMINE,
// 		  7, PETSC_NULL, 7, PETSC_NULL, &PSmat);
// 
// 
//   /* insert values into matrix object */
//   ios = mesh.dmn_part_ptr[mpi_rank];
//   for (ii=0; ii<m_mat; ++ii) {
//     i = ii + ios;
//     jj0 = mat_ia[ii];
//     jjn = mat_ia[ii+1];
//     for (jj=jj0; jj<jjn; ++jj) {
//       j = mat_ja[jj];
//       val = -1.; if (i == j) val = 4.;
// //       printf ("([%d] %d %d %g) ",mpi_rank,i,j,val);
//       MatSetValues (PSmat, 1, &i, 1, &j, &val, INSERT_VALUES);
//     }
//   }
// 
// 
//   /* assemble the matrix */
//   MatAssemblyBegin (PSmat, MAT_FINAL_ASSEMBLY);
//   MatAssemblyEnd (PSmat, MAT_FINAL_ASSEMBLY);
// 
// 
//   /* create an abstract krylov subspace object */
//   KSPCreate (PETSC_COMM_WORLD, &PSksp);
// 
// 
//   /* create solution and right-hand-side vector */
//   VecCreateMPI (PETSC_COMM_WORLD, m_mat, PETSC_DETERMINE, &PSrhs);
//   VecDuplicate (PSrhs, &PSsol);
// 
//   val = 1.;
//   VecSet (PSrhs, val);
// 
//   VecAssemblyBegin (PSrhs);
//   VecAssemblyEnd (PSrhs);
// 
//   val = 0.;
//   VecSet (PSsol, val);
// 
//   VecAssemblyBegin (PSsol);
//   VecAssemblyEnd (PSsol);
// 
// 
//   /* solve linear systems with the krylov subspace method selected */
//   KSPSetOperators (PSksp, PSmat, PSmat, SAME_NONZERO_PATTERN);
//   KSPSetType (PSksp, KSPFGMRES); 
//   KSPGMRESSetRestart (PSksp, 20);
//   KSPSetTolerances (PSksp, 1.0e-6, 1.0e-20, 1.0e6, 1000);
// 
// 
//   /* Extract PC context from the KSP context */
//   pts_err = KSPGetPC (PSksp, &PSpc); CHKERRQ (pts_err);
//   pts_err = PCSetType (PSpc, PCBJACOBI); CHKERRQ (pts_err);
//   KSPSetFromOptions (PSksp);
// 
//   KSPSolve (PSksp, PSrhs, PSsol);
// //   VecView (PSsol, PETSC_VIEWER_STDOUT_WORLD);
// 
// 
//   /* check */
//   val = 0.;
//   VecSet (PSrhs, val);
//   MatMult (PSmat, PSsol, PSrhs);
//   VecView (PSrhs, PETSC_VIEWER_STDOUT_WORLD);
// 
// 
//   VecGetValues(PSsol, m_mat, imat, phi);


  //------------------------------------------------------------//
  // pARMS METHOD
  //------------------------------------------------------------//


  // Set ownership range for each process and ownership of the elements
  int *vtxdist, *part;

  vtxdist = malloc ((mpi_size+1) * sizeof(int));
  part = malloc (m_mat * sizeof(int));  

  MPI_Allgather (&m_mat, 1, MPI_INT, &vtxdist[1], 1, MPI_INT, PETSC_COMM_WORLD);

  vtxdist[0] = 0;
  for (i=0; i<mpi_size; ++i) {
    vtxdist[i+1] += vtxdist[i];
  }

  for (i=0; i<m_mat; ++i) part[i] = mpi_rank;


  // Now create pARMS map object
  parms_MapCreateFromDist (&PAmap, vtxdist, part, PETSC_COMM_WORLD, 0, 1, NONINTERLACED);  


  // Free arrays no longer needed
  free (part);
  free (vtxdist);


  // Create a distributed matrix based on the parms_Map created above.
  parms_MatCreate (&PAmat, PAmap);


  // Insert entries into the matrix
  parms_MatSetValues (PAmat, m_mat, imat, mat_ia, mat_ja, mat_a, INSERT);


  // After calling the following statement, no entries can be inserted
  // into the matrix
  parms_MatSetup (PAmat);


  // Create a preconditioner object based on the matrix
  parms_PCCreate (&PApc, PAmat);


  // Setup the preconditioner type (eg. RAS)
  parms_PCSetType (PApc, PCRAS);


  // Setup the type of local ILU preconditioner (eg. ARMS)
  parms_PCSetILUType (PApc, PCARMS);


  // Setup the preconditioning matrix
  parms_PCSetup (PApc);


  // Create a solver object based on the matrix and the preconditioner
  parms_SolverCreate (&PAsol, PAmat, PApc);


  // set vector values
  for (i=0; i<m_mat; ++i) rhs[i] = 1.;
  for (i=0; i<m_mat; ++i) phi[i] = 0.;


  // Solve the linear system of equations
  parms_SolverApply (PAsol, rhs, phi);


  //------------------------------------------------------------//
  // END OF METHOD
  //------------------------------------------------------------//


  // write solution field to file
  char *pch;
  strcpy (fieldfile, argv[1]);
  pch = strrchr (fieldfile,'/');
  pch[0] = '\0';
  strcat (fieldfile, "/field.msh_p");
  sprintf (buf, "%d", mpi_rank);
  strcat (fieldfile, buf);

//   printf ("OUTPUT %s\n", fieldfile);

  write_gmsh_field (&mesh, fieldfile, mesh.n_elm-mesh.n_bfce, phi);


  pts_err = PetscFinalize (); CHKERRQ (pts_err);
  return (0);
}
