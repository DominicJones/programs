#include "partition.h"
#include "utils.h"

#include "mesh.h"
#include "mshf.h"
#include "map.h"
#include "graph.h"
#include "ghst_comm.h"

#ifdef USE_MPI
#include "mpi.h"
#include "parmetis.h"
#endif

#include "json/reader.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <fstream>
#include <streambuf>


int partition (struct mshf_t * mshf, struct mesh_t * mesh, char * in_file)
{
#ifdef USE_MPI
  int err = 0;
  char buf[256], out_file[256], mshffile[256];
  int n_dim, ty, dim, idx, idmn, ibnd, ielm, n_dmn, n_bnd, n_elm;
  int igbl, iloc;
  int *ord, *cnt, *blk, *ptr, *lst, *swap, *gswap;
  int ii, ii0, iin, i, jj0, jjn, jj, j, kk0, kkn, kk, k, kkz;
  int l, m, n, llz, mmz, nnz, p, op, np;
  int idx_offset = 0, min_fce_vrt;
  float *xlst;

  int key, value;
  struct Map *pair;
  struct Map *map = NULL;

  struct ghst_t *gdmn, *ldmn;
  struct ghst_t *gbnd, *lbnd;

  int mpi_size, mpi_rank;
  int mpi_tag = 100;
  MPI_Status* mpi_stat;
  MPI_Request* mpi_rqst;
  MPI_Comm mpi_comm_world;

  MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);

  mpi_comm_world = MPI_COMM_WORLD;
  mpi_stat = (MPI_Status*) malloc (mpi_size * sizeof(MPI_Status));
  mpi_rqst = (MPI_Request*) malloc (mpi_size * sizeof(MPI_Request));

  p = mesh->part = mpi_rank;
  np = mesh->n_part = mpi_size;

  n_dim = mesh->n_dim;


  // read mesh element data in a distributed manner
  read_gmsh_mesh (mesh, in_file, 'E', 'd');

  distributed_block (mesh->N_elm, np, &mesh->elm_blk);
  offset_to_count (np, mesh->elm_blk, &mesh->elm_cnt);


  // convert gmsh element ordering to metis ordering
  int vlst[8];
  for (i=0; i<mesh->n_elm; ++i) {
    ty = mesh->elm_tag_ptr[i];
    ty = mesh->elm_tag_lst[ty];
    idx = mshf->elm_idx[ty];
    dim = mshf->elm_dim[idx];

    jj0 = mesh->elm_vrt_ptr[i];
    jjn = mesh->elm_vrt_ptr[i+1];

    k = 0;

    for (jj=jj0; jj<jjn; ++jj) {
      vlst[k] = mesh->elm_vrt_lst[jj];
      k++;
    }

    // reverse 2D element vertices
    if (dim == 2) {
      for (jj=jj0; jj<jjn; ++jj) {
	k--;
        mesh->elm_vrt_lst[jj] = vlst[k];
      }
    }
  }


  // partition the distributed mesh
  swap = (int*) malloc (mesh->n_elm * sizeof(int));
  for (i=0; i<mesh->n_elm; ++i) swap[i]=INT_MAX;
  {
    int options[3] = {0, 1, 15};
    int ncon = 1;
    int* elmwgt = NULL;
    int wgtflag = 0;
    int edgecut;
    float ubvec = 1.05;
    float* tpwgts = (float*) malloc (np * sizeof(float));

    for (i=0; i<np; ++i) tpwgts[i] = 1./np;
    min_fce_vrt = n_dim;

    ParMETIS_V3_PartMeshKway (
      mesh->elm_blk, mesh->elm_vrt_ptr, mesh->elm_vrt_lst,
      elmwgt, &wgtflag, &idx_offset, &ncon, &min_fce_vrt,
      &mpi_size, tpwgts, &ubvec, options, &edgecut,
      swap, &mpi_comm_world);

    free (tpwgts);
  }


  // gather to each process the distributed partition lists
  mesh->elm_part = (int*) malloc (mesh->N_elm * sizeof(int));
  for (i=0; i<mesh->N_elm; ++i) mesh->elm_part[i]=INT_MAX;

  MPI_Allgatherv (
    swap, mesh->n_elm, MPI_INT,
    mesh->elm_part, mesh->elm_cnt, mesh->elm_blk, MPI_INT,
    MPI_COMM_WORLD);

  free (swap);


  // reset for next read
  free (mesh->elm_blk);
  free (mesh->elm_cnt);
  free (mesh->elm_idx);
  free (mesh->elm_tag_ptr);
  free (mesh->elm_tag_lst);
  free (mesh->elm_vrt_ptr);
  free (mesh->elm_vrt_lst);


  // read mesh element data in a partitioned manner
  read_gmsh_mesh (mesh, in_file, 'E', 'p');
  free (mesh->elm_idx);

  partitioned_block (mesh->N_elm, np, mesh->elm_part, &mesh->elm_blk);
  offset_to_count (np, mesh->elm_blk, &mesh->elm_cnt);
  free (mesh->elm_part);


  // store element "dimension": 1 for bnd, 0 for dmn
  m = mesh->n_elm;
  mesh->elm_dim = (int*) malloc (m * sizeof(int));

  for (i=0; i<m; ++i) {
    ty = mesh->elm_tag_ptr[i];
    ty = mesh->elm_tag_lst[ty];
    idx = mshf->elm_idx[ty];
    dim = mshf->elm_dim[idx];
    mesh->elm_dim[i] = n_dim - dim;
  }


  // split elements into boundary and domain
  split_graph (
    mesh->n_elm, mesh->elm_vrt_ptr, mesh->elm_vrt_lst,
    mesh->elm_dim,
    &mesh->n_bnd, &mesh->bnd_vrt_ptr, &mesh->bnd_vrt_lst,
    &mesh->n_dmn, &mesh->dmn_vrt_ptr, &mesh->dmn_vrt_lst);

  split_graph (
    mesh->n_elm, mesh->elm_tag_ptr, mesh->elm_tag_lst,
    mesh->elm_dim,
    &mesh->n_bnd, &mesh->bnd_tag_ptr, &mesh->bnd_tag_lst,
    &mesh->n_dmn, &mesh->dmn_tag_ptr, &mesh->dmn_tag_lst);


  free (mesh->elm_vrt_ptr);
  free (mesh->elm_vrt_lst);

  free (mesh->elm_tag_ptr);
  free (mesh->elm_tag_lst);


  mesh->bnd_blk = (int*) malloc ((np+1) * sizeof(int));
  MPI_Allgather (&mesh->n_bnd, 1, MPI_INT, &mesh->bnd_blk[1], 1, MPI_INT,
    MPI_COMM_WORLD);

  mesh->bnd_blk[0] = 0;
  for (i=0; i<np; ++i) {
    mesh->bnd_blk[i+1] += mesh->bnd_blk[i];
  }
  offset_to_count (np, mesh->bnd_blk, &mesh->bnd_cnt);


  mesh->dmn_blk = (int*) malloc ((np+1) * sizeof(int));
  MPI_Allgather (&mesh->n_dmn, 1, MPI_INT, &mesh->dmn_blk[1], 1, MPI_INT,
    MPI_COMM_WORLD);

  mesh->dmn_blk[0] = 0;
  for (i=0; i<np; ++i) {
    mesh->dmn_blk[i+1] += mesh->dmn_blk[i];
  }
  offset_to_count (np, mesh->dmn_blk, &mesh->dmn_cnt);


  mesh->bnd_idx = (int*) malloc (mesh->n_bnd * sizeof(int));
  for (i=0; i<mesh->n_bnd; ++i) {
    mesh->bnd_idx[i] = i + mesh->bnd_blk[p];
  }

  mesh->dmn_idx = (int*) malloc (mesh->n_dmn * sizeof(int));
  for (i=0; i<mesh->n_dmn; ++i) {
    mesh->dmn_idx[i] = i + mesh->dmn_blk[p];
  }


  // boundary ghost elements
  min_fce_vrt = n_dim - 1;

  ParMETIS_V3_Mesh2Dual (
    mesh->bnd_blk, mesh->bnd_vrt_ptr, mesh->bnd_vrt_lst,
    &idx_offset, &min_fce_vrt,
    &mesh->bnd_neig_ptr, &mesh->bnd_neig_lst,
    &mpi_comm_world);

  gbnd = (struct ghst_t*) malloc (sizeof(struct ghst_t));
  lbnd = (struct ghst_t*) malloc (sizeof(struct ghst_t));

  find_ghost_elements (gbnd, mesh->n_bnd, mesh->bnd_blk,
    mesh->bnd_neig_ptr, mesh->bnd_neig_lst);

  communicate_ghost_elements (
    gbnd, lbnd, mesh->bnd_blk,
    mesh->bnd_tag_ptr, mesh->bnd_tag_lst,
    mesh->bnd_vrt_ptr, mesh->bnd_vrt_lst);


  // domain ghost elements
  min_fce_vrt = n_dim;

  ParMETIS_V3_Mesh2Dual (
    mesh->dmn_blk, mesh->dmn_vrt_ptr, mesh->dmn_vrt_lst,
    &idx_offset, &min_fce_vrt,
    &mesh->dmn_neig_ptr, &mesh->dmn_neig_lst,
    &mpi_comm_world);

  gdmn = (struct ghst_t*) malloc (sizeof(struct ghst_t));
  ldmn = (struct ghst_t*) malloc (sizeof(struct ghst_t));

  find_ghost_elements (gdmn, mesh->n_dmn, mesh->dmn_blk,
    mesh->dmn_neig_ptr, mesh->dmn_neig_lst);

  communicate_ghost_elements (
    gdmn, ldmn, mesh->dmn_blk,
    mesh->dmn_tag_ptr, mesh->dmn_tag_lst,
    mesh->dmn_vrt_ptr, mesh->dmn_vrt_lst);


  // merge boundary and domain into elements
  int m_bnd, *bnd_ptr, *bnd_lst;
  int m_dmn, *dmn_ptr, *dmn_lst;

  gbnd->n_elm = gbnd->elm_blk[np];
  gdmn->n_elm = gdmn->elm_blk[np];


  merge_lists (
    mesh->n_bnd, mesh->bnd_idx,
    gbnd->n_elm, gbnd->elm_idx,
    &m_bnd, &bnd_lst);

  // add to dmn_idx the total num of bnds
  m = mesh->bnd_blk[np];
  for (i=0; i<mesh->n_dmn; ++i) mesh->dmn_idx[i] += m;
  for (i=0; i<gdmn->n_elm; ++i) gdmn->elm_idx[i] += m;

  merge_lists (
    mesh->n_dmn, mesh->dmn_idx,
    gdmn->n_elm, gdmn->elm_idx,
    &m_dmn, &dmn_lst);

  merge_lists (
    m_bnd, bnd_lst,
    m_dmn, dmn_lst,
    &mesh->n_elm, &mesh->elm_idx);

  free (bnd_lst);
  free (dmn_lst);


  merge_graphs (
    mesh->n_bnd, mesh->bnd_vrt_ptr, mesh->bnd_vrt_lst,
    gbnd->n_elm, gbnd->elm_vrt_ptr, gbnd->elm_vrt_lst,
    &m_bnd, &bnd_ptr, &bnd_lst);

  merge_graphs (
    mesh->n_dmn, mesh->dmn_vrt_ptr, mesh->dmn_vrt_lst,
    gdmn->n_elm, gdmn->elm_vrt_ptr, gdmn->elm_vrt_lst,
    &m_dmn, &dmn_ptr, &dmn_lst);

  merge_graphs (
    m_bnd, bnd_ptr, bnd_lst,
    m_dmn, dmn_ptr, dmn_lst,
    &mesh->n_elm, &mesh->elm_vrt_ptr, &mesh->elm_vrt_lst);

  free (bnd_ptr); free (bnd_lst);
  free (dmn_ptr); free (dmn_lst);


  merge_graphs (
    mesh->n_bnd, mesh->bnd_tag_ptr, mesh->bnd_tag_lst,
    gbnd->n_elm, gbnd->elm_tag_ptr, gbnd->elm_tag_lst,
    &m_bnd, &bnd_ptr, &bnd_lst);

  merge_graphs (
    mesh->n_dmn, mesh->dmn_tag_ptr, mesh->dmn_tag_lst,
    gdmn->n_elm, gdmn->elm_tag_ptr, gdmn->elm_tag_lst,
    &m_dmn, &dmn_ptr, &dmn_lst);

  merge_graphs (
    m_bnd, bnd_ptr, bnd_lst,
    m_dmn, dmn_ptr, dmn_lst,
    &mesh->n_elm, &mesh->elm_tag_ptr, &mesh->elm_tag_lst);

  free (bnd_ptr); free (bnd_lst);
  free (dmn_ptr); free (dmn_lst);


  // construct element partition index
  mesh->elm_part = (int*) malloc (mesh->n_elm * sizeof(int));

  blk = (int*) malloc (5 * sizeof(int));
  blk[0] = 0;
  blk[1] = mesh->n_bnd;
  blk[2] = gbnd->n_elm;
  blk[3] = mesh->n_dmn;
  blk[4] = gdmn->n_elm;
  for (i=1; i<5; ++i) blk[i] += blk[i-1];

  for (i=blk[0]; i<blk[4]; ++i) {
    j = mesh->elm_idx[i];
    if (i < blk[2]) {
      op = offset_lb (np, mesh->bnd_blk, j);
      if (op < 0) {err= 1; goto quit;}
    }
    else {
      j -= mesh->bnd_blk[np];
      op = offset_lb (np, mesh->dmn_blk, j);
      if (op < 0) {err= 1; goto quit;}
    }
    ++op;
    if (op != (p + 1)) op *= -1;
    mesh->elm_part[i] = op;
  }

  free (blk);


  // map element vertex indices
  m = mesh->n_elm;
  nnz = mesh->elm_vrt_ptr[m];
  swap = (int*) malloc (nnz * sizeof(int));
  memcpy (swap, mesh->elm_vrt_lst, nnz * sizeof(int));

  qsort (swap, nnz, sizeof(int), icmp);
  setcp (swap, nnz, INT_MAX);
  qsort (swap, nnz, sizeof(int), icmp);

  blk = (int*) malloc ((np+1) * sizeof(int));
  blk[0] = 0;
  for (i=0; i<nnz; ++i) {
    if (swap[i] < INT_MAX) ++blk[0];
  }
  nnz = blk[0];

  MPI_Allgather (&nnz, 1, MPI_INT, &blk[1], 1, MPI_INT,
    MPI_COMM_WORLD);

  blk[0] = 0;
  for (i=0; i<np; ++i) {
    blk[i+1] += blk[i];
  }

  offset_to_count (np, blk, &cnt);

  gswap = (int*) malloc (blk[np] * sizeof(int));

  MPI_Allgatherv (
    swap, nnz, MPI_INT,
    gswap, cnt, blk, MPI_INT,
    MPI_COMM_WORLD);

  qsort (gswap, blk[np], sizeof(int), icmp);
  setcp (gswap, blk[np], INT_MAX);
  qsort (gswap, blk[np], sizeof(int), icmp);

  j = 0;
  for (i=0; i<blk[np]; ++i) {
    if (gswap[i] < INT_MAX) {
      key = gswap[i];
      value = j++;
      map_add (&map, key, value);
    }
  }
  mesh->N_vrt = j;

  free (blk);
  free (gswap);

  mesh->n_vrt = 0;
  n = mesh->N_vrt;
  mesh->vrt_part = (int*) calloc (n, sizeof(int));

  for (i=0; i<cnt[p]; ++i) {
    key = swap[i];
    pair = map_find (&map, key);
    if (pair) {
      j = pair->value;
      mesh->vrt_part[j] = 1;
      ++mesh->n_vrt;
    }
  }

  free (swap);


  // read required mesh vertices
  read_gmsh_mesh (mesh, in_file, 'V', 'p');
  free (mesh->vrt_part);


  // increment so offset is one
  for (i=0; i<mesh->n_elm; ++i) ++mesh->elm_idx[i];


  // write mesh partition
  strcpy (out_file, in_file);
  strcat (out_file, "_p");
  sprintf (buf, "%d", p);
  strcat (out_file, buf);
  printf ("Writing GMSH file \"%s\"\n", out_file);
  write_gmsh_mesh (mesh, out_file, 'e');


  // clean up and return
  quit:
  // MPI_Finalize ();
  return (err);
#else
  return 0;
#endif
}
