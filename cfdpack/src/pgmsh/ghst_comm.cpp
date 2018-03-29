#include <stdlib.h>
#include "utils.h"
#include "graph.h"
#include "ghst_comm.h"
#include "mpi.h"


void find_ghost_elements (
  struct ghst_t* ghst,
  int n_elm, int* elm_blk, int* elm_neig_ptr, int* elm_neig_lst)
{
  int i, jj0, jjn, jj, j, k, n;
  int p, np, op;

  MPI_Comm_size (MPI_COMM_WORLD, &np);
  MPI_Comm_rank (MPI_COMM_WORLD, &p);


  ghst->elm_cnt = (int*) calloc (np, sizeof(int));

  for (i=0; i<n_elm; ++i) {
    jj0 = elm_neig_ptr[i];
    jjn = elm_neig_ptr[i+1];

    for (jj=jj0; jj<jjn; ++jj) {
      j = elm_neig_lst[jj];

      for (op=0; op<np; ++op) {
        if (j >= elm_blk[op] && j < elm_blk[op+1]) {
          if (op == p) break;
          ++ghst->elm_cnt[op];
        }
      }
    }
  }

  count_to_offset (np, ghst->elm_cnt, &ghst->elm_blk);
  memset (ghst->elm_cnt, 0, np * sizeof(int));

  n = ghst->elm_blk[np];
  ghst->elm_idx = (int*) malloc (n * sizeof(int));

  for (i=0; i<n_elm; ++i) {
    jj0 = elm_neig_ptr[i];
    jjn = elm_neig_ptr[i+1];

    for (jj=jj0; jj<jjn; ++jj) {
      j = elm_neig_lst[jj];

      for (op=0; op<np; ++op) {
        if (j >= elm_blk[op] && j < elm_blk[op+1]) {
          if (op == p) break;
          k = ghst->elm_blk[op] + ghst->elm_cnt[op];
          ghst->elm_idx[k] = j;
          ++ghst->elm_cnt[op];
        }
      }
    }
  }

  purge_graph (np, ghst->elm_blk, &ghst->elm_idx);

  free (ghst->elm_cnt);
  offset_to_count (np, ghst->elm_blk, &ghst->elm_cnt);
}


int communicate_ghost_elements (
  struct ghst_t* ghst, struct ghst_t* locl,
  int* elm_blk,
  int* elm_tag_ptr, int* elm_tag_lst,
  int* elm_vrt_ptr, int* elm_vrt_lst)
{
  int err;
  int ii, i, jj0, jjn, jj, j, kk0, kkn, kk, k, n, m, nnz;
  int p, np, op;

  int mpi_tag = 100;
  MPI_Request* mpi_rqst;

  MPI_Comm_size (MPI_COMM_WORLD, &np);
  MPI_Comm_rank (MPI_COMM_WORLD, &p);

  mpi_rqst = (MPI_Request*) malloc (np * sizeof(MPI_Request));


  // how many boundary definitions are required?
  locl->elm_cnt = (int*) calloc (np, sizeof(int));

  for (i=0; i<np; ++i) {
    if (i == p) continue;
    MPI_Irecv(&locl->elm_cnt[i], 1, MPI_INT, i, MPI_ANY_TAG,
      MPI_COMM_WORLD, &mpi_rqst[i]);
  }

  for (i=0; i<np; ++i) {
    if (i == p) continue;
    MPI_Send(&ghst->elm_cnt[i], 1, MPI_INT, i, mpi_tag,
      MPI_COMM_WORLD);
  }

  MPI_Barrier (MPI_COMM_WORLD);


  // which boundaries are they (what are their indices)?
  count_to_offset (np, locl->elm_cnt, &locl->elm_blk);
  m = locl->elm_blk[np];
  locl->elm_idx = (int*) malloc (m * sizeof(int));

  for (i=0; i<np; ++i) {
    if (i == p || locl->elm_cnt[i] == 0) continue;
    jj0 = locl->elm_blk[i];
    jjn = locl->elm_blk[i+1];
    MPI_Irecv(&locl->elm_idx[jj0], jjn-jj0, MPI_INT, i, MPI_ANY_TAG,
      MPI_COMM_WORLD, &mpi_rqst[i]);
  }

  for (i=0; i<np; ++i) {
    if (i == p || ghst->elm_cnt[i] == 0) continue;
    jj0 = ghst->elm_blk[i];
    jjn = ghst->elm_blk[i+1];
    MPI_Send(&ghst->elm_idx[jj0], jjn-jj0, MPI_INT, i, mpi_tag,
      MPI_COMM_WORLD);
  }

  MPI_Barrier (MPI_COMM_WORLD);


  // construct local boundary vertex graph
  m = locl->elm_blk[np];
  locl->elm_vrt_ptr = (int*) malloc ((m+1) * sizeof(int));

  locl->elm_vrt_ptr[0] = 0;
  for (ii=0; ii<m; ++ii) {
    i = locl->elm_idx[ii];
    op = offset_lb (np, elm_blk, i);
    if (op < 0) { err = 1; goto quit; }
    i -= elm_blk[op];
    jj0 = elm_vrt_ptr[i];
    jjn = elm_vrt_ptr[i+1];
    kk0 = locl->elm_vrt_ptr[ii];
    kkn = kk0 + (jjn - jj0);
    locl->elm_vrt_ptr[ii+1] = kkn;
  }

  offset_to_count (m, locl->elm_vrt_ptr, &locl->elm_vrt_cnt);

  nnz = locl->elm_vrt_ptr[m];
  locl->elm_vrt_lst = (int*) malloc (nnz * sizeof(int));

  for (ii=0; ii<m; ++ii) {
    i = locl->elm_idx[ii];
    op = offset_lb (np, elm_blk, i);
    i -= elm_blk[op];
    jj0 = elm_vrt_ptr[i];
    jjn = elm_vrt_ptr[i+1];
    kk0 = locl->elm_vrt_ptr[ii];
    for (jj=jj0, kk=kk0; jj<jjn; ++jj, ++kk) {
      locl->elm_vrt_lst[kk] = elm_vrt_lst[jj];
    }
  }


  // construct local boundary tag graph
  m = locl->elm_blk[np];
  locl->elm_tag_ptr = (int*) malloc ((m+1) * sizeof(int));

  locl->elm_tag_ptr[0] = 0;
  for (ii=0; ii<m; ++ii) {
    i = locl->elm_idx[ii];
    op = offset_lb (np, elm_blk, i);
    i -= elm_blk[op];
    jj0 = elm_tag_ptr[i];
    jjn = elm_tag_ptr[i+1];
    kk0 = locl->elm_tag_ptr[ii];
    kkn = kk0 + (jjn - jj0);
    locl->elm_tag_ptr[ii+1] = kkn;
  }

  offset_to_count (m, locl->elm_tag_ptr, &locl->elm_tag_cnt);

  nnz = locl->elm_tag_ptr[m];
  locl->elm_tag_lst = (int*) malloc (nnz * sizeof(int));

  for (ii=0; ii<m; ++ii) {
    i = locl->elm_idx[ii];
    op = offset_lb (np, elm_blk, i);
    i -= elm_blk[op];
    jj0 = elm_tag_ptr[i];
    jjn = elm_tag_ptr[i+1];
    kk0 = locl->elm_tag_ptr[ii];
    for (jj=jj0, kk=kk0; jj<jjn; ++jj, ++kk) {
      locl->elm_tag_lst[kk] = elm_tag_lst[jj];
    }
  }


  // send the requested boundary vertex count
  m = ghst->elm_blk[np];
  ghst->elm_vrt_cnt = (int*) malloc (m * sizeof(int));

  for (i=0; i<np; ++i) {
    if (i == p || ghst->elm_cnt[i] == 0) continue;
    ii = ghst->elm_blk[i];
    MPI_Irecv(&ghst->elm_vrt_cnt[ii], ghst->elm_cnt[i], MPI_INT, i, MPI_ANY_TAG,
      MPI_COMM_WORLD, &mpi_rqst[i]);
  }

  for (i=0; i<np; ++i) {
    if (i == p || locl->elm_cnt[i] == 0) continue;
    ii = locl->elm_blk[i];
    MPI_Send(&locl->elm_vrt_cnt[ii], locl->elm_cnt[i], MPI_INT, i, mpi_tag,
      MPI_COMM_WORLD);
  }

  MPI_Barrier (MPI_COMM_WORLD);


  // send the requested boundary vertex list
  m = ghst->elm_blk[np];
  count_to_offset (m, ghst->elm_vrt_cnt, &ghst->elm_vrt_ptr);
  nnz = ghst->elm_vrt_ptr[m];
  ghst->elm_vrt_lst = (int*) malloc (nnz * sizeof(int));

  for (i=0; i<np; ++i) {
    if (i == p || ghst->elm_cnt[i] == 0) continue;
    jj0 = ghst->elm_blk[i];
    jjn = ghst->elm_blk[i+1];
    kk0 = ghst->elm_vrt_ptr[jj0];
    kkn = ghst->elm_vrt_ptr[jjn];
    MPI_Irecv(&ghst->elm_vrt_lst[kk0], kkn-kk0, MPI_INT, i, MPI_ANY_TAG,
      MPI_COMM_WORLD, &mpi_rqst[i]);
  }

  for (i=0; i<np; ++i) {
    if (i == p || locl->elm_cnt[i] == 0) continue;
    jj0 = locl->elm_blk[i];
    jjn = locl->elm_blk[i+1];
    kk0 = locl->elm_vrt_ptr[jj0];
    kkn = locl->elm_vrt_ptr[jjn];
    MPI_Send(&locl->elm_vrt_lst[kk0], kkn-kk0, MPI_INT, i, mpi_tag,
      MPI_COMM_WORLD);
  }

  MPI_Barrier (MPI_COMM_WORLD);


  // send the requested boundary tag count
  m = ghst->elm_blk[np];
  ghst->elm_tag_cnt = (int*) malloc (m * sizeof(int));

  for (i=0; i<np; ++i) {
    if (i == p || ghst->elm_cnt[i] == 0) continue;
    ii = ghst->elm_blk[i];
    MPI_Irecv(&ghst->elm_tag_cnt[ii], ghst->elm_cnt[i], MPI_INT, i, MPI_ANY_TAG,
      MPI_COMM_WORLD, &mpi_rqst[i]);
  }

  for (i=0; i<np; ++i) {
    if (i == p || locl->elm_cnt[i] == 0) continue;
    ii = locl->elm_blk[i];
    MPI_Send(&locl->elm_tag_cnt[ii], locl->elm_cnt[i], MPI_INT, i, mpi_tag,
      MPI_COMM_WORLD);
  }

  MPI_Barrier (MPI_COMM_WORLD);


  // send the requested boundary tag list
  m = ghst->elm_blk[np];
  count_to_offset (m, ghst->elm_tag_cnt, &ghst->elm_tag_ptr);
  nnz = ghst->elm_tag_ptr[m];
  ghst->elm_tag_lst = (int*) malloc (nnz * sizeof(int));

  for (i=0; i<np; ++i) {
    if (i == p || ghst->elm_cnt[i] == 0) continue;
    jj0 = ghst->elm_blk[i];
    jjn = ghst->elm_blk[i+1];
    kk0 = ghst->elm_tag_ptr[jj0];
    kkn = ghst->elm_tag_ptr[jjn];
    MPI_Irecv(&ghst->elm_tag_lst[kk0], kkn-kk0, MPI_INT, i, MPI_ANY_TAG,
      MPI_COMM_WORLD, &mpi_rqst[i]);
  }

  for (i=0; i<np; ++i) {
    if (i == p || locl->elm_cnt[i] == 0) continue;
    jj0 = locl->elm_blk[i];
    jjn = locl->elm_blk[i+1];
    kk0 = locl->elm_tag_ptr[jj0];
    kkn = locl->elm_tag_ptr[jjn];
    MPI_Send(&locl->elm_tag_lst[kk0], kkn-kk0, MPI_INT, i, mpi_tag,
      MPI_COMM_WORLD);
  }

  MPI_Barrier (MPI_COMM_WORLD);


  // negate ghost boundary physical tag
  m = ghst->elm_blk[np];
  for (i=0; i<m; ++i) {
    jj0 = ghst->elm_tag_ptr[i];
    ghst->elm_tag_lst[jj0+1] *= -1;
  }


  free (mpi_rqst);


  quit:
  return (err);
}
