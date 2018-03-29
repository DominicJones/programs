#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "utils.h"
#include "settings.h"

#include "mesh.h"
#include "mshf.h"
#include "map.h"
#include "graph.h"
#include "xgrph.h"

#include "mpi.h"
#include "parmetis.h"


int mpi_size, mpi_rank;
int mpi_tag = 100;
MPI_Status* mpi_stat;
MPI_Request* mpi_rqst;
MPI_Comm mpi_comm_world;


extern "C" {
  void psolv_cfd_
  (
   int*, char*,
   int*, int*, int*, int*,
   int*, int*, int*, double*,
   int*, int*, int*, int*,
   int*, int*, int*,
   int*, int*, int*,
   int*, int*, int*,
   int*, int*, int*,
   int*, int*, int*,
   int*, int*, int*,
   int*, int*, int*,
   int*, int*, int*,
   int*, int*, int*,
   int*, int*, int*,
   int*, int*,
   int*, int*, int*,
   int*, int*,
   int*, int*,
   int*, int*);
}


int main (int argc, char* argv[])
{
  int err = 0;
  char buf[256], in_file[256], out_file[256], mshffile[256];
  int n_dim, ty, dim, idx, idmn, ibnd, ielm, n_dmn, n_bnd, n_elm;
  int igbl, iloc;
  int *ord, *cnt, *blk, *ptr, *lst, *swap, *gswap;
  int ii, ii0, iin, i, jj0, jjn, jj, j, kk0, kkn, kk, k, kkz;
  int l, m, n, llz, mmz, nnz, p, op, np, tt0, ph_ty;
  int idx_offset = 0, min_fce_vrt;
  float *xlst;

  int key, value;
  struct Map *pair;
  struct Map *map = NULL;

  struct mshf_t *mshf;
  struct mesh_t *mesh;
  struct ghst_t *gdmn, *ldmn;
  struct ghst_t *gbnd, *lbnd;

  FILE *fctrl;
  Settings *ctrl;


  // allocation and initialisation
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);

  mpi_comm_world = MPI_COMM_WORLD;
  mpi_stat = (MPI_Status*) malloc (mpi_size * sizeof(MPI_Status));
  mpi_rqst = (MPI_Request*) malloc (mpi_size * sizeof(MPI_Request));

  mshf = (struct mshf_t*) malloc (sizeof(struct mshf_t));
  mesh = (struct mesh_t*) malloc (sizeof(struct mesh_t));

  p = mesh->part = mpi_rank;
  np = mesh->n_part = mpi_size;


  // read control file
  fctrl = fopen(argv[1], "r");
  if (fctrl == NULL) {
    printf("error: no control file specified\n");
    err = 1; goto quit;
  }

  ctrl = settings_open(fctrl);
  fclose(fctrl);
  if (ctrl == NULL) {
    printf("error: problem closing the control file\n");
    err = 1; goto quit;
  }

  settings_get(ctrl, "Mesh", "File", in_file, sizeof(in_file));
  settings_get(ctrl, "Mesh", "Format", mshffile, sizeof(mshffile));
  mesh->n_dim = settings_get_int(ctrl, "Mesh", "Dimensions");

  settings_delete(ctrl);


  // read mesh format
  read_mesh_format (mshf, mshffile);


  // read mesh element data in a distributed manner
  if (np > 1) {
    strcat (in_file, "_p");
    sprintf (buf, "%d", p);
    strcat (in_file, buf);
  }
  read_gmsh_mesh (mesh, in_file, 'A', 's');


  // determine mesh partition connectivity
  // N.B. boundary faces are listed first
  connectivity (mshf, mesh);


  // split elements into boundary and domain
  swap = (int*) malloc (mesh->n_elm * sizeof(int));

  for (i=0; i<mesh->n_elm; ++i) {
    tt0 = mesh->elm_tag_ptr[i];
    ty = mesh->elm_tag_lst[tt0];
    idx = mshf->elm_idx[ty];
    dim = mshf->elm_dim[idx];
    swap[i] = mesh->n_dim - dim;
  }

  split_list (
    mesh->n_elm, mesh->elm_idx, swap,
    &mesh->n_bnd, &mesh->bnd_idx,
    &mesh->n_dmn, &mesh->dmn_idx);

  split_graph (
    mesh->n_elm, mesh->elm_vrt_ptr, mesh->elm_vrt_lst, swap,
    &mesh->n_bnd, &mesh->bnd_vrt_ptr, &mesh->bnd_vrt_lst,
    &mesh->n_dmn, &mesh->dmn_vrt_ptr, &mesh->dmn_vrt_lst);

  split_graph (
    mesh->n_elm, mesh->elm_tag_ptr, mesh->elm_tag_lst, swap,
    &mesh->n_bnd, &mesh->bnd_tag_ptr, &mesh->bnd_tag_lst,
    &mesh->n_dmn, &mesh->dmn_tag_ptr, &mesh->dmn_tag_lst);

  free (swap);


  // split boundary into real and ghost
  swap = (int*) malloc (mesh->n_bnd * sizeof(int));

  for (i=0; i<mesh->n_bnd; ++i) {
    tt0 = mesh->bnd_tag_ptr[i];
    ph_ty = mesh->bnd_tag_lst[tt0+1];
    if (ph_ty < 0) ph_ty = 0;
    swap[i] = ph_ty;
  }


  split_list (
    mesh->n_bnd, mesh->bnd_idx, swap,
    &m, &lst,
    &mesh->n_gbnd, &mesh->gbnd_idx);

  free (mesh->bnd_idx);

  mesh->bnd_idx = (int*) malloc (m * sizeof(int));
  memcpy (mesh->bnd_idx, lst, m * sizeof(int));

  free (lst);


  split_graph (
    mesh->n_bnd, mesh->bnd_vrt_ptr, mesh->bnd_vrt_lst, swap,
    &m, &ptr, &lst,
    &mesh->n_gbnd, &mesh->gbnd_vrt_ptr, &mesh->gbnd_vrt_lst);

  free (mesh->bnd_vrt_ptr);
  free (mesh->bnd_vrt_lst);

  mesh->bnd_vrt_ptr = (int*) malloc ((m+1) * sizeof(int));
  mesh->bnd_vrt_lst = (int*) malloc (ptr[m] * sizeof(int));
  memcpy (mesh->bnd_vrt_ptr, ptr, (m+1) * sizeof(int));
  memcpy (mesh->bnd_vrt_lst, lst, ptr[m] * sizeof(int));

  free (ptr);
  free (lst);


  split_graph (
    mesh->n_bnd, mesh->bnd_tag_ptr, mesh->bnd_tag_lst, swap,
    &m, &ptr, &lst,
    &mesh->n_gbnd, &mesh->gbnd_tag_ptr, &mesh->gbnd_tag_lst);

  free (mesh->bnd_tag_ptr);
  free (mesh->bnd_tag_lst);

  mesh->n_bnd = m;
  mesh->bnd_tag_ptr = (int*) malloc ((m+1) * sizeof(int));
  mesh->bnd_tag_lst = (int*) malloc (ptr[m] * sizeof(int));
  memcpy (mesh->bnd_tag_ptr, ptr, (m+1) * sizeof(int));
  memcpy (mesh->bnd_tag_lst, lst, ptr[m] * sizeof(int));

  free (ptr);
  free (lst);

  free (swap);
  mesh->n_bnd = m;


  // split domain into real and ghost
  swap = (int*) malloc (mesh->n_dmn * sizeof(int));

  for (i=0; i<mesh->n_dmn; ++i) {
    tt0 = mesh->dmn_tag_ptr[i];
    ph_ty = mesh->dmn_tag_lst[tt0+1];
    if (ph_ty < 0) ph_ty = 0;
    swap[i] = ph_ty;
  }


  split_list (
    mesh->n_dmn, mesh->dmn_idx, swap,
    &m, &lst,
    &mesh->n_gdmn, &mesh->gdmn_idx);

  free (mesh->dmn_idx);

  mesh->dmn_idx = (int*) malloc (m * sizeof(int));
  memcpy (mesh->dmn_idx, lst, m * sizeof(int));

  free (lst);


  split_graph (
    mesh->n_dmn, mesh->dmn_vrt_ptr, mesh->dmn_vrt_lst, swap,
    &m, &ptr, &lst,
    &mesh->n_gdmn, &mesh->gdmn_vrt_ptr, &mesh->gdmn_vrt_lst);

  free (mesh->dmn_vrt_ptr);
  free (mesh->dmn_vrt_lst);

  mesh->dmn_vrt_ptr = (int*) malloc ((m+1) * sizeof(int));
  mesh->dmn_vrt_lst = (int*) malloc (ptr[m] * sizeof(int));
  memcpy (mesh->dmn_vrt_ptr, ptr, (m+1) * sizeof(int));
  memcpy (mesh->dmn_vrt_lst, lst, ptr[m] * sizeof(int));

  free (ptr);
  free (lst);


  split_graph (
    mesh->n_dmn, mesh->dmn_tag_ptr, mesh->dmn_tag_lst, swap,
    &m, &ptr, &lst,
    &mesh->n_gdmn, &mesh->gdmn_tag_ptr, &mesh->gdmn_tag_lst);

  free (mesh->dmn_tag_ptr);
  free (mesh->dmn_tag_lst);

  mesh->dmn_tag_ptr = (int*) malloc ((m+1) * sizeof(int));
  mesh->dmn_tag_lst = (int*) malloc (ptr[m] * sizeof(int));
  memcpy (mesh->dmn_tag_ptr, ptr, (m+1) * sizeof(int));
  memcpy (mesh->dmn_tag_lst, lst, ptr[m] * sizeof(int));

  free (ptr);
  free (lst);

  free (swap);
  mesh->n_dmn = m;


  // set boundary block offsets
  mesh->bnd_blk = (int*) malloc ((np+1) * sizeof(int));
  MPI_Allgather (&mesh->n_bnd, 1, MPI_INT, &mesh->bnd_blk[1], 1, MPI_INT,
    MPI_COMM_WORLD);

  mesh->bnd_blk[0] = 0;
  for (i=0; i<np; ++i) {
    mesh->bnd_blk[i+1] += mesh->bnd_blk[i];
  }
  offset_to_count (np, mesh->bnd_blk, &mesh->bnd_cnt);


  // set domain block offsets
  mesh->dmn_blk = (int*) malloc ((np+1) * sizeof(int));
  MPI_Allgather (&mesh->n_dmn, 1, MPI_INT, &mesh->dmn_blk[1], 1, MPI_INT,
    MPI_COMM_WORLD);

  mesh->dmn_blk[0] = 0;
  for (i=0; i<np; ++i) {
    mesh->dmn_blk[i+1] += mesh->dmn_blk[i];
  }
  offset_to_count (np, mesh->dmn_blk, &mesh->dmn_cnt);


  // construct CSR boundary matrix
  min_fce_vrt = mesh->n_dim - 1;

  ParMETIS_V3_Mesh2Dual (
    mesh->bnd_blk, mesh->bnd_vrt_ptr, mesh->bnd_vrt_lst,
    &idx_offset, &min_fce_vrt,
    &mesh->bnd_neig_ptr, &mesh->bnd_neig_lst,
    &mpi_comm_world);

  n = mesh->bnd_blk[p];
  add_pole_entry (n, mesh->n_bnd,
    &mesh->bnd_neig_ptr, &mesh->bnd_neig_lst);


  // construct CSR domain matrix
  min_fce_vrt = mesh->n_dim;

  ParMETIS_V3_Mesh2Dual (
    mesh->dmn_blk, mesh->dmn_vrt_ptr, mesh->dmn_vrt_lst,
    &idx_offset, &min_fce_vrt,
    &mesh->dmn_neig_ptr, &mesh->dmn_neig_lst,
    &mpi_comm_world);

  add_pole_entry (mesh->dmn_blk[p], mesh->n_dmn,
    &mesh->dmn_neig_ptr, &mesh->dmn_neig_lst);


  // split connectivity into boundary and domain
  swap = (int*) malloc (mesh->n_xelm * sizeof(int));

  for (i=0; i<mesh->n_xelm; ++i) {
    jj0 = mesh->xelm_elm_ptr[i];
    j = mesh->xelm_elm_lst[jj0];
    --j;

    op = offset_lb (np, mesh->bnd_blk, j);
    swap[i] = 1;

    if (op < 0) {
      j -= mesh->bnd_blk[np];
      op = offset_lb (np, mesh->dmn_blk, j);
      if (op >= 0) {
        swap[i] = 0;
      } else {
        err = 1; goto quit;
      }
    }
  }

  split_graph (
    mesh->n_xelm, mesh->xelm_vrt_ptr, mesh->xelm_vrt_lst, swap,
    &mesh->n_xbnd, &mesh->xbnd_vrt_ptr, &mesh->xbnd_vrt_lst,
    &mesh->n_xdmn, &mesh->xdmn_vrt_ptr, &mesh->xdmn_vrt_lst);

  split_graph (
    mesh->n_xelm, mesh->xelm_elm_ptr, mesh->xelm_elm_lst, swap,
    &mesh->n_xbnd, &mesh->xbnd_elm_ptr, &mesh->xbnd_elm_lst,
    &mesh->n_xdmn, &mesh->xdmn_elm_ptr, &mesh->xdmn_elm_lst);

  free (swap);


  // pass control to CFD solver
  n = strlen (in_file);

  psolv_cfd_ (
    &n, in_file,
    &p, &np, mesh->bnd_blk, mesh->dmn_blk,

    &mesh->n_dim, &mesh->n_vrt, mesh->vrt_idx, mesh->x_vrt,

    &mshf->n_elm, mshf->elm_idx, mshf->elm_ty, mshf->elm_dim,
    mshf->elm_fce_ptr, mshf->fce_vrt_ptr, mshf->fce_vrt_lst,

    &mesh->n_bnd, mesh->bnd_idx, mesh->bnd_tag_ptr,
    mesh->bnd_tag_lst, mesh->bnd_vrt_ptr, mesh->bnd_vrt_lst,

    &mesh->n_gbnd, mesh->gbnd_idx, mesh->gbnd_tag_ptr,
    mesh->gbnd_tag_lst, mesh->gbnd_vrt_ptr, mesh->gbnd_vrt_lst,

    &mesh->n_dmn, mesh->dmn_idx, mesh->dmn_tag_ptr,
    mesh->dmn_tag_lst, mesh->dmn_vrt_ptr, mesh->dmn_vrt_lst,

    &mesh->n_gdmn, mesh->gdmn_idx, mesh->gdmn_tag_ptr,
    mesh->gdmn_tag_lst, mesh->gdmn_vrt_ptr, mesh->gdmn_vrt_lst,

    &mesh->n_xbnd, mesh->xbnd_vrt_ptr, mesh->xbnd_vrt_lst,
    mesh->xbnd_elm_ptr, mesh->xbnd_elm_lst,

    &mesh->n_xdmn, mesh->xdmn_elm_ptr, mesh->xdmn_elm_lst,
    mesh->xdmn_vrt_ptr, mesh->xdmn_vrt_lst,

    mesh->bnd_neig_ptr, mesh->bnd_neig_lst,
    mesh->dmn_neig_ptr, mesh->dmn_neig_lst
  );


  quit:
  MPI_Finalize ();
  if (err) fprintf (stderr, "An error occurred (%d)\n", err);
  return (err);
}
