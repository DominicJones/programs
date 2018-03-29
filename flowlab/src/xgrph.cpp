#include "mshf.h"
#include "mesh.h"
#include "graph.h"
#include "xgrph.h"
#ifdef USE_MPI
#include "mpi.h"
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <climits>


int connectivity (struct mshf_t* mshf, struct mesh_t* mesh)
{
  int err = 0;

  int i, j, k, l, m, n;
  int ii0, iin, ii, jj0, jjn, jj, kk0, kkn, kk, ll0, lln, ll;
  int tt0, ttn, tt;
  int n_dim, n_elm, n_fce, n_vrt, igl, ty, idx, dim, ph_ty;
  int *ord;
  int *fce_vrt_ptr, *fce_vrt_lst, *fce_elms;
  int *vrt_fce_ptr, *vrt_fce_lst;
  int *ghst_fce;

  int n_bnd_fce = 0, n_dmn_fce = 0;

  int np = 1, p = 0;

#ifdef USE_MPI
  MPI_Comm_size (MPI_COMM_WORLD, &np);
  MPI_Comm_rank (MPI_COMM_WORLD, &p);
#endif

  n_vrt = 0;
  n_fce = 0;

  for (i=0; i<mesh->n_elm; ++i) {
    tt0 = mesh->elm_tag_ptr[i];
    ty = mesh->elm_tag_lst[tt0];
    idx = mshf->elm_idx[ty];
    dim = mshf->elm_dim[idx];

    ph_ty = mesh->elm_tag_lst[tt0+1];

    jj0 = mshf->elm_fce_ptr[idx];
    jjn = mshf->elm_fce_ptr[idx+1];
    n_fce += (jjn - jj0);
    if (ph_ty > 0) n_dmn_fce += (jjn - jj0);
    for (jj=jj0; jj<jjn; ++jj) {
      kk0 = mshf->fce_vrt_ptr[jj];
      kkn = mshf->fce_vrt_ptr[jj+1];
      n_vrt += (kkn - kk0);
    }

    if (dim == mesh->n_dim - 1) {
      ++n_fce;
      if (ph_ty > 0) ++n_bnd_fce;
      jj0 = mshf->elm_fce_ptr[idx];
      jjn = mshf->elm_fce_ptr[idx+1];
      n_vrt += (jjn - jj0);
    }
  }


  fce_vrt_ptr = (int*) malloc ((n_fce+1) * sizeof(int));
  fce_vrt_ptr[0] = 0;
  fce_vrt_lst = (int*) malloc (n_vrt * sizeof(int));
  fce_elms = (int*) calloc ((n_fce*2), sizeof(int));
  ghst_fce = (int*) calloc ((n_fce), sizeof(int));


  // build the face vertex list and the face element list
  m = 0;

  for (i=0; i<mesh->n_elm; ++i) {
    tt0 = mesh->elm_tag_ptr[i];
    ty = mesh->elm_tag_lst[tt0];
    idx = mshf->elm_idx[ty];
    dim = mshf->elm_dim[idx];

    ph_ty = mesh->elm_tag_lst[tt0+1];

    igl= mesh->elm_idx[i];
    ii0 = mesh->elm_vrt_ptr[i];

    jj0 = mshf->elm_fce_ptr[idx];
    jjn = mshf->elm_fce_ptr[idx+1];

    for (jj=jj0; jj<jjn; ++jj) {
      kk0 = mshf->fce_vrt_ptr[jj];
      kkn = mshf->fce_vrt_ptr[jj+1];
      ll0 = fce_vrt_ptr[m];

      if (ph_ty < 0) ghst_fce[m] = 1;

      fce_elms[m*2] = (igl+0); // internal entered as [#, ]

      fce_vrt_ptr[m+1] = ll0 + (kkn - kk0);
      for (kk=kk0, ll=ll0; kk<kkn; ++kk, ++ll) {
        k = mshf->fce_vrt_lst[kk];
        fce_vrt_lst[ll] = mesh->elm_vrt_lst[ii0+k];
      }
      ++m;
    }

    if (dim == mesh->n_dim - 1) {
      jj0 = mshf->elm_fce_ptr[idx];
      jjn = mshf->elm_fce_ptr[idx+1];
      ll0 = fce_vrt_ptr[m];

      if (ph_ty < 0) ghst_fce[m] = 1;

      fce_elms[m*2+1] = -(igl+0); // boundary entered as [ ,-#]

      fce_vrt_ptr[m+1] = ll0 + (jjn - jj0);
      for (jj=jj0, ll=ll0, k=0; jj<jjn; ++jj, ++ll, ++k) {
        fce_vrt_lst[ll] = mesh->elm_vrt_lst[ii0+k];
      }
      ++m;
    }
  }


  int N_bnd_fce = 0, N_dmn_fce = 0;

#ifdef USE_MPI
  MPI_Allreduce (&n_bnd_fce, &N_bnd_fce, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (&n_dmn_fce, &N_dmn_fce, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  N_bnd_fce = n_bnd_fce;
  N_dmn_fce = n_dmn_fce;
#endif

  //   printf ("[%d] n_bnd_fce %d, N_bnd_fce %d\n", p, n_bnd_fce, N_bnd_fce);
  //   printf ("[%d] n_dmn_fce %d, N_dmn_fce %d\n", p, n_dmn_fce, N_dmn_fce);
  //   printf ("[%d] # global real faces %d\n", p, N_bnd_fce + N_dmn_fce);

  //   graph_to_dual(n_fce, fce_vrt_ptr, fce_vrt_lst,
  //                 &n_vrt, &vrt_fce_ptr, &vrt_fce_lst);

  int *fce_mtch, N_mtch;
  N_mtch = match_faces (n_fce, fce_vrt_ptr, fce_vrt_lst, fce_elms, ghst_fce, &fce_mtch);

  if (N_bnd_fce + N_dmn_fce != N_mtch) {
    printf ("[%d] Error in attempting to match element faces\n", p);
    exit;
  }


  mesh->n_xelm = 0;
  n_vrt = 0;
  for (i=0; i<n_fce; ++i) {
    if ((abs (fce_elms[i*2]) > 0 && abs (fce_elms[i*2+1]) > 0)
        && (fce_elms[i*2] > 0 || fce_elms[i*2+1] > 0)
        && (ghst_fce[i] == 0)) {

      jj0 = fce_vrt_ptr[i];
      jjn = fce_vrt_ptr[i+1];
      n_vrt += (jjn - jj0);
      ++mesh->n_xelm;
    }
  }


  mesh->xelm_elm_ptr = (int*) malloc ((mesh->n_xelm+1) * sizeof(int));
  mesh->xelm_elm_ptr[0] = 0;
  mesh->xelm_elm_lst = (int*) malloc ((mesh->n_xelm*2) * sizeof(int));

  mesh->xelm_vrt_ptr = (int*) malloc ((mesh->n_xelm+1) * sizeof(int));
  mesh->xelm_vrt_ptr[0] = 0;
  mesh->xelm_vrt_lst = (int*) malloc (n_vrt * sizeof(int));


  ii = 0;
  for (i=0; i<n_fce; ++i) {
    if ((abs (fce_elms[i*2]) > 0 && abs (fce_elms[i*2+1]) > 0)
        && (fce_elms[i*2] > 0 || fce_elms[i*2+1] > 0)
        && (ghst_fce[i] == 0)) {

      mesh->xelm_elm_ptr[ii+1] = mesh->xelm_elm_ptr[ii] + 2;

      jj0 = mesh->xelm_vrt_ptr[ii];
      jjn = jj0 + (fce_vrt_ptr[i+1] - fce_vrt_ptr[i]);
      mesh->xelm_vrt_ptr[ii+1] = jjn;

      mesh->xelm_elm_lst[ii*2] = fce_elms[i*2];
      mesh->xelm_elm_lst[ii*2+1] = fce_elms[i*2+1];

      for (jj=jj0, j=fce_vrt_ptr[i]; jj<jjn; ++jj, ++j) {
        mesh->xelm_vrt_lst[jj] = fce_vrt_lst[j];
      }
      ++ii;
    }
  }

  free (fce_elms);
  free (fce_vrt_ptr);
  free (fce_vrt_lst);
  free (ghst_fce);


  // "domain to boundary" faces listed first
  ord = (int*) malloc (mesh->n_xelm * sizeof(int));
  m = 0; n = mesh->n_xelm - 1;

  for (i=0; i<mesh->n_xelm; ++i) {
    jjn = mesh->xelm_elm_ptr[i+1];
    j = mesh->xelm_elm_lst[jjn-1];
    if (j < 0) {
      ord[i] = m; ++m;
    } else {
      ord[i] = n; --n;
    }
  }

  reorder_graph (mesh->n_xelm, mesh->xelm_elm_ptr, mesh->xelm_elm_lst, ord);
  reorder_graph (mesh->n_xelm, mesh->xelm_vrt_ptr, mesh->xelm_vrt_lst, ord);

  free (ord);

  // print_graph ("fce_elm", mesh->n_xelm, mesh->xelm_elm_ptr, mesh->xelm_elm_lst);
  // print_graph ("fce_vrt", mesh->n_xelm, mesh->xelm_vrt_ptr, mesh->xelm_vrt_lst);

 quit:
  return (err);
}
