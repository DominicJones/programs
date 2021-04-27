#include "mshf.h"
#include "mesh.h"
#include "graph.h"
#include "machine.h"

#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <climits>
#include <cmath>


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




int match_faces (int n_fce, int* fce_vrt_ptr, int* fce_vrt_lst, int* fce_elms, int* ghst_fce, int** fce_mtch)
{
  int n_mtch, N_mtch, n_vrt, intl, extl, ghst, bnd, neig, pair;
  int i, jj0, jjn, sjj, jj, j, iid0, iidn, iid, id, jj0a, jjna, sjja, sj, l;
  int *fce_vrt_slt;
  int *vrt_fce_ptr, *vrt_fce_lst;
  int *mtch;

  int icmb, n_mtch_cmb[8], N_mtch_cmb[8];
  int np, p;

#ifdef USE_MPI
  MPI_Comm_size (MPI_COMM_WORLD, &np);
  MPI_Comm_rank (MPI_COMM_WORLD, &p);
#else
  np = 1;
  p = 0;
#endif

  for (i=0; i<8; ++i) {n_mtch_cmb[i]=0; N_mtch_cmb[i]=0;}

  // make a sorted copy of the face vertex lists
  fce_vrt_slt = (int*) malloc (fce_vrt_ptr[n_fce] * sizeof(int));

  for (i=0; i<n_fce; ++i) {
    jj0 = fce_vrt_ptr[i];
    jjn = fce_vrt_ptr[i+1];
    for (jj=jj0; jj<jjn; ++jj) {
      fce_vrt_slt[jj] = fce_vrt_lst[jj];
    }
    qsort (&fce_vrt_slt[jj0], jjn-jj0, sizeof(int), icmp);
  }


  // construct the vertex face lists
  graph_to_dual (n_fce, fce_vrt_ptr, fce_vrt_slt, &n_vrt, &vrt_fce_ptr, &vrt_fce_lst);


  *fce_mtch = (int*) calloc (n_fce, sizeof(int));
  mtch = *fce_mtch;


  // find matching faces
  n_mtch = 0;

  for (i=0; i<n_fce; ++i) {

    // find a vertex of this face which has
    // the least number of associated faces
    jj0 = fce_vrt_ptr[i];
    jjn = fce_vrt_ptr[i+1];

    jj = shortest_dual_row (jjn-jj0, &fce_vrt_slt[jj0], n_vrt, vrt_fce_ptr, vrt_fce_lst);
    j = fce_vrt_slt[jj0+jj];

    // loop over the listed faces associated
    // to the sought vertex
    iid0 = vrt_fce_ptr[j];
    iidn = vrt_fce_ptr[j+1];

    for (iid=iid0; iid<iidn; ++iid) {
      id = vrt_fce_lst[iid];

      // ignore identical faces
      if (id - i == 0) continue;

      jj0a = fce_vrt_ptr[id];
      jjna = fce_vrt_ptr[id+1];

      // check that they have the same number of vertices
      if (jjn - jj0 == jjna - jj0a) {

	// check that this current face is not already matched
	if (mtch[id] == 0) {

	  // compare vertex lists of the two faces
	  sj = 0;
	  for (l=0; l<(jjn-jj0); ++l) {
	    sj += abs (fce_vrt_slt[jj0+l] - fce_vrt_slt[jj0a+l]);
	  }
	  if (sj == 0) {
	    if (fce_elms[i*2] > 0) {
	      intl = i; extl = id;
	    } else {
	      intl = id; extl = i;
	    }


            neig = 0;
	    if (fce_elms[extl*2] > 0) {
	      neig = 1; // neighbour is internal
	    }
	    else if (fce_elms[extl*2+1] < 0) {
	      neig = 2; // neighbour is a boundary
	    }


            pair = 0;
	    if (ghst_fce[intl] == 0 && ghst_fce[extl] == 0) {
	      pair = 1; // real - real
	    }
	    else if (ghst_fce[intl] == 0 && ghst_fce[extl] != 0) {
	      pair = 2; // real - ghost
	    }
	    else if (ghst_fce[intl] != 0 && ghst_fce[extl] == 0) {
	      pair = 3; // ghost - real
	    }
	    else if (ghst_fce[intl] != 0 && ghst_fce[extl] != 0) {
	      pair = 4; // ghost - ghost
	    }


        switch (neig) {
	    case 1:
	      fce_elms[intl*2+1] = fce_elms[extl*2];
	      break;
	    case 2:
	      fce_elms[intl*2+1] = fce_elms[extl*2+1];
	      break;
	    }

        icmb = 4 * (neig - 1) + (pair - 1);
        ++n_mtch_cmb[icmb];

	    mtch[i] = (id+1);
	    mtch[id] = -(i+1);
	    ++n_mtch;
	    break;
	  }
	}
      }
    }
  }

#ifdef USE_MPI
  MPI_Allreduce (&n_mtch_cmb, &N_mtch_cmb, 8, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  for (i=0; i<8; ++i) { N_mtch_cmb[i] = n_mtch_cmb[i]; }
#endif

  n_mtch = (n_mtch_cmb[0] + n_mtch_cmb[4]) * 2 + n_mtch_cmb[1];
  N_mtch = (N_mtch_cmb[0] + N_mtch_cmb[4]) * 2 + N_mtch_cmb[1];


  if (n_mtch_cmb[5] > 0) {
    printf ("[%d] error: real domain cell matched with ghost boundary cell\n",p);
    exit;
  }

  if (n_mtch_cmb[6] > 0) {
    printf ("[%d] error: ghost domain cell matched with real boundary cell\n",p);
    exit;
  }

  free (fce_vrt_slt);
  return (N_mtch);
}


void print_range (char const * name, int mrow, int* row)
{
  int i, jj0, jjn, jj, j;

  printf ("\n\"%s\" graph:\n", name);
  for (i=0; i<mrow; ++i) {
    jj0 = row[i];
    jjn = row[i+1];
    printf ("row %d: [%d %d)\n",i,jj0,jjn);
  }
  printf ("\n");
}


void print_graph (char const * name, int mrow, int* row, int* col)
{
  int i, jj0, jjn, jj, j;

  printf ("\n\"%s\" graph:\n", name);
  for (i=0; i<mrow; ++i) {
    jj0 = row[i];
    jjn = row[i+1];
    printf ("row %d: [%d %d), col[", i, jj0, jjn);
    for (jj=jj0; jj<jjn; ++jj) {
      j = col[jj];
      printf (" %d",j);
    }
    printf ("]\n");
  }
  printf ("\n");
}


void print_igraph (char const * name, int mrow, int* idx, int* row, int* col)
{
  int i, jj0, jjn, jj, j;

  printf ("\n\"%s\" graph:\n", name);
  for (i=0; i<mrow; ++i) {
    jj0 = row[i];
    jjn = row[i+1];
    printf ("row %d -> %d: [%d %d), col[", i, idx[i], jj0, jjn);
    for (jj=jj0; jj<jjn; ++jj) {
      j = col[jj];
      printf (" %d",j);
    }
    printf ("]\n");
  }
  printf ("\n");
}


void offset_to_count (int m, int* row, int** nrow)
{
  int i, *pnrow;
  pnrow = *nrow = (int*) malloc (m * sizeof(int));

  for (i=0; i<m; ++i) {
    pnrow[i] = row[i+1] - row[i];
  }
}


void count_to_offset (int m, int* nrow, int** row)
{
  int i, *prow;

  prow = *row = (int*) malloc ((m+1) * sizeof(int));

  prow[0] = 0;
  for (i=0; i<m; ++i) {
    prow[i+1] = prow[i] + nrow[i];
  }
}


void add_pole_entry (int os, int m, int** row, int** col)
{
  int i, j, jj0, jjn, jj, kk0, kkn, kk, nnz;
  int *row0, *col0, *row1, *col1;
  int *mem;

  row1 = *row;
  col1 = *col;

  row0 = (int*) malloc ((m+1) * sizeof(int));
  memcpy (row0, row1, (m+1) * sizeof(int));

  nnz = row1[m];

  col0 = (int*) malloc (nnz * sizeof(int));
  memcpy (col0, col1, nnz * sizeof(int));

  free (row1);
  free (col1);
  row1 = *row = (int*) calloc ((m+1) , sizeof(int));
  col1 = *col = (int*) calloc ((nnz+m) , sizeof(int));

  row1[0] = 0;

  for (i=0; i<m; ++i) {
    jj0 = row0[i];
    jjn = row0[i+1];

    kk0 = row1[i];
    kkn = kk0 + (jjn - jj0) + 1;
    row1[i+1] = kkn;

    col1[kk0] = os + i;

    for (jj=jj0, kk=kk0 + 1; jj<jjn; ++jj, ++kk) {
      col1[kk] = col0[jj];
    }

    qsort (&col1[kk0], kkn-kk0, sizeof(int), icmp);
  }

  free (row0);
  free (col0);
}


void purge_graph (int m, int* row, int** col)
{
  int i, j, jj0, jjn, jj;
  int nnz;
  int *cnt, *swap, *pcol;

  pcol = *col;
  nnz = row[m];
  swap = (int*) malloc (nnz * sizeof(int));
  memcpy (swap, pcol, nnz * sizeof(int));
  free (pcol);

  cnt = (int*) calloc (m, sizeof(int));

  nnz = 0;
  for (i=0; i<m; ++i) {
    jj0 = row[i];
    jjn = row[i+1];

    qsort (&swap[jj0], jjn-jj0, sizeof(int), icmp);
    setcp (&swap[jj0], jjn-jj0, INT_MAX);
    qsort (&swap[jj0], jjn-jj0, sizeof(int), icmp);

    for (jj=jj0; jj<jjn; ++jj) {
      if (swap[jj] < INT_MAX) {
        ++cnt[i];
        ++nnz;
      }
    }
  }

  pcol = *col = (int*) malloc (nnz * sizeof(int));

  j = 0;
  for (i=0; i<m; ++i) {
    jj0 = row[i];
    jjn = jj0 + cnt[i];
    for (jj=jj0; jj<jjn; ++jj) {
      pcol[j++] = swap[jj];
    }
  }
  free (swap);

  row[0] = 0;
  for (i=0; i<m; ++i) {
    row[i+1] = row[i] + cnt[i];
  }
  free (cnt);
}


void split_graph (
  int m_grph, int* grph_ia, int* grph_ja, int* flg,
  int* m_prt0, int** prt0_ia, int** prt0_ja,
  int* m_prt1, int** prt1_ia, int** prt1_ja)
{
  int i, j, jj0, jjn, jj, kk0, kk;
  int m, nnz, m0, nnz0, m1, nnz1;
  int *ia0, *ja0, *ia1, *ja1;

  m0 = nnz0 = m1 = nnz1 = 0;
  for (i=0; i<m_grph; ++i) {
    jj0 = grph_ia[i];
    jjn = grph_ia[i+1];
    if (flg[i]) {
      ++m0;
      nnz0 += (jjn - jj0);
    } else {
      ++m1;
      nnz1 += (jjn - jj0);
    }
  }

  *m_prt0 = m0;
  ia0 = *prt0_ia = (int*) malloc ((m0+1) * sizeof(int));
  ja0 = *prt0_ja = (int*) malloc ((nnz0) * sizeof(int));

  *m_prt1 = m1;
  ia1 = *prt1_ia = (int*) malloc ((m1+1) * sizeof(int));
  ja1 = *prt1_ja = (int*) malloc ((nnz1) * sizeof(int));

  m0 = nnz0 = m1 = nnz1 = 0;
  for (i=0; i<m_grph; ++i) {
    jj0 = grph_ia[i];
    jjn = grph_ia[i+1];
    if (flg[i]) {
      ia0[1+m0++] = jjn - jj0;
      for (jj=jj0; jj<jjn; ++jj) {
	ja0[nnz0++] = grph_ja[jj];
      }
    } else {
      ia1[1+m1++] = jjn - jj0;
      for (jj=jj0; jj<jjn; ++jj) {
	ja1[nnz1++] = grph_ja[jj];
      }
    }
  }

  ia0[0] = 0;
  for (i=1; i<m0+1; ++i) {
    ia0[i] += ia0[i-1];
  }

  ia1[0] = 0;
  for (i=1; i<m1+1; ++i) {
    ia1[i] += ia1[i-1];
  }
}


void merge_graphs (
  int m_prt0, int* prt0_ia, int* prt0_ja,
  int m_prt1, int* prt1_ia, int* prt1_ja,
  int* m_grph, int** grph_ia, int** grph_ja)
{
  int i, j, jj0, jjn, jj, kk0, kk;
  int m, nnz, m0, nnz0, m1, nnz1;
  int *ia, *ja;

  m0 = nnz0 = m1 = nnz1 = 0;

  for (i=0; i<m_prt0; ++i) {
    jj0 = prt0_ia[i];
    jjn = prt0_ia[i+1];
    ++m0;
    nnz0 += (jjn - jj0);
  }

  for (i=0; i<m_prt1; ++i) {
    jj0 = prt1_ia[i];
    jjn = prt1_ia[i+1];
    ++m1;
    nnz1 += (jjn - jj0);
  }

  m = m0 + m1;
  nnz = nnz0 + nnz1;

  *m_grph = m;
  ia = *grph_ia = (int*) malloc ((m+1) * sizeof(int));
  ja = *grph_ja = (int*) malloc ((nnz) * sizeof(int));

  memcpy (ja, prt0_ja, nnz0 * sizeof(int));
  memcpy (&ja[nnz0], prt1_ja, nnz1 * sizeof(int));

  j = 0;
  ia[j++] = 0;

  for (i=0; i<m_prt0; ++i) {
    jj0 = prt0_ia[i];
    jjn = prt0_ia[i+1];
    ia[j++] = jjn - jj0;
  }

  for (i=0; i<m_prt1; ++i) {
    jj0 = prt1_ia[i];
    jjn = prt1_ia[i+1];
    ia[j++] = jjn - jj0;
  }

  for (i=1; i<m+1; ++i) {
    ia[i] += ia[i-1];
  }
}


void reorder_graph (int m_grph, int* grph_ia, int* grph_ja, int* ord)
{
  int m, nnz;
  int i, j, jj0, jjn, jj, kk0, kk;
  int *swap_ia, *swap_ja;

  m = m_grph;
  swap_ia = (int*) malloc ((m+1) * sizeof(int));
  swap_ia[0] = 0;

  nnz = grph_ia[m];
  swap_ja = (int*) malloc (nnz * sizeof(int));

  for (i=0; i<m; ++i) {
    jj0 = grph_ia[i];
    jjn = grph_ia[i+1];
    swap_ia[ord[i]+1] = jjn - jj0;
  }

  for (i=0; i<m; ++i) {
    swap_ia[i+1] += swap_ia[i];
  }

  for (i=0; i<m; ++i) {
    jj0 = grph_ia[i];
    jjn = grph_ia[i+1];
    kk0 = swap_ia[ord[i]];
    for (jj=jj0, kk=kk0; jj<jjn; ++jj, ++kk) {
      swap_ja[kk] = grph_ja[jj];
    }
  }

  for (i=0; i<m+1; ++i) {
    grph_ia[i] = swap_ia[i];
  }

  for (i=0; i<nnz; ++i) {
    grph_ja[i] = swap_ja[i];
  }

  free (swap_ia);
  free (swap_ja);
}


void tuple_to_graph (int m_ord, int m_tup, int* tuple,
		     int* m_grph, int** grph_ia, int** grph_ja)
{
  int m, nnz, *ia, *ja, *chk;
  int i, ii, j, jj0, jjn, jj, k, kk0, kk;


  chk = (int*) calloc (m_tup, sizeof(int));


  // filter out tuples containing negative values
  m = 0;
  for (i=0; i<m_tup; ++i) {
    jj0 = m_ord*i;
    jjn = m_ord*(i+1);
    for (jj=jj0; jj<jjn; ++jj) {
      j = tuple[jj];
      if (jj == jjn - 1) {
	++m;
	chk[i] = 1;
      }
    }
  }

  if (m < 1) return;

  *m_grph = m;
  ia = *grph_ia = (int*) malloc ((m+1) * sizeof(int));
  ia[0] = 0;

  nnz = m_ord * m;
  ja = *grph_ja = (int*) malloc (nnz * sizeof(int));

  k = 0;
  for (i=0; i<m_tup; ++i) {
    if (chk[i] == 0) continue;
    jj0 = m_ord*i;
    jjn = m_ord*(i+1);
    ia[k+1] = ia[k] + m_ord;
    for (jj=jj0, kk=ia[k]; jj<jjn; ++jj, ++kk) {
      ja[kk] = tuple[jj];
    }
    ++k;
  }

  free (chk);
}


void graph_to_matrix (int m_grph, int* grph_ia, int* grph_ja,
		     int* m_mat, int** mat_ia, int** mat_ja)
{
  int m_dual, *dual_ia, *dual_ja;
  int m, nnz, *ia, *ja;
  int i, ii, j, jj0, jjn, jj, k, kk0, kk, diag, neig;


  graph_to_dual (m_grph, grph_ia, grph_ja, &m_dual, &dual_ia, &dual_ja);


  m = *m_mat = m_dual;
  ia = *mat_ia = (int*) malloc ((m+1) * sizeof(int));
  ia[0] = 0;


  nnz = m_dual + dual_ia[m_dual];
  ja = *mat_ja = (int*) malloc (nnz * sizeof(int));


  for (k=0; k<m_dual; ++k) {
    jj0 = dual_ia[k];
    jjn = dual_ia[k+1];
    diag = k;

    kk0 = ia[k];
    ia[k+1] = kk0 + (1 + (jjn - jj0));
    ja[kk0] = diag;

    for (jj=jj0, kk=kk0; jj<jjn; ++jj, ++kk) {
      j = dual_ja[jj];

      ii = grph_ia[j];
      neig = grph_ja[ii];

      if (diag == neig) {
	ii = grph_ia[j] + 1;
        neig = grph_ja[ii];
      }

      ja[kk+1] = neig;
    }
  }

  for (i=0; i<m_dual; ++i) {
    jj0 = ia[i];
    jjn = ia[i+1];
    qsort (&ja[jj0], jjn-jj0, sizeof(int), icmp);
  }

  free (dual_ia);
  free (dual_ja);
}


void graph_to_dual (int m_grph, int* grph_ia, int* grph_ja,
		   int* m_dual, int** dual_ia, int** dual_ja)
{
  int i, jj0, jjn, jj, j, jjd;
  int nnz_grph, nnz_dual;
  int* count;
  int *ia, *ja;

  nnz_grph = grph_ia[m_grph];

  *m_dual = 0;
  for (i=0; i<m_grph; ++i) {
    jj0 = grph_ia[i];
    jjn = grph_ia[i+1];
    for (jj=jj0; jj<jjn; ++jj) {
      j = grph_ja[jj];
      *m_dual = (j > *m_dual ? j : *m_dual);
    }
  }
  ++(*m_dual);

  count = (int*) calloc (*m_dual, sizeof(int));

  for (i=0; i<m_grph; ++i) {
    jj0 = grph_ia[i];
    jjn = grph_ia[i+1];
    for (jj=jj0; jj<jjn; ++jj) {
      j = grph_ja[jj];
      if (j >= 0) ++count[j];
    }
  }

  *dual_ia = (int*) malloc ((*m_dual+1) * sizeof(int));
  ia = *dual_ia;

  *dual_ja = (int*) malloc (nnz_grph * sizeof(int));
  ja = *dual_ja;

  ia[0] = 0;
  for (i=0; i<(*m_dual); ++i) {
    ia[i+1] = ia[i] + count[i];
    count[i] = 0;
  }

  for (i=0; i<m_grph; ++i) {
    jj0 = grph_ia[i];
    jjn = grph_ia[i+1];
    for (jj=jj0; jj<jjn; ++jj) {
      j = grph_ja[jj];
      if (j >= 0) {
        jjd = ia[j] + count[j];
        ja[jjd] = i;
        ++count[j];
      }
    }
  }

  free (count);
}


int shortest_dual_row (int n_row, int* row, int m_dual, int* dual_ia, int* dual_ja)
{
  int i, n, jj, j, jjd0, jjdn, jn;

  n = INT_MAX;
  i = 0;

  for (jj=0; jj<n_row; ++jj) {
    j = row[jj];

    jjd0 = dual_ia[j];
    jjdn = dual_ia[j+1];
    jn = jjdn - jjd0;

    if (jn < n) {
      n = jn;
      i = jj;
    }
  }

  return (i);
}


int dual_column_length(int * grph_ja, int djjn, int * dual_ia, int * dual_ja)
{
  int length = 0;
  for (int jj=0; jj<djjn; ++jj) {
    int i = grph_ja[jj];
    int jj0 = dual_ia[i];
    int jjn = dual_ia[i+1];
    length += jjn-jj0;
  }
  return length;
}


void dual_column_fill(int * grph_ja, int djjn, int * dual_ia, int * dual_ja, int ** grph_neig)
{
  int idx = 0;
  int * neig = *grph_neig;
  for (int jj=0; jj<djjn; ++jj) {
    int i = grph_ja[jj];
    int jj0 = dual_ia[i];
    int jjn = dual_ia[i+1];
    for (int k=jj0; k<jjn; ++k) {
      neig[idx] = dual_ja[k];
      ++idx;
    }
  }
}


// void mesh_to_dual_api(int m_msh, int* msh_ia, int* msh_ja, int os, int min_mtch, int** dual_ia, int** dual_ja)
// {
// #ifdef USE_MPI
//   MPI_Comm mpi_comm_world;
//   mpi_comm_world = MPI_COMM_WORLD;

//   ParMETIS_V3_Mesh2Dual (
//     mesh->dmn_blk, mesh->dmn_vrt_ptr, mesh->dmn_vrt_lst,
//     &idx_offset, &min_fce_vrt,
//     &mesh->dmn_neig_ptr, &mesh->dmn_neig_lst,
//     &mpi_comm_world);
// #else
//   mesh_to_dual (mesh->dmn_blk[np], mesh->dmn_vrt_ptr, mesh->dmn_vrt_lst,
//                 idx_offset, min_fce_vrt,
//                 &mesh->dmn_neig_ptr, &mesh->dmn_neig_lst);
// #endif
// }


void mesh_to_dual (int m_msh, int* msh_ia, int* msh_ja, int os, int min_mtch, int** dual_ia, int** dual_ja)
{
  int * msh_ja_s;
  sort_column_entries(m_msh, msh_ia, msh_ja, &msh_ja_s);

  int m_mshdual;
  int * mshdual_ia;
  int * mshdual_ja;
  graph_to_dual (m_msh, msh_ia, msh_ja_s, &m_mshdual, &mshdual_ia, &mshdual_ja);

  int nnz_max = 0;
  for (int i=0; i<m_msh; ++i) {
    int jj0 = msh_ia[i];
    int jjn = msh_ia[i+1];
    nnz_max += dual_column_length(&msh_ja_s[jj0], jjn-jj0, mshdual_ia, mshdual_ja);
  }

  int * ia = *dual_ia = (int*) malloc ((m_msh+1) * sizeof(int));
  int * ja = *dual_ja = (int*) malloc ((nnz_max) * sizeof(int));
  ia[0] = 0;

  std::vector<int> neig;
  for (int i=0; i<m_msh; ++i) {
    int jj0 = msh_ia[i];
    int jjn = msh_ia[i+1];
    int djj = dual_column_length(&msh_ja_s[jj0], jjn-jj0, mshdual_ia, mshdual_ja);
    if (neig.size() < djj) {
      neig.resize(djj);
    }
    int * neig_ptr = &neig[0];
    dual_column_fill(&msh_ja_s[jj0], jjn-jj0, mshdual_ia, mshdual_ja, &neig_ptr);
    std::sort(neig.begin(), neig.end());
    neig.erase(std::remove(neig.begin(), neig.end(), i), neig.end());

    if (min_mtch > 1) {
      // remove any unique values
      setuq (neig_ptr, neig.size(), INT_MAX);
      neig.erase(std::remove(neig.begin(), neig.end(), INT_MAX), neig.end());

      // remove duplicates iteratively
      for (int imtch=1; imtch<min_mtch; ++imtch) {
        set_last_cp (neig_ptr, neig.size(), INT_MAX);
        neig.erase(std::remove(neig.begin(), neig.end(), INT_MAX), neig.end());
      }
    }

    // remove duplicates above the minimum matches
    neig.erase(std::unique(neig.begin(), neig.end()), neig.end());

    ia[i+1] = ia[i] + neig.size();
    int k = 0;
    for (int jj=ia[i]; jj<ia[i+1]; ++jj) {
      ja[jj] = neig[k];
      ++k;
    }
  }

  free (msh_ja_s);
  free (mshdual_ia);
  free (mshdual_ja);
}


void sort_column_entries(int m_grph, int* grph_ia, int* grph_ja, int** grph_ja_s)
{
  *grph_ja_s = (int*) malloc (grph_ia[m_grph] * sizeof(int));
  int * ja_s = *grph_ja_s;

  for (int i=0; i<m_grph; ++i) {
    int jj0 = grph_ia[i];
    int jjn = grph_ia[i+1];
    for (int jj=jj0; jj<jjn; ++jj) {
      ja_s[jj] = grph_ja[jj];
    }
    qsort (&ja_s[jj0], jjn-jj0, sizeof(int), icmp);
  }
}




void find_ghost_elements (
  struct ghst_t* ghst,
  int n_elm, int* elm_blk, int* elm_neig_ptr, int* elm_neig_lst)
{
  int i, jj0, jjn, jj, j, k, n;
  int p = 0, np = 1, op;

#ifdef USE_MPI
  MPI_Comm_size (MPI_COMM_WORLD, &np);
  MPI_Comm_rank (MPI_COMM_WORLD, &p);
#endif

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
  int p = 0, np = 1, op;

#ifdef USE_MPI
  int mpi_tag = 100;
  MPI_Request* mpi_rqst;
  MPI_Comm_size (MPI_COMM_WORLD, &np);
  MPI_Comm_rank (MPI_COMM_WORLD, &p);
  mpi_rqst = (MPI_Request*) malloc (np * sizeof(MPI_Request));
#endif

  // how many boundary definitions are required?
  locl->elm_cnt = (int*) calloc (np, sizeof(int));

#ifdef USE_MPI
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
#endif

  // which boundaries are they (what are their indices)?
  count_to_offset (np, locl->elm_cnt, &locl->elm_blk);
  m = locl->elm_blk[np];
  locl->elm_idx = (int*) malloc (m * sizeof(int));

#ifdef USE_MPI
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
#endif


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

#ifdef USE_MPI
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
#endif

  // send the requested boundary vertex list
  m = ghst->elm_blk[np];
  count_to_offset (m, ghst->elm_vrt_cnt, &ghst->elm_vrt_ptr);
  nnz = ghst->elm_vrt_ptr[m];
  ghst->elm_vrt_lst = (int*) malloc (nnz * sizeof(int));

#ifdef USE_MPI
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
#endif

  // send the requested boundary tag count
  m = ghst->elm_blk[np];
  ghst->elm_tag_cnt = (int*) malloc (m * sizeof(int));

#ifdef USE_MPI
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
#endif

  // send the requested boundary tag list
  m = ghst->elm_blk[np];
  count_to_offset (m, ghst->elm_tag_cnt, &ghst->elm_tag_ptr);
  nnz = ghst->elm_tag_ptr[m];
  ghst->elm_tag_lst = (int*) malloc (nnz * sizeof(int));

#ifdef USE_MPI
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
#endif

  // negate ghost boundary physical tag
  m = ghst->elm_blk[np];
  for (i=0; i<m; ++i) {
    jj0 = ghst->elm_tag_ptr[i];
    ghst->elm_tag_lst[jj0+1] *= -1;
  }

#ifdef USE_MPI
  free (mpi_rqst);
#endif

  quit:
  return (err);
}




int icmp(const void *a, const void *b)
{
  const int *ia = (const int *)a;
  const int *ib = (const int *)b;
  return *ia  - *ib;
}


// SET REPEATED VALUES IN LIST OF INTEGERS
void setcp (int* lst, int n_lst, int iset)
{
  int i, j, k;

  i = 0;
  while (1) {

    for (j=i+1; j<n_lst; ++j) {
      if (lst[j] != lst[i]) {
        break;
      }
    }

    for (k=i+1; k<j; ++k) {
      lst[k] = iset;
    }

    i = j;
    if (i >= n_lst-1) break;
  }
}


// SET LAST REPEATED VALUE IN LIST OF INTEGERS
void set_last_cp (int* lst, int n_lst, int iset)
{
  int i, j, k;

  i = 0;
  while (1) {

    for (j=i+1; j<n_lst; ++j) {
      if (lst[j] != lst[i]) {
        break;
      }
    }

    lst[j-1] = iset;

    i = j;
    if (i >= n_lst-1) break;
  }
}


// SET UNIQUE VALUES IN LIST OF INTEGERS
void setuq (int* lst, int n_lst, int iset)
{
  int i;

  if (lst[0] != lst[1]) {
    lst[0] = iset;
  }

  for (i=1; i<n_lst-1; ++i) {
    if (lst[i] != lst[i-1] && lst[i] != lst[i+1]) {
      lst[i] = iset;
    }
  }

  if (lst[n_lst-1] != lst[n_lst-2]) {
    lst[n_lst-1] = iset;
  }
}


void print_list (char const * name, int m, int* lst)
{
  int i;

  printf ("\n\"%s\" list[0..%d]:", name, m);
  for (i=0; i<m; ++i) {
    printf (" %d", lst[i]);
  }
  printf ("\n");
}


int offset_lb (int m, int* blk, int val)
{
  int i, idx;

  idx = -1;
  for (i=0; i<m; ++i) {
    if (val >= blk[i] && val < blk[i+1]) {
      idx = i;
      break;
    }
  }

  return (idx);
}


void split_list (int m, int* lst, int* flg, int* m0, int** lst0, int* m1, int** lst1)
{
  int i, n0, n1;
  int *plst0, *plst1;

  n0 = n1 = 0;
  for (i=0; i<m; ++i) {
    if (flg[i]) {
      ++n0;
    } else {
      ++n1;
    }
  }

  plst0 = *lst0 = (int*) malloc (n0 * sizeof(int));
  memcpy (plst0, lst, n0  * sizeof(int));
  *m0 = n0;

  plst1 = *lst1 = (int*) malloc (n1 * sizeof(int));
  memcpy (plst1, &lst[n0], n1  * sizeof(int));
  *m1 = n1;
}


void merge_lists (int m0, int* lst0, int m1, int* lst1, int* m, int** lst)
{
  int *plst;
  plst = *lst = (int*) malloc ((m0+m1) * sizeof(int));
  memcpy (plst, lst0, m0  * sizeof(int));
  memcpy (&plst[m0], lst1, m1  * sizeof(int));
  *m = m0 + m1;
}


void uniq (int m0, int* lst0, int* m1, int** lst1)
{
  int i, pm1, *plst1;

  plst1 = *lst1;
  plst1 = (int*) malloc (m0 * sizeof(int));
  memcpy (plst1, lst0, m0 * sizeof(int));

  qsort (plst1, m0, sizeof(int), icmp);
  setcp (plst1, m0, INT_MAX);
  qsort (plst1, m0, sizeof(int), icmp);

  pm1 = 0;
  for (i = 0; i < m0; ++i) {
    if (plst1[i] < INT_MAX) {
      ++pm1;
    }
  }

  *m1 = pm1;
}



/*
subroutine init_ghost_update (bnd, dmn)
  use hash_m

  type(mesh_t),intent(inout)::bnd, dmn

  integer :: err = 0
  integer :: i, ii, j, jj1, jjn, k, m, n
  integer :: blk, n_blk, oblk
  integer, dimension(:), allocatable :: mpi_rqst
  integer, dimension(:,:), allocatable :: mpi_stat

  call MPI_Comm_size (mpi_comm_world, n_blk, err)
  call MPI_Comm_rank (mpi_comm_world, blk, err)

  allocate (mpi_rqst(n_blk)); mpi_rqst = 0
  allocate (mpi_stat(mpi_status_size, n_blk)); mpi_stat = 0


  !! setup ghost element index block offsets (ascending ghst_idx assumed)
  allocate (bnd%ghst_blk(n_blk + 1))
  bnd%ghst_blk = 0

  do i = 1, bnd%n_ghst
    call hash_get_value (bnd%gi2b, bnd%ghst_idx(i), oblk)
    bnd%ghst_blk(oblk + 1) = bnd%ghst_blk(oblk + 1) + 1
  end do

  bnd%ghst_blk(1) = 1
  do i = 2, n_blk + 1
    bnd%ghst_blk(i) = bnd%ghst_blk(i) + bnd%ghst_blk(i - 1)
  end do


  allocate (dmn%ghst_blk(n_blk + 1))
  dmn%ghst_blk = 0

  do i = 1, dmn%n_ghst
    call hash_get_value (dmn%gi2b, dmn%ghst_idx(i), oblk)
    dmn%ghst_blk(oblk + 1) = dmn%ghst_blk(oblk + 1) + 1
  end do

  dmn%ghst_blk(1) = 1
  do i = 2, n_blk + 1
    dmn%ghst_blk(i) = dmn%ghst_blk(i) + dmn%ghst_blk(i - 1)
  end do


  !! determine the ghosts possessed by this [partition]
  allocate (bnd%locl_blk(n_blk + 1))
  bnd%locl_blk = 0

  do i = 0, n_blk - 1
    if (i == blk) cycle
    call MPI_Irecv (bnd%locl_blk(i+2), 1, MPI_INTEGER, i, &
      mpi_any_tag, mpi_comm_world, mpi_rqst(i+1), err)
  end do

  do i = 0, n_blk - 1
    if (i == blk) cycle
    n = bnd%ghst_blk(i + 2) - bnd%ghst_blk(i + 1)
    call MPI_Send (n, 1, MPI_INTEGER, i, &
      100, mpi_comm_world, err)
  end do

  call MPI_Waitall (n_blk, mpi_rqst, mpi_stat, err)

  bnd%locl_blk(1) = 1
  do i = 2, n_blk + 1
    bnd%locl_blk(i) = bnd%locl_blk(i) + bnd%locl_blk(i - 1)
  end do


  allocate (dmn%locl_blk(n_blk + 1))
  dmn%locl_blk = 0

  do i = 0, n_blk - 1
    if (i == blk) cycle
    call MPI_Irecv (dmn%locl_blk(i+2), 1, MPI_INTEGER, i, &
      mpi_any_tag, mpi_comm_world, mpi_rqst(i+1), err)
  end do

  do i = 0, n_blk - 1
    if (i == blk) cycle
    n = dmn%ghst_blk(i + 2) - dmn%ghst_blk(i + 1)
    call MPI_Send (n, 1, MPI_INTEGER, i, &
      100, mpi_comm_world, err)
  end do

  call MPI_Waitall (n_blk, mpi_rqst, mpi_stat, err)

  dmn%locl_blk(1) = 1
  do i = 2, n_blk + 1
    dmn%locl_blk(i) = dmn%locl_blk(i) + dmn%locl_blk(i - 1)
  end do


  !! which elements are they (what are their indices)?
  n = bnd%locl_blk(n_blk + 1) - 1
  bnd%n_locl = n

  allocate (bnd%locl_idx(n))
  bnd%locl_idx = 0

  do i = 0, n_blk - 1
    jj1 = bnd%locl_blk(i + 1)
    jjn = bnd%locl_blk(i + 2) - 1
    n = jjn - jj1 + 1
    if (i == blk .or. n == 0) cycle
    call MPI_Irecv (bnd%locl_idx(jj1:jjn), n, MPI_INTEGER, i, &
      mpi_any_tag, mpi_comm_world, mpi_rqst(i+1), err)
  end do

  do i = 0, n_blk - 1
    jj1 = bnd%ghst_blk(i + 1)
    jjn = bnd%ghst_blk(i + 2) - 1
    n = jjn - jj1 + 1
    if (i == blk .or. n == 0) cycle
    call MPI_Send (bnd%ghst_idx(jj1:jjn), n, MPI_INTEGER, i, &
      100, mpi_comm_world, err)
  end do

  call MPI_Waitall (n_blk, mpi_rqst, mpi_stat, err)

  n = dmn%locl_blk(n_blk + 1) - 1
  dmn%n_locl = n

  allocate (dmn%locl_idx(n))
  dmn%locl_idx = 0

  do i = 0, n_blk - 1
    jj1 = dmn%locl_blk(i + 1)
    jjn = dmn%locl_blk(i + 2) - 1
    n = jjn - jj1 + 1
    if (i == blk .or. n == 0) cycle
    call MPI_Irecv (dmn%locl_idx(jj1:jjn), n, MPI_INTEGER, i, &
      mpi_any_tag, mpi_comm_world, mpi_rqst(i+1), err)
  end do

  do i = 0, n_blk - 1
    jj1 = dmn%ghst_blk(i + 1)
    jjn = dmn%ghst_blk(i + 2) - 1
    n = jjn - jj1 + 1
    if (i == blk .or. n == 0) cycle
    call MPI_Send (dmn%ghst_idx(jj1:jjn), n, MPI_INTEGER, i, &
      100, mpi_comm_world, err)
  end do

  call MPI_Waitall (n_blk, mpi_rqst, mpi_stat, err)

!   call print_graph ("boundary ghost indices", &
!      n_blk, bnd%ghst_blk, bnd%ghst_idx)
!
!   call print_graph ("boundary local indices", &
!      n_blk, bnd%locl_blk, bnd%locl_idx)
end subroutine
*/


/*
subroutine ghost_update_s (bnd, dmn, phi, ubnd, udmn)
  use mpi
  use hash_m

  type(mesh_t),intent(in)::bnd, dmn
  type(sfld_t),intent(inout)::phi
  integer,intent(in)::ubnd, udmn

  integer :: err = 0
  integer :: i, ii, j, jj1, jjn, k, m, n
  integer :: blk, n_blk, oblk
  integer, dimension(:), allocatable :: mpi_rqst
  integer, dimension(:,:), allocatable :: mpi_stat


  call MPI_Comm_size (mpi_comm_world, n_blk, err)
  call MPI_Comm_rank (mpi_comm_world, blk, err)

  allocate (mpi_rqst(n_blk)); mpi_rqst = 0
  allocate (mpi_stat(mpi_status_size, n_blk)); mpi_stat = 0


  !! send the requested elmement values
  if (ubnd /= 0) then
    n = bnd%locl_blk(n_blk + 1) - 1
    do i = 1, n
      j = bnd%locl_idx(i)
      k = j - bnd%elm_blk(blk+1) + 1
      phi%lbnd(i) = phi%bnd(k)
    end do


    do i = 0, n_blk - 1
      jj1 = bnd%ghst_blk(i + 1)
      jjn = bnd%ghst_blk(i + 2) - 1
      n = jjn - jj1 + 1
      if (i == blk .or. n == 0) cycle
      call MPI_Irecv (phi%gbnd(jj1:jjn), n, mpi_rp, i, &
        mpi_any_tag, mpi_comm_world, mpi_rqst(i+1), err)
    end do

    do i = 0, n_blk - 1
      jj1 = bnd%locl_blk(i + 1)
      jjn = bnd%locl_blk(i + 2) - 1
      n = jjn - jj1 + 1
      if (i == blk .or. n == 0) cycle
      call MPI_Send (phi%lbnd(jj1:jjn), n, mpi_rp, i, &
        100, mpi_comm_world, err)
    end do

    call MPI_Waitall (n_blk, mpi_rqst, mpi_stat, err)
  end if


  if (udmn /= 0) then
    n = dmn%locl_blk(n_blk + 1) - 1
    do i = 1, n
      j = dmn%locl_idx(i)
      k = j - dmn%elm_blk(blk+1) + 1
      phi%ldmn(i) = phi%dmn(k)
    end do


    do i = 0, n_blk - 1
      jj1 = dmn%ghst_blk(i + 1)
      jjn = dmn%ghst_blk(i + 2) - 1
      n = jjn - jj1 + 1
      if (i == blk .or. n == 0) cycle
      call MPI_Irecv (phi%gdmn(jj1:jjn), n, mpi_rp, i, &
        mpi_any_tag, mpi_comm_world, mpi_rqst(i+1), err)
    end do

    do i = 0, n_blk - 1
      jj1 = dmn%locl_blk(i + 1)
      jjn = dmn%locl_blk(i + 2) - 1
      n = jjn - jj1 + 1
      if (i == blk .or. n == 0) cycle
      call MPI_Send (phi%ldmn(jj1:jjn), n, mpi_rp, i, &
        100, mpi_comm_world, err)
    end do

    call MPI_Waitall (n_blk, mpi_rqst, mpi_stat, err)
  end if


!   do i = 1, size(phi%gbnd)
!     j = bnd%ghst_idx(i)
!     call hash_get_value (bnd%gi2b, j, k)
!     print *, blk+1, "phi_ghst",i,j,k,phi%gbnd(i)
!   end do
end subroutine
*/
