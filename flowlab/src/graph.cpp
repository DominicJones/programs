#include "utils.h"
#include "graph.h"
#ifdef USE_MPI
#include "mpi.h"
#endif

#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <climits>
#include <cmath>


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
