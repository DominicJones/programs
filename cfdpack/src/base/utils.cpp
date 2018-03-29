#include "utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>


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


void print_list (char* name, int m, int* lst)
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
