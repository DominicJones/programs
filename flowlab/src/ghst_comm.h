#ifndef ghst_comm_h
#define ghst_comm_h

#include "mesh.h"

void find_ghost_elements (
  struct ghst_t* ghst,
  int n_elm, int* elm_blk, int* elm_neig_ptr, int* elm_neig_lst);

int communicate_ghost_elements (
  struct ghst_t* ghst, struct ghst_t* locl,
  int* elm_blk,
  int* elm_tag_ptr, int* elm_tag_lst,
  int* elm_vrt_ptr, int* elm_vrt_lst);

#endif
