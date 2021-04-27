// -*- C++ -*-
#pragma once

struct mshf_t {
  int n_elm;

  int* elm_idx;
  int* elm_ty;
  int* elm_dim;

  int* elm_fce_ptr;
  int* fce_vrt_ptr;
  int* fce_vrt_lst;
};

void read_mesh_format (struct mshf_t* mshf, char* filename);
