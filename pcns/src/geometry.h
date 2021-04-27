// -*- C++ -*-
#pragma once

#include "vector.h"
#include "tensor.h"

#include <vector>
#include <unordered_map>


struct grphblk_t {
  int mblk=0;
  int mrow=0;
  int * nblk = nullptr;
  int * blk = nullptr;
  int * nrow = nullptr;
  int * row = nullptr;
  int * col = nullptr;
  int * ival = nullptr;
  int * rval = nullptr;
};

struct mshblk_t {
  int blk = 0;
  int n_blk = 0;
  int n_dim = 0;
  int n_vrt = 0;
  int n_elm = 0;
  int n_ghst = 0;
  int n_locl = 0;
  int n_fce = 0;
  int * vrt_idx = nullptr;
  int * elm_blk = nullptr;
  int * elm_idx = nullptr;
  int * ghst_blk = nullptr;
  int * ghst_idx = nullptr;
  int * locl_blk = nullptr;
  int * locl_idx = nullptr;
  grphblk_t elm_vrt, elm_tag, elm_neig;
  grphblk_t ghst_vrt, ghst_tag;
  grphblk_t fce_vrt, fce_elm;
  std::unordered_map<int, int> vi2i;
  std::unordered_map<int, int> gi2i, gi2b;
};

struct sfld_t {
  int id = 0;
  std::vector<double> bnd;
  std::vector<double> gbnd;
  std::vector<double> lbnd;
  std::vector<double> dmn;
  std::vector<double> gdmn;
  std::vector<double> ldmn;
};

struct vfld_t {
  int id = 0;
  std::vector<vector_t<3, double> > bnd;
  std::vector<vector_t<3, double> > gbnd;
  std::vector<vector_t<3, double> > lbnd;
  std::vector<vector_t<3, double> > dmn;
  std::vector<vector_t<3, double> > gdmn;
  std::vector<vector_t<3, double> > ldmn;
};

struct tfld_t {
  int id = 0;
  std::vector<tensor_t<3, double> > bnd;
  std::vector<tensor_t<3, double> > gbnd;
  std::vector<tensor_t<3, double> > lbnd;
  std::vector<tensor_t<3, double> > dmn;
  std::vector<tensor_t<3, double> > gdmn;
  std::vector<tensor_t<3, double> > ldmn;
};

struct geom_t {
  double vol_min = 0;
  double vol_max = 0;
  double vol_sum = 0;
  sfld_t vol;
  vfld_t x_vc, gl_min, gl_max, gl;
  tfld_t ls, gg;
  std::vector<vector_t<3, double> > x_vrt;
  std::vector<vector_t<3, double> > x_fc;
  std::vector<vector_t<3, double> > x_pc;
  std::vector<vector_t<3, double> > x_nc;
  std::vector<vector_t<3, double> > norm;
  std::vector<double> l_pn;
  std::vector<double> w_fp;
  std::vector<double> area;
  std::vector<double> volr;
};


void geometry(mshblk_t const &bnd, mshblk_t const &dmn, geom_t &xbnd, geom_t &xdmn);
