#ifndef mesh_h
#define mesh_h

struct mesh_t {
  int n_dim;
  int part, n_part;

  int N_vrt;
  int n_vrt;
  int *vrt_idx, *vrt_part;
  double *x_vrt;

  int N_elm;
  int *elm_dim, *elm_part;

  int n_gbnd, *gbnd_blk, *gbnd_cnt, *gbnd_idx, *gbnd_vrt_ptr, *gbnd_vrt_lst, *gbnd_tag_ptr, *gbnd_tag_lst;
  int n_gdmn, *gdmn_blk, *gdmn_cnt, *gdmn_idx, *gdmn_vrt_ptr, *gdmn_vrt_lst, *gdmn_tag_ptr, *gdmn_tag_lst;

  int n_bnd, *bnd_blk, *bnd_cnt, *bnd_idx, *bnd_vrt_ptr, *bnd_vrt_lst, *bnd_neig_ptr, *bnd_neig_lst, *bnd_tag_ptr, *bnd_tag_lst;
  int n_dmn, *dmn_blk, *dmn_cnt, *dmn_idx, *dmn_vrt_ptr, *dmn_vrt_lst, *dmn_neig_ptr, *dmn_neig_lst, *dmn_tag_ptr, *dmn_tag_lst;
  int n_elm, *elm_blk, *elm_cnt, *elm_idx, *elm_vrt_ptr, *elm_vrt_lst, *elm_neig_ptr, *elm_neig_lst, *elm_tag_ptr, *elm_tag_lst;

  int n_xbnd, *xbnd_elm_ptr, *xbnd_elm_lst, *xbnd_vrt_ptr, *xbnd_vrt_lst;
  int n_xdmn, *xdmn_elm_ptr, *xdmn_elm_lst, *xdmn_vrt_ptr, *xdmn_vrt_lst;
  int n_xelm, *xelm_elm_ptr, *xelm_elm_lst, *xelm_vrt_ptr, *xelm_vrt_lst;

  int n_pb_xbnd, *pb_xbnd_ptr, *pb_xbnd_lst;
  int n_pb_xdmn, *pb_xdmn_ptr, *pb_xdmn_lst;
};


struct ghst_t {
  int n_elm;
  int *elm_cnt, *elm_blk;
  int *elm_idx;
  int *elm_vrt_cnt, *elm_vrt_ptr, *elm_vrt_lst;
  int *elm_tag_cnt, *elm_tag_ptr, *elm_tag_lst;
};


void distributed_block (int n, int np, int** blk);
void partitioned_block (int n, int np, int* idx, int** blk);

void read_gmsh_mesh (struct mesh_t* mesh, char* filename, char section, char read_mode);
void write_gmsh_mesh (struct mesh_t* mesh, char* filename, char region);

extern "C" {
int write_gmsh_field_ (int* pflen, char* fname, int* pn, int* elm_idx, double* phi);
int write_gmsh_vector_field_ (int* pflen, char* fname, int* pn, int* dn, int* elm_idx, double* phi);
}

#endif
