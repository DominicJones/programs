#include "mesh.h"
#include "mshf.h"
#include "graph.h"
#include "args_parser.h"
#include "geometry.h"
#include "simulation.h"
#include "vector.h"
#include "machine.h"

#include "json/reader.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <string>
#include <fstream>
#include <streambuf>
#include <unordered_map>
#include <algorithm>
#include <iostream>


int mpi_size = 1, mpi_rank = 0;
int mpi_tag = 100;

#ifdef USE_MPI
MPI_Status* mpi_stat;
MPI_Request* mpi_rqst;
MPI_Comm mpi_comm_world;
#endif

extern "C" {
  void navierstokes_
  (
   int*, const char*,
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
  int dmn_blk_os;
  float *xlst;
  int key, value;

  struct mshf_t *mshf;
  struct mesh_t *mesh;
  struct ghst_t *gdmn, *ldmn;
  struct ghst_t *gbnd, *lbnd;

  args_parser_t args(argc, argv);

  if(args.has_option("-h") || args.size() < 2)
  {
    std::cout << "usage: " << args.arg(0) << " <json file>" << std::endl;
    return 0;
  }

  // machine_t::initialize(argc, argv);
  // mpi_rank = machine_t::rank();
  // mpi_size = machine_t::size();

  // allocation and initialisation
#ifdef USE_MPI
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);

  mpi_comm_world = MPI_COMM_WORLD;
  mpi_stat = (MPI_Status*) malloc (mpi_size * sizeof(MPI_Status));
  mpi_rqst = (MPI_Request*) malloc (mpi_size * sizeof(MPI_Request));
#endif

  mshf = (struct mshf_t*) malloc (sizeof(struct mshf_t));
  mesh = (struct mesh_t*) malloc (sizeof(struct mesh_t));

  std::string jsonfile{args.arg(1)};
  {
    std::ifstream jsonstream(jsonfile.c_str());
    std::string jsonfiletext((std::istreambuf_iterator<char>(jsonstream)), std::istreambuf_iterator<char>());

    Json::Value root;
    Json::Reader reader;
    bool success = reader.parse(jsonfiletext, root);
    if (!success) {
      std::cout  << "Failed to parse JSON file\n"
                 << reader.getFormatedErrorMessages();
      return 0;
    }

    std::string path(jsonfile.substr(0, jsonfile.find_last_of("/") + 1));
    strcpy(in_file, std::string(path + root["Mesh"]["File"].asString()).c_str());
    mesh->n_dim = root["Mesh"]["Dimensions"].asInt();
  }


  // read mesh format
  {
    std::string file("../etc/config/config.json");
    std::ifstream stream(file.c_str());
    std::string text((std::istreambuf_iterator<char>(stream)), std::istreambuf_iterator<char>());

    Json::Value root;
    Json::Reader reader;
    bool success = reader.parse(text, root);
    if (!success) {
      std::cout  << "Failed to parse JSON file\n"
                 << reader.getFormatedErrorMessages();
      return 0;
    }

    std::string path(file.substr(0, file.find_last_of("/") + 1));
    strcpy(mshffile, std::string(path + root["Mesh"]["Format"].asString()).c_str());

    read_mesh_format(mshf, mshffile);
  }


  // machine rank
  p = mesh->part = mpi_rank;
  np = mesh->n_part = mpi_size;


  // if (argc >= 3)
  {
    partition(mshf, mesh, in_file);
  }
  // else {
    // read mesh element data in a distributed manner
    if (np > 1) {
      strcat (in_file, "_p");
      sprintf (buf, "%d", p);
      strcat (in_file, buf);
    }
    read_gmsh_mesh (mesh, in_file, 'A', 's');
  // }


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
  mesh->bnd_blk = (int*) calloc ((np+1), sizeof(int));
#ifdef USE_MPI
  MPI_Allgather (&mesh->n_bnd, 1, MPI_INT, &mesh->bnd_blk[1], 1, MPI_INT,
    MPI_COMM_WORLD);
#else
  mesh->bnd_blk[1] = mesh->n_bnd;
#endif

  mesh->bnd_blk[0] = 0;
  for (i=0; i<np; ++i) {
    mesh->bnd_blk[i+1] += mesh->bnd_blk[i];
  }
  offset_to_count (np, mesh->bnd_blk, &mesh->bnd_cnt);


  // set domain block offsets
  mesh->dmn_blk = (int*) calloc ((np+1), sizeof(int));
#ifdef USE_MPI
  MPI_Allgather (&mesh->n_dmn, 1, MPI_INT, &mesh->dmn_blk[1], 1, MPI_INT,
    MPI_COMM_WORLD);
#else
  mesh->dmn_blk[1] = mesh->n_dmn;
#endif

  mesh->dmn_blk[0] = 0;
  for (i=0; i<np; ++i) {
    mesh->dmn_blk[i+1] += mesh->dmn_blk[i];
  }
  offset_to_count (np, mesh->dmn_blk, &mesh->dmn_cnt);


  // construct CSR boundary matrix
  min_fce_vrt = mesh->n_dim - 1;


#ifdef USE_MPI
  ParMETIS_V3_Mesh2Dual (
    mesh->bnd_blk, mesh->bnd_vrt_ptr, mesh->bnd_vrt_lst,
    &idx_offset, &min_fce_vrt,
    &mesh->bnd_neig_ptr, &mesh->bnd_neig_lst,
    &mpi_comm_world);
#else
  mesh_to_dual (mesh->bnd_blk[np], mesh->bnd_vrt_ptr, mesh->bnd_vrt_lst,
                idx_offset, min_fce_vrt,
                &mesh->bnd_neig_ptr, &mesh->bnd_neig_lst);
#endif


  // construct CSR domain matrix
  min_fce_vrt = mesh->n_dim;


#ifdef USE_MPI
  ParMETIS_V3_Mesh2Dual (
    mesh->dmn_blk, mesh->dmn_vrt_ptr, mesh->dmn_vrt_lst,
    &idx_offset, &min_fce_vrt,
    &mesh->dmn_neig_ptr, &mesh->dmn_neig_lst,
    &mpi_comm_world);
#else
  mesh_to_dual (mesh->dmn_blk[np], mesh->dmn_vrt_ptr, mesh->dmn_vrt_lst,
                idx_offset, min_fce_vrt,
                &mesh->dmn_neig_ptr, &mesh->dmn_neig_lst);
#endif


  // add diagonal entries to matrices
  n = mesh->bnd_blk[p];

  add_pole_entry (mesh->bnd_blk[p], mesh->n_bnd,
    &mesh->bnd_neig_ptr, &mesh->bnd_neig_lst);

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
        std::cout << "Error: split connectivity into boundary and domain" << std::endl;
        return 0;
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



  // if (false)
  {
    mshblk_t bnd;
    mshblk_t dmn;

    bnd.blk = p;
    bnd.n_blk = np;
    bnd.elm_blk = mesh->bnd_blk;
    bnd.n_elm = mesh->n_bnd;
    bnd.elm_idx = mesh->bnd_idx;
    bnd.elm_tag.row = mesh->bnd_tag_ptr;
    bnd.elm_tag.col = mesh->bnd_tag_lst;
    bnd.elm_vrt.row = mesh->bnd_vrt_ptr;
    bnd.elm_vrt.col = mesh->bnd_vrt_lst;
    bnd.elm_neig.row = mesh->bnd_neig_ptr;
    bnd.elm_neig.col = mesh->bnd_neig_lst;

    bnd.n_ghst = mesh->n_gbnd;
    bnd.ghst_idx = mesh->gbnd_idx;
    bnd.ghst_tag.row = mesh->gbnd_tag_ptr;
    bnd.ghst_tag.col = mesh->gbnd_tag_lst;
    bnd.ghst_vrt.row = mesh->gbnd_vrt_ptr;
    bnd.ghst_vrt.col = mesh->gbnd_vrt_lst;

    bnd.n_fce = mesh->n_xbnd;
    bnd.fce_elm.row = mesh->xbnd_elm_ptr;
    bnd.fce_elm.col = mesh->xbnd_elm_lst;
    bnd.fce_vrt.row = mesh->xbnd_vrt_ptr;
    bnd.fce_vrt.col = mesh->xbnd_vrt_lst;

    dmn.blk = p;
    dmn.n_blk = np;
    dmn.elm_blk = mesh->dmn_blk;
    dmn.n_elm = mesh->n_dmn;
    dmn.elm_idx = mesh->dmn_idx;
    dmn.elm_tag.row = mesh->dmn_tag_ptr;
    dmn.elm_tag.col = mesh->dmn_tag_lst;
    dmn.elm_vrt.row = mesh->dmn_vrt_ptr;
    dmn.elm_vrt.col = mesh->dmn_vrt_lst;
    dmn.elm_neig.row = mesh->dmn_neig_ptr;
    dmn.elm_neig.col = mesh->dmn_neig_lst;

    dmn.n_ghst = mesh->n_gdmn;
    dmn.ghst_idx = mesh->gdmn_idx;
    dmn.ghst_tag.row = mesh->gdmn_tag_ptr;
    dmn.ghst_tag.col = mesh->gdmn_tag_lst;
    dmn.ghst_vrt.row = mesh->gdmn_vrt_ptr;
    dmn.ghst_vrt.col = mesh->gdmn_vrt_lst;

    dmn.n_fce = mesh->n_xdmn;
    dmn.fce_elm.row = mesh->xdmn_elm_ptr;
    dmn.fce_elm.col = mesh->xdmn_elm_lst;
    dmn.fce_vrt.row = mesh->xdmn_vrt_ptr;
    dmn.fce_vrt.col = mesh->xdmn_vrt_lst;

    dmn.n_dim = mesh->n_dim;
    dmn.n_vrt = mesh->n_vrt;
    dmn.vrt_idx = mesh->vrt_idx;


    // offset domain block
    dmn_blk_os = bnd.elm_blk[bnd.n_blk];
    for (i = 0; i != (dmn.n_blk+1); ++i) {
      dmn.elm_blk[i] += dmn_blk_os;
    }

    // setup vertex index mapping
    dmn.vi2i.reserve(dmn.n_vrt);
    for (i = 0; i != dmn.n_vrt; ++i) {
      dmn.vi2i[dmn.vrt_idx[i]] = i;
    }

    // setup ghost index and block mappings
    bnd.gi2i.reserve(bnd.n_ghst);
    bnd.gi2b.reserve(bnd.n_ghst);
    for (i = 0; i != bnd.n_ghst; ++i) {
      jjn = bnd.ghst_tag.row[i+1] - 1;
      auto oblk = (-1) * bnd.ghst_tag.col[jjn];
      bnd.gi2b[bnd.ghst_idx[i]] = oblk;
      bnd.gi2i[bnd.ghst_idx[i]] = i;
    }

    dmn.gi2i.reserve(dmn.n_ghst);
    dmn.gi2b.reserve(dmn.n_ghst);
    for (i = 0; i != dmn.n_ghst; ++i) {
      jjn = dmn.ghst_tag.row[i+1] - 1;
      auto oblk = (-1) * dmn.ghst_tag.col[jjn];
      dmn.gi2b[dmn.ghst_idx[i]] = oblk;
      dmn.gi2i[dmn.ghst_idx[i]] = i;
    }

#ifdef USE_MPI
    // init_ghost_update(bnd, dmn);
#endif

    geom_t xbnd;
    geom_t xdmn;

    xdmn.x_vrt.resize(dmn.n_vrt);
    for (i = 0; i != dmn.n_vrt; ++i) {
      vector_t<3, double> xi{mesh->x_vrt[i*3+0], mesh->x_vrt[i*3+1], mesh->x_vrt[i*3+2]};
      xdmn.x_vrt[i] = xi;
    }

    // initialize domain geometry
    xdmn.x_vc.dmn.resize(dmn.n_elm);
    xdmn.x_vc.ldmn.resize(dmn.n_locl);
    xdmn.x_vc.gdmn.resize(dmn.n_ghst);

    xdmn.vol.dmn.resize(dmn.n_elm);
    xdmn.vol.ldmn.resize(dmn.n_locl);
    xdmn.vol.gdmn.resize(dmn.n_ghst);

    xdmn.gl_min.dmn.resize(dmn.n_elm);
    xdmn.gl_min.ldmn.resize(dmn.n_locl);
    xdmn.gl_min.gdmn.resize(dmn.n_ghst);
    xdmn.gl_max.dmn.resize(dmn.n_elm);
    xdmn.gl_max.ldmn.resize(dmn.n_locl);
    xdmn.gl_max.gdmn.resize(dmn.n_ghst);
    xdmn.gl.dmn.resize(dmn.n_elm);
    xdmn.gl.ldmn.resize(dmn.n_locl);
    xdmn.gl.gdmn.resize(dmn.n_ghst);

    xdmn.ls.dmn.resize(dmn.n_elm);
    xdmn.ls.ldmn.resize(dmn.n_locl);
    xdmn.ls.gdmn.resize(dmn.n_ghst);
    xdmn.gg.dmn.resize(dmn.n_elm);
    xdmn.gg.ldmn.resize(dmn.n_locl);
    xdmn.gg.gdmn.resize(dmn.n_ghst);

    xdmn.x_fc.resize(dmn.n_fce);
    xdmn.x_pc.resize(dmn.n_fce);
    xdmn.x_nc.resize(dmn.n_fce);
    xdmn.norm.resize(dmn.n_fce);
    xdmn.volr.resize(dmn.n_elm);
    xdmn.area.resize(dmn.n_fce);
    xdmn.l_pn.resize(dmn.n_fce);
    xdmn.w_fp.resize(dmn.n_fce);

    // print_range("bnd.elm_blk [cpp]", bnd.n_blk, bnd.elm_blk);
    // print_range("dmn.elm_blk [cpp]", dmn.n_blk, dmn.elm_blk);
    // print_graph("dmn.fce_elm [cpp]", dmn.n_fce, dmn.fce_elm.row, dmn.fce_elm.col);

    geometry(bnd, dmn, xbnd, xdmn);

    if (p == 0) std::cout << "Volume " << xdmn.vol_sum << std::endl;
  }


  // pass control to CFD solver
  n = jsonfile.length();

  // re-adjust domain offets
  for (i = 0; i != (np+1); ++i) {
    mesh->dmn_blk[i] -= dmn_blk_os;
  }

  navierstokes_ (
    &n, jsonfile.c_str(),
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

  machine_t::finalize();
}
