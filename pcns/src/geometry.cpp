#include "geometry.h"
#include "mesh.h"
#include "mshf.h"
#include "graph.h"
#include "vector.h"
#include "tensor.h"
#include "machine.h"

#include <iostream>


void geometry(mshblk_t const &bnd, mshblk_t const &dmn, geom_t &xbnd, geom_t &xdmn)
{
  int n_dim,bc,bnd1,bndl,bndu,dmn1,dmnl,dmnu;
  int i,jj1,jjn,jj,j,j1,jn,jb,jv,kk1,kkn,skk,kk,k,k1,k2,kv,k1v;
  int j11,j1l,jn1,jnl,jb1,jbl;
  double fvol,area1,l_pn1,area_sum,vol_sum;
  std::vector<double> bc_area;
  vector_t<3, double> edge1,edge,cprod,unit_z,vec0;
  vector_t<3, double> dx,x1,x2,x_fc1,x_pc1,x_nc1,norm1;

  n_dim = dmn.n_dim;
  unit_z = {double(0), double(0), double(1)};
  vec0 = {double(0), double(0), double(0)};

  bnd1 = bnd.elm_blk[0];
  bndl = bnd.elm_blk[bnd.blk];
  bndu = bnd.elm_blk[bnd.blk+1];

  dmn1 = dmn.elm_blk[0];
  dmnl = dmn.elm_blk[dmn.blk];
  dmnu = dmn.elm_blk[dmn.blk+1];

  for (i = 0; i != dmn.n_elm; ++i) {
    jj1 = dmn.elm_vrt.row[i];
    jjn = dmn.elm_vrt.row[i+1];

    x1 = vec0;
    for (jj = jj1; jj != jjn; ++jj) {
      j = dmn.elm_vrt.col[jj];
      jv = dmn.vi2i.find(j)->second;
      x1 = x1 + xdmn.x_vrt[jv];
    }

    x1 = x1 / double(jjn - (jj1 + 1));
    xdmn.x_vc.dmn[i] = x1;
    xdmn.vol.dmn[i] = 0;
  }

#ifdef USE_MPI
  // ghost_update_v(bnd, dmn, xdmn.x_vc, 0, 1);
#endif

  bc_area.resize(100);
  std::fill(bc_area.begin(), bc_area.end(), 0);

  area_sum = 0;
  vol_sum = 0;

  for (i = 0; i != dmn.n_fce; ++i) {
    jj1 = dmn.fce_elm.row[i];
    jjn = dmn.fce_elm.row[i+1];

    j1 = dmn.fce_elm.col[jj1];
    j11 = (j1 - dmn1);
    j1l = (j1 - dmnl);
    j1 = (j1 - dmnl);
    j1 = j1 - 1; // [0..n) indices

    jn = dmn.fce_elm.col[jjn-1];

    kk1 = dmn.fce_vrt.row[i];
    kkn = dmn.fce_vrt.row[i+1];

    // face centre
    x1 = vec0;
    for (kk = kk1; kk != kkn; ++kk) {
      k = dmn.fce_vrt.col[kk];
      kv = dmn.vi2i.find(k)->second;
      x1 = x1 + xdmn.x_vrt[kv];
    }

    x1 = x1 / double(kkn - kk1);
    x_fc1 = x1;
    xdmn.x_fc[i] = x1;

    // face normal and area
    k1 = dmn.fce_vrt.col[kk1];
    k1v = dmn.vi2i.find(k1)->second;
    x1 = vec0;

    for (kk = kk1 + 1; kk != kkn; ++kk) {
      k = dmn.fce_vrt.col[kk];
      kv = dmn.vi2i.find(k)->second;
      edge = vec0;
      edge = xdmn.x_vrt[kv] - xdmn.x_vrt[k1v];

      auto const kkn_case{kkn - kk1};
      if (kkn_case == 2) {
        cprod = cross(edge, unit_z);
        x1 = cprod;
      }
      else if (kkn_case > 2) {
        auto const kk_case{kk - (kk1+1)};
        if (kk_case > 0) {
          cprod = cross(edge1, edge);
          x1 = x1 + cprod / double(2);
        }
        edge1 = edge;
      }
    }


    area1 = mag(x1);
    norm1 = x1 / area1;
    xdmn.area[i] = area1;
    xdmn.norm[i] = norm1;


    // volume
    auto const ivol_cmp = 0;
    fvol = x_fc1[ivol_cmp] * norm1[ivol_cmp] * area1;
    xdmn.vol.dmn[j1] += fvol;


    // aux. pole position
    dx = x_fc1 - xdmn.x_vc.dmn[j1];
    x_pc1 = x_fc1 - norm1 * dot(dx, norm1);
    xdmn.x_pc[i] = x_pc1;

    // internal and internal-ghost faces:
    if (jn >= 1) {

      // internal face:
      if (jn >= (dmnl+1) && jn < (dmnu+1)) {
        jb = (jn - (dmnl+1)) + 1;
        jb = jb - 1; // [0..n) indices

        // volume
        xdmn.vol.dmn[jb] -= fvol;

        // aux. neig. position
        dx = x_fc1 - xdmn.x_vc.dmn[jb];
        x_nc1 = x_fc1 - norm1 * dot(dx, norm1);
        xdmn.x_nc[i] = x_nc1;

        // pole to neig. length
        dx = x_nc1 - x_pc1;
        l_pn1 = mag(dx);

        // interpolation factor
        dx = x_fc1 - x_pc1;
        xdmn.w_fp[i] = mag(dx) / l_pn1;

      // internal-ghost face:
      } else {
        jb = dmn.gi2i.find(jn)->second;
        jb = jb - 1; // [0..n) indices

        // aux. neig. position
        dx = x_fc1 - xdmn.x_vc.gdmn[jb];
        x_nc1 = x_fc1 - norm1 * dot(dx, norm1);
        xdmn.x_nc[i] = x_nc1;

        // pole to neig. length
        dx = x_nc1 - x_pc1;
        l_pn1 = mag(dx);

        // interpolation factor
        xdmn.w_fp[i] = double(0.5);
      }

    // boundary faces:
    } else {
      jb = ((-1)*jn - (bndl+1)) + 1;
      jb = jb - 1; // [0..n) indices

      // checks
      bc = bnd.elm_tag.row[jb];
      bc = bnd.elm_tag.col[bc+1];
      bc_area[bc] += area1;

      area_sum += area1;
      vol_sum += fvol;

      // pole to face length
      dx = x_fc1 - x_pc1;
      l_pn1 = mag(dx);
    }

    xdmn.l_pn[i] = l_pn1;
  }


#ifdef USE_MPI
  // ghost_update_s(bnd, dmn, xdmn.vol, 0, 1);
#endif

  xdmn.vol_min = 1.0e+20;
  xdmn.vol_max = 0;
  xdmn.vol_sum = 0;

  for (i = 0; i != dmn.n_elm; ++i) {
    xdmn.volr[i] = 1 / (xdmn.vol.dmn[i] + 1.0e-20);
    xdmn.vol_min = std::min(xdmn.vol_min, xdmn.vol.dmn[i]);
    xdmn.vol_max = std::max(xdmn.vol_max, xdmn.vol.dmn[i]);
    xdmn.vol_sum += xdmn.vol.dmn[i];
  }

  if(xdmn.vol_min < 1.0e-20) {
    std::cout << "invalid volume(s) found" << std::endl;
  }
}
