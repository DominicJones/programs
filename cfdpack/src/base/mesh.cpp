#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mesh.h"


float gmsh_ver[] = {2.2, 0., 8.};
const int type_to_nvrt[] = {0,2,3,4,4,8,6,5};


void distributed_block (int n, int np, int** blk)
{
  int i, m, *pblk;
  pblk = *blk = (int*) malloc ((np+1) * sizeof(int));
  pblk[0] = 0;
  m = n / np;
  for (i=1; i<np; ++i) {
    pblk[i] = pblk[i-1] + m;
  }
  pblk[np] = n;
}


void partitioned_block (int n, int np, int* idx, int** blk)
{
  int i, j, *pblk;
  pblk = *blk = (int*) calloc ((np+1), sizeof(int));
  for (i=0; i<n; ++i) {
    j = idx[i];
    ++pblk[j+1];
  }
  for (i=1; i<np+1; ++i) {
    pblk[i] += pblk[i-1];
  }
}


void read_gmsh_mesh (struct mesh_t* mesh, char* filename, char section, char read_mode)
{
  int ierr = 0;
  FILE* fp;
  char buf[256];
  int i, j, k, l, m, n, p, np, idx, tag[12], n_tag, vrt[8], n_vrt;
  int jj0, jjn, jj, k1, k2, k3, *blk;
  double x_vrt[3];

  // open mesh file
  fp = fopen (filename, "r");

  if (fp == NULL) {
    fprintf (stderr, "Cannot open input file \"%s\".\n", filename);
    exit(1);
  }

  // read mesh file
  while (fscanf (fp, "%s", &buf[0]) == 1) {

    // Vertex data
    if ((section == 'V' || section == 'A') &&
      strcmp (buf, "$Nodes") == 0) {

      ierr = fscanf (fp, "%d", &n);

      switch (read_mode) {
      case 's':
        mesh->N_vrt = mesh->n_vrt = m = n;
        mesh->vrt_idx = (int*) malloc ((n) * sizeof(int));
        mesh->x_vrt = (double*) malloc ((n*3) * sizeof(double));

        for (i=0; i<n; ++i) {
          ierr = fscanf (fp, "%d", &mesh->vrt_idx[i]);
          ierr = fscanf (fp, "%lg", &mesh->x_vrt[i*3+0]);
          ierr = fscanf (fp, "%lg", &mesh->x_vrt[i*3+1]);
          ierr = fscanf (fp, "%lg", &mesh->x_vrt[i*3+2]);
        }
        break;
      case 'p':
        mesh->N_vrt = n;
        m = mesh->n_vrt;
        mesh->vrt_idx = (int*) malloc ((m) * sizeof(int));
        mesh->x_vrt = (double*) malloc ((m*3) * sizeof(double));

        j = 0;
        for (i=0; i<n; ++i) {
          ierr = fscanf (fp, "%d", &idx);
          ierr = fscanf (fp, "%lg", &x_vrt[0]);
          ierr = fscanf (fp, "%lg", &x_vrt[1]);
          ierr = fscanf (fp, "%lg", &x_vrt[2]);

          if (mesh->vrt_part[i] == 1) {
            mesh->vrt_idx[j] = idx;
            mesh->x_vrt[j*3+0] = x_vrt[0];
            mesh->x_vrt[j*3+1] = x_vrt[1];
            mesh->x_vrt[j*3+2] = x_vrt[2];
            ++j;
          }
        }
        break;
      }
    }

    // Element data
    else if ((section == 'E' || section == 'A') &&
      strcmp (buf, "$Elements") == 0) {

      ierr = fscanf (fp, "%d", &n);

      p = mesh->part;
      np = mesh->n_part;

      switch (read_mode) {
      case 's':
        mesh->N_elm = mesh->n_elm = m = n;
        break;
      case 'd':
        distributed_block (n, np, &blk);
        mesh->N_elm = n;
        mesh->n_elm = m = blk[p+1] - blk[p];
        break;
      case 'p':
        partitioned_block (n, np, mesh->elm_part, &blk);
        mesh->N_elm = n;
        mesh->n_elm = m = blk[p+1] - blk[p];
        break;
      }

      mesh->elm_idx = (int*) malloc (m * sizeof(int));
      mesh->elm_tag_ptr = (int*) malloc ((m+1) * sizeof(int));
      mesh->elm_tag_lst = (int*) malloc ((m*8) * sizeof(int));
      mesh->elm_vrt_ptr = (int*) malloc ((m+1) * sizeof(int));
      mesh->elm_vrt_lst = (int*) malloc ((m*8) * sizeof(int));

      mesh->elm_tag_ptr[0] = 0;
      mesh->elm_vrt_ptr[0] = 0;

      k1 = 0; k2 = 0; k3 = 0;
      for (i=0; i<n; ++i) {
        ierr = fscanf (fp, "%d", &idx);

        ierr = fscanf (fp, "%d", &tag[0]);
        ierr = fscanf (fp, "%d", &n_tag);

        n_tag = n_tag + 1;
        n_vrt = type_to_nvrt[tag[0]];

        for (j=1; j<n_tag; ++j) {
          ierr = fscanf (fp, "%d", &tag[j]);
        }
        if (n_tag >= 4 && tag[3] < 1) n_tag = 3;

        for (j=0; j<n_vrt; ++j) {
          ierr = fscanf (fp, "%d", &vrt[j]);
        }

        if (read_mode == 's' ||
          (read_mode == 'd' && i >= blk[p] && i < blk[p+1]) ||
          (read_mode == 'p' && mesh->elm_part[i] == p)) {

          mesh->elm_idx[k1] = idx;
          mesh->elm_tag_ptr[k1+1] = n_tag;
          mesh->elm_vrt_ptr[k1+1] = n_vrt;
          ++k1;

          for (j=0; j<n_tag; ++j) {
            mesh->elm_tag_lst[k2++] = tag[j];
          }

          for (j=0; j<n_vrt; ++j) {
            mesh->elm_vrt_lst[k3++] = vrt[j];
          }
        }
      }

      for (i=1; i<m+1; ++i) {
        mesh->elm_tag_ptr[i] += mesh->elm_tag_ptr[i-1];
        mesh->elm_vrt_ptr[i] += mesh->elm_vrt_ptr[i-1];
      }

      switch (read_mode) {
      case 'd': case 'p':
        free (blk);
      }
    }
  }

  fclose (fp);
}


void write_gmsh_mesh (struct mesh_t* mesh, char* filename, char region)
{
  FILE* fp;
  char buf[256];
  int i, jj0, jjn, jj, j, kk0, kkn, kk, k, n;
  int p, np, dim, ph_ty, n_tag, part;
  int n_elm, *elm_idx, *elm_tag_ptr, *elm_tag_lst, *elm_vrt_ptr, *elm_vrt_lst;

  switch (region) {
  case 'e':
    n_elm = mesh->n_elm;
    elm_idx = mesh->elm_idx;
    elm_tag_ptr = mesh->elm_tag_ptr;
    elm_tag_lst = mesh->elm_tag_lst;
    elm_vrt_ptr = mesh->elm_vrt_ptr;
    elm_vrt_lst = mesh->elm_vrt_lst;
    break;
  case 'b':
    n_elm = mesh->n_bnd;
    elm_idx = mesh->bnd_idx;
    elm_tag_ptr = mesh->bnd_tag_ptr;
    elm_tag_lst = mesh->bnd_tag_lst;
    elm_vrt_ptr = mesh->bnd_vrt_ptr;
    elm_vrt_lst = mesh->bnd_vrt_lst;
    break;
  case 'd':
    n_elm = mesh->n_dmn;
    elm_idx = mesh->dmn_idx;
    elm_tag_ptr = mesh->dmn_tag_ptr;
    elm_tag_lst = mesh->dmn_tag_lst;
    elm_vrt_ptr = mesh->dmn_vrt_ptr;
    elm_vrt_lst = mesh->dmn_vrt_lst;
    break;
  }

  p = mesh->part;
  np = mesh->n_part;

  part = 0;
  if (np > 1 && (p >= 0 && p < np) && region == 'e') part = 1;

  // open mesh file
  fp = fopen (filename, "w");

  if (fp == NULL) {
    fprintf (stderr, "Cannot create output file \"%s\".\n", filename);
    exit(1);
  }

  // Header data
  fprintf (fp, "$MeshFormat\n");
  fprintf (fp, "%g %g %g\n", gmsh_ver[0], gmsh_ver[1], gmsh_ver[2]);
  fprintf (fp, "$EndMeshFormat\n");

  // Vertex data
  fprintf (fp, "$Nodes\n%d\n", mesh->n_vrt);
  for (i=0; i<mesh->n_vrt; ++i) {
    fprintf (fp, "%d %.10lg %.10lg %.10lg\n",
      mesh->vrt_idx[i],
      mesh->x_vrt[i*3+0],
      mesh->x_vrt[i*3+1],
      mesh->x_vrt[i*3+2]);
  }
  fprintf (fp, "$EndNodes\n");

  // Element data
  fprintf (fp, "$Elements\n%d\n", n_elm);

  for (i=0; i<n_elm; ++i) {
    jj0 = elm_tag_ptr[i];
    jjn = elm_tag_ptr[i+1];

    kk0 = elm_vrt_ptr[i];
    kkn = elm_vrt_ptr[i+1];

    n_tag = jjn - jj0 - 1;
    if (part) n_tag += 2;


    fprintf (fp, "%d", elm_idx[i]);


    fprintf (fp, " %d", elm_tag_lst[jj0]);
    fprintf (fp, " %d", n_tag);

    for (jj=jj0+1; jj<jjn; ++jj) {
      fprintf (fp, " %d", elm_tag_lst[jj]);
    }

    if (part) {
      fprintf (fp, " %d", 1);
      fprintf (fp, " %d", mesh->elm_part[i]);
    }

    for (kk=kk0; kk<kkn; ++kk) {
      fprintf (fp, " %d", elm_vrt_lst[kk]);
    }

    fprintf (fp, "\n");
  }

  fprintf (fp, "$EndElements\n");
  fclose (fp);
}


extern "C" {
int write_gmsh_field_ (int* pflen, char* fname, int* pn, int* elm_idx, double* phi)
{
  int err = 0;
  char filename[256];
  FILE *fp;
  int i, n_phi, flen;

  flen = *pflen;
  n_phi = *pn;
  strncpy (filename, fname, flen);
  filename[flen] = '\0';

  fp = fopen (filename, "w");

  if (fp == NULL) {
    fprintf (stderr, "Cannot create output file \"%s\".\n", filename);
    err = 1; goto quit;
  }

  printf ("Writing %s\n", filename);

  // Header data
  fprintf (fp, "$MeshFormat\n");
  fprintf (fp, "%g %g %g\n", gmsh_ver[0], gmsh_ver[1], gmsh_ver[2]);
  fprintf (fp, "$EndMeshFormat\n");

  // Element data
  fprintf (fp, "$ElementData\n");

  // string tags: name of view
  fprintf (fp, "%d\n", 1);
  fprintf (fp, "%s\n", filename);

  // real tags: time
  fprintf (fp, "%d\n", 1);
  fprintf (fp, "%g\n", 0.);

  // int tags: time step index, number of field components,
  //   number of elements, partition index
  fprintf (fp, "%d\n", 4);
  fprintf (fp, "%d\n%d\n%d\n%d\n", 0, 1, n_phi, 0);

  // element centred data
  for (i=0; i<n_phi; ++i) {
    fprintf (fp, "%d %lg\n", elm_idx[i], phi[i]);
  }

  fprintf (fp, "$EndElementData\n");
  fclose (fp);

  quit:
  return (err);
}
}
