#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mshf.h"


void read_mesh_format (struct mshf_t* mshf, char* filename)
{
  int ierr = 0;
  FILE* fp;
  char buf[256];
  int i, j, k, n;


  fp = fopen (filename, "r");

  if (fp == NULL) {
    fprintf (stderr, "[E] Cannot open input file.\n");
    exit(1);
  }


  while (fscanf (fp, "%s", &buf[0]) == 1) {


    if (strcmp (buf, "Element.Types") == 0) {
      ierr = fscanf (fp, "%d", &mshf->n_elm);

      mshf->elm_idx = (int*) malloc ((20) * sizeof(int)); // magic number 20
      mshf->elm_ty = (int*) malloc (mshf->n_elm * sizeof(int));

      for (i=0; i<mshf->n_elm; ++i) {
	ierr = fscanf (fp, "%d", &j);
	mshf->elm_ty[i] = j;
	mshf->elm_idx[j] = i;
      }
    }


    else if (strcmp (buf, "Element.Dimensions") == 0) {
      ierr = fscanf (fp, "%d", &mshf->n_elm);

      mshf->elm_dim = (int*) malloc (mshf->n_elm * sizeof(int));

      for (i=0; i<mshf->n_elm; ++i) {
	ierr = fscanf (fp, "%d", &mshf->elm_dim[i]);
      }
    }


    else if (strcmp (buf, "Element.Face_Offsets") == 0) {
      ierr = fscanf (fp, "%d", &mshf->n_elm);
      --mshf->n_elm;

      mshf->elm_fce_ptr = (int*) malloc ((mshf->n_elm+1) * sizeof(int));

      for (i=0; i<mshf->n_elm+1; ++i) {
	ierr = fscanf (fp, "%d", &mshf->elm_fce_ptr[i]);
      }
    }


    else if (strcmp (buf, "Face.Vertex_Offsets") == 0) {
      ierr = fscanf (fp, "%d", &n);

      mshf->fce_vrt_ptr = (int*) malloc (n * sizeof(int));

      for (i=0; i<n; ++i) {
	ierr = fscanf (fp, "%d", &mshf->fce_vrt_ptr[i]);
      }
    }


    else if (strcmp (buf, "Face.Vertex_Lists") == 0) {
      ierr = fscanf (fp, "%d", &n);

      mshf->fce_vrt_lst = (int*) malloc (n * sizeof(int));

      for (i=0; i<n; ++i) {
	ierr = fscanf (fp, "%d", &mshf->fce_vrt_lst[i]);
      }
    }
  }


  fclose (fp);
}
