#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "settings.h"


int main(int argc, char** argv){
  int ierr = 0;
  int result;

  FILE *f;
  Settings *settings;

  char title[256];
  char gd_mthd[256];
  double u_wall[3];

  char sw[2];
  int iw = 1;
  char wall_id[256];


  /* Initialize */
  f = fopen(argv[1], "r");
  if (f == NULL) {
    printf("error: no file specified\n"); ierr = 1;
    goto quit;
  }

  settings = settings_open(f);
  fclose(f);
  if (settings == NULL) {
    printf("error: problem closing the file\n"); ierr = 2;
    goto quit;
  }


  /* Retrieve */
  result = settings_get(settings, "Case", "Title", title, sizeof(title));
  if (result != 0) { printf("%s\n", title); }

  sprintf(sw, "%d", iw);
  printf("[%s] {%d}\n",sw, iw);

  strcpy(wall_id, "Wall ");
  strcat(wall_id, sw);

  result = settings_get_double_tuple(settings, "Boundary Conditions", wall_id, u_wall, sizeof(u_wall));
  if (result != 0) {
    printf("Wall 1: [%g %g %g]\n", u_wall[0], u_wall[1], u_wall[2]);
  }


  result = settings_get(settings, "Gradient", "Method", gd_mthd, sizeof(gd_mthd));
  if (result != 0) { printf("%s\n", gd_mthd); }
  if (strcmp(gd_mthd, "Least Squares") == 0) { printf("High quality mesh required using %s\n", gd_mthd); }


  /* Finalize */
  settings_delete(settings);


 quit:
  return ierr;
}

