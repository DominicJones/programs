#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main (int argc, char* argv[])
{
  int err = 0;
  char f1[256], f2[256], f3[256];
  int lf1, lf2, lf3;

  if (argc < 3) {
    printf ("%s: expected <meshfile> <mshffile> <casedir>\n", argv[0]);
    err = 1; goto quit;
  }

  strcpy (f1, argv[1]);
  lf1 = strlen (f1);

  strcpy (f2, argv[2]);
  lf2 = strlen (f2);

  strcpy (f3, argv[3]);
  lf3 = strlen (f3);


  gpde_main_ (&lf1, f1, &lf2, f2, &lf3, f3);


  quit:
  if (err) fprintf (stderr, "%s: error %d\n", argv[0], err);
  return err;
}
