#include <stdio.h>
#include "settings.h"


static void iter(const char *key, const char *value, const void *obj)
{
  printf("menu item: %s menu position: %s\n", key, value);
}


int main(){
  FILE *f;
  Settings *settings;
  char buf[255];
  int result;

  f = fopen("settings.txt", "r");
  if (f == NULL) {
    /* Handle error... */
  }
  settings = settings_open(f);
  fclose(f);
  if (settings == NULL) {
    /* Handle read error... */
  }

  /* Insert a new key-value pair */
  settings_set(settings, "Application", "Has Started", "True");

  /* Retrieve a value */
  result = settings_get(settings, "Application", "Version", buf, sizeof(buf));
  if (result == 0) {
    /* Handle value not found... */
  }
  printf("version of this application: %s\n", buf)
    /* Insert a key-value pair that will create a new section */
    settings_set(settings, "Window Properties", "X", "0");

  /* Save the settings object to disk in textual form */
  f = fopen("settings.txt", "w");
  if (f == NULL) {
    /* Handle error... */
  }
  result = settings_save(settings, f);
  fclose(f);
  if (result == 0) {
    /* Handle write error... */
  }


  // reading a tuple
  {
    int argb[4];
    int result;

    result = settings_get_int_tuple(settings, "Foreground",
                                    "Text Color", argb, sizeof(argb));
    if (result == 0) {
      /* Handle value not found... */
    }
    printf("alpha: %d red: %d green: %d blue: %d\n",
           argb[0], argb[1], argb[2], argb[3]);
  }


  // iterate over key-value pair
  settings_enum(settings, "Edit Menu Options", iter, NULL);


  /* When done, destroy the settings object */
  settings_delete(settings);
}

