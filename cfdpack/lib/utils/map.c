#include "map.h"


struct Map *map_handle = NULL;


void map_add (struct Map **map, int key, int value) {
  struct Map *s;
  s = (struct Map*) malloc (sizeof(struct Map));
  s->key = key;
  s->value = value;
  map_handle = *map;
  HASH_ADD_INT (map_handle, key, s);
  *map = map_handle;
}


struct Map *map_find (struct Map **map, int key) {
  struct Map *s;
  map_handle = *map;
  HASH_FIND_INT (map_handle, &key, s);
  *map = map_handle;
  return s;
}


void map_delete (struct Map **map) {
  struct Map *curr, *tmp;
  map_handle = *map;
  HASH_ITER (hh, map_handle, curr, tmp) {
    HASH_DEL (map_handle, curr);
    free (curr);
  }
  *map = map_handle;
}