#ifdef __cplusplus
extern "C" {
#endif

#ifndef MAP_H
#define MAP_H

#include <stdlib.h>
#include "uthash.h"

struct Map {
  int key;
  int value;
  UT_hash_handle hh;
};

void map_add (struct Map **map, int key, int value);

struct Map *map_find (struct Map **map, int key);

void map_delete (struct Map **map);

#endif

#ifdef __cplusplus
}
#endif
