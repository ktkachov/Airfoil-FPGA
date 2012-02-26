#ifndef HASH_MAP_H
#define HASH_MAP_H

#include <stdint.h>
#include <limits.h>

#define EMPTY_ENTRY UINT_MAX

/*TODO: Generalise for any value type*/
struct entry {
  uint32_t k, v;
  struct entry* next;
};

typedef struct hash_map_struct {
  uint32_t size;
  struct entry* arr;
} hash_map;

uint32_t hash(uint32_t);

hash_map* createHashMap(uint32_t size) {
  hash_map* res = malloc(sizeof(*res));
  res->size = size;
  res->arr = malloc(size * sizeof(*res->arr));
  for (uint32_t i = 0; i < size; ++i) {
    res->arr[i].k = EMPTY_ENTRY;
    res->arr[i].v = EMPTY_ENTRY;
    res->arr[i].next = NULL;
  }
  return res;
}

inline uint32_t hash(uint32_t a) {
   a = (a+0x7ed55d16) + (a<<12);
   a = (a^0xc761c23c) ^ (a>>19);
   a = (a+0x165667b1) + (a<<5);
   a = (a+0xd3a2646c) ^ (a<<9);
   a = (a+0xfd7046c5) + (a<<3);
   a = (a^0xb55a4f09) ^ (a>>16);
   return a;
}

int addToHashMap(hash_map* m, uint32_t k, uint32_t v) {
  uint32_t hk = hash(k) % m->size;
  if (m->arr[hk].k == EMPTY_ENTRY) {
    m->arr[hk].k = k;
    m->arr[hk].v = v;
    return 1;
  }
  struct entry* e = &(m->arr[hk]);
  while (e->next != NULL) {
    if (e->k == k) {
      return 0;
    }
    e = e->next;
  }
  if (e->k == k) {
    return 0;
  }
  struct entry* ne = malloc(sizeof(*ne));
  ne->k = k;
  ne->v = v;
  ne->next = NULL;
  e->next = ne;
  return 1;
} 

uint32_t contains(hash_map* m, uint32_t k) {
  uint32_t hk = hash(k) % m->size;
  if (m->arr[hk].k == EMPTY_ENTRY) {
    return 0;
  }
  struct entry* e = &(m->arr[hk]);
  while (e != NULL) {
    if (e->k == k) {
      return 1;
    }
    e = e->next;
  }
  return 0;
}

uint32_t getValue(hash_map* m, uint32_t k) {
  uint32_t hk = hash(k) % m->size;
  if (m->arr[hk].k == k) {
    return m->arr[hk].v;
  }
  if (m->arr[hk].k == EMPTY_ENTRY) {
    return EMPTY_ENTRY;
  }
  struct entry* e = &(m->arr[hk]);
  while (e != NULL) {
    if (e->k == k) {
      return e->v;
    }
    e = e->next;
  }
  return EMPTY_ENTRY;
}

void freeEntryList(struct entry* e) {
  if (e == NULL) {
    return;
  }
  struct entry* n = e->next;
  free(e);
  freeEntryList(n);
}

void destroyHashMap(hash_map* m) {
  for (uint32_t i = 0; i < m->size; ++i) {
     freeEntryList(m->arr[i].next);
  }
  free(m->arr);
}

#endif
