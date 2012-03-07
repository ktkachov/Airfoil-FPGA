#ifndef AIR_UTILS_H
#define AIR_UTILS_H

#include <stdint.h>
#include <limits.h>

#define EMPTY_ENTRY UINT_MAX

/*TODO: Generalise for any value type*/

/*Semi dynamic array, if this was C++, I would have used an std::vector instead*/
typedef struct array_struct {
  uint32_t* arr;
  uint32_t len;
  uint32_t maxLen;
} arr_t;

/*Helper functions to manipulate arr_t structures*/
void addToArr(arr_t* a, uint32_t e) {
  if (a->len == a->maxLen - 1) {
    a->maxLen += 64;
    a->arr = realloc(a->arr, a->maxLen * sizeof(*(a->arr)));
  }
  a->arr[a->len++] = e;
}

void initArr(arr_t* a, uint32_t size) {
  a->len = 0;
  a->maxLen = size;
  a->arr = malloc(a->maxLen * sizeof(*(a->arr)));
}

void destroyArr(arr_t* a) {
  free(a->arr);
  free(a);
}

/*Hash map*/
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

arr_t* hash_map_key_set(hash_map* m) {
  arr_t* res = malloc(sizeof(*res));
  initArr(res, m->size);
  for (uint32_t i = 0; i < m->size; ++i) {
    struct entry *p = &m->arr[i];
    while (p != NULL && p->k != EMPTY_ENTRY) {
      addToArr(res, p->k);
      p = p->next;
    }
  }
  res->maxLen = res->len;
  res->arr = realloc(res->arr, res->len * sizeof(*res->arr));
  return res;
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
  free(m);
}

typedef hash_map hash_set;
hash_set* createHashSet(uint32_t);


int addToHashSet(hash_set* m, uint32_t k) {
  uint32_t hk = k % m->size;
  uint32_t v = 0;
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

hash_set* setDiff(hash_set* s1, hash_set* s2) {
  hash_set* res = createHashSet(s1->size);
  for (uint32_t i = 0; i < s1->size; ++i) {
    struct entry* e = &s1->arr[i];
    while (e->next != NULL && e->k != EMPTY_ENTRY) {
      if (!contains(s2, e->k)) {
        addToHashSet(res, e->k);
      }
      e = e->next;
    }
  }
  return res;
}

hash_set* setIntersection(hash_set* s1, hash_set* s2) {
  hash_set* res = createHashSet(s1->size);
  for (uint32_t i = 0; i < s1->size; ++i) {
    struct entry* e = &s1->arr[i];
    while (e->next != NULL && e->k != EMPTY_ENTRY) {
      if (contains(s2, e->k)) {
        addToHashSet(res, e->k);
      }
      e = e->next;
    }
  }
  return res;
}

hash_set* setUnion(hash_set* s1, hash_set* s2) {
  hash_set* res = createHashSet(s1->size);
  for (uint32_t i = 0; i < s1->size; ++i) {
    struct entry* e = &s1->arr[i];
    while (e->next != NULL && e->k != EMPTY_ENTRY) {
      addToHashSet(res, e->k);
      e = e->next;
    }
  }
  for (uint32_t i = 0; i < s2->size; ++i) {
    struct entry* e = &s2->arr[i];
    while (e->next != NULL && e->k != EMPTY_ENTRY) {
      addToHashSet(res, e->k);
      e = e->next;
    }
  }
  return res;
}

arr_t* toArr(hash_set* s) {
  return hash_map_key_set(s);
}

hash_set* fromArr(arr_t* a) {
  hash_set* res = createHashSet(a->len / 2);
  for (uint32_t i = 0; i < a->len; ++i) {
    addToHashSet(res, a->arr[i]);
  }
  return res;
}

hash_set* createHashSet(uint32_t size) {
  return createHashMap(size);
}

void destroyHashSet(hash_set* s) {
  destroyHashMap(s);
}


#endif
