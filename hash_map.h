#ifndef HASH_MAP_H
#define HASH_MAP_H

#define EMPTY_ENTRY -1

/*TODO: Generalise for any value type*/
struct entry {
  int k, v;
  struct entry* next;
};

typedef struct hash_map_struct {
  int size;
  struct entry* arr;
} hash_map;

int hash(int, hash_map*);

hash_map* createHashMap(int size) {
  hash_map* res = malloc(sizeof(*res));
  res->size = size;
  res->arr = malloc(size * sizeof(*res->arr));
  for (int i = 0; i < size; ++i) {
    res->arr[i].k = EMPTY_ENTRY;
    res->arr[i].v = EMPTY_ENTRY;
    res->arr[i].next = NULL;
  }
  return res;
}

inline int hash(int k, hash_map* m) {
  return k % m->size;
}

void addToHashMap(hash_map* m, int k, int v) {
  int hk = hash(k, m);
  if (m->arr[hk].k == EMPTY_ENTRY) {
    m->arr[hk].k = k;
    m->arr[hk].v = v;
    return;
  }
  struct  entry* e = &(m->arr[hk]);
  while (e->next != NULL) {
    if (e->k == k) {
      return;
    }
    e = e->next;
  }
  if (e->k == k) {
    return;
  }
  struct  entry* ne = malloc(sizeof(*ne));
  ne->k = k;
  ne->v = v;
  ne->next = NULL;
  e->next = ne;
} 

int contains(hash_map* m, int k) {
  int hk = hash(k, m);
  if (m->arr[hk].k == EMPTY_ENTRY) {
    return 0;
  }
  struct  entry* e = &(m->arr[hk]);
  while (e != NULL) {
    if (e->k == k) {
      return 1;
    }
    e = e->next;
  }
  return 0;
}

int getValue(hash_map* m, int k) {
  int hk = hash(k, m);
  if (m->arr[hk].k == k) {
    return m->arr[hk].v;
  }
  if (m->arr[hk].k == EMPTY_ENTRY) {
    return EMPTY_ENTRY;
  }
  struct  entry* e = &(m->arr[hk]);
  while (e != NULL) {
    if (e->k == k) {
      return e->v;
    }
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
  for (int i = 0; i < m->size; ++i) {
     freeEntryList(m->arr[i].next);
  }
  free(m->arr);
}

#endif
