/*
 * Copyright (C) 2008-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */
  
#ifndef HEAPC_H_
#define HEAPC_H_
  
#include "fastCMaths.h"

//#include "macros.h"
  
//BEGIN_DECL

//#define heap_empty(h) (!(h)->num_entries)

struct heap {
        uint32_t max_entries;
        uint32_t num_entries;
        int64_t *entries; 
};

extern int32_t heap_empty(struct heap *heap);
extern struct heap *heap_create(uint32_t max_entries);
extern void heap_destroy(struct heap *heap);
extern void heapify(struct heap *heap, int64_t i);
extern int32_t heap_expand(struct heap *heap);
extern int32_t heap_insert(struct heap *heap, int64_t entry);
extern int64_t heap_extract(struct heap *heap);
extern int64_t heap_peek(struct heap *heap); 
extern void heap_clean(struct heap *heap);

#define HEADROOM 100

//END_DECL
 
#endif
