/*
 * Copyright (C) 2008-2011 by Benedict Paten (benedictpaten@gmail.com)
 *
 * Released under the MIT license, see LICENSE.txt
 */
/*
 * This is adapted from GNU code as stated below
 */

#include "heapC.h"
#include "fastCMaths.h"  
#include "commonC.h"
  
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h> 
#include <limits.h>

struct heap *heap_create(uint64_t max_entries)
{
    struct heap *heap;
    uint64_t i;
    
    heap = (struct heap *) st_malloc(sizeof(struct heap));
    if (!heap)
        return 0;
    memset(heap, 0, sizeof(struct heap));
    heap->max_entries = max_entries;
    heap->num_entries = 0;
    heap->entries = (int64_t *) st_malloc(sizeof(int64_t) * max_entries);
    if (!heap->entries) {
        free(heap);
        assert(0);
        return 0;
    }
    //memset(heap->entries, LONG_64_MIN, sizeof(LONG_64) * max_entries);
    for(i=0; i<max_entries; i++) {
        heap->entries[i] = INT64_MIN;
    }
    return heap;
}

void heap_destroy(struct heap *heap)
{   
    if (heap->entries)
        free(heap->entries);
    free(heap);
}

void heapify(struct heap *heap, int64_t i)
{
    uint64_t left;
    uint64_t right;
    int64_t largest;

    left = 2 * i;
    right = left + 1;

    if (left < heap->num_entries &&
        heap->entries[left] > heap->entries[i])
        largest = left;
    else
        largest = i;

    if (right < heap->num_entries &&
        heap->entries[right] > heap->entries[largest])
        largest = right;

    if (largest != i) {
        int64_t tmp = heap->entries[i];
        heap->entries[i] = heap->entries[largest];
        heap->entries[largest] = tmp;
        heapify(heap, largest);
    }
}

int64_t heap_expand(struct heap *heap)
{
    int64_t *tmp;
    tmp =
        (int64_t *) realloc(heap->entries,
                 heap->max_entries * 2 * sizeof(int64_t));
    if (!tmp) {
        assert(0);
        return -1;
    }
    heap->entries = tmp;
    heap->max_entries *= 2;
    return 0;
}

int64_t heap_insert(struct heap *heap, int64_t entry)
{
    int64_t i;

    /* Expand the heap if necessary */
    if (heap->num_entries == heap->max_entries) {
        if (heap_expand(heap) < 0) {
            assert(0);
            return -1;
        }
    }

    /* Put the new entry at the bottom of the heap */
    i = heap->num_entries++;

    /* Percolate the new entry up to where it belongs */
    while (i > 0 && heap->entries[i / 2] < entry) {
        heap->entries[i] = heap->entries[i / 2];
        i /= 2;
    }

    heap->entries[i] = entry; 

    return 0;
}

int64_t heap_extract(struct heap *heap)
{
    int64_t max;

    if (!heap->num_entries)
        return INT64_MIN;

    max = heap->entries[0];
    heap->entries[0] = heap->entries[heap->num_entries - 1];
    heap->num_entries--;
    heapify(heap, 0);

    return max;
}

int64_t heap_peek(struct heap *heap) {
    if (!heap->num_entries)
        return INT64_MIN;
        
    return heap->entries[0];
}

int64_t heap_empty(struct heap *heap) { 
    return !(heap->num_entries);
}


void heap_clean(struct heap *heap) 
{
    uint64_t i;
    
    //memset(heap->entries, LONG_64_MIN, heap->num_entries * sizeof(LONG_64));
    for(i=0; i<heap->num_entries; i++) {
        heap->entries[i] = INT64_MIN;
    }
    heap->num_entries = 0;
}

