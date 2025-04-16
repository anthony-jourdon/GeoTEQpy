#ifndef __max_heap_h__
#define __max_heap_h__

typedef struct {
  double   distance;   // distance to the point
  long int global_idx; // index in the original data set
} heap_node;

typedef struct {
  int       npoints; // total number of points
  int       cnt;     // current number of points in the heap
  heap_node *nodes;  // array of heap nodes
} maxheap;

/* prototypes */
void MaxHeapCreate(long int size, maxheap **heap);
void MaxHeapDestroy(maxheap *heap);
void MaxHeapReset(maxheap *heap);
void MaxHeapInsert(maxheap *heap, heap_node *insert_node);
void MaxHeapReplaceRoot(maxheap *heap, heap_node *replace_node);
void MaxHeapView(maxheap *heap);
#endif