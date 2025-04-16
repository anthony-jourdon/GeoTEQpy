#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "max_heap.h"

static int GetParentIdx(int idx)     { return (idx - 1) / 2; }
static int GetLeftChildIdx(int idx)  { return (2 * idx) + 1; }
static int GetRightChildIdx(int idx) { return (2 * idx) + 2; }

static void MaxHeapNodeInit(heap_node *node)
{
  node->distance   = 1e32;
  node->global_idx = -1;
  return;
}

void MaxHeapCreate(long int size, maxheap **heap)
{
  maxheap *_heap;
  int     i;

  _heap = (maxheap *)malloc(sizeof(maxheap));
  _heap->npoints = size;
  _heap->cnt     = 0;
  _heap->nodes   = (heap_node *)malloc(sizeof(heap_node) * size);
  for (i=0; i<size; i++) {
    MaxHeapNodeInit(&_heap->nodes[i]);
  }
  *heap = _heap;
  return;
}

void MaxHeapDestroy(maxheap *heap)
{
  if (!heap) return;
  if (heap->nodes) { free(heap->nodes); }
  heap->nodes = NULL;
  free(heap);
  heap = NULL;
  return;
}

void MaxHeapReset(maxheap *heap)
{
  int i;
  if (!heap) return;
  heap->cnt = 0;
  for (i=0; i<heap->npoints; i++) {
    MaxHeapNodeInit(&heap->nodes[i]);
  }
  return;
}

static void MaxHeapSwap(heap_node *a, heap_node *b)
{
  heap_node tmp;
  tmp = *a;
  *a = *b;
  *b = tmp;
  return;
}

static void MaxHeapMoveUp(maxheap *heap, int idx)
{
  int parent_idx;
  heap_node *node, *parent;

  node       = &heap->nodes[idx];        // node to move up
  parent_idx = GetParentIdx(idx);        // parent index
  parent     = &heap->nodes[parent_idx]; // parent node

  // while the distance of the node is greater than its parent
  while (idx > 0 && node->distance > parent->distance) {
    MaxHeapSwap(node, parent);             // swap the nodes
    idx        = parent_idx;               // update the index
    parent_idx = GetParentIdx(idx);        // get the new parent index
    node       = &heap->nodes[idx];        // move up the node
    parent     = &heap->nodes[parent_idx]; // get the new parent node
  }
}

static void MaxHeapMoveDown(maxheap *heap, int idx)
{
  int greater_idx, left_idx, right_idx;

  greater_idx = idx;
  left_idx    = GetLeftChildIdx(idx);
  right_idx   = GetRightChildIdx(idx);

  // if the left child is greater than the current node set it as the greater
  if (left_idx < heap->cnt && heap->nodes[left_idx].distance > heap->nodes[greater_idx].distance) {
    greater_idx = left_idx;
  }
  // if the right child is greater than the previously set greater node, set it as the greater
  if (right_idx < heap->cnt && heap->nodes[right_idx].distance > heap->nodes[greater_idx].distance) {
    greater_idx = right_idx;
  }
  // if the greater is not the current node, swap and move down
  if (greater_idx != idx) {
    MaxHeapSwap(&heap->nodes[idx], &heap->nodes[greater_idx]);
    MaxHeapMoveDown(heap, greater_idx);
  }
  return;
}

void MaxHeapInsert(maxheap *heap, heap_node *insert_node)
{
  heap_node *node;
  if (!heap) {
    fprintf(stderr, "[max_heap error] MaxHeapInsert() called with NULL heap.\n");
    return;
  }
  if (heap->cnt >= heap->npoints) {
    fprintf(stderr, "[max_heap error] MaxHeapInsert() called with full heap.\n");
    return;
  }
  if (!insert_node) {
    fprintf(stderr, "[max_heap error] MaxHeapInsert() called with NULL node.\n");
    return;
  }

  node = &heap->nodes[heap->cnt];
  memcpy(node, insert_node, sizeof(heap_node)); // copy the node to insert
  MaxHeapMoveUp(heap, heap->cnt); // move the node up in the heap
  heap->cnt++;                // increment the count of nodes in the heap

  return;
}

void MaxHeapReplaceRoot(maxheap *heap, heap_node *replace_node)
{
  heap_node *node;
  if (!heap) {
    fprintf(stderr, "[max_heap error] MaxHeapReplace() called with NULL heap.\n");
    return;
  }
  if (!replace_node) {
    fprintf(stderr, "[max_heap error] MaxHeapReplace() called with NULL node.\n");
    return;
  }

  node = &heap->nodes[0]; // root node
  memcpy(node, replace_node, sizeof(heap_node)); // copy the node to insert
  MaxHeapMoveDown(heap, 0); // move the node down in the heap

  return;
}

void MaxHeapView(maxheap *heap)
{
  int i;
  if (!heap) {
    fprintf(stderr, "[max_heap error] MaxHeapView() called with NULL heap.\n");
    return;
  }
  printf("MaxHeapView\n");
  printf(" npoints: %d\n", heap->npoints);
  printf(" cnt:     %d\n", heap->cnt);
  for (i=0; i<heap->cnt; i++) {
    heap_node *node = &heap->nodes[i];
    printf(" [%d] distance: %g global_idx: %ld\n", i, node->distance, node->global_idx);
  }
  return;
}