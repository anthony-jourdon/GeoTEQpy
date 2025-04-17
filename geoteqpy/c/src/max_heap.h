/*
====================================================================================================
  Copyright (c) 2024, 
  Anthony Jourdon, 

  project:  GeoTEQpy
  filename: max_heap.h

  This file is part of GeoTEQpy.

  GeoTEQpy is free software: you can redistribute it and/or modify it under the terms 
  of the GNU General Public License as published by the Free Software Foundation, either 
  version 3 of the License, or any later version.

  GeoTEQpy is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with GeoTEQpy. 
  If not, see <https://www.gnu.org/licenses/>.
====================================================================================================
*/

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