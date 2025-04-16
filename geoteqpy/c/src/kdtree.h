/*
 
 K-d tree
 Copyright (C) 2012 RosettaCode User:Ledrug.
 Permission is granted to copy, distribute and/or modify this document
 under the terms of the GNU Free Documentation License, Version 1.2
 or any later version published by the Free Software Foundation;
 with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.
 A copy of the license is included in the section entitled "GNU
 Free Documentation License".
 
 C source from:
 https://rosettacode.org/wiki/K-d_tree
 
 Page revision information:
 14:20, 8 April 2017 Trizen (Talk | contribs) m . . (73,308 bytes) (0) . . (->{{header|Sidef}}:  updated code) (undo)
 
 Notes regarding usage can be found under "C Entry":
 https://rosettacode.org/wiki/Talk:K-d_tree#C_Entry
 
 Caution / implementation limitation:
 - The method does not behave correctly if two input nodes (different pointer but identical x[] values)
 are placed within the kdtree.

*/
 
#ifndef __kdtree_h__
#define __kdtree_h__

#include "max_heap.h"

#define KDTR_MAX_DIM 3

typedef struct _p_kd_node_t* kd_node;
struct _p_kd_node_t {
  double  x[KDTR_MAX_DIM];
  kd_node left,right;
  int     index;
};

typedef struct _p_KDTree *KDTree;
struct _p_KDTree {
  int     npoints,cnt;
  kd_node root,point;
  int     dim;
  int     visited;
  int     setup;
};

/* prototypes */
void kdtr_node_init(kd_node n);
void KDTreeCreate(int dim,KDTree *_k);
void KDTreeDestroy(KDTree *_k);
void KDTreeReset(KDTree kt);
void KDTreeView(KDTree kt);
void KDTreeSetPoints(KDTree k,int np);
void KDTreeGetPoints(KDTree k,int *n,kd_node *nodes);
void KDTreeInsertPoint(KDTree k,double coor[]);
void KDTreeSetup(KDTree kt);
void KDTreeFindNearest(KDTree k,double coor[],kd_node *nearest,double *sep);
void KDTreeFindKNearest(KDTree kt, double coor[], maxheap *heap);
#endif
