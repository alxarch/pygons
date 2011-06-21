//      split.h
//      
//      Copyright 2011 Alexandros Sigalas <alxarch@gmail.com>
//      
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//      
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//      
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.

#ifndef _SPLIT
#define _SPLIT

#include <stdlib.h>
#include "math.h"
#include "polygon.h"

typedef struct spl_node{
	vec co;
	int r, flags;
	Polygon * pl;
	struct spl_node * nxt;
	struct spl_node * prv;
} SplitNode;

typedef struct spl_graph{
	SplitNode * nodes;
	uint last;
	Plane *pln;
	Polygon *pl;
	int r;
} SplitGraph;

SplitGraph * build_split_graph(Polygon * pl, Plane * pln);

void kill_split_graph(SplitGraph ** sg);

void split_graph(SplitGraph *sg, PolyList * fparts, PolyList * bparts);

#endif
