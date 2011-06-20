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

struct SplitNode{
	struct Vec co;
	int r, flags;
	struct Polygon * pl;
	struct SplitNode * nxt;
	struct SplitNode * prv;
};

struct SplitGraph{
	struct SplitNode * nodes;
	uint last;
	struct Plane *pln;
	struct Polygon *pl;
	int r;
};

struct SplitGraph * build_split_graph(struct Polygon * pl, struct Plane * pln);

void kill_split_graph(struct SplitGraph ** sg);

struct Polygon * new_split(struct SplitNode * nd);

void split_graph(struct SplitGraph *sg, struct PL_List * fparts, struct PL_List * bparts);

#endif
