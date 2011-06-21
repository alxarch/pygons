//      boolean.h
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


#ifndef BOOL
#define BOOL

#include <stdlib.h>
#include "math.h"
#include "polygon.h"

enum BL_Op{ SUBTRACT_A_B, SUBTRACT_B_A, UNION, INTERSECT};
enum BL_Set{ A=0, B };
enum itype{ ENEX, EXEN, EN, EX, NOTYPE=0, INVALID };

typedef struct bpoly{
	Polygon * pl;
	int rel, hole;
	enum BL_Set set;
	int has_entry;
} BL_Poly;

typedef struct bnode{
	vec co;
	int flags;
	struct bnode *nxt, *prv, *isx;
	BL_Poly *pl;
	enum itype t;
	int entry;
} BL_Node;

typedef struct bgraph{
	BL_Node *nodes;
	BL_Poly *polys;
	BL_Poly *pl_a;
	BL_Poly *pl_b;
	uint last, size, pl_num, inum;
} BL_Graph;

BL_Graph * bl_build_graph(Polygon * pl_a, Polygon * pl_b);

PolyList * bl_operation(BL_Graph * g, const enum BL_Op op);

void bl_kill_graph(BL_Graph * g);

#endif
