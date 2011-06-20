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

#include "math.h"
#include "polygon.h"

enum BL_Op{ SUBTRACT_A_B, SUBTRACT_B_A, UNION, INTERSECT};
enum BL_Set{ A=0, B };
enum itype{ ENEX, EXEN, EN, EX, NOTYPE=0, INVALID };

typedef struct BL_Poly{
	struct Polygon * pl;
	int rel;
	int hole;
	enum BL_Set set;
	int has_entry;
} bpoly;

typedef struct BL_Node{
	struct Vec co;
	int flags;
	struct BL_Node *nxt, *prv, *isx;
	struct BL_Poly *pl;
	enum itype t;
	int entry;
} bnode;

typedef struct BL_Graph{
	struct BL_Node *nodes;
	unsigned int last, size;
	struct BL_Poly *polys;
	struct BL_Poly *pl_a;
	struct BL_Poly *pl_b;
	unsigned int pl_num;
	unsigned int inum;
} bgraph;

bgraph * bl_build_graph(polygon * pl_a, polygon * pl_b);

polylist * bl_operation(struct BL_Graph * g, const enum BL_Op op);

void bl_kill_graph(bgraph * g);

#endif
