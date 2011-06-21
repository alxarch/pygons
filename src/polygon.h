//      polygon.h
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


#ifndef POLYGON
#define POLYGON

#include <stdlib.h>
#include "globals.h"
#include "mathutils.h"
#include "bounds.h"

enum PL_Order{ CW=0, CCW };

enum PL_Flags{ 
	HIDDEN	=1, 
	SECTION	=2, 
	BACKF	=4, 
	CLIPPED =8, 
	HOLE    =16,
	SPLIT   =32 
};

struct polygon;

typedef struct PolyList{
	struct polygon ** polys;
	uint last, size;
} PolyList;

typedef struct polyvert{
	vec co;
	uint flags;
} PolyVert;

typedef struct polygon{
	struct polyvert * points;
	uint size, last;
	struct bounds bb;
	struct PolyList *holes;
	uint flags;
} Polygon;

Polygon * pl_init(const uint count);

Polygon * pl_copy(const Polygon * pl, void (*copy_method)(vec dest, const vec src), const int holes);

void pl_kill(Polygon ** pl);

PolyVert * pl_append(Polygon ** pl, const vec co, const int flag);
PolyVert * pl_append_vertex(Polygon ** pl, const PolyVert * vx);

PolyVert *pl_insert(Polygon * pl, const uint pos, const vec co, const int flags);
PolyVert *pl_insert_point(Polygon * pl, const uint pos, const vec co);
PolyVert *pl_insert_vertex(Polygon * pl, const uint pos, const PolyVert * vx);

PolyVert * pl_pop(Polygon * pl);
PolyVert * pl_remove(Polygon * pl, const uint idx);
PolyVert * pl_remove_vertex(Polygon * pl, const PolyVert * vx);

void pl_extend_points(Polygon * pl, const vec * points, const uint count);
void pl_extend_verts(Polygon * pl, const PolyVert * verts, const uint count);

int pl_add_hole(Polygon * pl, Polygon * hole);
int pl_self_intersects(const Polygon *pl, AxisPair ax);
void pl_reverse(Polygon * pl);

int pl_rel_point_2d(const Polygon * pl, const vec p, const AxisPair ax, const int holes);
int pl_rel_plane(const Polygon * pl, const Plane * pn, int *rmap);

double pl_area_2d(const Polygon * pl, const AxisPair ax, const int holes);
double pl_area_3d(const Polygon * pl, const int holes);

AxisPair pl_best_axis_pair(Polygon * pl);

void pl_normal(const Polygon * pl, vec normal);

uint pl_simplify_2d(Polygon * pl, AxisPair ax);
uint pl_simplify_3d(Polygon * pl);
uint pl_rm_doubles_2d(Polygon * pl, AxisPair ax);
uint pl_rm_doubles_3d(Polygon * pl);

enum PL_Order pl_order(const Polygon * pl);
enum PL_Order pl_order_2d(const Polygon *  pl, const AxisPair ax);

void pl_project_to_plane( Polygon * pl,
                          const Plane * const pn,
                          const vec direction);

int pl_split_by_plane(Polygon * pl,
					  Plane * pln,
					  PolyList * fparts, 
					  PolyList * bparts);

// PolyList ____________________________________________________________

PolyList * pll_init(const uint count);

void pll_kill(PolyList ** pll);
void pll_killall(PolyList ** pll);

PolyList * pll_copy(PolyList * pll, void (*copy_method)(vec dest, const vec src));

void pll_append(PolyList ** pll, Polygon * pl);
Polygon * pll_insert(PolyList * pll, const uint pos, const Polygon * pl);
Polygon * pll_pop(PolyList * pll);

Polygon * pll_remove(PolyList * pll, const uint idx);

int pll_remove_pl(PolyList * pll, const Polygon * pl);

Polygon * pll_replace(PolyList * pll, const uint idx, Polygon * pl);

int pll_replace_pl(PolyList * pll, const Polygon * to_replace, Polygon * pl);

void pll_extend(PolyList ** pll, Polygon ** polys, const uint count);

void pll_sort_by_area(PolyList * pll);
uint pll_count_points(PolyList * pll);

#endif
