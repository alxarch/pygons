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

typedef enum PL_Order{ CW=0, CCW } order;
enum PL_Flags{ 
	HIDDEN=1, 
	SECTION=2, 
	BACKF=4, 
	CLIPPED=8, 
	HOLE=16,
	SPLIT=32 
};

struct Polygon;

typedef struct PL_List{
	struct Polygon ** polys;
	uint last, size;
} polylist;

typedef struct PL_Vertex{
	vec co;
	int flags;
} polyvert;

typedef struct Polygon{
	polyvert * points;
	uint size, last;
	struct Bounds bb;
	struct PL_List *holes;
	int flags;
} polygon;

polygon * pl_init(const uint count);

polygon * pl_copy(const polygon * pl, void (*copy_method)(vec * dest, const vec * src), const int holes);

void pl_kill(polygon ** pl);

polyvert * pl_append(polygon ** pl, const vec * co, const int flag);
polyvert * pl_append_vertex(polygon ** pl, const polyvert * vx);

polyvert *pl_insert(polygon *pl, const uint pos, const vec * co, const int flag);
polyvert *pl_insert_point(polygon * pl, const uint pos, const vec *co);
polyvert *pl_insert_vertex(polygon * pl, const uint pos, const polyvert * vx);

polyvert * pl_pop(polygon * pl);
polyvert * pl_remove(polygon * pl, const uint idx);
polyvert * pl_remove_vertex(polygon * pl, const polyvert * vx);

void pl_extend_points(polygon * pl, const vec * points, const uint count);
void pl_extend_verts(polygon * pl, const polyvert * verts, const uint count);

int pl_add_hole(polygon * pl, polygon * hole);
int pl_self_intersects(const polygon *pl, enum axis_pair ax);
void pl_reverse(polygon * pl);

int pl_rel_point_2d(const polygon * pl, const vec * p, const enum axis_pair ax, const int holes);
int pl_rel_plane(const polygon * pl, const plane * pn, int *rmap);

double pl_area_2d(const polygon * pl, const enum axis_pair ax, const int holes);
double pl_area_3d(const polygon * pl, const int holes);

enum axis_pair pl_best_axis_pair(struct Polygon * pl);

void pl_normal(const polygon * pl, vec * normal);

uint pl_simplify_2d(polygon * pl, enum axis_pair ax);
uint pl_simplify_3d(polygon * pl);
uint pl_rm_doubles_2d(polygon * pl, enum axis_pair ax);
uint pl_rm_doubles_3d(polygon * pl);

enum PL_Order pl_order(const polygon * pl);
enum PL_Order pl_order_2d(const polygon *  pl, const enum axis_pair ax);

void pl_project_to_plane( struct Polygon * pl,
                          const plane * const pn,
                          const vec * const direction);

int pl_split_by_plane(struct Polygon * pl,
					  const struct Plane * const pln,
					  struct PL_List * fparts, 
					  struct PL_List * bparts);

// polylist ____________________________________________________________

polylist * pll_init(const uint count);

void pll_kill(polylist ** pll);
void pll_killall(polylist ** pll);

polylist * pll_copy(polylist * pll, void (*copy_method)(vec * dest, const vec * src));

void pll_append(polylist ** pll, polygon * pl);
polygon * pll_insert(polylist * pll, const uint pos, const polygon * pl);
polygon * pll_pop(polylist * pll);

polygon * pll_remove(polylist * pll, const uint idx);

int pll_remove_pl(polylist * pll, const polygon * pl);

polygon * pll_replace(polylist * pll, const uint idx, polygon * pl);

int pll_replace_pl(polylist * pll, const polygon * to_replace, polygon * pl);

void pll_extend(polylist ** pll, polygon ** polys, const uint count);

void pll_sort_by_area(polylist * pll);
uint pll_count_points(polylist * pll);

#endif
