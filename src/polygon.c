//      polygon.c
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


#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "globals.h"
#include "mathutils.h"
#include "split.h"
#include "polygon.h"

// Polygon _____________________________________________________________

Polygon * 
pl_init(const uint count)
{
    Polygon * pl= (Polygon *)calloc(1, sizeof(Polygon));
    assert(pl);
    
    pl->points= (PolyVert *)malloc(count * sizeof(PolyVert));
    assert(pl->points);
    
    pl->size= count;
    
    return pl;
}

/*
 * Copies a polygon struct.
 * 
 * Also allows to create transformed copies of the struct Polygon by using a 
 * function that copies the points.
 * 
 * @param Polygon * 	pl		Pointer to the polygon to be copied.
 * @param void * 		   	cpy 	A pointer to a function to transform the coords.
 * @param int   			holes 	Whether or not to also copy the holes.
 * 
 * @return Polygon *  A pointer to the copied polygon.
 * 
 */
Polygon * 
pl_copy(const Polygon * pl, 
		void (*cpy)(vec dest, const vec src), 
		const int holes)
{
    assert(pl);
    
    Polygon * cp= (Polygon *)malloc(sizeof(Polygon));
    assert(cp);
    
    memcpy(cp, pl, sizeof(Polygon));
    
    cp->points= (PolyVert *)malloc(pl->size * sizeof(PolyVert));
    assert(cp->points);
    
    uint i;
    
    if (cpy){
        for (i= 0; i < pl->last; i++){
            cpy(cp->points[i].co, pl->points[i].co);
            cp->points[i].flags = pl->points[i].flags;
        }
        cpy(cp->bb.min, pl->bb.min);
        cpy(cp->bb.max, pl->bb.max);
    }
    else{
        memcpy(cp->points, pl->points, pl->last * sizeof(PolyVert));
        // min & max have already been copied by memcpy.
    }
    
    cp->holes = (holes) ? pll_copy(pl->holes, cpy) : NULL;
    
    return cp;
}

/*
 * Frees memory of struct Polygon and nullifies the pointer.
 */
void pl_kill(Polygon ** pl_pt)
{
	assert(pl_pt);
	
	Polygon * pl= *pl_pt;
	
	if(pl){
		if (pl->points) free(pl->points);
		if (pl->holes) pll_killall(&pl->holes);
		free(pl);
	}
}

static void 
pl_check_bounds(Polygon * pl, const vec p)
{
    if (pl->last > 1){
		check_bounds(&pl->bb, p);
    }
    else{
        v_cpy(pl->bb.min, p);
        v_cpy(pl->bb.max, p);
    }
}

AxisPair 
pl_best_axis_pair(Polygon * pl)
{
	vec n;
	pl_normal(pl, n);
	return fabs(n[2]) > ZEROLENGTH ? XY : (fabs(n[0]) > fabs(n[1]) ? XZ : YZ);
}

static void 
pl_recalc_bounds(Polygon * pl)
{
    assert(pl);
    
    v_cpy(pl->bb.min, pl->points[0].co);
    v_cpy(pl->bb.max, pl->points[0].co);
    
    uint i;
    for (i=1; i < pl->last; i++){
        check_bounds(&pl->bb, pl->points[i].co);
	}
}

static Polygon * 
pl_check_pointer(Polygon ** pl_pt){
	assert(pl_pt);
	
	Polygon * result= *pl_pt;
	
	if (!result){
		result= *pl_pt= pl_init(10);
	}
	
	return result;
}

void pl_expand(Polygon * pl, const uint extra_space)
{
    assert(pl);
    
    pl->points= (PolyVert *)realloc(pl->points, pl->size += extra_space);
    assert(pl->points);
}

PolyVert * 
pl_append(Polygon **pl_pt, const vec co, int flags)
{
    assert(co);

	Polygon * pl= pl_check_pointer(pl_pt);
    
    if (pl->last > pl->size) pl_expand(pl, pl->size);
	
    PolyVert * vx= pl->points + pl->last++;
    v_cpy(vx->co, co);
    vx->flags= flags;
    
    pl_check_bounds(pl, co);
    
    return vx;
}

PolyVert * 
pl_append_vertex(Polygon ** pl_pt, const PolyVert * const pt)
{
    assert(pt);
    
    Polygon * pl= pl_check_pointer(pl_pt);
    
    if (pl->last > pl->size) pl_expand(pl, pl->size);
    
    PolyVert * vx= pl->points + pl->last++;
    *vx= *pt;
    
    pl_check_bounds(pl, pt->co);
    
    return vx;
}

void pl_check_size(Polygon * pl)
{
	if(pl->last < pl->size) {
		pl_expand(pl, pl->size);
	}
}

PolyVert * 
pl_insert_point(Polygon * pl, uint pos, const vec p)
{
    assert(pl);
    assert(p);
    pl_check_size(pl);
    
    uint i= pl->last++;
    uint j= pl->last;
    for (; i >= pos; j=i--){
        pl->points[j]= pl->points[i];
	}
	
    PolyVert *vx= pl->points + i;
    v_cpy(vx->co, p);
    vx->flags= 0;
    
    pl_check_bounds(pl, p);
    
    return vx;
}

PolyVert * pl_insert_vertex(Polygon * pl, uint pos, const PolyVert * p)
{
    assert(pl);
    assert(p);
    
    pl_check_size(pl);
    
    uint i= pl->last++;
    uint j= pl->last;
    for (; i >= pos; j=i--){
        pl->points[j]= pl->points[i];
	}
	
    PolyVert * vx= pl->points + i;
    *vx= *p;
    
    pl_check_bounds(pl, p->co);
    
    return vx;
}

PolyVert * pl_pop(Polygon * pl)
{
    assert(pl);
    
    if (pl->last){
        PolyVert * result= (PolyVert *)malloc(sizeof(PolyVert));
        pl->last--;
        PolyVert * last= pl->points + pl->last;
        v_cpy(result->co, last->co);
        result->flags= last->flags;
        
        return result;
    }
    
    return NULL;
}

PolyVert * pl_remove(Polygon * pl, const uint idx)
{
    assert(pl);
    
    if (idx < pl->last){
        PolyVert *result= (PolyVert*)malloc(sizeof(PolyVert));
        memcpy(result, pl->points+idx, sizeof(PolyVert));
        uint i=idx;
        PolyVert *vx= pl->points + i++;
        for(;i < pl->last; i++)
                *vx++=pl->points[i];
        pl->last--;
        pl_recalc_bounds(pl);
    }
    return NULL;
}

PolyVert * 
pl_remove_vertex(Polygon * pl, const PolyVert * const vx)
{
    assert(pl);
    assert(vx);
    
    uint i;
    for (i=0; i < pl->last; i++){
        if (pl->points + i == vx) return pl_remove(pl, i);
    }
    
    return NULL;
}

void 
pl_extend_verts(Polygon * pl, const PolyVert * const points, uint count)
{
    assert(pl);
    assert(points);
    
    pl_check_size(pl);

    uint i;
    for (i=0; i < count; i++){
        pl->points[pl->last++]= points[i];
	}
}

void 
pl_extend_points(Polygon * pl, const vec * points, uint count)
{
    assert(pl);
    assert(points);
    
	PolyVert *tmp;
	
	pl_check_size(pl);
	
    uint i;
    for (i=0 ; i < count; i++){
		tmp= pl->points + pl->last++;
		tmp->flags= 0;
		v_cpy(tmp->co, points[i]);
	}
}

int 
pl_add_hole(Polygon * pl, Polygon * hole)
{
    //TODO:fix this
    pll_append(&pl->holes, hole);
    return 1;
}

int 
pl_rel_point_2d(const Polygon * pl, 
				const vec p, 
				AxisPair ax, 
				const int do_holes)
{
    assert(pl);
    assert(p);
    
    if (!point_in_bounds(&pl->bb, p)) return OUT;
    
    uint i;
    if (do_holes && pl->holes){
        int rel;
        for (i=0; i < pl->holes->last; i++){
            if ((rel= pl_rel_point_2d(pl->holes->polys[i], p, ax, 0)) != OUT){
                return (rel == IN) ? OUT : ON;
            }
        }
    }
    uint j, wn= 0;
    double x, y, x0, y0, x1, y1, ux, uy;
    int _x= (ax == YZ) ? 1 : 0;
    int _y= (ax == XY) ? 1 : 2;
    double *v0= (double*)p, *v1;
    x= *(v0 + _x);
    y= *(v0 + _y);
    for (i=(pl->last - 1), j=0; j < pl->last; i=j++){
        v0= (double *)&pl->points[i].co;
        v1= (double *)&pl->points[j].co;
        x0= *(v0 + _x);
        y0= *(v0 + _y);
        x1= *(v1 + _x);
        y1= *(v1 + _y);
        ux= x1 - x0;
        uy= y1 - y0;
        if (fabs((x - x0) * uy - (y - y0) * ux) < E){
            if (fabs(ux) > fabs(uy)){
                if ((x0 < x && x < x1) || (x1 < x && x < x0)){
                    return ON;
                }
            }
            else if ((y0 < y && y < y1) || (y1 < y && y < y0)){
                    return ON;
            }
        }

        if (y0 <= y){
            if ((y1 > y) && ((x - x0) * uy - (y - y0) * ux > 0.0)){
                wn += 1;
            }
        }else if((y1 <= y) && ((x - x0) * uy - (y - y0) * ux < 0.0)){
            wn -= 1;
        }
    }
    return (wn) ? IN : OUT;
}

double pl_signed_area_2d(const Polygon * pl, AxisPair ax)
{
    assert(pl);
    
    double area= 0.0;
    
    if (pl->last < 3) return area;
    
    uint i, j;
    for (i=(pl->last - 1), j=0; j < pl->last; i=j++){
        area += v_perp(pl->points[i].co, pl->points[j].co, ax);
	}
    
    area/= 2.0;
    
    return area;
}

void 
pl_area_vec(const Polygon * pl, vec dest)
{
    assert(pl);
    
    v_zero(dest);
    
    if (pl->last < 3) return;
    
    uint i, j;
    vec tmp;
    for (i=pl->last - 1, j=0; j < pl->last; i=j++){
		v_cross(tmp, pl->points[i].co, pl->points[j].co);
		v_add(dest, dest, tmp);
	}
}

double 
pl_area_3d(const Polygon * pl, const int do_holes)
{
    assert(pl);
    
    double area= 0.0;
    
    if (pl->last < 3) return area;
    
    vec tmp;
    pl_area_vec(pl, tmp);
    area= v_length(tmp);
    
    if (do_holes && pl->holes){
        uint i;
        for (i=0; i < pl->holes->last; i++){
            pl_area_vec(pl->holes->polys[i], tmp);
            area-= v_length(tmp);
        }
    }
    
    return area;
}

double 
pl_area_2d(const Polygon * pl, AxisPair ax, const int do_holes)
{
    assert(pl);
    
    double area= 0.0;
    if (pl->last < 3) return area;
    
    area= fabs(pl_signed_area_2d(pl, ax));
    
    if (do_holes && pl->holes){
        uint i;
        for (i=0; i < pl->holes->last; i++){
            area-= fabs(pl_signed_area_2d(pl->holes->polys[i], ax));
		}
    }
    
    return area;
}

enum PL_Order 
pl_order_2d(const Polygon * pl, AxisPair ax)
{
    return (pl_signed_area_2d(pl, ax) > 0.0) ? CW : CCW;
}

enum PL_Order 
pl_order(const Polygon * pl)
{
    vec n;
    
    pl_normal(pl, n);
    
    if (fabs(n[2]) < ZEROLENGTH){
        return (fabs(n[0]) > fabs(n[1])) ? 
			pl_order_2d(pl, XZ) : 
			pl_order_2d(pl, YZ);
    }
    
    return pl_order_2d(pl, XY);
	
}
    
void 
pl_normal(const Polygon * pl, vec normal)
{
    assert(pl);
    
    if (pl->last < 3) v_zero(normal);
    else{
		pl_area_vec(pl, normal);
		v_normal(normal, normal);
	}
}

int 
pl_rel_plane(const Polygon * pl, const Plane * const pn, int * rmap)
{
	assert(pl);
	assert(pn);
	
    int R = 0;
    uint i;
    if (rmap){
        for (i=0; i < pl->last; i++)
            R |= *rmap++= rel_point_plane(pl->points[i].co, pn);
        
    }
    else{
        for (i=0; i < pl->last; i++)
            R |= rel_point_plane(pl->points[i].co, pn);
    }
    return R;
}

int 
pl_self_intersects(const Polygon * pl, AxisPair ax)
{
    vec diffs[pl->last];
    uint i, j;
    for(i=pl->last - 1, j=0; j < pl->last; i=j++){
		v_sub(diffs[i], pl->points[j].co, pl->points[i].co);
	}
	
    uint k, m;
    vec w;
    double d, ratio;
    for(i=pl->last - 1, j=0; j < pl->last; i=j++){
        for (k=(pl->last - 1), m=0; j < pl->last; k=m++){
			
            if (i == k || i == m || j == k || j == m) continue;
            
            d = v_perp(diffs[i], diffs[k], ax);
            
            if (fabs(d) < E) continue;
            
            v_sub(w, pl->points[i].co, pl->points[k].co);
            
            ratio = v_perp(diffs[i], w, ax) / d;
            
            if (0.0 < ratio && ratio < 1.0){
                ratio = v_perp(diffs[k], w, ax) / d;
                if (0.0 < ratio && ratio < 1.0) return 1;
            }
        }
    }
    
    return 0;
}

uint 
pl_simplify_2d(Polygon * pl, AxisPair ax)
{
    assert(pl);
    
    if (pl->last < 3) return pl->last;
    
    uint count= 0;
    uint i, j;
    vec diffs[pl->last];
    PolyVert *new_points= (PolyVert *)malloc(pl->last * sizeof(PolyVert));
    assert(new_points);
    
    for (i=(pl->last - 1), j=0; j < pl->last; i=j++){
        v_sub(diffs[i], pl->points[i].co, pl->points[j].co);
	}
    for (i=(pl->last - 1), j= 0; j < pl->last; i=j++){
        if (fabs(v_perp(diffs[i], diffs[j], ax)) > ZEROAREA)
            new_points[count++]= pl->points[j];
    }
    
    free(pl->points);
    pl->points= new_points;
    pl->last= count;
    
    if (pl->holes){
        for (i=0; i < pl->holes->last; i++){
            pl_simplify_2d(pl->holes->polys[i], ax);
		}
    }
    
    return count;
};

uint 
pl_simplify_3d(Polygon * pl)
{
    assert(pl);
    if (pl->last < 3) return pl->last;

	vec n;
    pl_normal(pl, n);
	AxisPair ax;
	if (fabs(n[2]) < ZEROLENGTH){
		ax= (fabs(n[1]) > fabs(n[0])) ? XZ : YZ;
	}
	else{
		ax= XY;
	}
	
	return pl_simplify_2d(pl, ax);
};

uint 
pl_rm_doubles_2d(Polygon * pl, AxisPair ax)
{
    PolyVert *new_points= (PolyVert*)malloc(pl->last * sizeof(PolyVert));
    assert(new_points);
    
    uint i, j, count= 0;
    for (i= (pl->last - 1), j= 0; j < pl->last; i=j++){
        if (!same_2d(pl->points[i].co, pl->points[j].co, ax)){
            new_points[count++]= pl->points[j];
        }
    }
    
    if (count < pl->last){
        free(pl->points);
        pl->points= new_points;
        pl->last= count;
    }
    
    pl_recalc_bounds(pl);
    
    if (pl->holes){
        for (i=0; i < pl->holes->last; i++){
            pl_rm_doubles_2d(pl->holes->polys[i], ax);
		}
    }
    
    return count;
}

uint pl_rm_doubles_3d(Polygon * pl)
{
    PolyVert * new_points = (PolyVert *)malloc(pl->last * sizeof(PolyVert));
    assert(new_points);
    
    uint i, j, count=0;
    for (i= (pl->last - 1), j= 0; j < pl->last; i=j++){
        if (!same_3d(pl->points[i].co, pl->points[j].co)){
            new_points[count++]= pl->points[j];
        }
    }
    
    if (count < pl->last){
        free(pl->points);
        pl->points= new_points;
        pl->last= count;
    }
    
    pl_recalc_bounds(pl);
    
    if (pl->holes){
        for (i=0; i < pl->holes->last; i++){
            pl_rm_doubles_3d(pl->holes->polys[i]);
		}
    }
    
    return count;
}

// keeps first point same
void pl_reverse(Polygon * pl)
{
    assert(pl);
    
    if (pl->last < 3) return; //nothing to reverse
    
    PolyVert *rpoints = (PolyVert *)calloc(pl->size, sizeof(PolyVert));
    PolyVert *last = pl->points + pl->last - 1;
    
    uint i;
    rpoints[0]= pl->points[0];
    rpoints[0].flags= last->flags;
    for (i=1; i < pl->last; i++){
        rpoints[i]= *last--;
        rpoints[i].flags= last->flags;
    }
    
    free(pl->points);
    pl->points= rpoints;
    
    if (pl->holes){
        for (i=0; i < pl->holes->last; i++)
            pl_reverse(pl->holes->polys[i]);
    }
}

void pl_project_to_plane(Polygon * pl, 
                         const Plane * const pn, 
                         const vec direction)
{
    assert(pl);
    assert(pn);
    
    double d= v_dot(pn->normal, direction);
    
    // otherwise projection is impossible
    assert(fabs(d) > E);
    
    vec tmp;
    double * pt;
    uint i;
    double N;
    for (i=0; i < pl->last; i++){
		pt= pl->points[i].co;
        v_sub(tmp, pt, pn->point);
        N= -v_dot(pn->normal, tmp);
        v_scale(tmp, direction, N/d);
        v_add(pt, pt, tmp);
    }
    v_sub(tmp, pl->bb.min, pn->point);
    N= -v_dot(pn->normal, tmp);
    v_scale(tmp, direction, N/d);
    v_add(pl->bb.min, pl->bb.min, tmp);
    
    v_sub(tmp, pl->bb.max, pn->point);
    N= -v_dot(pn->normal, tmp);
    v_scale(tmp, direction, N/d);
    v_add(pl->bb.max, pl->bb.max, tmp);
}
			
int 
pl_split_by_plane(Polygon * pl, 
				  Plane * pln, 
				  PolyList * fparts, 
				  PolyList * bparts)
{
	SplitGraph * sg = build_split_graph(pl, pln);
	
	switch(sg->r){
		case ON:
			break;
		case BK:
		case ONBK:
			pll_append(&bparts, pl);
			break;
		case FR:
		case ONFR:
			pll_append(&fparts, pl);
			break;
		default:
			split_graph(sg, fparts, bparts);
	}
	
	return sg->r;
}


// PL_List _____________________________________________________________

/*
 * Initializes a PL_List struct.
 */
 
PolyList * 
pll_init(const uint count)
{
    PolyList * pll= (PolyList *)calloc(1, sizeof(PolyList));
    assert(pll);
    
    pll->polys= (Polygon**)malloc(count * sizeof(Polygon**));
    assert(pll->polys);
    
    pll->size= count;
    
    return pll;
}

/*
 * Ensures a PL_List struct is present at the given pointer.
 */

PolyList * 
pll_check_pointer(PolyList ** pll_pt)
{
    assert(pll_pt);
    
    PolyList * result= *pll_pt;
    if (!result)
        result= *pll_pt= pll_init(10);
        
    return result;
}

/*
 * Increases the size of the PL_List.
 * 
 */
 
void 
pll_expand(PolyList * pll, const uint count)
{
    assert(pll);
    pll->polys= (Polygon **)realloc(pll->polys, pll->size+= count);
    
    assert(pll->polys);
}


/*
 * Checks whether the PL_List is full and doubles its size if so.
 * 
 */
 
void 
pll_check_size(PolyList * pll)
{
	 if(pll->last > pll->size){
        pll_expand(pll, 2 * pll->size);
    }
}

/*
 * Appends a polygon to a PL_List struct
 * 
 */
void pll_append(PolyList ** pll_pt, Polygon * const pl)
{

    PolyList * pll= pll_check_pointer(pll_pt);
    pll_check_size(pll);
   
    pll->polys[pll->last++]= pl;
}


Polygon * 
pll_pop(PolyList * pll)
{
    assert(pll);
    
    return (pll->last) ? pll->polys[--pll->last] : NULL;
}

PolyList * 
pll_copy(PolyList * pll, void (*copy_method)(vec dest, const vec src))
{
    if (!pll) return NULL;
    
    PolyList * cp = pll_init(pll->size);
    
    uint i;
    for (i=0; i < pll->last; i++){
		pll_append(&cp, pl_copy(pll->polys[i], copy_method, 1));
    }
    
    return cp;
}

void pll_kill(PolyList ** pll_pt)
{
    assert(pll_pt);
    
    PolyList * pll= *pll_pt;
    assert(pll);
    
    if (pll->polys) free(pll->polys);
    
    free(pll);
    
    *pll_pt= NULL;
}

void pll_killall(PolyList ** pll_pt)
{
    assert(pll_pt);
    
    PolyList * pll= *pll_pt;
    assert(pll);
    
    if (pll->polys){
        uint i;
        for (i=0; i < pll->last; i++)
			pl_kill(&pll->polys[i]);
        
        free(pll->polys);
    }
    
    free(pll);
    
    *pll_pt= NULL;
}

int pll_remove_pl(PolyList * pll, const Polygon * pl)
{
    assert(pll);
    
    uint i;
    for (i=0; i < pll->last; i++){
        if (pll->polys[i] == pl) break;
    }
    
    if (i == pll->last)
        return 0;
    else
        pll_remove(pll, i);
    return 1;
}

Polygon * pll_remove(PolyList * pll, const uint idx)
{
    assert(pll);
    if (idx < pll->last){
        Polygon * result= pl_copy(pll->polys[idx], NULL, 1);
        uint i=idx;
        Polygon ** pl_pt= pll->polys + i++;
        for (; i < pll->last; i++)
            *pl_pt++= pll->polys[i];
        pll->last--;
        return result;
    }
    return NULL;
}

Polygon * pll_replace(PolyList * pll, const uint idx, Polygon * pl)
{
    assert(pll);
    assert(pl);
    Polygon * result= pll->polys[idx];
    pll->polys[idx]= pl;
    return result;
}

int pll_replace_pl(PolyList * pll, const Polygon * to_replace, Polygon * pl)
{
    assert(pll);
    
    uint i;
    for (i=0; i < pll->last; i++){
        if (pll->polys[i] == to_replace){
            pll_replace(pll, i, pl);
            return 1;
		}
    }
    
    return 0;
}

void pll_extend(PolyList ** pll_pt, Polygon ** polys, const uint count)
{
    assert(polys);
    
    PolyList * pll= pll_check_pointer(pll_pt);
    
    if(pll->last + count > pll->size){
        pll_expand(pll, (pll->size * 2) + count);
    }
    
    uint i;
    for (i=0; i < count; i++){
        pll->polys[pll->last++]= polys[i];
	}
}
 
uint pll_count_points(PolyList * pll)
{
    if(!pll) return 0;
    
    uint i, count=0;
    for (i=0; i < pll->last; i++)
        count+= pll->polys[i]->last;
    
    return count; 
}

int pl_cmp_area(const void *pl_pt_0, const void *pl_pt_1)
{
	double area_0= pl_area_3d(*((Polygon **)pl_pt_0), 1);
	double area_1= pl_area_3d(*((Polygon **)pl_pt_1), 1);
	if (area_0 < area_1)
		return 1;
	else if(area_0 > area_1)
		return 0;
	else
		return -1;
}
	
void pll_sort_by_area(PolyList * pll)
{
	qsort(pll->polys, pll->last, sizeof(struct Polygon*), pl_cmp_area);
}

