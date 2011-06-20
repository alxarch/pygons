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

struct Polygon * 
pl_init(const uint count)
{
    struct Polygon * pl= (struct Polygon *)calloc(1, sizeof(struct Polygon));
    assert(pl);
    
    pl->points= (struct PL_Vertex *)malloc(count * sizeof(struct PL_Vertex));
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
 * @param struct Polygon * 	pl		Pointer to the polygon to be copied.
 * @param void * 		   	cpy 	A pointer to a function to transform the coords.
 * @param int   			holes 	Whether or not to also copy the holes.
 * 
 * @return struct Polygon *  A pointer to the copied polygon.
 * 
 */
struct Polygon * 
pl_copy(const struct Polygon * pl, 
		void (*cpy)(vec * dest, const vec * const src), 
		const int holes)
{
    assert(pl);
    
    struct Polygon * cp= (struct Polygon *)malloc(sizeof(struct Polygon));
    assert(cp);
    
    memcpy(cp, pl, sizeof(struct Polygon));
    
    cp->points= (struct PL_Vertex *)malloc(pl->size * sizeof(struct PL_Vertex));
    assert(cp->points);
    
    uint i;
    
    if (cpy){
        for (i= 0; i < pl->last; i++){
            cpy(&cp->points[i].co, &pl->points[i].co);
            cp->points[i].flags = pl->points[i].flags;
        }
        cpy(&cp->bb.min, &pl->bb.min);
        cpy(&cp->bb.max, &pl->bb.max);
    }
    else{
        memcpy(cp->points, pl->points, pl->last * sizeof(struct PL_Vertex));
        // min & max have already been copied by memcpy.
    }
    
    cp->holes = (holes) ? pll_copy(pl->holes, cpy) : NULL;
    
    return cp;
}

/*
 * Frees memory of struct Polygon and nullifies the pointer.
 */
void pl_kill(struct Polygon ** pl_pt)
{
	assert(pl_pt);
	
	struct Polygon * pl= *pl_pt;
	
	if(pl){
		if (pl->points) free(pl->points);
		if (pl->holes) pll_killall(&pl->holes);
		free(pl);
	}
}

static void 
pl_check_bounds(struct Polygon * pl, const vec *p)
{
    if (pl->last > 1){
		check_bounds(&pl->bb, p);
    }
    else{
        v_copy(&pl->bb.min, p);
        v_copy(&pl->bb.max, p);
    }
}

enum axis_pair 
pl_best_axis_pair(struct Polygon * pl)
{
	vec n;
	pl_normal(pl, &n);
	return fabs(n.z) > ZEROLENGTH ? XY : (fabs(n.x) > fabs(n.y) ? XZ : YZ);
}

static void pl_recalc_bounds(struct Polygon * pl)
{
    assert(pl);
    
    v_copy(&pl->bb.min, &pl->points[0].co);
    v_copy(&pl->bb.max, &pl->points[0].co);
    
    uint i;
    for (i=1; i < pl->last; i++){
        check_bounds(&pl->bb, &pl->points[i].co);
	}
}

static struct Polygon * pl_check_pointer(struct Polygon ** pl_pt){
	assert(pl_pt);
	
	struct Polygon * result= *pl_pt;
	
	if (!result){
		result= *pl_pt= pl_init(10);
	}
	
	return result;
}

void pl_expand(struct Polygon * pl, const uint extra_space)
{
    assert(pl);
    
    pl->points= (struct PL_Vertex *)realloc(pl->points, pl->size += extra_space);
    assert(pl->points);
}

struct PL_Vertex * 
pl_append(struct Polygon **pl_pt, const vec * const co, int flags)
{
    assert(co);

	struct Polygon * pl= pl_check_pointer(pl_pt);
    
    if (pl->last > pl->size) pl_expand(pl, pl->size);
	
    struct PL_Vertex * vx= pl->points + pl->last++;
    vx->co= *co;
    vx->flags= flags;
    
    pl_check_bounds(pl, co);
    
    return vx;
}

struct PL_Vertex * 
pl_append_vertex(struct Polygon ** pl_pt, const struct PL_Vertex * const pt)
{
    assert(pt);
    
    struct Polygon * pl= pl_check_pointer(pl_pt);
    
    if (pl->last > pl->size) pl_expand(pl, pl->size);
    
    polyvert * vx= pl->points + pl->last++;
    *vx= *pt;
    
    pl_check_bounds(pl, &pt->co);
    
    return vx;
}

void pl_check_size(struct Polygon * pl)
{
	if(pl->last < pl->size) {
		pl_expand(pl, pl->size);
	}
}

polyvert * pl_insert_point(struct Polygon * pl, uint pos, const vec * p)
{
    assert(pl);
    assert(p);
    pl_check_size(pl);
    
    uint i= pl->last++;
    uint j= pl->last;
    for (; i >= pos; j=i--){
        pl->points[j]= pl->points[i];
	}
	
    polyvert *vx= pl->points + i;
    vx->co= *p;
    vx->flags= 0;
    
    pl_check_bounds(pl, p);
    
    return vx;
}

polyvert * pl_insert_vertex(struct Polygon * pl, uint pos, const polyvert * p)
{
    assert(pl);
    assert(p);
    
    pl_check_size(pl);
    
    uint i= pl->last++;
    uint j= pl->last;
    for (; i >= pos; j=i--){
        pl->points[j]= pl->points[i];
	}
	
    polyvert * vx= pl->points + i;
    *vx= *p;
    
    pl_check_bounds(pl, &p->co);
    
    return vx;
}

polyvert * pl_pop(struct Polygon * pl)
{
    assert(pl);
    
    if (pl->last){
        polyvert * result= (polyvert *)malloc(sizeof(polyvert));
        pl->last--;
        polyvert * last= pl->points + pl->last;
        v_copy(&result->co, &last->co);
        result->flags= last->flags;
        
        return result;
    }
    
    return NULL;
}

polyvert * pl_remove(struct Polygon * pl, const uint idx)
{
    assert(pl);
    
    if (idx < pl->last){
        polyvert *result= (polyvert*)malloc(sizeof(polyvert));
        memcpy(result, pl->points+idx, sizeof(polyvert));
        uint i=idx;
        polyvert *vx= pl->points + i++;
        for(;i < pl->last; i++)
                *vx++=pl->points[i];
        pl->last--;
        pl_recalc_bounds(pl);
    }
    return NULL;
}

polyvert * pl_remove_vertex(struct Polygon * pl, const polyvert * const vx)
{
    assert(pl);
    assert(vx);
    
    uint i;
    for (i=0; i < pl->last; i++){
        if (pl->points + i == vx) return pl_remove(pl, i);
    }
    
    return NULL;
}

void pl_extend_verts(struct Polygon * pl, const polyvert * const points, uint count)
{
    assert(pl);
    assert(points);
    
    pl_check_size(pl);

    uint i;
    for (i=0; i < count; i++){
        pl->points[pl->last++]= points[i];
	}
}

void pl_extend_points(struct Polygon * pl, const vec * const points, uint count)
{
    assert(pl);
    assert(points);
    
	polyvert *tmp;
	
	pl_check_size(pl);
	
    uint i;
    for (i=0 ; i < count; i++){
		tmp= pl->points + pl->last++;
		tmp->flags= 0;
		tmp->co= points[i];
	}
}

int pl_add_hole(struct Polygon * pl, struct Polygon * hole)
{
    //TODO:fix this
    pll_append(&pl->holes, hole);
    return 1;
}

int pl_rel_point_2d(const struct Polygon * pl, const vec * const p, enum axis_pair ax, const int do_holes)
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

double pl_signed_area_2d(const struct Polygon * pl, enum axis_pair ax)
{
    assert(pl);
    
    double area= 0.0;
    
    if (pl->last < 3) return area;
    
    uint i, j;
    for (i=(pl->last - 1), j=0; j < pl->last; i=j++){
        area+= v_perp(&pl->points[i].co, &pl->points[j].co, ax);
	}
    
    area/= 2.0;
    
    return area;
}

void pl_area_vec(const struct Polygon * pl, vec * dest)
{
    assert(pl);
    
    v_zero(dest);
    
    if (pl->last < 3) return;
    
    uint i, j;
    vec tmp;
    for (i=pl->last - 1, j=0; j < pl->last; i=j++){
		v_cross(&tmp, &pl->points[i].co, &pl->points[j].co);
		v_add(dest, dest, &tmp);
	}
}

double pl_area_3d(const struct Polygon * pl, const int do_holes)
{
    assert(pl);
    
    double area= 0.0;
    
    if (pl->last < 3) return area;
    
    vec tmp;
    pl_area_vec(pl, &tmp);
    area= v_length(&tmp);
    
    if (do_holes && pl->holes){
        uint i;
        for (i=0; i < pl->holes->last; i++){
            pl_area_vec(pl->holes->polys[i], &tmp);
            area-= v_length(&tmp);
        }
    }
    
    return area;
}

double pl_area_2d(const struct Polygon * pl, enum axis_pair ax, const int do_holes)
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

enum PL_Order pl_order_2d(const struct Polygon * pl, enum axis_pair ax)
{
    return (pl_signed_area_2d(pl, ax) > 0.0) ? CW : CCW;
}

enum PL_Order pl_order(const struct Polygon * pl)
{
    vec n;
    
    pl_normal(pl, &n);
    
    if (fabs(n.z) < ZEROLENGTH){
        return (fabs(n.x) > fabs(n.y)) ? 
			pl_order_2d(pl, XZ) : 
			pl_order_2d(pl, YZ);
    }
    else{
        return pl_order_2d(pl, XY);
	}
}
    
void pl_normal(const struct Polygon * pl, vec *normal)
{
    assert(pl);
    
    if (pl->last < 3) v_zero(normal);
    else{
		pl_area_vec(pl, normal);
		v_normal(normal, normal);
	}
}

int pl_rel_plane(const struct Polygon * pl, const plane * const pn, int* rmap)
{
	assert(pl);
	assert(pn);
	
    int R= 0;
    uint i;
    if (rmap){
        for (i=0; i < pl->last; i++)
            R|= *rmap++= rel_point_plane(&pl->points[i].co, pn);
        
    }
    else{
        for (i=0; i < pl->last; i++)
            R|= rel_point_plane(&pl->points[i].co, pn);
    }
    return R;
}

int pl_self_intersects(const struct Polygon * pl, enum axis_pair ax)
{
    vec diffs[pl->last];
    uint i, j;
    for(i=pl->last - 1, j=0; j < pl->last; i=j++){
		v_sub(diffs + i, &pl->points[j].co, &pl->points[i].co);
	}
	
    uint k, m;
    vec w;
    double d, ratio;
    for(i=pl->last - 1, j=0; j < pl->last; i=j++){
        for (k=(pl->last - 1), m=0; j < pl->last; k=m++){
			
            if (i == k || i == m || j == k || j == m) continue;
            
            d= v_perp(diffs + i, diffs + k, ax);
            
            if (fabs(d) < E) continue;
            
            v_sub(&w, &pl->points[i].co, &pl->points[k].co);
            
            ratio= v_perp(diffs + i, &w, ax) / d;
            
            if (0.0 < ratio && ratio < 1.0){
                ratio= v_perp(diffs + k, &w, ax) / d;
                if (0.0 < ratio && ratio < 1.0) return 1;
            }
        }
    }
    
    return 0;
}

uint pl_simplify_2d(struct Polygon * pl, enum axis_pair ax)
{
    assert(pl);
    
    if (pl->last < 3) return pl->last;
    
    uint count= 0;
    uint i, j;
    vec diffs[pl->last];
    polyvert *new_points= (polyvert*)malloc(pl->last * sizeof(polyvert));
    assert(new_points);
    
    for (i=(pl->last - 1), j=0; j < pl->last; i=j++){
        v_sub(diffs + i, &pl->points[i].co, &pl->points[j].co);
	}
    for (i=(pl->last - 1), j= 0; j < pl->last; i=j++){
        if (fabs(v_perp(diffs + i, diffs + j, ax)) > ZEROAREA)
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

uint pl_simplify_3d(struct Polygon * pl)
{
    assert(pl);
    if (pl->last < 3) return pl->last;

	vec n;
    pl_normal(pl, &n);
	enum axis_pair ax;
	if (fabs(n.z) < ZEROLENGTH){
		ax= (fabs(n.y) > fabs(n.x)) ? XZ : YZ;
	}
	else{
		ax= XY;
	}
	
	return pl_simplify_2d(pl, ax);
};

uint pl_rm_doubles_2d(struct Polygon * pl, enum axis_pair ax)
{
    polyvert *new_points= (polyvert*)malloc(pl->last * sizeof(polyvert));
    assert(new_points);
    
    uint i, j, count= 0;
    for (i= (pl->last - 1), j= 0; j < pl->last; i=j++){
        if (!same_2d(&pl->points[i].co, &pl->points[j].co, ax)){
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

uint pl_rm_doubles_3d(struct Polygon * pl)
{
    polyvert *new_points= (polyvert*)malloc(pl->last * sizeof(polyvert));
    assert(new_points);
    
    uint i, j, count= 0;
    for (i= (pl->last - 1), j= 0; j < pl->last; i=j++){
        if (!same_3d(&pl->points[i].co, &pl->points[j].co)){
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
void pl_reverse(struct Polygon * pl)
{
    assert(pl);
    
    if (pl->last < 3) return; //nothing to reverse
    
    polyvert *rpoints= (polyvert *)calloc(pl->size, sizeof(polyvert));
    polyvert *last= pl->points + pl->last - 1;
    
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

void pl_project_to_plane(struct Polygon * pl, 
                         const plane * const pn, 
                         const vec * const direction)
{
    assert(pl);
    assert(pn);
    assert(direction);
    
    double d= v_dot(&pn->normal, direction);
    assert(fabs(d) > E);// otherwise projection is impossible
    vec tmp, *pt;
    uint i;
    double N;
    for (i=0; i < pl->last; i++){
		pt= &pl->points[i].co;
        v_sub(&tmp, pt, &pn->point);
        N= -v_dot(&pn->normal, &tmp);
        v_scale(&tmp, direction, N/d);
        v_add(pt, pt, &tmp);
    }
    v_sub(&tmp, &pl->bb.min, &pn->point);
    N= -v_dot(&pn->normal, &tmp);
    v_scale(&tmp, direction, N/d);
    v_add(&pl->bb.min, &pl->bb.min, &tmp);
    
    v_sub(&tmp, &pl->bb.max, &pn->point);
    N= -v_dot(&pn->normal, &tmp);
    v_scale(&tmp, direction, N/d);
    v_add(&pl->bb.max, &pl->bb.max, &tmp);
}
			
int pl_split_by_plane(struct Polygon * pl, 
					  const struct Plane * pln, 
					  struct PL_List * fparts, 
					  struct PL_List * bparts
					 )
{
	struct SplitGraph * sg = build_split_graph(pl, pln);
	
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
struct PL_List * 
pll_init(const uint count)
{
    struct PL_List * pll= (struct PL_List *)calloc(1, sizeof(struct PL_List));
    assert(pll);
    
    pll->polys= (struct Polygon**)malloc(count * sizeof(struct Polygon**));
    assert(pll->polys);
    
    pll->size= count;
    
    return pll;
}

/*
 * Ensures a PL_List struct is present at the given pointer.
 */
struct PL_List * 
pll_check_pointer(struct PL_List ** pll_pt)
{
    assert(pll_pt);
    
    struct PL_List * result= *pll_pt;
    if (!result)
        result= *pll_pt= pll_init(10);
        
    return result;
}

/*
 * Increases the size of the PL_List.
 * 
 */
void pll_expand(polylist * pll, const uint count)
{
    assert(pll);
    pll->polys= (struct Polygon **)realloc(pll->polys, pll->size+= count);
    
    assert(pll->polys);
}


/*
 * Checks whether the PL_List is full and doubles its size if so.
 * 
 */
void pll_check_size(struct PL_List * pll){
	 if(pll->last > pll->size){
        pll_expand(pll, 2 * pll->size);
    }
}

/*
 * Appends a polygon to a PL_List struct
 * 
 */
void pll_append(struct PL_List ** pll_pt, struct Polygon * const pl)
{

    struct PL_List * pll= pll_check_pointer(pll_pt);
    pll_check_size(pll);
   
    pll->polys[pll->last++]= pl;
}


struct Polygon * pll_pop(struct PL_List * pll)
{
    assert(pll);
    
    return (pll->last) ? pll->polys[--pll->last] : NULL;
}

polylist * pll_copy(polylist * pll, void (*copy_method)(vec * dest, const vec * const src))
{
    if (!pll) return NULL;
    
    polylist * cp= pll_init(pll->size);
    
    uint i;
    for (i=0; i < pll->last; i++){
		pll_append(&cp, pl_copy(pll->polys[i], copy_method, 1));
    }
    
    return cp;
}

void pll_kill(polylist ** pll_pt)
{
    assert(pll_pt);
    
    polylist * pll= *pll_pt;
    assert(pll);
    
    if (pll->polys) free(pll->polys);
    
    free(pll);
    
    *pll_pt= NULL;
}

void pll_killall(polylist ** pll_pt)
{
    assert(pll_pt);
    
    polylist * pll= *pll_pt;
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

int pll_remove_pl(polylist * pll, const struct Polygon * pl)
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

struct Polygon * pll_remove(polylist * pll, const uint idx)
{
    assert(pll);
    if (idx < pll->last){
        struct Polygon * result= pl_copy(pll->polys[idx], NULL, 1);
        uint i=idx;
        struct Polygon ** pl_pt= pll->polys + i++;
        for (; i < pll->last; i++)
            *pl_pt++= pll->polys[i];
        pll->last--;
        return result;
    }
    return NULL;
}

struct Polygon * pll_replace(polylist * pll, const uint idx, struct Polygon * pl)
{
    assert(pll);
    assert(pl);
    struct Polygon * result= pll->polys[idx];
    pll->polys[idx]= pl;
    return result;
}

int pll_replace_pl(polylist * pll, const struct Polygon * to_replace, struct Polygon * pl)
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

void pll_extend(polylist ** pll_pt, struct Polygon ** polys, const uint count)
{
    assert(polys);
    
    polylist * pll= pll_check_pointer(pll_pt);
    
    if(pll->last + count > pll->size){
        pll_expand(pll, (pll->size * 2) + count);
    }
    
    uint i;
    for (i=0; i < count; i++){
        pll->polys[pll->last++]= polys[i];
	}
}
 
uint pll_count_points(polylist * pll)
{
    if(!pll) return 0;
    
    uint i, count=0;
    for (i=0; i < pll->last; i++)
        count+= pll->polys[i]->last;
    
    return count; 
}

int pl_cmp_area(const void *pl_pt_0, const void *pl_pt_1)
{
	double area_0= pl_area_3d(*((struct Polygon **)pl_pt_0), 1);
	double area_1= pl_area_3d(*((struct Polygon **)pl_pt_1), 1);
	if (area_0 < area_1)
		return 1;
	else if(area_0 > area_1)
		return 0;
	else
		return -1;
}
	
void pll_sort_by_area(polylist * pll)
{
	qsort(pll->polys, pll->last, sizeof(struct Polygon*), pl_cmp_area);
}

