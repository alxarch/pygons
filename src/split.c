//      split.c
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
#include "math.h"
#include "bounds.h"
#include "polygon.h"
#include "split.h"

SplitGraph * build_split_graph(Polygon * pl, Plane * pln)
{
	SplitGraph * sg = (SplitGraph *)malloc(sizeof(SplitGraph));
	assert(sg);
	
	sg->last = pl->last;
	
	sg->nodes = (SplitNode *)malloc(2 * sg->last * sizeof(SplitNode));
	assert(sg->nodes);
	sg->r = 0;
	sg->pl = pl;
	sg->pln = pln;
	
	SplitNode * nd = sg->nodes;
	PolyVert * vx;
	uint i, j;
	for(i=(sg->last-1), j=0; j < sg->last; i=j++){
		vx = pl->points + i;
		v_cpy(nd->co, vx->co);
		nd->flags = vx->flags;

		nd->nxt = sg->nodes + j;
		nd->prv = i ? sg->nodes + i - 1 : sg->nodes + sg->last - 1;
		nd->r = rel_point_plane(nd->co, pln);

		sg->r |= nd->r;

		nd++;
	}
	
	if(sg->r&CROSS){
		SplitNode * newnode;
		nd = sg->nodes;
		while(nd->nxt != sg->nodes){
			if((nd->r | nd->nxt->r) == CROSS){
				newnode = sg->nodes + sg->last++;
				line_intersect_plane(nd->co, nd->nxt->co, pln, newnode->co);
				newnode->flags = nd->flags;
				newnode->r = ON;
				newnode->nxt = nd->nxt;
				nd->nxt = newnode;
				newnode->prv = nd;
				newnode->nxt->prv = newnode;
				
				sg->r |= ON;
			}
		
			nd = nd->nxt;
		}
	}
	
	return sg;
}


	
static Polygon * 
new_split(SplitNode * nd)
{
	Polygon * split = pl_init(nd->pl->last);
	split->flags = nd->pl->flags | SPLIT;
	pl_append(&split, nd->co, nd->flags);
	
	return split;
}

/*
 * Frees all memory occupied by a SplitGraph struct.
 */
void 
kill_split_graph(SplitGraph ** sg_pt)
{
	assert(sg_pt);
	
	SplitGraph * sg = *sg_pt;
	
	if(sg){
		if(sg->nodes) free(sg->nodes);
		free(sg);
	}

}

/*
 * Helper struct used for sorting polygon points according to distance.
 */
struct Dist {
    double distance;
    uint idx;
};

/*
 * Helper function that compares Dist structs used for sorting.
 */
static int 
cmp_dists(const void* dst_a, const void* dst_b)
{
    double d_a = ((struct Dist *)dst_a)->distance;
    double d_b = ((struct Dist *)dst_b)->distance;
    if (d_a < d_b)
        return -1;
    else if (d_a > d_b)
        return 1;
    else
        return 0;
}
            
static void 
join_holes(Polygon * pl, PolyList * holes)
{
	
	vec ref;
	v_cpy(ref, pl->points[pl->last - 1].co);
	uint i,limit= holes->last;
	struct Dist dists[limit];
	Polygon * hole = NULL;
	for (i=0; i < limit; i++){
		hole= holes->polys[i];
		
		dists[i].idx= i;
		
		// Store the distance between the first point of the hole and the last point of the polygon.
		dists[i].distance= dist(hole->points[0].co, ref);
	}
	
	qsort(&dists, limit, sizeof(struct Dist), cmp_dists);

	// Add all holes points to the part in order of distance.
	for (i=0; i < limit; i++){
		hole= holes->polys[dists[i].idx];
		pl_extend_verts(pl, hole->points, hole->last);
	}
	
	// Points are now safely transfered to part
	pll_killall(&holes); 
}
                

/*
 * Proccesses all splits of a side of the plane.
 * 
 * Its purpose is to find out which of the splits are holes and to 
 * distribute them to the correct polygon.
 * 
 */
static void 
process_splits(PolyList **splits_pt, enum PL_Order orig, AxisPair ax)
{
    assert(splits_pt);
    
    PolyList * splits= *splits_pt;
    assert(splits);
    
    PolyList * append = NULL,
		     * holes = pll_init(splits->last),
		     * parts = pll_init(holes->size);
				      
    Polygon * split;
    uint i, j;
    int r;
    for (i=0; i < splits->last; i++){
        split= splits->polys[i];
		
		// Remove defective polygons.
        if (split->last < 3){
            pl_kill(&split);
            continue;
        }
        
        // If a polygon's order is opposite to the original order it's a hole
        append = (pl_order_2d(split, ax) == orig) ? parts : holes;
        pll_append(&append, split);
    }

	
    if (holes->last){
        PolyList *leftovers =  NULL, *holes2join = NULL;
        Polygon *pl = NULL, *hole = NULL;
		
		// Iterate over each part and find which holes belong to it.
        for (i=0; i < parts->last; i++){
			
            pl= parts->polys[i];
            
            holes2join= leftovers= NULL;
            
            for (j=0; j < holes->last; j++){
				
                hole = holes->polys[j];
                r = pl_rel_point_2d(pl, hole->points[0].co, ax, 1);
                if (r & OUT){
					// If any point of the hole is outside the part it 
					// belongs to another part.
                    pll_append(&leftovers, hole);
				}
                else if (hole->flags & SPLIT){
					// If the hole was the result of a split it needs 
					// to be joined to the outline of the part.
                    pll_append(&holes2join, hole);
				}
                else{
					// This is a hole of this part.
                    pl_add_hole(pl, hole);
				}
            }
            
            free(holes->polys);
            free(holes);
            holes= leftovers;
            
            if (holes2join) join_holes(pl, holes2join);
            
            if (!holes) break; // No more holes to distribute.
            
        }
    }
    
    pll_kill(splits_pt);
    
    *splits_pt= parts;
}
     

void 
split_graph(SplitGraph *sg, PolyList * fparts, PolyList * bparts)
{
	uint maxnodes = sg->last;
	Polygon * split = pl_init(maxnodes);
	split->flags = sg->pl->flags | SPLIT;
	
	
	SplitNode * nd = sg->nodes;
	int fcross, bcross, R;
	bcross = nd->r == BK || nd->r == ON;
	fcross = !bcross;
	
	PolyList * fsplits = pll_init(maxnodes),
		     * bsplits = pll_init(maxnodes),
		     * poplast = fcross ? fsplits : bsplits;
	
	
	while(nd->nxt != sg->nodes){
		pl_append(&split, nd->co, nd->flags);
		
		if (nd->r == ON){
			R = nd->nxt->r | nd->prv->r;
			switch(R){
				case ON:
					// Both previous, this and next points are on the plane.
					continue; 
					break;
				case BK:
				case ONBK:
					if (++bcross % 2){
						pll_append(&bsplits, split);
						split = new_split(nd);
						if (R == BK) bcross++;
					}
					break;
				case FR:
				case ONFR:
					if (++fcross % 2){
						pll_append(&fsplits, split);
						split= new_split(nd);
						if (R == FR) fcross++;
					}
					break;
			}
			
		}
		
		nd = nd->nxt;
	}
	
	if (fcross % 2){
		pll_append(&fsplits, split);
	}
	else{
		pll_append(&bsplits, split);
	}
	
	if (poplast->last > 1){
		Polygon * tojoin= poplast->polys[--poplast->last];
		
		pl_extend_verts(poplast->polys[0], tojoin->points, tojoin->last);
		
		pl_kill(&tojoin);
	}
	
	if (sg->pl->holes){
		uint i;
		for (i=0; i < sg->pl->holes->last; i++){
			pl_split_by_plane(sg->pl->holes->polys[i], sg->pln, fsplits, bsplits);
		}
	}
		
	AxisPair ax = pl_best_axis_pair(sg->pl);
	enum PL_Order ord = pl_order_2d(sg->pl, ax);
	
	process_splits(&fsplits, ord, ax);
	process_splits(&bsplits, ord, ax);
	
	// Put splits back to the appropriate parts lists
	pll_extend(&fparts, fsplits->polys, fsplits->last);
	pll_kill(&fsplits);
	
	pll_extend(&bparts, bsplits->polys, bsplits->last);
	pll_kill(&bsplits);
	
}
                
                
                

