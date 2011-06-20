//      boolean.c
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
//      boolean.c
//      

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "math.h"
#include "globals.h"
#include "polygon.h"
#include "boolean.h"

bpoly * bl_parse_poly(bgraph * g, polygon * pl, enum BL_Set set, int hole)
{
    assert(pl->last > 2);
    unsigned int i;
    
    if (!hole && pl->holes){
		for (i=0; i < pl->holes->last; i++){
          bl_parse_poly(g, pl->holes->polys[i], set, 1);
		}
    }
    bpoly * result= g->polys + g->pl_num++;
    bnode * new_node, * offset= g->nodes + g->last;
    polyvert * vx;
    for (i=0; i < pl->last; i++){
        vx= pl->points + i;
        new_node= offset + i;
        new_node->nxt= new_node + 1;
        new_node->prv= new_node - 1;
        //new_node->isx= NULL;
        //new_node->t= NOTYPE;
        new_node->co= vx->co;
        new_node->flags= vx->flags;
        new_node->pl= result;
    }
    g->last+= i;
    new_node->nxt= offset;
    offset->nxt= new_node;
    result->set= set;
    result->hole= hole;
    return result;
}

int bl_expand_graph(bgraph * g, unsigned int count)
{
    unsigned int limit= g->last;
    unsigned int nxt_diffs[limit];
    unsigned int prv_diffs[limit];
    unsigned int isx_diffs[limit];
    unsigned int i;
    bnode * nd, * old_base= g->nodes;
    
    for (i=0; i < limit; i++){
        nd= g->nodes + i;
        nxt_diffs[i]= nd->nxt - g->nodes;
        prv_diffs[i]= nd->prv - g->nodes;
        isx_diffs[i]= (nd->isx) ? nd->isx - g->nodes : 0;
    }
    
    g->size+= count;
    
    g->nodes= (bnode *)realloc(g->nodes, g->size * sizeof(bnode));
    assert(g->nodes);
    
    if (g->nodes != old_base){
        for(i=0; i < limit; i++){
            nd= g->nodes + i;
            nd->nxt= g->nodes + nxt_diffs[i];
            nd->prv= g->nodes + prv_diffs[i];
            nd->isx= (nd->isx) ? g->nodes + isx_diffs[i] : NULL;
        }
        return 1;
    }
    else
        return 0;
}

unsigned int bl_find_intersections(bgraph * g)
{
    bnode *a0, *a1, *b0, *b1, *isect_a=NULL, *isect_b=NULL;
    vec u, w, v, I;
    double ratio_a, ratio_b, d;
    unsigned int i, j;
    
    for(i=0; i < g->last; i++){
        
        a0= g->nodes + i;
        
        if (a0->pl->set == B) continue;
        
        a1= a0->nxt;
        
        v_sub(&u, &a1->co, &a0->co);
        
        for(j=0; j < g->last; j++){
            
            b0= g->nodes + j;
            
            if (b0->pl->set == A) continue;
            
            b1= b0->nxt;
            
            if (b0->isx == a0 || b0->isx == a1 ||
                b1->isx == a0 || b1->isx == a1)
                   continue;
            
            v_sub(&v, &b1->co, &b0->co);
            
            d= v_perp(&u, &v, XY);
            
            if (fabs(d) < E) continue;
            
            v_sub(&w, &a0->co, &b0->co);
            
            ratio_a= v_perp(&v, &w, XY) / d;
            
            if (-0.1 < ratio_a && ratio_a < 1.1){
				
                ratio_b= v_perp(&u, &w, XY) / d;
                
                if (-0.1 < ratio_b && ratio_b < 1.1){
                    v_scale(&I, &u, ratio_a);
                    v_add(&I, &I, &a0->co);
                }
                else{
					continue;
				}
            }
            else {
				continue;
			}
    
            if (g->last + 2 > g->size){
                if (bl_expand_graph(g, g->size * 2)){
                    a0= g->nodes + i;
                    a1= a0->nxt;
                    b0= g->nodes + j;
                    b1= b0->nxt;
                }
            }
            
            isect_a= isect_b= NULL;
            
            if (same_2d(&a0->co, &I, XY)){
                if (!a0->isx){
                    isect_a= a0;
				}
                else if (ratio_a < 0.0){
                    continue;
				}
            }
            else if (same_2d(&a1->co, &I, XY)){
                if (!a1->isx){
                    isect_a= a1;
				}
                else if (ratio_a > 1.0){
                    continue;
				}
            }
            else if (!point_in_segment_2d(&I, &a0->co, &a1->co, XY)){
                continue;
			}
            
            if (same_2d(&b0->co, &I, XY)){
                if (!b0->isx){
                    isect_b= b0;
                }
                else if (ratio_b < 0.0){
					continue;
				}
            }
            else if (same_2d(&b1->co, &I, XY)){
                if (!b1->isx)
                    isect_b= b1;
                else if (ratio_b > 1.0)
                    continue;
            }
            else if (!point_in_segment_2d(&I, &b0->co, &b1->co, XY))
                continue;
            
            if (!isect_a){
                isect_a= g->nodes + g->last++;
                isect_a->nxt= a1;
                isect_a->prv= a0;
                a0->nxt= isect_a;
                a1->prv= isect_a;
                isect_a->co= I;
                isect_a->flags= a0->flags;
                isect_a->pl= a0->pl;
                //isect_a->t= NOTYPE;
            
                v_sub(&u, &I, &a0->co);
                a1= isect_a;
            }
            if (!isect_b){
                isect_b= g->nodes + g->last++;
                isect_b->nxt= b1;
                isect_b->prv= b0;
                isect_b->co= I;
                isect_b->flags= b0->flags;
                isect_b->pl= b0->pl;
                b0->nxt= isect_b;
                b1->prv= isect_b;
                //isect_b->t= NOTYPE;
            }
            
            isect_a->isx= isect_b;
            isect_b->isx= isect_a;
            g->inum++;
        }
    }    
    return g->inum;
}

void bl_classify_intersections(bgraph * g)
{
    unsigned int i;
    bnode * nd;
    vec mid;
    int rnext, rprev, *R;
    bpoly * other;
        
    for (i=0; i < g->last; i++){
        nd= g->nodes + i;
        
        if (!nd->isx) continue;
        
        other= nd->isx->pl;
        if (nd->prv->isx &&
            ((nd->prv->isx == nd->isx->prv) ||
             (nd->prv->isx == nd->isx->nxt)))
                rprev= ON;
        else{
            midpoint(&mid, &nd->co, &nd->prv->co);
            rprev= pl_rel_point_2d(other->pl, &mid, XY, 0);
            if (other->hole){
                if (rprev == OUT)
                    rprev = IN;
                else if(rprev == IN)
                    rprev = OUT;
            }
        }
        
        if (nd->nxt->isx &&
            ((nd->nxt->isx == nd->isx->prv) ||
             (nd->nxt->isx == nd->isx->nxt)))
                rnext= ON;
        else{
            midpoint(&mid, &nd->co, &nd->nxt->co);
            rnext= pl_rel_point_2d(other->pl, &mid, XY, 0);
            if (other->hole){
                if (rnext == OUT)
                    rnext = IN;
                else if(rnext == IN)
                    rnext = OUT;
            }
        }
        
        R= &nd->pl->rel;
        *R|= rprev | rnext;
        if (nd->t!= NOTYPE)
            continue;
        else if (rprev > rnext){
            nd->t= EX;
            nd->isx->t= EN;
        }
        else if (rprev < rnext){
            nd->t= EN;
            nd->isx->t= EX;
        }
        else if (rprev == ON){
            nd->t= INVALID;
            nd->isx->t= INVALID;
            g->inum--;
        }
        else if (rprev == IN){
            nd->t= EXEN;
            nd->isx->t= ENEX;
        }
        else{
            nd->t= ENEX;
        }
    }
}

bgraph * bl_build_graph(polygon * pl_a, polygon * pl_b)
{
    bgraph * g= (bgraph*)calloc(1, sizeof(bgraph));
    assert(g);
    
    // Allocate enough nodes so no expanding is needed in most cases.
    g->size= (pl_a->last + pl_b->last + pll_count_points(pl_a->holes) + pll_count_points(pl_b->holes)) * 4;
    g->nodes= (struct BL_Node*)calloc(g->size, sizeof(struct BL_Node));
    assert(g->nodes);

	
    unsigned int total_polys= 2;
    total_polys+= (pl_a->holes) ? pl_a->holes->last : 0;
    total_polys+= (pl_b->holes) ? pl_b->holes->last : 0;
    g->polys= (bpoly *)calloc(total_polys, sizeof(bpoly)); 
	assert(g->polys);

    g->pl_a= bl_parse_poly(g, pl_a, A, 0);
    g->pl_b= bl_parse_poly(g, pl_b, B, 0);
    if (bl_find_intersections(g)){
        bl_classify_intersections(g);
    }
    else{
        g->pl_a->rel= pl_rel_point_2d(pl_b, &pl_a->points[0].co, XY, 0);
        g->pl_a->rel= pl_rel_point_2d(pl_a, &pl_b->points[0].co, XY, 0);
    }
    
    return g;
}

polylist * bl_handle_out_in(polygon * pl_out, polygon * pl_in){
    return pll_init(10);
}

void bl_special_case(bgraph * g, enum BL_Op op, polylist ** parts) 
{
    int a= g->pl_a->rel;
    int b= g->pl_b->rel;
    if ((a == OUT) && (b == OUT)){
        switch (op){
            case UNION:
                pll_append(parts, pl_copy(g->pl_a->pl, NULL, 1));
                pll_append(parts, pl_copy(g->pl_b->pl, NULL, 1));
                break;
            case INTERSECT:
                *parts= pll_init(2);
                break;
            case SUBTRACT_B_A:
                pll_append(parts, pl_copy(g->pl_b->pl, NULL, 1));
                break;
            case SUBTRACT_A_B:
                pll_append(parts, pl_copy(g->pl_a->pl, NULL, 1));
                break;
        }
    }
    else if ((A == ONOUT) && (B == ONOUT)){
        switch (op){
            case UNION:
                *parts= pll_init(2);
                break;
            case INTERSECT:
                *parts= pll_init(2);
                break;
            case SUBTRACT_B_A:
                pll_append(parts, pl_copy(g->pl_b->pl, NULL, 1));
                break;
            case SUBTRACT_A_B:
                pll_append(parts, pl_copy(g->pl_a->pl, NULL, 1));
                break;
        }
    }
    else if (A == OUT && B == IN){
        switch (op){
            case UNION:
                pll_append(parts, pl_copy(g->pl_a->pl, NULL, 1));
                break;
            case INTERSECT:
                pll_append(parts, pl_copy(g->pl_b->pl, NULL, 1));
                break;
            case SUBTRACT_A_B:
                *parts= bl_handle_out_in(g->pl_a->pl, g->pl_b->pl);
                break;
            case SUBTRACT_B_A:
                *parts= pll_init(2);
                break;
        }
    }
    else if ((A == IN) && (B == OUT)){
        switch (op){
            case UNION:
                pll_append(parts, pl_copy(g->pl_a->pl, NULL, 1));
                break;
            case INTERSECT:
                pll_append(parts, pl_copy(g->pl_b->pl, NULL, 1));
                break;
            case SUBTRACT_A_B:
                *parts= pll_init(2);
                break;
            case SUBTRACT_B_A:
                *parts= bl_handle_out_in(g->pl_b->pl, g->pl_a->pl);
                break;
        }
    }
    else if ((A == ONIN) && (B == ONOUT)){
        switch (op){
            case UNION:
                pll_append(parts, pl_copy(g->pl_b->pl, NULL, 1));
                break;
            case INTERSECT:
                pll_append(parts, pl_copy(g->pl_a->pl, NULL, 1));
                break;
            case SUBTRACT_A_B:
                *parts= pll_init(2);
                break;
            case SUBTRACT_B_A:
                break;
        }
    }
    else if ((A == ONOUT) && (B == ONIN)){
        switch (op){
            case UNION:
                pll_append(parts, pl_copy(g->pl_a->pl, NULL, 1));
                break;
            case INTERSECT:
                pll_append(parts, pl_copy(g->pl_b->pl, NULL, 1));
                break;
            case SUBTRACT_A_B:
                break;
            case SUBTRACT_B_A:
                *parts= pll_init(2);
                break;
        }
    }
}

/*
 * Checks if a node can be used as an entry point for a given boolean 
 * operation.
 * 
 * @param bnode nd The node to check.
 * @param BL_Op op The operation for which to check if node is entry point.
 * 
 * @return int
 */
int bl_entry_point(bnode * nd, enum BL_Op op)
{
    if (!nd->isx)
        return 0;
    else if ((nd->pl->set == A) || (op == SUBTRACT_B_A)){
        if (op == INTERSECT){
            if (nd->t== EN)
                nd->pl->has_entry= nd->entry= 1;
            else if (nd->t== ENEX && nd->isx->t== ENEX)
                nd->pl->has_entry= nd->entry= 2;
            else
                return 0;
        }
        else if (op == UNION){
            if (nd->t== EX && 
                nd->nxt->t!= EX && 
                nd->nxt->t!= INVALID)
                    nd->pl->has_entry= nd->entry= 1;
            else if (nd->t== ENEX && nd->isx->t== ENEX)
                    nd->pl->has_entry= nd->entry= 2;
            else
                return 0;
        }
        else{
            if (nd->t== EX &&
                nd->nxt->t!= EX &&
                nd->nxt->t!= INVALID)
                    nd->pl->has_entry= nd->entry= 1;
            else if (((nd->t== ENEX && !nd->pl->hole) ||
                      (nd->t== EXEN && nd->pl->hole)) &&
                     ((nd->isx->t== EXEN && !nd->isx->pl->hole) ||
                      (nd->isx->t== ENEX && nd->isx->pl->hole)))
                nd->pl->has_entry= nd->entry= 2;
            else
                return 0;
        }
    }
    else{
        return 0;
    }
    return 1;
}


/*
 * Executes boolean operation op.
 * 
 * @param BL_Graph* 	g  	The graph to conduct the operation on.
 * @param enum BL_Op    op 	The operation to be performed. 
 * 
 * @return Polylist* Pointer to a polylist struct of resulting polys.
 */
polylist * bl_operation(bgraph * g, enum BL_Op op)
{
    assert(g);
    
    polylist * parts= NULL;
    
    bl_special_case(g, op, &parts);
    
    //check special case
    if (parts) return parts;
    
    polylist * holes= NULL;
    polyvert * vx;
    bnode * nd;
    polygon * part;
    unsigned int entry_points= 0;
    unsigned int limit=g->last;
    unsigned int i, j;
    //reset entry points
    for (i=0; i < limit; i++)
        entry_points+= bl_entry_point(g->nodes + i, op);
    //printf("entrypoints found: %d\n", entry_points);

    if (entry_points){
        int lt;
        bnode * start;
        parts= pll_init(entry_points);
        holes= pll_init(entry_points);
        
        for (i=0; i < limit; i++){
            nd= g->nodes + i;
            if (nd->entry > 0){
                // printf("walking... start is: %d\n", i);
                part= pl_init(limit);
                start= nd;
                j= 0;
                lt= NOTYPE;
                do{
                    if (++j > limit){
                        // sth went wrong, endless walk
                        // printf("endless walk, num of steps: %d\n", j);
                        break;
                    }
                    vx= pl_append(&part, &nd->co, nd->flags);

                    if (nd->isx){
                        nd->entry= --nd->isx->entry;
                        if ((op == UNION && (nd->t== EN || (nd->t== ENEX && nd->isx->t== ENEX)))
                                ||
                            (op == INTERSECT && nd->t== EX)
                                ||
                            (nd->t!= lt && !((nd->t== EXEN) && ((nd->isx->t== ENEX) ||
                                                                ((nd->pl->set==A) ? (op==SUBTRACT_B_A) : (op==SUBTRACT_A_B))
                                                               )
                                            )
                            )
                           )
                            nd= nd->isx;
                        lt= nd->t;
                    }
                    if ((nd->pl->set == A) ? (op == SUBTRACT_B_A) : (op == SUBTRACT_A_B)){
                        nd= nd->prv;
                        vx->flags |= HIDDEN;
                    }
                    else{
                        nd= nd->nxt;
					}
                }while (nd != start && nd != start->isx);
                
                // printf("part has %d points\n", part->last);
                
                if (part->last < 3){
                    //invalid result
                    // printf("Invalid part, part has %d points\n", part->last);
                    pl_kill(&part);
                }
                else{
                    // printf("valid part, part has %d points\n", part->last);
                    if (pl_order(part) == CW){
                        pll_append(&parts, part);
                        printf("part is normal\n");
                    }
                    else{
                        printf("part is hole\n");
                        pll_append(&holes, part);
                    }
                }
            }
        }
        // printf("found %d parts\n", parts->last);
        // printf("found %d holes\n", holes->last);
    }
    //handle polygons that had no entry points
    bpoly * bpl;
    for (i=0; i < g->pl_num; i++){
        bpl= g->polys + i;
        if (!bpl->has_entry){
            part= NULL;
            switch (op){
                case UNION:
                    if (bpl->rel == OUT)
                        part= bpl->pl;
                    break;
                case SUBTRACT_A_B:
                    if (((bpl->set == B) && (bpl->rel == IN)) ||
                        ((bpl->set == A) && (bpl->rel == OUT)))
                        part= bpl->pl;
                    break;
                case SUBTRACT_B_A:
                    if (((bpl->set == A) && (bpl->rel == IN)) ||
                        ((bpl->set == B) && (bpl->rel == OUT)))
                        part= bpl->pl;
                    break;
                case INTERSECT:
                    if (bpl->rel == IN)
                        part= bpl->pl;
                    break;
            }
            if (part){
                if (bpl->hole)
                    pll_append(&holes, part);
                else
                    pll_append(&parts, part);
            }
        }
    }
    
    polygon * hole, * parthole;
    int c;
    unsigned int k;
    for (i=0; i < holes->last; i++){
        hole= holes->polys[i];
        for (j=0; j < parts->last; j++){
            part= parts->polys[j];
            //get some classification other than on if possible;
            c= ON;
            for (k=0; c == ON && k < hole->last; k++)
                c= pl_rel_point_2d(part, &hole->points[k].co, XY, 1); 
            if (c==IN){
                polylist * newholes= pll_init(part->holes->size + holes->size);
                for (k=0; k < part->holes->last; k++){
                    parthole= part->holes->polys[k];
                    if (pl_rel_point_2d(hole, &parthole->points[0].co, XY, 0) != IN)
                        pll_append(&newholes, parthole);
                    else
                        pl_kill(&parthole);
                }
                pll_append(&newholes, hole);
                pll_kill(&part->holes);
                part->holes= newholes;
                holes->polys[i]= NULL;
                break;
            }
        }
    }
    // any holes remaining need to be converted to "real" polys
    for (i=0; i < holes->last; i++){
        hole= holes->polys[i];
        if (hole){
            pl_reverse(hole);
            pll_append(&parts, hole);
        }
    }
    pll_kill(&holes);
    return parts;
}

/*
 * Frees graph allocated memory.
 */
void bl_kill_graph(bgraph * g)
{
    free(g->nodes);
    free(g->polys);
    free(g);
}
