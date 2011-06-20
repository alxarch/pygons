//      bounds.c
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


#include <assert.h>
#include "globals.h"
#include "math.h"
#include "bounds.h"

void check_bounds(bounds * b, const vec * const p)
{
    assert(b);
    assert(p);
    if(b->min.x > p->x)
        b->min.x= p->x;
    else if(b->max.x < p->x)
        b->max.x= p->x;

    if(b->min.y > p->y)
        b->min.y= p->y;
    else if(b->max.y < p->y)
        b->max.y= p->y;

    if(b->min.z > p->z)
        b->min.z= p->z;
    else if(b->max.z < p->z)
        b->max.x= p->z;
}

int disjoint_2d(const bounds * b0, const bounds * b1, enum axis_pair ax)
{
    switch (ax){
        case XY:
            return (b0->min.x > b1->max.x || b0->max.x < b1->min.x ||
                    b0->min.y > b1->max.y || b0->max.y < b1->min.y);
            break;
        case XZ:
            return (b0->min.x > b1->max.x || b0->max.x < b1->min.x ||
                    b0->min.z > b1->max.z || b0->max.z < b1->min.z);
            break;
        case YZ:
            return (b0->min.y > b1->max.y || b0->max.y < b1->min.y ||
                    b0->min.z > b1->max.z || b0->max.z < b1->min.z);
            break;
    }
    return 0;
}

int disjoint_3d(const bounds * b0, const bounds * b1)
{
    assert(b0);
    assert(b1);
    return !(((b0->min.x < b1->min.x && b1->min.x < b0->max.x) ||
             (b0->min.x < b1->max.x && b1->max.x < b0->max.x) ||
             (b1->min.x < b0->min.x && b0->min.x < b1->max.x) ||
             (b1->min.x < b0->max.x && b0->max.x < b1->max.x)) &&
            ((b0->min.y < b1->min.y && b1->min.y < b0->max.y) ||
             (b0->min.y < b1->max.y && b1->max.y < b0->max.y) ||
             (b1->min.y < b0->min.y && b0->min.y < b1->max.y) ||
             (b1->min.y < b0->max.y && b0->max.y < b1->max.y)) &&
            ((b0->min.z < b1->min.z && b1->min.z < b0->max.z) ||
             (b0->min.z < b1->max.z && b1->max.z < b0->max.z) ||
             (b1->min.z < b0->min.z && b0->min.z < b1->max.z) ||
             (b1->min.z < b0->max.z && b0->max.z < b1->max.z)));
}

int point_in_bounds(const bounds * b, const vec * const p)
{
    return p->x > b->min.x - ZEROLENGTH &&
           p->x < b->max.x + ZEROLENGTH &&
           p->y > b->min.y - ZEROLENGTH &&
           p->y < b->max.y + ZEROLENGTH &&
           p->z > b->min.z - ZEROLENGTH &&
           p->z < b->max.z + ZEROLENGTH;
}

int point_in_bounds_2d(const bounds * b, const vec * const p, enum axis_pair ax)
{
    switch (ax){
        case XY:
            return !(p->x < b->min.x - ZEROLENGTH ||
                     p->x > b->max.x + ZEROLENGTH ||
                     p->y < b->min.y - ZEROLENGTH ||
                     p->y > b->max.y + ZEROLENGTH);
            break;
        case XZ:
            return !(p->x < b->min.x - ZEROLENGTH ||
                     p->x > b->max.x + ZEROLENGTH ||
                     p->z < b->min.z - ZEROLENGTH ||
                     p->z > b->max.z + ZEROLENGTH);
            break;
        case YZ:
            return !(p->y < b->min.y - ZEROLENGTH ||
                     p->y > b->max.y + ZEROLENGTH ||
                     p->z < b->min.z - ZEROLENGTH ||
                     p->z > b->max.z + ZEROLENGTH);
            break;
        default:
            return !(p->x < b->min.x - ZEROLENGTH ||
                     p->x > b->max.x + ZEROLENGTH ||
                     p->y < b->min.y - ZEROLENGTH ||
                     p->y > b->max.y + ZEROLENGTH);
            break;
    }
}

