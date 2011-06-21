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

void check_bounds(Bounds * b, const vec p)
{
    assert(b);
    assert(p);
    if(b->min[0] > p[0])
        b->min[0]= p[0];
    else if(b->max[0] < p[0])
        b->max[0]= p[0];

    if(b->min[1] > p[1])
        b->min[1]= p[1];
    else if(b->max[1] < p[1])
        b->max[1]= p[1];

    if(b->min[2] > p[2])
        b->min[2]= p[2];
    else if(b->max[2] < p[2])
        b->max[0]= p[2];
}

int disjoint_2d(const Bounds * b0, const Bounds * b1, AxisPair ax)
{
    switch (ax){
        case XY:
            return (b0->min[0] > b1->max[0] || b0->max[0] < b1->min[0] ||
                    b0->min[1] > b1->max[1] || b0->max[1] < b1->min[1]);
            break;
        case XZ:
            return (b0->min[0] > b1->max[0] || b0->max[0] < b1->min[0] ||
                    b0->min[2] > b1->max[2] || b0->max[2] < b1->min[2]);
            break;
        case YZ:
            return (b0->min[1] > b1->max[1] || b0->max[1] < b1->min[1] ||
                    b0->min[2] > b1->max[2] || b0->max[2] < b1->min[2]);
            break;
    }
    return 0;
}

int disjoint_3d(const Bounds * b0, const Bounds * b1)
{
    assert(b0);
    assert(b1);
    return !(((b0->min[0] < b1->min[0] && b1->min[0] < b0->max[0]) ||
              (b0->min[0] < b1->max[0] && b1->max[0] < b0->max[0]) ||
              (b1->min[0] < b0->min[0] && b0->min[0] < b1->max[0]) ||
              (b1->min[0] < b0->max[0] && b0->max[0] < b1->max[0])) &&
             ((b0->min[1] < b1->min[1] && b1->min[1] < b0->max[1]) ||
              (b0->min[1] < b1->max[1] && b1->max[1] < b0->max[1]) ||
              (b1->min[1] < b0->min[1] && b0->min[1] < b1->max[1]) ||
              (b1->min[1] < b0->max[1] && b0->max[1] < b1->max[1])) &&
             ((b0->min[2] < b1->min[2] && b1->min[2] < b0->max[2]) ||
              (b0->min[2] < b1->max[2] && b1->max[2] < b0->max[2]) ||
              (b1->min[2] < b0->min[2] && b0->min[2] < b1->max[2]) ||
              (b1->min[2] < b0->max[2] && b0->max[2] < b1->max[2])));
}

int point_in_bounds(const Bounds * b, const vec p)
{
    return p[0] > b->min[0] - ZEROLENGTH &&
           p[0] < b->max[0] + ZEROLENGTH &&
           p[1] > b->min[1] - ZEROLENGTH &&
           p[1] < b->max[1] + ZEROLENGTH &&
           p[2] > b->min[2] - ZEROLENGTH &&
           p[2] < b->max[2] + ZEROLENGTH;
}

int point_in_bounds_2d(const Bounds * b, const vec p, AxisPair ax)
{
    switch (ax){
        case XY:
            return !(p[0] < b->min[0] - ZEROLENGTH ||
                     p[0] > b->max[0] + ZEROLENGTH ||
                     p[1] < b->min[1] - ZEROLENGTH ||
                     p[1] > b->max[1] + ZEROLENGTH);
            break;
        case XZ:
            return !(p[0] < b->min[0] - ZEROLENGTH ||
                     p[0] > b->max[0] + ZEROLENGTH ||
                     p[2] < b->min[2] - ZEROLENGTH ||
                     p[2] > b->max[2] + ZEROLENGTH);
            break;
        case YZ:
            return !(p[1] < b->min[1] - ZEROLENGTH ||
                     p[1] > b->max[1] + ZEROLENGTH ||
                     p[2] < b->min[2] - ZEROLENGTH ||
                     p[2] > b->max[2] + ZEROLENGTH);
            break;
        default:
            return !(p[0] < b->min[0] - ZEROLENGTH ||
                     p[0] > b->max[0] + ZEROLENGTH ||
                     p[1] < b->min[1] - ZEROLENGTH ||
                     p[1] > b->max[1] + ZEROLENGTH);
            break;
    }
}

