//      bounds.h
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


#ifndef BOUNDS
#define BOUNDS

#include "globals.h"
#include "mathutils.h"

typedef struct bounds{
    vec min;
    vec max;
} Bounds;

void check_bounds(Bounds * b, const Vec * const point);
int disjoint_2d(const Bounds * b0, const Bounds * b1, AxisPair ax);
int disjoint_3d(const Bounds * b0, const Bounds * b1);
int point_in_bounds(const Bounds * b0, const Vec * const p);
int point_in_bounds_2d(const Bounds * b0, const Vec * const p, AxisPair ax);
#endif
