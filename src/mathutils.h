//      mathutils.h
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


#ifndef MATHUTILS
#define MATHUTILS

#include <math.h>
#include "globals.h"

#define RAD2DEG 0.0174532925199432954743716805978692L
#define DEG2RAD 57.295779513082322864647721871733665L

typedef struct Vec4d{
	double x,y,z,w;
} v4d;

typedef struct Vec{
	double x, y, z;
} vec;

typedef struct Plane{
	struct Vec point;
	struct Vec normal;
} plane;

typedef double matrix[4][4];

// matrix ______________________________________________________________
	
void mx_copy(matrix dest, matrix src);

void mx_set_rotation( matrix dest, const float rx, const float ry, const float rz);

int mx_invert(matrix mx);

void mx_transpose(matrix mx);

void mx_set_identity(matrix dest);

void mx_set_translation(matrix dest, const double dx, const double dy, const double dz);

void mx_set_rot_x(matrix mx, const float rx);

void mx_set_rot_y(matrix mx, const float ry);

void mx_set_rot_z(matrix mx, const float rz);

void mx_mul_scalar(matrix mx, matrix src, const double s);

void mx_mul(matrix dest, matrix mx_a, matrix mx_b);

// vector ______________________________________________________________

void v_zero(vec * v);

void v_add(vec * dest, const vec * const v0, const vec * const v1);
void v_sub(vec * dest, const vec * v0, const vec * const v1);
void v_scale(vec * dest, const vec * const v0, const double s);
void v_cross(vec * dest, const vec * const v0, const vec * const v1);
double v_dot(const vec * const v0, const vec * const v1);
void v_mul_mx(vec * dest, const vec * v, matrix mx);

void v_normal(const vec * const v, vec * dest);

double v_length(const vec * const v);
double dist(const vec * const v0, const vec * const v1);

double v_perp(const vec * const v0, const vec * const v1, enum axis_pair ax);

void midpoint(vec * dest, const vec * const v0, const vec * const v1);

int same_2d( const vec * const v0, const vec * const v1, enum axis_pair ax);
int same_3d( const vec * const v0, const vec * const v1);

// vector copy methods _________________________________________________

void v_copy(vec * dest, const vec * const src);

void v_swp_xy(vec * dest, const vec * const src);
                    
void v_swp_xz(vec * dest, const vec * const src);
                    
void v_swp_yz(vec * dest, const vec * const src);
                    
void v_flat_x(vec * dest, const vec * const src);
                    
void v_flat_y(vec * dest, const vec * const src);
                    
void v_flat_z(vec * dest, const vec * const src);

void v_flat_persp(vec * dest, const vec * const src);

void v_print(vec v);

float v_angle(const vec * const v0, const vec * const v1);
float v_angle_r(const vec * const v0, const vec * const v1);

// other _______________________________________________________________ 

int
rel_point_plane(const vec * const v, const plane * const pn);

int
point_in_segment_2d(const vec * const p,const vec * const p0, const vec * const p1,  const enum axis_pair ax); 

int
line_intersect_plane(const vec * const v0, const vec * const v1, const plane * const pn, vec * dest);

#endif
