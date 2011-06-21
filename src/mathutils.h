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

typedef struct vec4d{
	double x,y,z,w;
} V4d;

typedef struct vec{
	double x, y, z;
} Vec;

typedef struct Plane{
	struct Vec point;
	struct Vec normal;
} plane;

typedef double matrix[4][4];

// matrix ______________________________________________________________
	
void mx_copy(matrix dest, matrix src);

void mx_set_rotation( Matrix dest, const float rx, const float ry, const float rz);

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

void v_zero(Vec * v);

void v_add(Vec * dest, const Vec * const v0, const Vec * const v1);
void v_sub(Vec * dest, const Vec * v0, const Vec * const v1);
void v_scale(Vec * dest, const Vec * const v0, const double s);
void v_cross(Vec * dest, const Vec * const v0, const Vec * const v1);
double v_dot(const Vec * const v0, const Vec * const v1);
void v_mul_mx(Vec * dest, const Vec * v, matrix mx);

void v_normal(const Vec * const v, Vec * dest);

double v_length(const Vec * const v);
double dist(const Vec * const v0, const Vec * const v1);

double v_perp(const Vec * const v0, const Vec * const v1, enum axis_pair ax);

void midpoint(Vec * dest, const Vec * const v0, const Vec * const v1);

int same_2d( const Vec * const v0, const Vec * const v1, enum axis_pair ax);
int same_3d( const Vec * const v0, const Vec * const v1);

// vector copy methods _________________________________________________

void v_copy(Vec * dest, const Vec * const src);

void v_swp_xy(Vec * dest, const Vec * const src);
                    
void v_swp_xz(Vec * dest, const Vec * const src);
                    
void v_swp_yz(Vec * dest, const Vec * const src);
                    
void v_flat_x(Vec * dest, const Vec * const src);
                    
void v_flat_y(Vec * dest, const Vec * const src);
                    
void v_flat_z(Vec * dest, const Vec * const src);

void v_flat_persp(Vec * dest, const Vec * const src);

void v_print(Vec v);

float v_angle(const Vec * const v0, const Vec * const v1);
float v_angle_r(const Vec * const v0, const Vec * const v1);

// other _______________________________________________________________ 

int
rel_point_plane(const Vec * const v, const plane * const pn);

int
point_in_segment_2d(const Vec * const p,const Vec * const p0, const Vec * const p1,  const enum axis_pair ax); 

int
line_intersect_plane(const Vec * const v0, const Vec * const v1, const plane * const pn, Vec * dest);

#endif
