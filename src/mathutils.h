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

typedef double vec[3];
typedef double matrix[4][4];

typedef struct plane{
	vec point, normal;
} Plane;

typedef enum axis_pair{ XY, XZ, YZ} AxisPair;

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

void v_zero(vec v);

void v_add(vec dest, const vec v0, const vec v1);
void v_sub(vec dest, const vec v0, const vec v1);
void v_scale(vec dest, const vec v0, const double s);
void v_cross(vec dest, const vec v0, const vec v1);
double v_dot(const vec v0, const vec v1);
void v_mul_mx(vec dest, const vec v, matrix mx);

void v_normal(const vec v, vec dest);

double v_length(const vec v);
double dist(const vec v0, const vec v1);

double v_perp(const vec v0, const vec v1, AxisPair ax);

void midpoint(vec dest, const vec v0, const vec v1);

int same_2d( const vec v0, const vec v1, AxisPair ax);
int same_3d( const vec v0, const vec v1);

// vector copy methods _________________________________________________

void v_cpy(vec dest, const vec src);

void v_swp_xy(vec dest, const vec src);
                    
void v_swp_xz(vec dest, const vec src);
                    
void v_swp_yz(vec dest, const vec src);
                    
void v_flat_x(vec dest, const vec src);
                    
void v_flat_y(vec dest, const vec src);
                    
void v_flat_z(vec dest, const vec src);

void v_flat_persp(vec dest, const vec src);

void v_print(vec v);

float v_angle(const vec v0, const vec v1);
float v_angle_r(const vec v0, const vec v1);

// other _______________________________________________________________ 

int
rel_point_plane(const vec v, const Plane * pn);

int
point_in_segment_2d(const vec p,const vec p0, const vec p1,  const AxisPair ax); 

int
line_intersect_plane(const vec v0, const vec v1, const Plane * const pn, vec dest);

#endif
