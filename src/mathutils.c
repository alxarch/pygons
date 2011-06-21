//      mathutils.c
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


#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "globals.h"
#include "mathutils.h"

// matrix ______________________________________________________________  


void mx_copy(matrix dst_mx, matrix src_mx)
{
    memcpy(dst_mx, src_mx, sizeof(matrix));
}

void mx_set_identity(matrix mx){

                mx[0][1]=   mx[0][2]=   mx[0][3] = 
    mx[1][0]=               mx[1][2]=   mx[1][3] = 
    mx[2][0]=   mx[2][1]=               mx[2][3] = 
    mx[3][0]=   mx[3][1]=   mx[3][2]=               0.0;

    mx[0][0]=
                mx[1][1]=   
                            mx[2][2]=
                                        mx[3][3]=   1.0;
}

void mx_set_rot_x(matrix mx, const float rot_x)
{
    double sin_d= sin(rot_x);
    double cos_d= cos(rot_x);
    mx_set_identity(mx);
    mx[1][2]= sin_d;
    mx[2][2]= mx[1][1]= cos_d;
    mx[2][1]= -sin_d;
}

/*
 *                            | cos  0 -sin  0 |
 *  Matrix - Set Y Rotation   |  0   1   0   0 |
 *                            | sin  0  cos  0 |
 *                            |  0   0   0   1 |
 */
void mx_set_rot_y(matrix mx, const float rot_y)
{
    double sin_d= sin(rot_y);
    double cos_d= cos(rot_y);
    mx_set_identity(mx);
    mx[2][0]= sin_d;
    mx[0][0]= mx[2][2]= cos_d;
    mx[0][2]= -sin_d;
}   

/*
 *                            |  cos sin  0   0 |
 *  Matrix - Set Z Rotation   | -sin cos  0   0 |
 *                            |   0   0   1   0 |
 *                            |   0   0   0   1 |
 */
void mx_set_rot_z(matrix mx, const float rot_z)
{
    double sin_d= sin(rot_z);
    double cos_d= cos(rot_z);
    mx_set_identity(mx);
    mx[0][1]= sin_d;
    mx[0][0]= mx[1][1]= cos_d;
    mx[1][0]= -sin_d;
}

void mx_set_rotation(matrix mx, const float rot_x, const float rot_y, const float rot_z)
{
    matrix tmp;
    mx_set_rot_x(mx, rot_x);
    mx_set_rot_y(tmp, rot_y);
    mx_mul(mx, mx, tmp);
    mx_set_rot_z(tmp, rot_z);
    mx_mul(mx, mx, tmp);
}

void mx_set_translation(matrix mx, const double dx, const double dy, const double dz)
{
    mx_set_identity(mx);
    mx[3][0]= dx;
    mx[3][1]= dx;
    mx[3][2]= dx;
}


/*
 * |a0 a1 a2 a3|    |a0 b0 c0 d0|
 * |b0 b1 b2 b3| => |a1 b1 c1 d1| 
 * |c0 c1 c2 c3|    |a2 b2 c2 d2|
 * |d0 d1 d2 d3|    |a3 b3 c3 d3|
 * 
 * 
 */
void mx_transpose(matrix mx)
{
    
    // Use a temp buffer
    double tmp[6];
    
    // Diagonal
    //dmx[0][0]=amx[0][0]; dmx[1][1]=amx[1][1]; dmx[2][2]=amx[2][2]; dmx[3][3]=amx[3][3];

    tmp[0]= mx[0][1];   tmp[1]= mx[0][2];   tmp[2]= mx[0][3];
                        tmp[3]= mx[1][2];   tmp[4]= mx[1][3];
                                            tmp[5]= mx[2][3];
                                                
    mx[0][1]=mx[1][0]; mx[1][0]=tmp[0];
    mx[0][2]=mx[2][0]; mx[2][0]=tmp[1];
    mx[0][3]=mx[3][0]; mx[3][0]=tmp[2];  
    mx[1][2]=mx[2][1]; mx[2][1]=tmp[3];
    mx[1][3]=mx[3][1]; mx[3][1]=tmp[4]; 
    mx[2][3]=mx[3][2]; mx[3][2]=tmp[5];
}

int mx_invert(matrix mx)
{
    matrix s;
    mx_copy(s, mx);
    mx_transpose(s);
    v4d t0, t1, t2;
    double d;
    int result;

    t0.x= s[2][2] * s[3][3];       //  = A2z*A3w   A2z*A3w
    t0.y= s[2][3] * s[3][2];       //  = A2w*A3z   A3z*A2w
    t0.z= s[2][1] * s[3][3];       //  = A2y*A3w   A1z*A3w
    t0.w= s[2][2] * s[3][1];       //  = A2z*A3y   A2z*A1w

    t1.x= s[2][3] * s[3][1];       // += A2w*A3y   A3z*A1w
    t1.y= s[2][0] * s[3][3];       // += A2x*A3w   A0z*A3w
    t1.z= s[2][3] * s[3][0];       // += A2w*A3x   A3z*A0w
    t1.w= s[2][0] * s[3][2];       // += A2x*A3z   A0z*Azw

    t2.x= s[2][1] * s[3][2];       // += A2y*A3z   A1z*Azw
    t2.y= s[2][2] * s[3][0];       // += A2z*A3x   A2z*A0w
    t2.z= s[2][0] * s[3][1];       // += A2x*A3y   A0z*A1w
    t2.w= s[2][1] * s[3][0];       // += A2y*A3x   A1z*A0w

    mx[0][0]  = t0.x*s[1][1] + t1.x*s[1][2] + t2.x*s[1][3];
    mx[0][0] -= t0.y*s[1][1] + t0.z*s[1][2] + t0.w*s[1][3];

    mx[0][1]  = t0.y*s[1][0] + t1.y*s[1][2] + t2.y*s[1][3];
    mx[0][1] -= t0.x*s[1][0] + t1.z*s[1][2] + t1.w*s[1][3];

    mx[0][2]  = t0.z*s[1][0] + t1.z*s[1][1] + t2.z*s[1][3];
    mx[0][2] -= t1.x*s[1][0] + t1.y*s[1][1] + t2.w*s[1][3];

    mx[0][3]  = t0.w*s[1][0] + t1.w*s[1][1] + t2.w*s[1][2];
    mx[0][3] -= t2.x*s[1][0] + t2.y*s[1][1] + t2.z*s[1][2];

    mx[1][0]  = t0.y*s[0][1] + t0.z*s[0][2] + t0.w*s[0][3];
    mx[1][0] -= t0.x*s[0][1] + t1.x*s[0][2] + t2.x*s[0][3];

    mx[1][1]  = t0.x*s[0][0] + t1.z*s[0][2] + t1.w*s[0][3];
    mx[1][1] -= t0.y*s[0][0] + t1.y*s[0][2] + t2.y*s[0][3];

    mx[1][2]  = t1.x*s[0][0] + t1.y*s[0][1] + t2.w*s[0][3];
    mx[1][2] -= t0.z*s[0][0] + t1.z*s[0][1] + t2.z*s[0][3];

    mx[1][3]  = t2.x*s[0][0] + t2.y*s[0][1] + t2.z*s[0][2];
    mx[1][3] -= t0.w*s[0][0] + t1.w*s[0][1] + t2.w*s[0][2];

                                    //     A^T        A
    t0.x = s[0][2]*s[1][3];         //  = A0y*A1w   A2x*A3y
    t0.y = s[0][3]*s[1][2];         //  = A0w*A1z   A2x*A2y
    t0.z = s[0][1]*s[1][3];         //  = A0y*A1w   A1x*A3y
    t0.w = s[0][2]*s[1][1];         //  = A0z*A1y   A2x*A1y

    t1.x = s[0][3]*s[1][1];         // += A0w*A1y   A3x*A1y
    t1.y = s[0][0]*s[1][3];         // += A0x*A1w   A0x*A3y
    t1.z = s[0][3]*s[1][0];         // += A0w*A1x   A3x*A0y
    t1.w = s[0][0]*s[1][2];         // += A0x*A1z   A0x*A2y

    t2.x = s[0][1]*s[1][2];         // += A0y*A1z   A1x*A2y
    t2.y = s[0][2]*s[1][0];         // += A0z*A1x   A2x*A0y
    t2.z = s[0][0]*s[1][1];         // += A0x*A1y   A0x*A1y
    t2.w = s[0][1]*s[1][0];         // += A0y*A1x   A1x*A0y

    mx[2][0]  = t0.x*s[3][1] + t1.x*s[3][2] + t2.x*s[3][3];
    mx[2][0] -= t0.y*s[3][1] + t0.z*s[3][2] + t0.w*s[3][3];

    mx[2][1]  = t0.y*s[3][0] + t1.y*s[3][2] + t2.y*s[3][3];
    mx[2][1] -= t0.x*s[3][0] + t1.z*s[3][2] + t1.w*s[3][3];

    mx[2][2]  = t0.z*s[3][0] + t1.z*s[3][1] + t2.z*s[3][3];
    mx[2][2] -= t1.x*s[3][0] + t1.y*s[3][1] + t2.w*s[3][3];

    mx[2][3]  = t0.w*s[3][0] + t1.w*s[3][1] + t2.w*s[3][2];
    mx[2][3] -= t2.x*s[3][0] + t2.y*s[3][1] + t2.z*s[3][2];

    mx[3][0]  = t0.z*s[2][2] + t0.w*s[2][3] + t0.y*s[2][1];
    mx[3][0] -= t1.x*s[2][2] + t2.x*s[2][3] + t0.x*s[2][1];

    mx[3][1]  = t1.w*s[2][3] + t0.x*s[2][0] + t1.z*s[2][2];
    mx[3][1] -= t2.y*s[2][3] + t0.y*s[2][0] + t1.y*s[2][2];

    mx[3][2]  = t1.y*s[2][1] + t2.w*s[2][3] + t1.x*s[2][0];
    mx[3][2] -= t1.z*s[2][1] + t2.z*s[2][3] + t0.z*s[2][0];

    mx[3][3]  = t2.z*s[2][2] + t2.x*s[2][0] + t2.y*s[2][1];
    mx[3][3] -= t2.w*s[2][2] + t0.w*s[2][0] + t1.w*s[2][1];

    // calculate determinant
    d= s[0][0]*mx[0][0] + s[0][1]*mx[0][1] + s[0][2]*mx[0][2] + s[0][3]*mx[0][3];
    result = (d != 0.0);

    if(result){
        d= 1.0 / d;
        mx_mul_scalar(mx, mx, d);
    }

    return result;
}

void mx_mul(matrix dest, matrix mx_a, matrix mx_b)
{
    matrix tmp;
    unsigned int i;
    double *pt_a, *pt_b, *pt_d;
    
    pt_b= (double *)mx_b;
    pt_a= (double *)mx_a;
    pt_d= (double *)tmp;
/*
//persumably outputs better assemply code
    for (i=0; i < 4; i++ ){
        *(pt_d+0)= ((*(pt_b+0)) * (*(pt_a+0))) + ((*(pt_b+4)) * (*(pt_a+1))) + ((*(pt_b+8)) * (*(pt_a+2))) + ((*(pt_b+12)) * (*(pt_a+3)));
        *(pt_d+1)= ((*(pt_b+1)) * (*(pt_a+0))) + ((*(pt_b+5)) * (*(pt_a+1))) + ((*(pt_b+9)) * (*(pt_a+2))) + ((*(pt_b+13)) * (*(pt_a+3)));
        *(pt_d+2)= ((*(pt_b+2)) * (*(pt_a+0))) + ((*(pt_b+6)) * (*(pt_a+1))) + ((*(pt_b+10))* (*(pt_a+2))) + ((*(pt_b+14)) * (*(pt_a+3)));
        *(pt_d+3)= ((*(pt_b+3)) * (*(pt_a+0))) + ((*(pt_b+7)) * (*(pt_a+1))) + ((*(pt_b+11))* (*(pt_a+2))) + ((*(pt_b+15)) * (*(pt_a+3)));
        pt_a+=4;
        pt_d+=4;
    }
*/
    for (i=0; i < 4; i++ ){
        pt_d[0]= (pt_b[0] * pt_a[0]) + (pt_b[4] * pt_a[1]) + (pt_b[8]  * pt_a[2]) + (pt_b[12] * pt_a[3]);
        pt_d[1]= (pt_b[1] * pt_a[0]) + (pt_b[5] * pt_a[1]) + (pt_b[9]  * pt_a[2]) + (pt_b[13] * pt_a[3]);
        pt_d[2]= (pt_b[2] * pt_a[0]) + (pt_b[6] * pt_a[1]) + (pt_b[10] * pt_a[2]) + (pt_b[14] * pt_a[3]);
        pt_d[3]= (pt_b[3] * pt_a[0]) + (pt_b[7] * pt_a[1]) + (pt_b[11] * pt_a[2]) + (pt_b[15] * pt_a[3]);
        pt_a+=4;
        pt_d+=4;
    }
    mx_copy(dest, tmp);
}

void mx_mul_scalar(matrix dst, matrix src, const double s)
{
    unsigned int i;
    double *dst_pt= (double*)dst;
    double *src_pt= (double*)src;
    for (i=0; i < 16; i++)
        *dst_pt++= (*src_pt++) * s;
}

// Vector ______________________________________________________________

void v_mul_mx(Vec * dest, const Vec * const v, matrix mx)
{
    vec tmp;
    tmp.x= (mx[0][0] * v->x) + (mx[1][0] * v->y) + (mx[2][0] * v->z) + (mx[3][0]);
    tmp.y= (mx[0][1] * v->x) + (mx[1][1] * v->y) + (mx[2][1] * v->z) + (mx[3][1]);
    tmp.z= (mx[0][2] * v->x) + (mx[1][2] * v->y) + (mx[2][2] * v->z) + (mx[3][2]);
    v_copy(dest, &tmp);
}

void v_zero(Vec * v)
{
    v->x= 0.0;
    v->y= 0.0;
    v->z= 0.0;
}

void v_add(Vec * dest, const Vec * const v0, const Vec * const v1)
{
    dest->x= v0->x + v1->x;
    dest->y= v0->y + v1->y;
    dest->z= v0->z + v1->z;
}

void v_sub(Vec * dest, const Vec * v0, const Vec * const v1)
{
    dest->x= v0->x - v1->x;
    dest->y= v0->y - v1->y;
    dest->z= v0->z - v1->z;
}

void v_scale(Vec * dest, const Vec * const v, const double s)
{
    dest->x= v->x * s;
    dest->y= v->y * s;
    dest->z= v->z * s;
}


double v_dot(const Vec * const v0, const Vec * const v1)
{
    return v0->x * v1->x + v0->y * v1->y + v0->z * v1->z;
}

void v_cross(Vec * dest, const Vec * const v0, const Vec * const v1)
{
    Vec c;
    c.x = (v0->y - v1->y) * (v0->z + v1->z);
    c.y = (v0->z - v1->z) * (v0->x + v1->x);
    c.z = (v0->x - v1->x) * (v0->y + v1->y);
    v_copy(dest, &c);
}

int same_2d( const Vec * const v0, const Vec * const v1, enum axis_pair ax)
{
    int r;
    switch (ax){
        case XY:
            r= (fabs(v0->x - v1->x) < ZEROLENGTH) &&
               (fabs(v0->y - v1->y) < ZEROLENGTH);
            break;
        case XZ:
            r= (fabs(v0->x - v1->x) < ZEROLENGTH) &&
               (fabs(v0->z - v1->z) < ZEROLENGTH);
        case YZ:
            r= (fabs(v0->y - v1->y) < ZEROLENGTH) &&
               (fabs(v0->z - v1->z) < ZEROLENGTH);
            break;        
        default:
            r= (fabs(v0->x - v1->x) < ZEROLENGTH) &&
               (fabs(v0->y - v1->y) < ZEROLENGTH);
            break;
    }
    return r;
}

int same_3d( const Vec * const v0, const Vec * const v1)
{
    return (fabs(v0->x - v1->x) < ZEROLENGTH &&
            fabs(v0->y - v1->y) < ZEROLENGTH &&
            fabs(v0->z - v1->z) < ZEROLENGTH);
}

double v_perp(const Vec * const v0,  const Vec * const v1, enum axis_pair ax)
{
    switch (ax){
        case XY:
            return v0->x * v1->y - v0->y * v1->x;
        case XZ:
            return v0->x * v1->z - v0->z * v1->x;
        case YZ:
            return v0->y * v1->z - v0->z * v1->y;
        default:
            return v0->x * v1->y - v0->y * v1->x;
    }
}

int rel_point_plane(const Vec * const p, const plane * const pn)
{
    Vec tmp;
    v_sub(&tmp, p, &pn->point);
    double d= v_dot(&tmp, &pn->normal);
    if (d > ZEROLENGTH)
        return FR;
    else if( d < -ZEROLENGTH)
        return BK;
    else
        return ON;
}

double v_length(const Vec * const v)
{
    return sqrt(v_dot(v, v));
}

void v_normal(const Vec * const v, Vec *dest)
{
    double L= v_length(v);
    if (fabs(L) > E)
        v_scale(dest, v, 1.0 / L);
    else
        v_zero(dest);
}

void midpoint(Vec *dest, const Vec * const v0, const Vec * const v1)
{
    Vec tmp;
    v_sub(&tmp, v1, v0);
    v_scale(&tmp, &tmp, 0.5);
    v_add(dest, v0, &tmp);
}

double dist( const Vec * const v0, const Vec * const v1)
{
    Vec tmp;
    v_sub(&tmp, v1, v0);
    return v_length(&tmp);
}

// Vector copy methods _________________________________________________

void v_copy(Vec *dest, const Vec * const src)
{
    dest->x= src->x;
    dest->y= src->y;
    dest->z= src->z;
}

void v_swp_xy(Vec *dest, const Vec * const src)
{
    double tmp= src->x;
    dest->x= src->y;
    dest->y= tmp;
    dest->z= src->z;
}

void v_swp_xz(Vec *dest, const Vec * const src)
{
    double tmp= src->x;
    dest->x= src->z;
    dest->y= src->y;
    dest->z= tmp;
}

void v_swp_yz(Vec *dest, const Vec * const src)
{
    double tmp= src->z;
    dest->x= src->x;
    dest->y= src->z;
    dest->z= tmp;
}

void v_flat_x(Vec *dest, const Vec * const src)
{
    dest->x= 0.0;
    dest->y= src->y;
    dest->z= src->z;
}

void v_flat_y(Vec *dest, const Vec * const src)
{
    dest->x= src->x;
    dest->y= 0.0;
    dest->z= src->z;
}

void v_flat_z(Vec *dest, const Vec * const src)
{
    dest->x= src->x;
    dest->y= src->y;
    dest->z= 0.0;
}

void v_flat_persp(Vec *dest, const Vec * const src)
{
    double z= src->z; 
    if (fabs(z) < E)
        v_copy(dest, src);
    else{
        dest->x= src->x / z;
        dest->x= src->x / z;
        dest->z= 0.0;
    }
}

void v_print(Vec v)
{
    printf("Vec: %.4f, %.4f, %.4f\n", v.x, v.y, v.z);
}


// other methods _______________________________________________________

int point_in_segment_2d( const Vec * const p, 
						 const Vec * const p0,  
						 const Vec * const p1, 
						 const enum axis_pair ax)
{
    Vec diff;
    v_sub(&diff, p1, p0);
    if (ax == XZ){
        if (fabs(diff.x) > fabs(diff.z))
            return ((p0->x < p->x && p->x < p1->x) ||
                    (p0->x > p->x && p->x > p1->x));
        else
            return ((p0->z < p->z && p->z < p1->z) ||
                    (p0->z > p->z && p->z > p1->z));
    }
    else if (ax == YZ){
        if (fabs(diff.y) > fabs(diff.z))
            return ((p0->y < p->y && p->y < p1->y) ||
                    (p0->y > p->y && p->y > p1->y));
        else
            return ((p0->z < p->z && p->z < p1->z) ||
                    (p0->z > p->z && p->z > p1->z));
          }
    else{
        if (fabs(diff.x) > fabs(diff.y))
            return ((p0->x < p->x && p->x < p1->x) ||
                    (p0->x > p->x && p->x > p1->x));
        else
            return ((p0->y < p->y && p->y < p1->y) ||
                    (p0->y > p->y && p->y > p1->y));
    }
}

int line_intersect_plane( const Vec * const p0, 
						  const Vec * const p1,
                          const plane * const pn, 
                          Vec * isx)
{
    Vec u, w;
    v_sub(&u, p1, p0);
    v_sub(&w, p0, &pn->point);
    double d= v_dot(&pn->normal, &w);
    double N= -v_dot(&pn->normal, &u);
    if (fabs(d) < E){
        if (fabs(N) < E){
            v_copy(isx, p0);
            return 1;
        }
        return 0;
    }
    v_scale(&u, &u, N/d);
    v_add(isx, p0, &u);

    return 1;
}

float v_angle(const Vec * const v0, const Vec * const v1)
{
	Vec n0, n1;
	v_normal(v0, &n0);
	v_normal(v1, &n1);
	float angle= acos(v_dot(&n0, &n1));
	return angle * RAD2DEG;
}

float v_angle_r(const Vec * const v0, const Vec * const v1)
{
	Vec n0, n1;
	v_normal(v0, &n0);
	v_normal(v1, &n1);
	float angle= acos(v_dot(&n0, &n1));
	return angle;
}
