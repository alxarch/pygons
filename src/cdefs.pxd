#       cdefs.pxd
#       
#       Copyright 2011 Alexandros Sigalas <alxarch@gmail.com>
#       
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.


cdef extern from "math.h":
    double fabs(double d)nogil
    double sqrtf(double d)nogil

cdef extern from "stdlib.h":
    ctypedef unsigned long size_t
    void free(void *ptr) nogil
    void *malloc(size_t size) nogil
    void *realloc(void *ptr, size_t size) nogil
    void *calloc(size_t nelem, size_t elsize) nogil

cdef extern from "globals.h":
    ctypedef unsigned int uint
    int IN, ON, ONIN, OUT, ONOUT, CR, ONCR, FR, ONFR, BK, ONBK
    enum axis_pair:
        XY, XZ, YZ   

cdef extern from "mathutils.h":
    ctypedef double matrix[4][4]
    struct Vec:
        double x,y,z

    struct Plane:
        Vec point
        Vec normal
       
    void v_mul_mx (Vec * dest, Vec * v,  matrix mx) nogil
    void v_add(Vec * dest, Vec * v0, Vec * v1) nogil
    void v_sub(Vec * dest, Vec * v0, Vec * v1) nogil
    void v_scale (Vec * dest, Vec * v, double s) nogil
    void v_cross (Vec * dest, Vec * v0,  Vec *v1) nogil
    double v_dot (Vec * v0, Vec * v1) nogil
    double v_perp (Vec * v0,  Vec * v1, axis_pair ax) nogil
    int same_2d(Vec * v0,  Vec * v1, axis_pair ax) nogil
    int same_3d(Vec * v0,  Vec * v1) nogil
    
    double v_length (Vec * v) nogil
    void v_normal (Vec * dest, Vec * v) nogil
    void v_zero (Vec * v) nogil

    void v_copy(Vec * dest, Vec * src) nogil
    void v_swp_xy (Vec *dest, Vec *src) nogil
    void v_swp_yz (Vec *dest, Vec *src) nogil
    void v_swp_xz (Vec *dest, Vec *src) nogil
    void v_flat_x (Vec *dest, Vec *src) nogil
    void v_flat_y (Vec *dest, Vec *src) nogil
    void v_flat_z (Vec *dest, Vec *src) nogil
    void v_flat_persp (Vec *dest, Vec *src) nogil
    
    void mx_copy(matrix dst_mx,  matrix src_mx)nogil
    void mx_identity "mx_set_identity"(matrix mx)nogil
    void mx_rotation "mx_set_rotation" (matrix mx,  float rot_x,  float rot_y,  float rot_z)nogil
    void mx_translation "mx_set_translation" (matrix mx,  double dx,  double dy,  double dz)nogil
    void mx_transpose "mx_transpose" (matrix mx) nogil
    int mx_invert "mx_invert" (matrix mx) nogil
    void mx_mul_mx "mx_mul" (matrix mx,  matrix mx_a,  matrix mx_b) nogil
    void mx_mul_scalar "mx_mul_scalar" (matrix dst, matrix mx,  double s) nogil

    int rel_point_plane( Vec *point,  Plane * pn) nogil
    int line_intersect_plane( Vec *v0,  Vec *v1,  Plane *pn, Vec *isx) nogil
    int point_in_segment_2d(Vec *p,  Vec *p0,  Vec *p1, axis_pair ax) nogil
    void midpoint(Vec * dest, Vec *v0, Vec * v1) nogil
    double dist(Vec * v0,  Vec * v1) nogil

cdef extern from "bounds.h":
    struct Bounds:
        Vec min, max

    void check_bounds(Bounds * b, Vec * point)nogil
    int disjoint_2d(Bounds * b0, Bounds * b1, axis_pair ax)nogil
    int disjoint_3d(Bounds * b0, Bounds * b1)nogil
    int point_in_bounds_2d (Bounds * b0, Vec *p, axis_pair ax)nogil
    int point_in_bounds (Bounds * b, Vec * p)nogil

cdef extern from "polygon.h":
    enum PL_Order:
        CW=0
        CCW

    enum PL_Attr:
        HIDDEN, HOLE, SECTION, BACKF, CLIPPED

    struct PL_Vertex:
        Vec co
        int flags

    struct Polygon
    struct PL_List:
        Polygon **polys
        unsigned int size, last

    struct Polygon:
        PL_Vertex *points
        uint size, last
        int flags
        Bounds bb
        PL_List *holes
        
    Polygon * pl_init (uint count) nogil
    Polygon * pl_copy (Polygon * pl, void (*cpy)(Vec *dest, Vec *src), int do_holes) nogil
    int pl_add_hole (Polygon * pl, Polygon * hole) nogil

    void pl_kill (Polygon **pl) nogil

    PL_Vertex* pl_append (Polygon **pl, Vec *co, int flags) nogil
    PL_Vertex* pl_append_vertex (Polygon **pl, PL_Vertex *pt) nogil

    PL_Vertex* pl_insert_point (Polygon *pl, Vec *p, uint pos) nogil
    PL_Vertex* pl_insert_vertex (Polygon *pl, PL_Vertex *p, uint pos) nogil

    void pl_extend_verts (Polygon *pl, PL_Vertex *points, uint count) nogil
    void pl_extend_points (Polygon *pl, Vec * points, uint count) nogil
    void pl_reverse (Polygon *pl) nogil
    void pl_project_to_plane (Polygon *pl,  Plane *pn, Vec *direction) nogil
    
    void pl_normal (Polygon *pl, Vec *normal) nogil
    
    double pl_area_3d(Polygon *pl, int do_holes) nogil
    double pl_area_2d(Polygon *pl, axis_pair ax, int do_holes) nogil
    
    PL_Order pl_order(Polygon * pl) nogil
    PL_Order pl_order_2d(Polygon * pl) nogil

    int pl_rel_point_2d (Polygon * pl, Vec * p, axis_pair ax, int do_holes) nogil
    int pl_rel_plane (Polygon * pl, Plane * pn, int *rmap) nogil
    int pl_self_intersects (Polygon * pl, axis_pair ax) nogil

    unsigned int pl_rm_doubles_2d (Polygon * pl, axis_pair ax) nogil
    unsigned int pl_rm_doubles_3d (Polygon * pl) nogil
    unsigned int pl_simplify_3d (Polygon * pl) nogil
    unsigned int pl_simplify_2d (Polygon * pl, axis_pair ax) nogil

    void pl_split_by_plane (Polygon * pl,\
                            Plane *pn,\
                            PL_List *fparts,\
                            PL_List *bparts) nogil


    PL_List * pll_init (unsigned int count) nogil
    void pll_kill (PL_List ** pll) nogil
    void pll_killall (PL_List ** pll) nogil

    void pll_append(PL_List * *pll, Polygon * pl)nogil

    Polygon * pll_pop (PL_List * pll)nogil
    
    PL_List * pll_copy (PL_List * pll, void (*copy_method)(Vec *dest, Vec *src))nogil
    
    int pll_remove_pl (PL_List * pll, Polygon * pl) nogil
    Polygon * pll_rem_poly_idx (PL_List * pll, unsigned int idx) nogil
    unsigned int pll_count_points (PL_List * pll) nogil
    void pll_extend (PL_List ** pll, Polygon ** polys, unsigned int count) nogil
    Polygon * pll_replace (PL_List * pll, unsigned int i, Polygon * pl) nogil
    void pll_replace_pl (PL_List * pll, unsigned int i, Polygon * pl) nogil

    void pll_sort_by_area (PL_List * pll) nogil

cdef extern from "boolean.h":
    enum itype:
        NOTYPE= 0, EN, EX, ENEX, EXEN, INVALID

    enum BL_Set:
        A=0, B
        
    enum BL_Op:
        UNION, SUBTRACT_A_B, SUBTRACT_B_A, INTERSECT

    struct BL_Node:
        pass
    
    struct BL_Poly:
        Polygon * pl
        int hole
        int has_entry
        BL_Set set
        int rel

    struct BL_Graph:
        BL_Node* nodes
        unsigned int last, size, inum
        BL_Poly* polys
        BL_Poly* pl_a
        BL_Poly* pl_b
        unsigned int pl_num

    BL_Graph* bl_build_graph(Polygon * pl_a, Polygon * pl_b)nogil
    void bl_kill_graph(BL_Graph *g)nogil
    PL_List * bl_operation(BL_Graph *g, BL_Op op)nogil




