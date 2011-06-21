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



cdef extern from "globals.h":
    ctypedef unsigned int uint
    int IN, ON, ONIN, OUT, ONOUT, CR, ONCR, FR, ONFR, BK, ONBK
    enum axis_pair:
        XY, XZ, YZ
    ctypedef axis_pair AxisPair

cdef extern from "mathutils.h":
    ctypedef double matrix[4][4]

    ctypedef double vec[3]

    struct plane:
        vec point
        vec normal
    ctypedef plane Plane

    void v_mul_mx (vec dest, vec v,  matrix mx) nogil
    void v_add(vec dest, vec v0, vec v1) nogil
    void v_sub(vec dest, vec v0, vec v1) nogil
    void v_scale (vec dest, vec v, double s) nogil
    void v_cross (vec dest, vec v0,  vec v1) nogil
    double v_dot (vec v0, vec v1) nogil
    double v_perp (vec v0,  vec v1, AxisPair ax) nogil
    int same_2d(vec v0,  vec v1, AxisPair ax) nogil
    int same_3d(vec v0,  vec v1) nogil

    double v_length (vec v) nogil
    void v_normal (vec dest, vec v) nogil
    void v_zero (vec v) nogil

    void v_cpy(vec dest, vec src) nogil
    void v_swp_xy (vec dest, vec src) nogil
    void v_swp_yz (vec dest, vec src) nogil
    void v_swp_xz (vec dest, vec src) nogil
    void v_flat_x (vec dest, vec src) nogil
    void v_flat_y (vec dest, vec src) nogil
    void v_flat_z (vec dest, vec src) nogil
    void v_flat_persp (vec dest, vec src) nogil

    void mx_copy(matrix dst_mx,  matrix src_mx)nogil
    void mx_identity "mx_set_identity"(matrix mx)nogil
    void mx_rotation "mx_set_rotation" (matrix mx,  float rot_x,  float rot_y,  float rot_z)nogil
    void mx_translation "mx_set_translation" (matrix mx,  double dx,  double dy,  double dz)nogil
    void mx_transpose "mx_transpose" (matrix mx) nogil
    int mx_invert "mx_invert" (matrix mx) nogil
    void mx_mul_mx "mx_mul" (matrix mx,  matrix mx_a,  matrix mx_b) nogil
    void mx_mul_scalar "mx_mul_scalar" (matrix dst, matrix mx,  double s) nogil

    int rel_point_plane( vec point,  Plane *pn) nogil
    int line_intersect_plane( vec v0,  vec v1,  Plane *pn, vec isx) nogil
    int point_in_segment_2d(vec p,  vec p0,  vec p1, AxisPair ax) nogil
    void midpoint(vec dest, vec v0, vec v1) nogil
    double dist(vec v0,  vec v1) nogil

cdef extern from "bounds.h":
    struct bounds:
        vec min, max

    ctypedef bounds Bounds

    void check_bounds(Bounds *b, vec point)nogil
    int disjoint_2d(Bounds *b0, Bounds *b1, AxisPair ax)nogil
    int disjoint_3d(Bounds *b0, Bounds *b1)nogil
    int point_in_bounds_2d (Bounds *b0, vec p, AxisPair ax)nogil
    int point_in_bounds (Bounds *b, vec p)nogil

cdef extern from "polygon.h":
    enum PL_Order:
        CW=0
        CCW

    enum PL_Attr:
        HIDDEN, HOLE, SECTION, BACKF, CLIPPED

    struct polyvert:
        vec co
        int flags

    ctypedef polyvert PolyVert

    struct polygon
    struct polylist:
        polygon **polys
        uint size, last
    ctypedef polylist PolyList

    struct polygon:
        PolyVert *points
        uint size, last, flags
        Bounds bb
        PolyList *holes

    ctypedef polygon Polygon

    Polygon *pl_init(uint count) nogil
    Polygon *pl_copy(Polygon *pl, void (*cpy)(vec dest, vec src), int holes) nogil
    int pl_add_hole(Polygon *pl, Polygon *hole) nogil

    void pl_kill (Polygon **pl) nogil

    PolyVert *pl_append (Polygon **pl, vec co, int flags) nogil
    PolyVert *pl_append_vertex (Polygon **pl, PolyVert *pt) nogil

    PolyVert *pl_insert_point (Polygon *pl, vec p, uint pos) nogil
    PolyVert *pl_insert_vertex (Polygon *pl, PolyVert *p, uint pos) nogil

    void pl_extend_verts (Polygon *pl, PolyVert *points, uint count) nogil
    void pl_extend_points (Polygon *pl, vec *points, uint count) nogil
    void pl_reverse (Polygon *pl) nogil
    void pl_project_to_plane (Polygon *pl,  Plane *pn, vec direction) nogil

    void pl_normal (Polygon *pl, vec normal) nogil

    double pl_area_3d(Polygon *pl, int do_holes) nogil
    double pl_area_2d(Polygon *pl, AxisPair ax, int holes) nogil

    PL_Order pl_order(Polygon *pl) nogil
    PL_Order pl_order_2d(Polygon *pl) nogil

    int pl_rel_point_2d (Polygon *pl, vec p, AxisPair ax, int holes) nogil
    int pl_rel_plane (Polygon *pl, Plane pn, int *rmap) nogil
    int pl_self_intersects (Polygon *pl, AxisPair ax) nogil

    unsigned int pl_rm_doubles_2d (Polygon *pl, AxisPair ax) nogil
    unsigned int pl_rm_doubles_3d (Polygon *pl) nogil
    unsigned int pl_simplify_3d (Polygon *pl) nogil
    unsigned int pl_simplify_2d (Polygon *pl, AxisPair ax) nogil

    void pl_split_by_plane (Polygon *pl,\
                            Plane *pn,\
                            PolyList *fparts,\
                            PolyList *bparts) nogil


    PolyList * pll_init (unsigned int count) nogil
    void pll_kill (PolyList ** pll) nogil
    void pll_killall (PolyList ** pll) nogil

    void pll_append(PolyList * *pll, Polygon * pl)nogil

    Polygon * pll_pop (PolyList * pll)nogil

    PolyList * pll_copy (PolyList * pll, void (*copy_method)(vec*dest, vec*src))nogil

    int pll_remove_pl (PolyList * pll, Polygon * pl) nogil
    Polygon * pll_rem_poly_idx (PolyList * pll, unsigned int idx) nogil
    unsigned int pll_count_points (PolyList * pll) nogil
    void pll_extend (PolyList ** pll, Polygon ** polys, unsigned int count) nogil
    Polygon * pll_replace (PolyList * pll, unsigned int i, Polygon * pl) nogil
    void pll_replace_pl (PolyList * pll, unsigned int i, Polygon * pl) nogil

    void pll_sort_by_area (PolyList * pll) nogil

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
    PolyList * bl_operation(BL_Graph *g, BL_Op op)nogil




