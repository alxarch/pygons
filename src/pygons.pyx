#       pygons.pyx
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


cimport cdefs as c

# Vector _______________________________________________________________

cdef c.Vec asvec(object ob):
    cdef c.Vec v
    if isinstance(ob, Vector):
        return Vector.value(ob)
    try:
        v.x = ob[0]
        v.y = ob[1]
        v.z = ob[2]
        
    except:
        raise TypeError("object cannot become vector")
    return v

cdef void asmatrix(c.matrix m, object ob):
    
    m[0][0]= ob[0][0]
    m[0][1]= ob[0][1]
    m[0][2]= ob[0][2]
    m[0][3]= ob[0][3]
    
    m[1][0]= ob[1][0]
    m[1][1]= ob[1][1]
    m[1][2]= ob[1][2]
    m[1][3]= ob[1][3]
    
    m[2][0]= ob[2][0]
    m[2][1]= ob[2][1]
    m[2][2]= ob[2][2]
    m[2][3]= ob[2][3]
    
    m[3][0]= ob[3][0]
    m[3][1]= ob[3][1]
    m[3][2]= ob[3][2]
    m[3][3]= ob[3][3]

cdef class Plane

cdef class Vector:
    cdef:
        c.Vec *v
        readonly bint wrapped
        
    property x:
        def __get__(self):
            return self.v.x
        def __set__(self, double value):
            self.v.x = value

    property y:
        def __get__(self):
            return self.v.y
        def __set__(self, double value):
            self.v.y = value

    property z:
        def __get__(self):
            return self.v.z
        def __set__(self, double value):
            self.v.z = value

    property length:
        def __get__(self):
            return c.v_length(self.v)

    property normal:
        def __get__(self):
            cdef Vector normal= Vector()
            c.v_normal(normal.v, self.v)
            return normal

    def __cinit__(self, object co=[0.0,0.0,0.0]):
        self.wrapped = 0
        if co is None:
            self.v = NULL
        else:
            self.v = <c.Vec *>c.malloc(sizeof(c.Vec))
            self.v.x = <double>(co[0])
            self.v.y = <double>(co[1])
            self.v.z = <double>(co[2])
            

    def __dealloc__(self):
        if self.wrapped == 0:
            c.free(self.v)

    def __str__(self):
        cdef object s
        s = "[%.4f, %.4f, %.4f]"%(self.x, self.y, self.z)
        if self.wrapped:
            s += "(wrapped data)"
        s += "\n"
        return s

    def __add__(Vector x, Vector y):
        cdef Vector result= Vector()
        c.v_add(result.v, x.v, y.v);
        return result

    def __sub__(Vector x, Vector y):
        cdef Vector result= Vector()
        c.v_sub(result.v, x.v, y.v);
        return result

    def __mul__(x , y):
        cdef Vector result= Vector()
        if isinstance(x, Vector):
            if isinstance(y, Vector):
                return c.v_dot(Vector.pointer(x), Vector.pointer(y))
            else:
                c.v_scale(result.v, Vector.pointer(x), y)
        else:
            c.v_scale(result.v, Vector.pointer(y), x)
        return result
        
    def __iadd__(self, Vector other):
        c.v_add(self.v, self.v, other.v)
        return self

    def __isub__(self, Vector other):
        c.v_sub(self.v, self.v, other.v)
        return self

    def __imul__(self, double other):
        c.v_scale(self.v, self.v, other)
        return self

    def __getitem__(self, int i):
        cdef double *pt=<double *>self.v
        if 0 <= i < 3:
            return pt[i]
        else:
            raise IndexError("0 <= index < 3")

    def __setitem__(self, unsigned int i, double value):
        cdef double *pt= <double *>self.v
        if i < 3:
            pt[i] = value
        else:
            raise IndexError("0 <= i < 3")

    def dot(self, Vector other):
        return c.v_dot(self.v, other.v)

    def copy(self):
        if self.v == NULL:
            return Vector()
        return Vector([self.v.x, self.v.y, self.v.z])

    def same3D(self, object other):
        cdef c.Vec tmp= asvec(other)
        return c.same_3d(self.v, &tmp)

    def same2D(self, object other):
        cdef c.Vec tmp=  asvec(other)
        return c.same_2d(self.v, &tmp, c.XY)

    def perp(self, object other):
        cdef c.Vec tmp=  asvec(other)
        return c.v_perp(self.v, &tmp, c.XY)

    def midTo(self, object other):
        cdef:
            c.Vec tmp = asvec(other)
            Vector mid = Vector()
        c.midpoint(mid.v, self.v, &tmp)
        return mid

    def apply_matrix(self, object obj):
        cdef c.matrix mx
        asmatrix(mx, obj)
        c.v_mul_mx(self.v, self.v, mx)

    def rel2Plane(self, Plane pl):
        return c.rel_point_plane(self.v, pl.pn)

    cdef void wrap(self, c.Vec * v):
        if self.wrapped == 0:
            c.free(self.v)
            self.wrapped = 1
        self.v = v

    cdef void grab(self, c.Vec* v):
        if self.wrapped:
            self.wrapped= 0
        else:
            c.free(self.v)
        self.v = v

    cdef c.Vec value(self):
        return self.v[0]

    cdef c.Vec* pointer(self):
        return self.v

# Plane ________________________________________________________________

cdef class Plane:
    cdef:
        c.Plane * pn
        readonly bint wrapped

    property point:
        def __get__(self):
            cdef Vector result = Vector()
            result.wrap(&self.pn.point)
            return result

        def __set__(self, value):
            self.pn.point = asvec(value)

    property normal:
        def __get__(self):
            cdef Vector result= Vector()
            c.v_normal(result.v, &self.pn.normal)
            return result
        def __set__(self, value):
            self.pn.normal = asvec(value)

    def __cinit__(self, object normal=[0.0,0.0,1.0], object point=[0.0,0.0,0.0]):
        self.wrapped = 0
        if normal is None:
            self.pn = NULL
        else:
            self.pn = <c.Plane *>c.malloc(sizeof(c.Plane))
            self.pn.normal = asvec(normal)
            self.pn.point = asvec(point)

    def __dealloc__(self):
        if self.wrapped == 0:
            c.free(self.pn)

    def __str__(self):
        cdef object s
        
        s = "Plane(%s, %s)"%(self.normal, self.point)
        if self.wrapped:
            s += "(wrapped data)"
        s += "\n"
        return s

    cdef wrap(self, c.Plane * pn):
        if(self.wrapped == 0):
            c.free(self.pn)
        self.pn = pn
        self.wrapped = 1

# Bounds _______________________________________________________________

cdef class Bounds:
    cdef c.Bounds * b

    property xmin:
        def __get__(self):
            return self.b.min.x

    property ymin:
        def __get__(self):
            return self.b.min.y

    property zmin:
        def __get__(self):
            return self.b.min.z

    property xmax:
        def __get__(self):
            return self.b.max.x

    property ymax:
        def __get__(self):
            return self.b.max.y

    property zmax:
        def __get__(self):
            return self.b.max.z
    
    def __cinit__(self, object ob):
        cdef:
            unsigned int i
            c.Vec tmp
        if ob is None:
            self.b= NULL;
        else:
            self.b= <c.Bounds *>c.malloc(sizeof(c.Bounds))
            tmp= asvec(ob[0])
            c.v_copy(&self.b.min, &tmp)
            c.v_copy(&self.b.max, &tmp)
            for i from 0 < i < len(ob):
                tmp= asvec(ob[i])
                c.check_bounds(self.b, &tmp)
    
    def __dealloc__(self):
        c.free(self.b)
    
    def disjoint2D(self, object other):
        cdef Bounds tmp
        if other is None:
            raise ValueError("Object is None")
        elif isinstance(other, Bounds):
            return c.disjoint_2d(self.b, Bounds.pointer(other), c.XY)
        else:
            tmp= Bounds(other)
            return c.disjoint_2d(self.b, tmp.b, c.XY)
        
    def disjoint3D(self, object other):
        cdef Bounds tmp
        if other is None:
            raise ValueError("Object is None")
        elif isinstance(other, Bounds):
            return c.disjoint_3d(self.b, Bounds.pointer(other))
        else:
            tmp= Bounds(other)
            return c.disjoint_3d(self.b, tmp.b)
    
    def pointInside3D(self, object other):
        cdef c.Vec tmp= asvec(other)
        return c.point_in_bounds(self.b, &tmp)
    
    def pointInside2D(self, object other):
        cdef c.Vec tmp = asvec(other)
        return c.point_in_bounds_2d(self.b, &tmp, c.XY)
    
    cdef c.Bounds * pointer(self):
        return self.b

cdef Bounds b2B(c.Bounds * bb):
    cdef Bounds result= Bounds(None)
    result.b= bb
    return result

# Polygon ______________________________________________________________

cdef class Polyvert:
    cdef:
        c.PL_Vertex *vx
        readonly bint wrapped

    property co:
        def __get__(self):
            cdef Vector result= Vector(None)
            result.wrap(&self.vx.co)
            return result
        
        def __set__(self, object value):
            cdef c.Vec tmp= asvec(value)
            c.v_copy(&self.vx.co, &tmp)
    
    property flags:
        def __get__(self):
            return self.vx.flags
        def __set__(self, int value):
            self.vx.flags= value
    
    def __cinit__(self, object co, int flags=0):
        cdef c.Vec tmp
        self.wrapped= 0
        if co is None:
            self.vx= NULL
        else:
            self.vx= <c.PL_Vertex *>c.malloc(sizeof(c.PL_Vertex))
            tmp= asvec(co)
            c.v_copy(&self.vx.co, &tmp)
            self.vx.flags= flags

    def __str__(self):
        return "Vertex Oject" #TODO: Vertex string output for debugging
    
    cdef void wrap(self, c.PL_Vertex *vx):
        if not self.wrapped:
            c.free(self.vx)
            self.wrapped= 1
        self.vx= vx

    def __dealloc__(self):
        if not self.wrapped:
            c.free(self.vx)

cdef class Poly

cdef class PolyList:
    cdef:
        c.PL_List **pll_pt
        unsigned int index
        readonly bint wrapped
        
    def __cinit__(self, list polygons=list()):
        self.wrapped= 0
        self.index= 0
        self.pll_pt= <c.PL_List **>c.malloc(sizeof(c.PL_List*))

        cdef unsigned int i, limit
        if polygons is None:
            self.pll_pt[0]= NULL
        if not polygons:
            self.pll_pt[0]= c.pll_init(10)
        else:
            limit= len(polygons)
            self.pll_pt[0]= c.pll_init(limit)
            for i from 0 <= i < limit:
                self.append(polygons[i])

    def __dealloc__(self):
        if not self.wrapped:
            c.pll_killall(self.pll_pt)
        c.free(self.pll_pt)
    
    def __getitem__(self, unsigned int i):
        if (self.pll_pt[0] == NULL):
            raise IndexError("PolyList holds no data")
        if (self.pll_pt[0].last <= i):
            raise IndexError("Index out of bounds")
        cdef Poly result= Poly(coords=None)
        result.wrap(self.pll_pt[0].polys[i])
        return result

    def __setitem__(self, unsigned int i, Poly data):
        cdef c.Polygon * tmp
        if self.pll_pt[0] == NULL:
            raise IndexError("PolyList holds no data")
        if (i < self.pll_pt[0].last):
            tmp= c.pll_replace(self.pll_pt[0], i, data.pl)
            c.pl_kill(&tmp)
        else:
            raise IndexError("Index out of bounds")

    def __next__(self):
        cdef Poly result= Poly(coords=None)
        if self.index == self.pll_pt[0].last:
            self.index= 0
            raise StopIteration
        else:
            result.wrap(self.pll_pt[0].polys[self.index])
            self.index += 1
            return result

    cpdef append(self, Poly pl):
        c.pll_append(self.pll_pt, pl.pl)

    cdef wrap(self, c.PL_List ** data):
        if not self.wrapped:
            c.pll_killall(self.pll_pt)
            self.wrapped= 1
        self.pll_pt= data;

    cdef grab(self, c.PL_List * data):
        if self.wrapped:
            self.wrapped= 0
        else:
            c.pll_killall(self.pll_pt)
        self.pll_pt[0]= data

cdef PolyList pll2PLL(c.PL_List * pll):
    cdef PolyList result= PolyList(polygons=None)
    result.grab(pll)
    return result


cdef class Poly:
    cdef:
        c.Polygon * pl
        unsigned int index
        readonly bint wrapped
        
    # properties ################
    property bounds:
        def __get__(self):
            return b2B(&self.pl.bb)
            
    property area_3d:
        def __get__(self):
            return c.pl_area_3d(self.pl, 1)

    property area_2d:
        def __get__(self):
            return c.pl_area_2d(self.pl, c.XY, 1)

    property normal:
        def __get__(self):
            cdef Vector result= Vector()
            c.pl_normal(self.pl, result.v)
            return result

    property plane:
        def __get__(self):
            cdef Plane result = Plane()
            c.pl_normal(self.pl, &result.pn.normal) 
            c.v_copy(&self.pl.points[0].co, &result.pn.point)
            return result
        
    
    property flags:
        def __get__(self):
            return self.pl.flags
        def __set__(self, int value):
            self.pl.flags= value

    property order:
        def __get__(self):
            return c.pl_order(self.pl)

    property holes:
        def __get__(self):
            cdef PolyList result= PolyList(polygons=None)
            result.wrap(&self.pl.holes)
            return result
    
    # methods ################
    def __cinit__(self, object coords, int flags=0):

        self.index= 0
        self.wrapped=0
        cdef unsigned int i, limit, size
        cdef c.Vec co
        if coords is None:
            self.pl= NULL
        else:
            limit= len(coords)
            size= 100 if limit < 10 else 2 * limit
            self.pl= c.pl_init(size)
            self.pl.flags = flags
            for i from 0 <= i < limit:
                co= asvec(coords[i])
                c.pl_append(&self.pl, &co, 0)
        
    cdef c.Polygon * pointer(self):
        return self.pl

    def __dealloc__(self):
        c.pl_kill(&self.pl)

    def __next__(self):
        # returns wrapped Vector object
        cdef:
            Polyvert vx= Polyvert(co=None)
            c.PL_Vertex *v= self.pl.points + self.index

        if self.index == self.pl.last:
            self.index = 0
            raise StopIteration
        else:
            vx.wrap(v)
            self.index += 1
            return vx

    def __str__(self):
        cdef object s
        s = "Polygon:\n"
        for p in self:
            s += "\t" + str(p) + "\n"
        return s

    def __len__(self):
        return self.pl.last

    def __getitem__(self, unsigned int i):
        # returns wrapped Vector object
        cdef:
            Polyvert v= Polyvert(None)
            c.PL_Vertex *vx

        if self.pl == NULL:
            raise IndexError("Poly holds no data")
        elif i < self.pl.last:
            vx= self.pl.points + i
            v.wrap(vx)
            return v
        else:
            raise IndexError("index out of bounds")

    def __setitem__(self, int i, object p):
        cdef c.Vec tmp
        if i < self.pl.last:
            tmp= asvec(p)
            c.v_copy(&self.pl.points[i].co, &tmp)
        else:
            raise IndexError()

    def copy(self, bint holes=True):
        return pl2PL(c.pl_copy(self.pl, NULL, holes))

    def append(self, object co, int flags=0):
        cdef c.Vec pt= asvec(co)
        c.pl_append(&self.pl, &pt, flags);

    def addHole(self, Poly hole):
        return c.pl_add_hole(self.pl, hole.pl)

    def reverse(self):
        c.pl_reverse(self.pl)

    def relPoint2D(self, Vector p, bint holes=False):
        return c.pl_rel_point_2d(self.pl, p.v, c.XY, holes)

    def disjoint2D(self, Poly other):
        return c.disjoint_2d(&self.pl.bb, &other.pl.bb, c.XY)

    def disjoint3D(self, Poly other):
        return c.disjoint_3d(&self.pl.bb, &other.pl.bb)

    def swap2D(self, int axp=c.XY):
        cdef c.Polygon * old_pl= self.pl
        if axp == c.XZ:
            self.pl= c.pl_copy(old_pl, c.v_swp_xz, 1)
        elif axp == c.YZ:
            self.pl= c.pl_copy(old_pl, c.v_swp_yz, 1)
        else:
            self.pl= c.pl_copy(old_pl, c.v_swp_xy, 1)
        c.pl_kill(&old_pl)

    def project2Plane(self, Plane plane, Vector direction):
        c.pl_project_to_plane(self.pl, plane.pn, direction.v)

    def project2Poly(self, Poly other, Vector direction):
        self.project2Plane(other.plane, direction)

    def split2Plane(self, Plane plane):
        cdef c.PL_List * fparts= NULL
        cdef c.PL_List * bparts= NULL
        c.pl_split_by_plane(self.pl, plane.pn, fparts, bparts)
        return pll2PLL(fparts), pll2PLL(bparts)
        
    def removeDoubles2D(self, c.axis_pair axp=c.XY):
        return c.pl_rm_doubles_2d(self.pl, axp)

    def removeDoubles3D(self):
        return c.pl_rm_doubles_3d(self.pl)

    def simplify2D(self, c.axis_pair axp=c.XY):
        return c.pl_simplify_2d(self.pl, axp)

    def simplify3D(self):
        return c.pl_simplify_3d(self.pl)

    def flattenZ(self):
        cdef c.Polygon * old_pl= self.pl
        self.pl= c.pl_copy(self.pl, c.v_flat_z, 1)
        c.pl_kill(&old_pl)
    
    def flattenPersp(self):
        cdef c.Polygon * old_pl= self.pl
        self.pl= c.pl_copy(self.pl, c.v_flat_persp, 1)
        c.pl_kill(&old_pl)

    cdef void grab(self, c.Polygon * data):
        if self.wrapped:
            self.wrapped= 0
        else:
            c.pl_kill(&self.pl)
        self.pl= data
        
    cdef void wrap(self, c.Polygon * data):
        if self.wrapped == 0:
            self.wrapped = 1
            c.pl_kill(&self.pl)
        self.pl= data

cdef Poly pl2PL(c.Polygon * pl):
    cdef Poly result= Poly(None)
    result.grab(pl)
    return result

# Boolean ______________________________________________________

bool_op= {"A-B":c.SUBTRACT_A_B,\
          "B-A":c.SUBTRACT_B_A,\
          "A+B":c.UNION,\
          "A|B":c.INTERSECT}

cdef class BooleanGraph:
    cdef c.BL_Graph *g
    property rel_a:
        def __get__(self):
            return self.g.pl_a.rel

    property rel_b:
        def __get__(self):
            return self.g.pl_b.rel
        
    def __cinit__(self, Poly pl_a, Poly pl_b):
        self.g= c.bl_build_graph(pl_a.pl, pl_b.pl)
    
    def __str__(self):
        #TODO:BooleanGraph string output for debugging
        return "BooleanGraph" 
    
    def __call__(self, op):
        if op != c.SUBTRACT_B_A and\
           op != c.SUBTRACT_A_B and\
           op != c.UNION and\
           op != c.INTERSECT:
               raise ValueError("Unsupported operation")
        cdef PolyList result= PolyList(polygons=None)
        result.grab(c.bl_operation(self.g, op))
        return result
        
