from netgen.occ import *
from ngsolve import Mesh

import sys
sys.path.append("../build/rigid_body")
from rigid_body import *
import rigid_body.bla as bla

import pythreejs as p3

from IPython.display import display
import ipywidgets as widgets
# set up OCC CAD model
box = Box(Pnt(-0.5,-0.5,-0.5), Pnt(0.5,0.5,0.5))
geo = OCCGeometry(box)
mesh = Mesh(geo.GenerateMesh(maxh=0.5))
# important: move the center of mass into the origin
box = box.Move((-box.center[0], -box.center[1], -box.center[2]))
def mass_matrix_from_solid(obj):
    "extracts the mass matrix of the TopAbs_ShapeEnum.SOLID obj, using the figures computed by netgen"
    
    # copy the inertia matrix from netgen
    inertia_matrix = bla.Matrix(3,3)
    for i in range(3):
        for j in range(3):
            inertia_matrix[i, j] = obj.inertia[i, j]

    # copy the center of mass from netgen
    center_of_mass = bla.Vector(3)
    for i in range(3): center_of_mass[i] = obj.center[i]

    # rearrange it in C++ to make the mass matrix (the elegant way, using MatrixView)
    mass_matrix = mass_matrix_from_inertia(inertia_matrix, center_of_mass, obj.mass)
    return mass_matrix

mass_matrix = mass_matrix_from_solid(box)
# print(mass_matrix)
#for i in range(18):
#    for j in range(18):
#        print(F"{round(mass_matrix[i, j], 5)}, ", end="")
#    print("\\")
def extract_vertices(mesh: Mesh):
    "extracts a p3js compatible vertex list from a netgen.occ Mesh"
    combined_vertices = [] # flat array of vertices of all faces
    all_vertices = [vert.point for vert in mesh.vertices] # all available vertices

    # throw all vertices of all faces into combined_vertices
    for face in mesh.faces:
        for vertex in face.vertices:
            point = all_vertices[vertex.nr]
            combined_vertices.append(point)
            
    return combined_vertices
rbs = RBSystem()
# set up physics simulation object
r1 = RigidBody()
r1.setMass(mass_matrix)

p1 = Transformation()
p1.setTranslation(0,0,0)
p1.setRotation(0,0,1)
p1.setRotation(1,1,1)
p1.setRotation(2,2,1)
r1.q=p1

t1 = Transformation()
t1.setTranslation(0.001,0.0,0.0)
#t1.setRotation(0,2,0.001)
#t1.setRotation(2,0,-0.001)
r1.dq=t1



r2 = RigidBody()
r2.setMass(mass_matrix)

p2 = Transformation()
p2.setTranslation(3,3,3)
p2.setRotation(0,0,1)
p2.setRotation(1,1,1)
p2.setRotation(2,2,1)
r2.q=p2

t2 = Transformation()
#t2.setTranslation(0.000,0.1,0)
#t2.setRotation(0,2,0.02)
#t2.setRotation(2,0,-0.02)
r2.dq=t2

#Add body to system
c1 = rbs.addBody(r1)
c1.pos=(1,1,1)
c2 = rbs.addBody(r2)
c2.pos=(-1,-1,-1)
#b = Beam(c1,c2,1.732);
#rbs.addBeam(b)

s = Spring(c1,c2, 1.7320508075688772 ,10)
rbs.addSpring(s)

# save a restore point for the reset button
rbs.saveState()
rbs.simulate(3,100)