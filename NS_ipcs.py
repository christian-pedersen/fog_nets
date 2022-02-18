from dolfin import *
import numpy as np
from LagrangianParticles import LagrangianParticles

# constants

L = 3.5
H = 2.
rho = Constant(1.0)
nu = Constant(0.001)
mu = rho*nu
U_m = 1.5
U_0 = 1.0
x0, y0 = 0.2, 0.2
#cylinder_radius = 0.05
U_ = 2.0*U_m / 3
cylinder_radius = 0.25
xa = np.array([0.15, 0.2])
xe = np.array([0.25, 0.2])
dt = 0.00005
end_time = 20

folder = '_results/U=%g/'%U_0

# create mesh and boundaries
mesh = Mesh('mesh/test_mesh_copy.xml')
#bc1 = MeshFunction("size_t", mesh, "mesh/test_mesh_physical_region.xml")

class Periodic_sides(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0) and on_boundary #

    def map(self, x, y):
        y[0] = x[0]
        y[1] = x[1] - H

class Top(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], H) and on_boundary

class Inlet(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], 0) and on_boundary

class Outlet(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[0], L) and on_boundary

class Cylinders(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > 0 and x[0] < (L)  and x[1] > 0 and x[1] < H and on_boundary

boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
boundaries.set_all(0)

cylinders = Cylinders()
cylinders.mark(boundaries, 4)
top = Top()
top.mark(boundaries, 0)
inlet = Inlet()
inlet.mark(boundaries, 1)
outlet = Outlet()
outlet.mark(boundaries, 2)
periodic_sides = Periodic_sides()
periodic_sides.mark(boundaries, 3)

domain = File('domain.pvd')
domain << boundaries
#plot(boundaries, interactive=True)



# function spaces and functions

V = VectorFunctionSpace(mesh, 'CR', 1, constrained_domain=Periodic_sides())
Q = FunctionSpace(mesh, 'DG', 0, constrained_domain=Periodic_sides())

u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

u0 = Function(V)
p0 = Function(Q)
ustar = Function(V)
u1 = Function(V)
p1 = Function(Q)

uhalf = 0.5*(u0+u)

def epsilon(u):
    return 0.5*(grad(u) + grad(u).T)
def sigma(u, p):
    """
        Cauchy stress tensor
    """
    return 2*mu*epsilon(u) - p*Identity(len(u))

n = FacetNormal(mesh)

# boundary conditions

inlet_velocity = Constant((U_0, 0.))
outlet_pressure = Constant(0)
solids_velocity = Constant((0, 0))

inletBC = DirichletBC(V, inlet_velocity, boundaries, 1)
outletBC = DirichletBC(Q, outlet_pressure, boundaries, 2)
#noslip_walls = DirichletBC(V, solids_velocity, walls)
noslip_cylinder = DirichletBC(V, solids_velocity, boundaries, 4)
#noslip_cylinder = DirichletBC(V, solids_velocity, bc1)

velocityBC = [noslip_cylinder, inletBC]
pressureBC = [outletBC]
"""
# equation for tentative velocity
eq1 = rho*dot((u-u0)/dt, v)*dx + rho*dot(dot(u0, grad(u0)), v)*dx \
    + inner(sigma(uhalf, p0), epsilon(v))*dx - mu*dot(grad(uhalf).T*n, v)*ds \
    + inner(p0*n, v)*ds
eq1 = rho*dot((u-u0)/dt, v)*dx + rho*dot(grad(u0)*u0, v)*dx \
     - mu*inner(grad(uhalf).T*n, v)*ds + dot(p0*n, v)*ds \
     + inner(sigma(uhalf, p0), epsilon(v))*dx
a1 = lhs(eq1)
L1 = rhs(eq1)
# equation for pressure correction
a2 = dt*dot(grad(p), grad(q))*dx
L2 = dt*dot(grad(p0), grad(q))*dx - rho*div(ustar)*q*dx
# equation for velocity correction
a3 = rho*inner(u1, v)*dx
L3 = rho*inner(ustar, v)*dx - dt*inner(grad(p1- p0), v)*dx
"""


# tentative velocity equation
eq1 = rho*dot((u-u0)/dt, v)*dx + rho*dot(grad(u0)*u0, v)*dx \
     - mu*inner(grad(uhalf).T*n, v)*ds + dot(p0*n, v)*ds \
     + inner(sigma(uhalf, p0), epsilon(v))*dx
a1 = lhs(eq1)
L1 = rhs(eq1)

# pressure correction equation
a2 = dt*dot(grad(p), grad(q))*dx
L2 = dt*dot(grad(p0), grad(q))*dx - rho*div(u1)*q*dx

# corrected velocity equation
a3 = rho*dot(u, v)*dx
L3 = rho*dot(u1, v)*dx + dot(dt*grad(p0-p1), v)*dx

A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

t = dt

u_ = XDMFFile(folder+'u.xdmf')
p_ = XDMFFile(folder+'p.xdmf')

u_.parameters['flush_output'] = True
p_.parameters['flush_output'] = True


# _u = File('results_coarse/velocity.pvd')
# _p = File('results_coarse/pressure.pvd')

pix = 0

#results = open('results_coarse.txt', 'w')
#results.write('t\tcd\tcl\tdp\n')

ds = Measure('ds', subdomain_data=boundaries)

while t < end_time +DOLFIN_EPS:
    print('calculation progress:',(t-dt)/end_time*100,'%')

    # calculate tentative velocity
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in velocityBC]
    solve(A1, u1.vector(), b1)
    #solve(A1 == b1, ustar, velocityBC)

    # correct pressure
    b2 = assemble(L2)
    [bc.apply(A2, b2) for bc in pressureBC]
    solve(A2, p1.vector(), b2)

    # correct velocity
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in velocityBC]
    solve(A3, u1.vector(), b3)


    # calculating drag/lift coefficients
    tau = -p1*Identity(2) + rho*nu*(grad(u1) + grad(u1).T)

    force = dot(tau, n)

    Drag = -assemble(force[0]*ds(4))
    Lift = -assemble(force[1]*ds(4))

    CD = 2*Drag / (1.0*(2.0/3.0*1.5)**2*0.1)#(rho*Umean**2*cylinder_diameter)
    CL = 2*Lift / (1.0*(2.0/3.0*1.5)**2*0.1) #(rho*Umean**2*cylinder_diameter)

    # calculating pressure diffence
    # between front and back of cylinder

    #delta_pressure = p1(xa) - p1(xe)




    if pix == 10:
        pix = 0
    if pix == 0:
        u1.rename('u', 'f')
        u_.write(u1, t)
        p1.rename('p', 'f')
        p_.write(p1, t)
    pix += 1
    t += dt

    #results.write('%s\t%s\t%s\t%s\n' % (str(t), str(CD), str(CL), str(delta_pressure)))

    u0.assign(u1)
    p0.assign(p1)

results.close()
#plot(u1, interactive=True)
