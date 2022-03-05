from dolfin import *
import numpy as np
from LagrangianParticles import LagrangianParticles
import matplotlib.pyplot as plt

# constants

radius = 0.15
x_dist = 0.2

x_center = (5+2)*radius + x_dist

rho = Constant(1.0)
nu = Constant(0.001)
mu = rho*nu
U_0 = 1.0
dt = 0.0005
end_time = 20

seed_chance = np.asarray([0,1]) # 50% chance of seeding per time step

folder = '_results/U=%g/'%U_0

# create mesh and boundaries
mesh = Mesh('mesh/mesh.xml')
x = mesh.coordinates()
L, H = max(x[:,0]), max(x[:,1])


class Periodic_sides(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] < 1e-8 and on_boundary #

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

class Interior_Cylinders(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > DOLFIN_EPS and x[0] < L and x[1] > DOLFIN_EPS and x[1] < (H  -0.1) and on_boundary

class Boundary_Cylinders(SubDomain):
    def inside(self, x, on_boundary):
        return (x[0]-x_center)*(x[0]-x_center) + x[1]*x[1] - radius*radius < 1e-08  or (x[0]-x_center)*(x[0]-x_center) + (x[1]-H)*(x[1]-H) - radius*radius < 1e-8  and on_boundary


boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
boundaries.set_all(0)


cylinders = Interior_Cylinders()
cylinders.mark(boundaries, 4)
cylinders = Boundary_Cylinders()
cylinders.mark(boundaries, 4)
periodic_sides = Periodic_sides()
periodic_sides.mark(boundaries, 3)
top = Top()
top.mark(boundaries, 0)
inlet = Inlet()
inlet.mark(boundaries, 1)
outlet = Outlet()
outlet.mark(boundaries, 2)

domain = File('domain.pvd')
domain << boundaries
#plot(boundaries, interactive=True)



# function spaces and functions

V = VectorFunctionSpace(mesh, 'CG', 2, constrained_domain=Periodic_sides())
Q = FunctionSpace(mesh, 'CG', 1, constrained_domain=Periodic_sides())

lp = LagrangianParticles(V)

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
noslip_cylinder = DirichletBC(V, solids_velocity, boundaries, 4)


velocityBC = [noslip_cylinder, inletBC]
pressureBC = [outletBC]


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


pix = 0

step = 0
fig = plt.figure()
plt.ion()

while t < end_time +DOLFIN_EPS:
    print('calculation progress:',(t-dt)/end_time*100,'%')

    # calculate tentative velocity
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in velocityBC]
    solve(A1, u1.vector(), b1)

    # correct pressure
    b2 = assemble(L2)
    [bc.apply(A2, b2) for bc in pressureBC]
    solve(A2, p1.vector(), b2)

    # correct velocity
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in velocityBC]
    solve(A3, u1.vector(), b3)

    
    lp.step(u1, dt=dt)
    if np.random.choice(seed_chance) != 0:
        lp.add_particles(np.asarray([[0, H*np.random.uniform(0.001, 0.999)]]))
    step += 1
    print(lp.total_number_of_particles(), '<<< particles')

    lp.scatter(fig)
    fig.suptitle('At step %d' % step)
    fig.canvas.draw()
    plt.pause(0.0001)

    fig.clf()
    

    if pix == 10:
        u1.rename('u', 'f')
        u_.write(u1, t)
        p1.rename('p', 'f')
        p_.write(p1, t)
        pix = 0
    pix += 1
    t += dt

    u0.assign(u1)
    p0.assign(p1)
