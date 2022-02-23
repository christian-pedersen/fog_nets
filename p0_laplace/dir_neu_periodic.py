from dolfin import *
from dir_neu_periodic_mixed import Top2Bottom, Right2Left
from periodicity import stitch_periodic
from dir_neu import CellCenterDistance
import numpy as np


def solve_poisson(boundaries, f, dir_data, neu_data, periodic_data):
    '''Poisson Hdiv x P0 solver'''
    periodic_data = tuple(sorted(periodic_data))
    assert periodic_data in ((3, 4), (1, 2))

    if periodic_data == (3, 4):
        mapping = (4, 3, lambda x: np.c_[x[:, 0], x[:, 1]-1])  
    else:
        raise ValueError

    boundaries = stitch_periodic(boundaries, mapping, tol=1E-13)
    mesh = boundaries.mesh()
    cell = mesh.ufl_cell()
    
    Velm = FiniteElement('Discontinuous Lagrange', cell, 0)        
    print(f'Unconstrained dim(W)) = {FunctionSpace(mesh, Velm).dim()}')
    V = FunctionSpace(mesh, Velm)
    print(f'Constrained dim(V) = {V.dim()}')
    
    u = TrialFunction(V)
    v = TestFunction(V)

    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

    # FIXME: this doesn't compute what we need on periodic
    hVof = CellCenterDistance(mesh)
    print(np.sort(hVof.vector().get_local()))
    
    a = (avg((1/hVof))*inner(jump(u), jump(v))*dS
         +sum((1/hVof)*inner(u, v)*ds(tag) for tag in dir_data))

    n = FacetNormal(mesh)
    L = (inner(f, v)*dx
         +sum(inner(v, dot(neu_data[tag], n))*ds(tag) for tag in neu_data))
    # DG0 we need to account for Dirichlet bcs which are enforced in the form
    L += sum((1/hVof)*inner(dir_data[tag], v)*ds(tag) for tag in dir_tags)    

    A, b = map(assemble, (a, L))

    uh = Function(V)
    solve(A, uh.vector(), b)

    return uh

# --------------------------------------------------------------------

if __name__ == '__main__':
    # Manufactured problem, just to show that with Dirichlet Neumann we
    # converge
    import ulfy

    mesh = UnitSquareMesh(2, 2)
    x, y = SpatialCoordinate(mesh)

    u = sin(pi*x*(x-1)*y*(y-1))  # Works with sine (only because of null Neumann then)
    f = -div(grad(u))
    grad_u = grad(u)

    as_expr = lambda v: ulfy.Expression(v, degree=4)
    u_true, f_true, grad_u_true = map(as_expr, (u, f, grad_u))
    
    # Geometry
    left = CompiledSubDomain('near(x[0], 0)')
    right = CompiledSubDomain('near(x[0], 1)')
    bottom = CompiledSubDomain('near(x[1], 0)')
    top = CompiledSubDomain('near(x[1], 1)')
    subdomains = [left, right, bottom, top]

    periodic_data = (3, 4)  # 
    
    dir_tags = (1, )
    dir_data = {tag: u_true for tag in dir_tags}    
    # Neumann is the rest
    neu_data = {tag: grad_u_true for tag in ({1, 2, 3, 4} - set(periodic_data + dir_tags))}

    for n in (4, 8, 16, 32, 64):
        mesh = UnitSquareMesh(n, n, 'crossed')

        boundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
        [subd.mark(boundaries, tag) for tag, subd in enumerate(subdomains, 1)]

        uh = solve_poisson(boundaries, f=f_true,
                           dir_data=dir_data, neu_data=neu_data, periodic_data=periodic_data)
        
        print(f'|u-uh|_0 = {errornorm(u_true, uh, "L2"):.4E}')
