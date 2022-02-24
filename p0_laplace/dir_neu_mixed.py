from dolfin import *


def solve_poisson_mixed(boundaries, f, dir_data, neu_data):
    '''Poisson Hdiv x P0 solver'''
    mesh = boundaries.mesh()
    cell = mesh.ufl_cell()
    
    Welm = MixedElement([FiniteElement('Raviart-Thomas', cell, 1),
                         FiniteElement('Discontinuous Lagrange', cell, 0)])
    W = FunctionSpace(mesh, Welm)
    u, p = TrialFunctions(W)
    v, q = TestFunctions(W)

    n = FacetNormal(mesh)
    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
    
    a = inner(u, v)*dx + inner(p, div(v))*dx + inner(q, div(u))*dx
    L = -inner(f, q)*dx + sum(inner(dir_data[tag], dot(v, n))*ds(tag) for tag in dir_tags)

    bcs = [DirichletBC(W.sub(0), value, boundaries, tag)
           for tag, value in neu_data.items()]

    A, b = assemble_system(a, L, bcs)

    wh = Function(W)
    solve(A, wh.vector(), b)

    _, ph = wh.split(deepcopy=True)

    return ph

# --------------------------------------------------------------------

if __name__ == '__main__':
    # Manufactured problem, just to show that with Dirichlet Neumann we
    # converge
    import ulfy

    mesh = UnitSquareMesh(2, 2)
    x, y = SpatialCoordinate(mesh)

    u = sin(pi*(x-y))
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

    dir_tags = (1, 2, 3)
    dir_data = {tag: u_true for tag in dir_tags}    
    # Neumann is the rest
    neu_data = {tag: grad_u_true for tag in (set((1, 2, 3, 4)) - set(dir_tags))}    
    for n in (4, 8, 16, 32, 64):
        mesh = UnitSquareMesh(n, n, 'crossed')

        boundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
        [subd.mark(boundaries, tag) for tag, subd in enumerate(subdomains, 1)]

        uh = solve_poisson_mixed(boundaries, 
                                 f=f_true, dir_data=dir_data, neu_data=neu_data)
        
        print(f'|u-uh|_0 = {errornorm(u_true, uh, "L2"):.4E}')
