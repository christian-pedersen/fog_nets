from dolfin import *


def CellCenterDistance(mesh, strict=True, debug=False):
    '''Discontinuous Lagrange Trace function that holds the cell-to-cell distance'''
    # Cell-cell distance for the interior facet is defined as a distance 
    # of midpoints of the cells that share the facet. For exterior facet
    # we take the distance of cell midpoint and the facet midpoint
    Q = FunctionSpace(mesh, 'DG', 0)
    V = FunctionSpace(mesh, 'CG', 1)
    L = FunctionSpace(mesh, 'Discontinuous Lagrange Trace', 0)
    
    cK, fK = CellVolume(mesh), FacetArea(mesh)
    q, l = TestFunction(Q), TestFunction(L)
    # The idea here to first assemble component by component the cell 
    # and (exterior) facet midpoint
    cell_centers, facet_centers = [], []
    for xi in SpatialCoordinate(mesh):
        qi = Function(Q)
        # Pretty much use the definition that a midpoint is int_{cell} x_i/vol(cell)
        # It's thanks to eval in q that we take values :)
        assemble((1/cK)*inner(xi, q)*dx, tensor=qi.vector())
        cell_centers.append(qi)
        # Same here but now our mean is over an edge
        li = Function(L)
        assemble((1/fK)*inner(xi, l)*ds, tensor=li.vector())
        facet_centers.append(li)
    # We build components to vectors
    cell_centers, facet_centers = map(as_vector, (cell_centers, facet_centers))
        
    distances = Function(L)
    # FIXME: This might not be necessary but it's better to be certain
    dS_, ds_ = dS(metadata={'quadrature_degree': 0}), ds(metadata={'quadrature_degree': 0})

    # Check whether vector connecting cell centers is orthognal to facet tangent,
    # i.e. colinear with
    if debug:
        n_interior = jump(cell_centers)/sqrt(dot(jump(cell_centers), jump(cell_centers)))
        n_exterior = (cell_centers-facet_centers)/sqrt(dot(cell_centers-facet_centers, cell_centers-facet_centers))
        n = FacetNormal(mesh)

        orthogonality = Function(L)
        assemble((1/fK('+'))*inner(dot(n_interior, n('+')), l('+'))*dS_+
                 (1/fK)*inner(dot(n_exterior, n), l)*ds_,
                 orthogonality.vector())

        og = orthogonality.vector().get_local()
        n_notog = sum(np.abs(np.abs(og) - 1) > 1E-10)
        print('There are {} facets not orthogonal to cell-cell vectors'.format(n_notog))
    
    # Finally we assemble magniture of the vector that is determined by the
    # two centers
    assemble(((1/fK('+'))*inner(sqrt(dot(jump(cell_centers), jump(cell_centers))), l('+'))*dS_+
              (1/fK)*inner(sqrt(dot(cell_centers-facet_centers, cell_centers-facet_centers)), l)*ds_),
             distances.vector())
        
    return distances    


def solve_poisson(boundaries, f, dir_data, neu_data):
    '''Poisson P0 solver'''
    mesh = boundaries.mesh()

    V = FunctionSpace(mesh, 'DG', 0)
    u = TrialFunction(V)
    v = TestFunction(V)

    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

    hVof = CellCenterDistance(mesh)    
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

        uh = solve_poisson(boundaries, 
                           f=f_true, dir_data=dir_data, neu_data=neu_data)
        
        print(f'|u-uh|_0 = {errornorm(u_true, uh, "L2"):.4E}')
