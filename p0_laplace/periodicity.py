from make_mesh_cpp import build_mesh
import dolfin as df
import numpy as np


def stitch_periodic(boundaries, mapping, tol=1E-13):
    '''Fake the mesh to be periodic'''
    mesh = boundaries.mesh()
    assert mesh.topology().dim()-1 == boundaries.dim()
    
    domain_tag, range_tag, mapping = mapping
    assert domain_tag != range_tag
    # These tags must be interior boundary tags ...
    dS = df.Measure('dS', domain=mesh, subdomain_data=boundaries)
    assert all(abs(df.assemble(df.Constant(1)*dS(tag))) < tol for tag in (domain_tag, range_tag))

    ds = df.Measure('ds', domain=mesh, subdomain_data=boundaries)    
    # ... that actually are in the mesh
    lens = [np.round(df.assemble(df.Constant(1)*ds(tag)), -int(np.log10(tol))) for tag in (domain_tag, range_tag)]
    val, = set(lens)
    assert val > 0

    # Now we want to build a mapping between vertices
    fdim = boundaries.dim()
    _, f2v = mesh.init(fdim, 0), mesh.topology()(fdim, 0)

    domain_vidx = np.unique(np.hstack([f2v(f.index()) for f in df.SubsetIterator(boundaries, domain_tag)]))
    range_vidx = np.unique(np.hstack([f2v(f.index()) for f in df.SubsetIterator(boundaries, range_tag)]))
    assert len(domain_vidx) == len(range_vidx)
    
    # Next we want to map the coordinates
    x = mesh.coordinates()
    domain_x, range_x = x[domain_vidx], x[range_vidx]
    range_y = mapping(domain_x)

    # Now we want range_vidx[mapping[i]] for domin_vidx[i]
    lookup = np.zeros(len(domain_vidx), dtype='uintp')
    for i, yi in enumerate(range_y):
        dist = np.linalg.norm(range_x - yi, 2, axis=1)
        j = np.argmin(dist)
        assert dist[j] < tol
        lookup[i] = j

    # And build a mesh -- in some sense we wrap around the mesh
    pidx = np.ones(len(x), dtype=bool)
    pidx[lookup] = False
    pidx, = np.where(pidx)
    new_coords = x[pidx]    # Ignore the range vertices
    
    vtx_map = {old: new for new, old in enumerate(pidx)}
    periodic = dict(zip(range_vidx, domain_vidx))    
    vtx_map.update({r: vtx_map[periodic[r]] for r in range_vidx})

    cells = mesh.cells()
    # Handle periodicity in old
    new_cells = (periodic.get(v, v) for v in cells.flat)
    # Now rewrite in new
    new_cells = np.fromiter((vtx_map[v] for v in new_cells),
                            dtype=cells.dtype).reshape(cells.shape)

    new_mesh = build_mesh(new_coords, new_cells, mesh.ufl_cell())
    dim = new_mesh.topology().dim()-1

    new_boundaries = type(boundaries)(new_mesh, dim, 0)
    new_boundaries_array = new_boundaries.array()
    
    # These facets as
    _, new_f2v = new_mesh.init(dim, 0), new_mesh.topology()(dim, 0)
    new_facets = {tuple(sorted(new_f2v(f))): f for f in range(new_mesh.num_entities(dim))}
    
    _, f2v = mesh.init(dim, 0), mesh.topology()(dim, 0)
    boundaries = boundaries.array()
    # Now we want to rebuild the facet function
    tags = np.unique(boundaries)
    for tag in tags:
        if tag == range_tag:
            continue
        
        tagged_facets, = np.where(boundaries == tag)
        # As vertex
        new_tagged_facets = [new_facets[tuple(sorted([vtx_map[v] for v in f2v(f)]))]
                             for f in tagged_facets]
        # Encode in new
        new_boundaries_array[new_tagged_facets] = tag

    return new_boundaries

# --------------------------------------------------------------------

if __name__ == '__main__':
    from dolfin import *
    
    left = CompiledSubDomain('near(x[0], 0)')
    right = CompiledSubDomain('near(x[0], 1)')
    bottom = CompiledSubDomain('near(x[1], 0)')
    top = CompiledSubDomain('near(x[1], 1)')
    subdomains = [left, right, bottom, top]

    mapping = (4, 3, lambda x: np.c_[x[:, 0], x[:, 1]-1])  
    
    mesh = UnitSquareMesh(1, 4)
    boundaries = MeshFunction('size_t', mesh, mesh.topology().dim()-1, 0)
    [subd.mark(boundaries, tag) for tag, subd in enumerate(subdomains, 1)]

    new_boundaries = stitch_periodic(boundaries, mapping, tol=1E-13)

    # Basically, our goal was that the mapping[0] tag is now view as an
    # interior facet so
    dS = Measure('dS', domain=new_boundaries.mesh(), subdomain_data=new_boundaries)
    print(assemble(Constant(1)*dS(4)), '>>>')

    dS = Measure('dS', domain=boundaries.mesh(), subdomain_data=boundaries)
    print(assemble(Constant(1)*dS(4)), '<<<')
    ds = Measure('ds', domain=boundaries.mesh(), subdomain_data=boundaries)
    print(assemble(Constant(1)*ds(4)), '<<<')
    
    
