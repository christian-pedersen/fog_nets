# Periodic bcs for pressure solver with P0 elements

The pressure solver should be a component in the seggregated N-S solver where
velocity should be in CR (as we are trying to avoid discretizations with dofs
tied to vertices).

- Basis idea of the P0 solver is shown in `dir_neu.py`; we can handle Dirichlet and Neuamann bcs but not
periodic (this is a FEniCS limitation)
- So the alternative is to consider a mixed formulation where we can enforce periodicity of the flux. The
solver works, however, we are solving for an additional vector field
- It should be possible to FEniCS such that we can solve only for pressure in P0

# TODO
-[ ] What is the CellCenterDistance on the `stiched` mesh
-[ ] Related to above; fix the periodic solver
-[ ] What is the mapping between the P0 spaces on stitched/un-stitched meshes?