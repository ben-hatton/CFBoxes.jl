module CFBoxes

# boundary topology exports
export Periodic, Bounded

# boundary condition exports
export BCType, DirichletBCType, NeumannBCType, PeriodicBCType
export bc_type
export AbstractBoundaryConditionSide, FunctionalBoundaryConditionSide, PeriodicBoundaryConditionSide
export NeumannBC, DirichletBC, PeriodicBC
export BoundaryCondition2D, DiscreteBoundaryCondition2D
export BoundaryConditions2D
export check_bcs, display_bcs

# box exports
export AbstractBox, AbstractBoxHalo, Box2D, Box3D
export xdim, ydim, zdim, Xdim, Ydim, Zdim
export xrange, yrange, zrange, Xrange, Yrange, Zrange
export xrange_interior, yrange_interior, zrange_interior
export Xrange_interior, Yrange_interior, Zrange_interior
export xgrid, ygrid, zgrid, Xgrid, Ygrid, Zgrid
export xgrid_interior, ygrid_interior, zgrid_interior
export Xgrid_interior, Ygrid_interior, Zgrid_interior
export dims, halo_size, boundary_topology
export xpoint, ypoint, zpoint, Xpoint, Ypoint, Zpoint
export alloc_xz, alloc_Xz, alloc_xZ, alloc_XZ
export interior, add_halo

# centered operator exports
export avg_x, avg_z, avg_X, avg_Z
export avg_xz, avg_xZ, avg_Xz, avg_XZ
export dif_x, dif_z, dif_X, dif_Z

# advection scheme exports
export AdvectionScheme, AdvEnergyCons, AdvEnstrophyCons

include("julia/boundary_topology.jl")
include("julia/boxes.jl")
include("julia/boundary_conditions.jl")
include("julia/centered_operators.jl")
include("julia/advection_schemes.jl")

end
