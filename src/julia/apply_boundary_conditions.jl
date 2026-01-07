"""
    apply_bc_xz!(model, bc::BoundaryCondition2D, data)

Apply boundary conditions to a scalar array positioned on cell centre (xz).
"""
function apply_bc_xz!(model, bc::BoundaryCondition2D, data::AbstractArray)
    apply_bc_xz_left!(model, model.domain, data, bc.left)
    apply_bc_xz_right!(model, model.domain, data, bc.right)
    apply_bc_xz_top!(model, model.domain, data, bc.top)
    apply_bc_xz_bottom!(model, model.domain, data, bc.bottom)
end

"""
    apply_bc_xZ!(model, bc::BoundaryCondition2D, data)

Apply boundary conditions to a scalar array positioned on lateral edge (Xz)
For example, for horizontal velocity `u` on C-grid.
"""
function apply_bc_Xz!(model, bc::BoundaryCondition2D, data::AbstractArray)
    apply_bc_Xz_left!(model, model.domain, data, bc.left)
    apply_bc_Xz_right!(model, model.domain, data, bc.right)
    apply_bc_Xz_top!(model, model.domain, data, bc.top)
    apply_bc_Xz_bottom!(model, model.domain, data, bc.bottom)
end

"""
    apply_bc_xZ!(model, bc::BoundaryCondition2D, data)

Apply boundary conditions to a scalar array positioned on vertical edge (xZ)
For example, for vertical velocity `w` on C-grid.
"""
function apply_bc_xZ!(model, bc::BoundaryCondition2D, data::AbstractArray)
    apply_bc_xZ_bottom!(model, model.domain, data, bc.bottom)
    apply_bc_xZ_top!(model, model.domain, data, bc.top)
    apply_bc_xZ_left!(model, model.domain, data, bc.left)
    apply_bc_xZ_right!(model, model.domain, data, bc.right)
end

# to implement: apply_bc_XZ!

"""
    apply_bc_xz_left!(model, domain, data, bc_left)
    apply_bc_xz_right!(model, domain, data, bc_right)
    apply_bc_xz_top!(model, domain, data, bc_top)
    apply_bc_xz_bottom!(model, domain, data, bc_bottom)

Boundary-conditions for cell-centre (xz) quantities on each side, using centred-scheme extrapolation where needed.

These methods are dispatched on the boundary-condition type (`DirichletBCType`, `NeumannBCType`, `PeriodicBCType`)

Model `model` is assumed to include 
- `mgr` for loop management
- `dx`, `dz` for grid spacing
"""
function apply_bc_xz_left!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_left::FunctionalBoundaryConditionSide{DirichletBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dz) = model
    @with mgr, let jrange = axes(data, 2)
        for j in jrange
            z = zpoint(domain, dz, j)
            data[Hx, j] = 2 * bc_left.bc_func(z) - data[Hx+1, j]
        end
    end
end
function apply_bc_xz_left!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_left::FunctionalBoundaryConditionSide{NeumannBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dx, dz) = model
    @with mgr, let jrange = axes(data, 2)
        for j in jrange
            z = zpoint(domain, dz, j)
            data[Hx, j] = - bc_left.bc_func(z) * dx + data[Hx+1, j]
        end
    end
end
function apply_bc_xz_left!(model, domain::Box2D, data, ::AbstractBoundaryConditionSide{PeriodicBCType})
    apply_periodic_left!(model, domain, data)
end

function apply_bc_xz_right!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_right::FunctionalBoundaryConditionSide{DirichletBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dz) = model
    @with mgr, let jrange = axes(data, 2)
        for j in jrange
            z = zpoint(domain, dz, j)
            data[Mx+Hx+1, j] = 2 * bc_right.bc_func(z) - data[Mx+Hx, j]
        end
    end
end
function apply_bc_xz_right!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_right::FunctionalBoundaryConditionSide{NeumannBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dx, dz) = model
    @with mgr, let jrange = axes(data, 2)
        for j in jrange
            z = zpoint(domain, dz, j)
            data[Mx+Hx+1, j] = - bc_right.bc_func(z) * dx + data[Mx+Hx, j]
        end
    end
end
function apply_bc_xz_right!(model, domain::Box2D, data, ::AbstractBoundaryConditionSide{PeriodicBCType})
    apply_periodic_right!(model, domain, data)
end 

function apply_bc_xz_bottom!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_bottom::FunctionalBoundaryConditionSide{DirichletBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dx) = model
    @with mgr, let irange = axes(data, 1)
        for i in irange
            x = xpoint(domain, dx, i)
            data[i, Hz] = 2 * bc_bottom.bc_func(x) - data[i, Hz+1]
        end
    end
end
function apply_bc_xz_bottom!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_bottom::FunctionalBoundaryConditionSide{NeumannBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dx, dz) = model
    @with mgr, let irange = axes(data, 1)
        for i in irange
            x = xpoint(domain, dx, i)
            data[i, Hz] = - bc_bottom.bc_func(x) * dz + data[i, Hz+1]
        end
    end
end
function apply_bc_xz_bottom!(model, domain::Box2D, data, ::AbstractBoundaryConditionSide{PeriodicBCType})
    apply_periodic_bottom!(model, domain, data)
end

function apply_bc_xz_top!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_top::FunctionalBoundaryConditionSide{DirichletBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dx) = model
    @with mgr, let irange = axes(data, 1)
        for i in irange
            x = xpoint(domain, dx, i)
            data[i, Mz+Hz+1] = 2 * bc_top.bc_func(x) - data[i, Mz+Hz]
        end
    end
end
function apply_bc_xz_top!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_top::FunctionalBoundaryConditionSide{NeumannBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dx, dz) = model
    @with mgr, let irange = axes(data, 1)
        for i in irange
            x = xpoint(domain, dx, i)
            data[i, Mz+Hz+1] = bc_top.bc_func(x) * dz + data[i, Mz+Hz]
        end
    end
end
function apply_bc_xz_top!(model, domain::Box2D, data, ::AbstractBoundaryConditionSide{PeriodicBCType})
    apply_periodic_top!(model, domain, data)
end

"""
    apply_bc_Xz_left!(...)
    apply_bc_Xz_right!(...)
    apply_bc_Xz_top!(...)
    apply_bc_Xz_bottom!(...)

Boundary-conditions for lateral edge quantities (Xz) on each side, using centred-scheme extrapolation where needed.

These methods are dispatched on the boundary-condition type (`DirichletBCType`, `NeumannBCType`, `PeriodicBCType`)

Model `model` is assumed to include
- `mgr` for loop management
- `dx`, `dz` for grid spacing
"""
function apply_bc_Xz_left!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_left::FunctionalBoundaryConditionSide{DirichletBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dz) = model
    @with mgr, let (irange, jrange) = (1:Hx+1, axes(data, 2))
        for i in irange, j in jrange
            z = zpoint(domain, dz, j)
            data[i, j] = bc_left.bc_func(z)
        end
    end
end
function apply_bc_Xz_left!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_left::FunctionalBoundaryConditionSide{NeumannBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dx, dz) = model
    @with mgr, let jrange = axes(data, 2)
        for j in jrange
            z = zpoint(domain, dz, j)
            data[Hx, j] = - 2 * bc_left.bc_func(z) * dx + data[Hx+2, j]
        end
    end
end
function apply_bc_Xz_left!(model, domain::Box2D, data, ::AbstractBoundaryConditionSide{PeriodicBCType})
    apply_periodic_left!(model, domain, data)
end

function apply_bc_Xz_right!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_right::FunctionalBoundaryConditionSide{DirichletBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dz) = model
    @with mgr, let (irange, jrange) = (1:Hx+1, axes(data, 2))
        for i in irange, j in jrange
            z = zpoint(domain, dz, j)
            data[Mx+Hx+i, j] = bc_right.bc_func(z)
        end
    end
end
function apply_bc_Xz_right!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_right::FunctionalBoundaryConditionSide{NeumannBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dx, dz) = model
    @with mgr, let jrange = axes(data, 2)
        for j in jrange
            z = zpoint(domain, dz, j)
            data[Mx+Hx+2, j] = 2 * bc_right.bc_func(z) * dx + data[Mx+Hx, j]
        end
    end
end
function apply_bc_Xz_right!(model, domain::Box2D, data, ::AbstractBoundaryConditionSide{PeriodicBCType})
    apply_periodic_right!(model, domain, data)
end

function apply_bc_Xz_bottom!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_bottom::FunctionalBoundaryConditionSide{DirichletBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dx) = model
    @with mgr, let irange = axes(data, 1)
        for I in irange
            X = Xpoint(domain, dx, I)
            data[I, Hz] = 2 * bc_bottom.bc_func(X) - data[I, Hz+1]
        end
    end
end
function apply_bc_Xz_bottom!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_bottom::FunctionalBoundaryConditionSide{NeumannBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dx, dz) = model
    @with mgr, let irange = axes(data, 1)
        for I in irange
            X = Xpoint(domain, dx, I)
            data[I, Hz] = - bc_bottom.bc_func(X) * dz + data[I, Hz+1]
        end
    end 
end
function apply_bc_Xz_bottom!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, ::AbstractBoundaryConditionSide{PeriodicBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    apply_periodic_bottom!(model, domain, data)
end

function apply_bc_Xz_top!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_top::FunctionalBoundaryConditionSide{DirichletBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dx) = model
    @with mgr, let irange = axes(data, 1)
        for I in irange
            X = Xpoint(domain, dx, I)
            data[I, Mz+Hz+1] = 2 * bc_top.bc_func(X) - data[I, Mz+Hz]
        end
    end
end
function apply_bc_Xz_top!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_top::FunctionalBoundaryConditionSide{NeumannBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dx, dz) = model
    @with mgr, let irange = axes(data, 1)
        for I in irange
            X = Xpoint(domain, dx, I)
            data[I, Mz+Hz+1] = - bc_top.bc_func(X) * dz + data[I, Mz+Hz]
        end
    end
end
function apply_bc_Xz_top!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, ::AbstractBoundaryConditionSide{PeriodicBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    apply_periodic_top!(model, domain, data)
end

"""
    apply_bc_xZ_left!(...)
    apply_bc_xZ_right!(...)
    apply_bc_xZ_top!(...)
    apply_bc_xZ_bottom!(...)

Boundary-conditions for `w` on each side, using a centred-scheme extrapolation where needed.

These methods are dispatched on the boundary-condition type (`DirichletBCType`, `NeumannBCType`, `PeriodicBCType`)
"""
function apply_bc_xZ_left!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_left::FunctionalBoundaryConditionSide{DirichletBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dz) = model
    @with mgr, let jrange = axes(data, 2)
        for J in jrange
            Z = zpoint(domain, dz, J)
            data[Hx, J] = 2 * bc_left.bc_func(Z) - data[Hx+1, J]
        end
    end
end
function apply_bc_xZ_left!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_left::FunctionalBoundaryConditionSide{NeumannBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dx, dz) = model
    @with mgr, let jrange = axes(data, 2)
        for J in jrange
            Z = Zpoint(domain, dz, J)
            data[Hx, J] = - bc_left.bc_func(Z) * dx + data[Hx+1, J]
        end
    end
end
function apply_bc_xZ_left!(model, domain::Box2D, data, ::AbstractBoundaryConditionSide{PeriodicBCType})
    apply_periodic_left!(model, domain, data)
end

function apply_bc_xZ_right!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_right::FunctionalBoundaryConditionSide{DirichletBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dz) = model
    @with mgr, let jrange = axes(data, 2)
        for J in jrange
            Z = zpoint(domain, dz, J)
            data[Mx+Hx+1, J] = 2 * bc_right.bc_func(Z) - data[Mx+Hx, J]
        end
    end
end
function apply_bc_xZ_right!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_right::FunctionalBoundaryConditionSide{NeumannBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dx, dz) = model
    @with mgr, let jrange = axes(data, 2)
        for J in jrange
            Z = Zpoint(domain, dz, J)
            data[Mx+Hx+1, J] = bc_right.bc_func(Z) * dx + data[Mx+Hx, J]
        end
    end
end
function apply_bc_xZ_right!(model, domain::Box2D, data, ::AbstractBoundaryConditionSide{PeriodicBCType})
    apply_periodic_right!(model, domain, data)
end

function apply_bc_xZ_bottom!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_bottom::FunctionalBoundaryConditionSide{DirichletBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dx) = model
    @with mgr, let (irange, jrange) = (axes(data, 1), 1:Hz+1)
        for i in irange, j in jrange
            x = xpoint(domain, dx, i)
            data[i, j] = bc_bottom.bc_func(x)
        end
    end
end
function apply_bc_xZ_bottom!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_bottom::FunctionalBoundaryConditionSide{NeumannBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dx, dz) = model
    @with mgr, let irange = axes(data, 1)
        for i in irange
            x = xpoint(domain, dx, i)
            data[i, Hz] = -2 * bc_bottom.bc_func(x) * dz + data[i, Hz+2]
        end
    end
end
function apply_bc_xZ_bottom!(model, domain::Box2D, data, ::AbstractBoundaryConditionSide{PeriodicBCType})
    apply_periodic_bottom!(model, domain, data)
end 

function apply_bc_xZ_top!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_top::FunctionalBoundaryConditionSide{DirichletBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dx) = model
    @with mgr, let (irange, jrange) = (axes(data, 1), 1:Hz+1)
        for i in irange, j in jrange
            x = xpoint(domain, dx, i)
            data[i, Mz+Hz+j] = bc_top.bc_func(x)
        end
    end
end
function apply_bc_xZ_top!(model, domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data, bc_top::FunctionalBoundaryConditionSide{NeumannBCType}) where {Mx, Mz, Hx, Hz, Bx, Bz}
    (; mgr, dx, dz) = model
    @with mgr, let irange = axes(data, 1)
        for i in irange
            x = xpoint(domain, dx, i)
            data[i, Mz+Hz+2] = 2 * bc_top.bc_func(x) * dz + data[i, Mz+Hz]
        end
    end
end
function apply_bc_xZ_top!(model, domain::Box2D, data, ::AbstractBoundaryConditionSide{PeriodicBCType})
    apply_periodic_top!(model, domain, data)
end

"""
    apply_periodic_left!(mgr, domain, data)
    apply_periodic_right!(mgr, domain, data)
    apply_periodic_bottom!(mgr, domain, data)
    apply_periodic_top!(mgr, domain, data)

Functions to apply periodic boundary conditions to state/dstate on each side of a 2D domain.

Inputs:
- `mgr`: LoopManager object containing parallelization info
- `domain`: Box domain object
- `data`: Data array to which the periodic BCs are applied
"""
function apply_periodic_left!(model, ::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data) where {Mx, Mz, Hx, Hz, Bx, Bz}
    @with model.mgr, let (irange, jrange) = (1:Hx, axes(data, 2))
        @vec for i in irange, j in jrange                      
            data[i, j] = data[i+Mx, j]
        end
    end
end
function apply_periodic_right!(model, ::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data) where {Mx, Mz, Hx, Hz, Bx, Bz}
    @with model.mgr, let (irange, jrange) = (1:Hx, axes(data, 2))
        @vec for i in irange, j in jrange                      
            data[i+Mx+Hx, j] = data[i+Hx, j]
        end
    end
end
function apply_periodic_bottom!(model, ::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data) where {Mx, Mz, Hx, Hz, Bx, Bz}
    @with model.mgr, let (irange, jrange) = (axes(data, 1), 1:Hz)
        @vec for i in irange, j in jrange                      
            data[i, j+Mz+Hz] = data[i, j+Hz]
        end
    end
end
function apply_periodic_top!(model, ::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, data) where {Mx, Mz, Hx, Hz, Bx, Bz}
    @with model.mgr, let (irange, jrange) = (axes(data, 1), 1:Hz)
        @vec for i in irange, j in jrange                      
            data[i, j] = data[i, j+Mz]
        end
    end
end

"""
    periodize!(model, data)
    periodize!(model, tuple::Tuple)
Periodization of general arrays according to the model's boundary topology.
"""
function periodize!(model, data)
    return periodize!(model.mgr, model.domain, boundary_topology(model.domain), data)
end

"""
    periodize(mgr, domain, boundary_topology, data)
Periodization of general arrays according to the specified boundary topology.
"""
function periodize!(mgr, ::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, boundary::Tuple{Bounded, Bounded}, data) where {Mx, Mz, Hx, Hz, Bx, Bz}
    return nothing
end
function periodize!(mgr, ::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, boundary::Tuple{Periodic, Bounded}, data) where {Mx, Mz, Hx, Hz, Bx, Bz}
    @with mgr, let (irange, jrange) = (1:Hx, axes(data, 2))
        @vec for i in irange, j in jrange                      
            data[i, j]      = data[i+Mx, j]
            data[i+Mx+Hx, j]  = data[i+Hx, j]
        end
    end
end
function periodize!(mgr, ::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}, boundary::Tuple{Bounded, Periodic}, data) where {Mx, Mz, Hx, Hz, Bx, Bz}
    @with mgr, let (irange, jrange) = (axes(data, 1), 1:Hz)
        @vec for i in irange, j in jrange                      
            data[i, j]      = data[i, j+Mz]
            data[i, j+Mz+Hz]  = data[i, j+Hz]
        end
    end
end
periodize!(model, tuple::Tuple) = unwrapper!(periodize!, model, tuple)

"""
    unwrapper!(func, model, tuple::Tuple)
    unwrapper!(func, model, tuple1::Tuple, tuple2::Tuple)

Utility functions to unwrap tuples and apply `func` to each element.
"""
function unwrapper!(func, model, tuple::Tuple)
    foreach(tuple) do m
        func(model, m)
    end
end
function unwrapper!(func, model, tuple1::Tuple, tuple2::Tuple)
    @assert length(tuple1) == length(tuple2)
    for i in eachindex(tuple1, tuple2)
        func(model, tuple1[i], tuple2[i])
    end
end