## boundary condition types and structures

# boundary condition categories
"""
    BCType

Abstract type for boundary condition category
"""
abstract type BCType end

"""
    DirichletBCType, NeumannBCType, PeriodicBCType

Categories of boundary condition
"""
struct DirichletBCType <: BCType end
struct NeumannBCType <: BCType end
struct PeriodicBCType <: BCType end


# single side boundary condition
""" 
    AbstractBoundaryConditionSide{BCType}

Abstract type for boundary condition of category BCType on a single side of a domain
"""
abstract type AbstractBoundaryConditionSide{BCType} end

"""
    PeriodicBoundaryConditionSide{PeriodicBCType}

Periodic boundary condition on a single side
"""
struct PeriodicBoundaryConditionSide <: AbstractBoundaryConditionSide{PeriodicBCType} end

"""
    FunctionalBoundaryConditionSide{BCType}

Boundary condition on a single side defined by a function
"""
struct FunctionalBoundaryConditionSide{BCType} <: AbstractBoundaryConditionSide{BCType}
    bc_func::Function
end

"""
    DiscreteBoundaryConditionSide{BCType}

Discrete boundary condition on a single side defined by a vector of values
"""
struct DiscreteBoundaryConditionSide{BCType} <: AbstractBoundaryConditionSide{BCType}
    bc_vec::Vector
end

"""
    NeumannBC(func::Function)
    NeumannBC(value::Number)
    NeumannBC(vec::Vector)

Construct Neumann boundary condition on a single side
"""
NeumannBC(func::Function)   = FunctionalBoundaryConditionSide{NeumannBCType}(func)
NeumannBC(value::Number)    = FunctionalBoundaryConditionSide{NeumannBCType}(x->value)
NeumannBC(vec::Vector)      = DiscreteBoundaryConditionSide{NeumannBCType}(vec)

"""
    DirichletBC(func::Function)
    DirichletBC(value::Number)
    DirichletBC(vec::Vector)

Construct Dirichlet boundary condition on a single side
"""
DirichletBC(func::Function) = FunctionalBoundaryConditionSide{DirichletBCType}(func)
DirichletBC(value::Number)  = FunctionalBoundaryConditionSide{DirichletBCType}(x->value)
DirichletBC(vec::Vector)    = DiscreteBoundaryConditionSide{DirichletBCType}(vec)

"""
    PeriodicBC()

Construct Periodic boundary condition on a single side
"""
PeriodicBC() = PeriodicBoundaryConditionSide()

"""
    bc_type(bc::AbstractBoundaryConditionSide{BCType}) where BCType

Returns the boundary condition category BCType of the input boundary condition side
"""
bc_type(::AbstractBoundaryConditionSide{BCType}) where BCType = BCType

# domain boundary condition
"""
    AbstractBoundaryConditionDomain

Abstract type for boundary condition on a domain
"""
abstract type AbstractBoundaryConditionDomain end

"""
    AbstractBoundaryCondition2D

Abstract types for boundary condition on 2D domain
"""
abstract type AbstractBoundaryCondition2D <: AbstractBoundaryConditionDomain end

"""
    AbstractBoundaryCondition3D

Abstract type for boundary condition on 3D domain
"""
abstract type AbstractBoundaryCondition3D <: AbstractBoundaryConditionDomain end

"""
    BoundaryCondition2D(left, right, bottom, top)

Boundary condition on 2D rectangular domain
"""
struct BoundaryCondition2D <: AbstractBoundaryCondition2D
    left::AbstractBoundaryConditionSide     # x-direction
    right::AbstractBoundaryConditionSide
    bottom::AbstractBoundaryConditionSide   # z-direction
    top::AbstractBoundaryConditionSide
end

"""
    DiscreteBoundaryCondition2D(left, right, bottom, top)

Discrete boundary condition on 2D rectangular domain
"""
struct DiscreteBoundaryCondition2D <: AbstractBoundaryCondition2D
    bottom::Union{PeriodicBoundaryConditionSide, DiscreteBoundaryConditionSide}
    top::Union{PeriodicBoundaryConditionSide, DiscreteBoundaryConditionSide}
    left::Union{PeriodicBoundaryConditionSide, DiscreteBoundaryConditionSide}
    right::Union{PeriodicBoundaryConditionSide, DiscreteBoundaryConditionSide}
end

"""
    BoundaryCondition3D(left, right, front, back, bottom, top)

Boundary condition on 3D rectangular domain
"""
struct BoundaryCondition3D <: AbstractBoundaryCondition3D
    left::AbstractBoundaryConditionSide     # x-direction
    right::AbstractBoundaryConditionSide
    front::AbstractBoundaryConditionSide    # y-direction
    back::AbstractBoundaryConditionSide
    bottom::AbstractBoundaryConditionSide   # z-direction
    top::AbstractBoundaryConditionSide
end

"""
    BoundaryCondition2D(domain, var_bcs::NamedTuple)
    BoundaryCondition2D(var_bcs::NamedTuple)
    BoundaryCondition2D(; left, right, bottom, top)
    BoundaryCondition2D(domain; left, right, bottom, top)

Constructs boundary condition on 2D rectangular domain from a NamedTuple of boundary conditions 
for each side (left, right, bottom, top) or from keyword arguments.

Periodic boundary conditions are automatically assigned for sides with periodic topology in the domain.
Error is raised if any required boundary condition is missing, or if incompatible BCs are provided.
"""
BoundaryCondition2D(domain, var_bcs::NamedTuple) = BoundaryCondition2D(domain; var_bcs...)
BoundaryCondition2D(; left, right, bottom, top) = BoundaryCondition2D(left, right, bottom, top)

function BoundaryCondition2D(var_bcs::NamedTuple) 
    for key in (:bottom, :top, :left, :right)
        if !haskey(var_bcs, key)
            error("Missing boundary condition for $key side. Call function with domain argument to set default BCs.")
        end
    end
    return BoundaryCondition2D(; var_bcs...)
end

function BoundaryCondition2D(domain::Box2D{Mx, Mz, Hx, Hz, Bx, Bz}; left=nothing, right=nothing, bottom=nothing, top=nothing, kwargs...) where {Mx, Mz, Hx, Hz, Bx, Bz}
    return BoundaryCondition2D(; 
            left   = Bx <: Periodic ? PeriodicBC() : ( (isnothing(left)   || bc_type(left)   <:PeriodicBCType ) ? error("Missing/incompatible boundary condition for left side.")   : left   ),
            right  = Bx <: Periodic ? PeriodicBC() : ( (isnothing(right)  || bc_type(right)  <:PeriodicBCType ) ? error("Missing/incompatible boundary condition for right side.")  : right  ),
            bottom = Bz <: Periodic ? PeriodicBC() : ( (isnothing(bottom) || bc_type(bottom) <:PeriodicBCType ) ? error("Missing/incompatible boundary condition for bottom side.") : bottom ),
            top    = Bz <: Periodic ? PeriodicBC() : ( (isnothing(top)    || bc_type(top)    <:PeriodicBCType ) ? error("Missing/incompatible boundary condition for top side.")    : top    )
    )
end


"""
    BoundaryCondition3D(domain, var_bcs::NamedTuple)
    BoundaryCondition3D(var_bcs::NamedTuple)
    BoundaryCondition3D(; left, right, bottom, top)
    BoundaryCondition3D(domain; left, right, bottom, top)

Constructs boundary condition on 3D rectangular domain from a NamedTuple of boundary conditions 
for each side (left, right, bottom, top) or from keyword arguments.

Periodic boundary conditions are automatically assigned for sides with periodic topology in the domain.
Error is raised if any required boundary condition is missing.
"""
BoundaryCondition3D(domain, var_bcs::NamedTuple) = BoundaryCondition3D(domain; var_bcs...)
BoundaryCondition3D(var_bcs::NamedTuple) = BoundaryCondition3D(; var_bcs...)
BoundaryCondition3D(; left, right, front, back, bottom, top) = BoundaryCondition3D(left, right, front, back, bottom, top)

function BoundaryCondition3D(domain::Box3D{Mx, My, Mz, Hx, Hy, Hz, Bx, By, Bz}; left=nothing, right=nothing, front=nothing, back=nothing, bottom=nothing, top=nothing, kwargs...) where {Mx, My, Mz, Hx, Hy, Hz, Bx, By, Bz}
    return BoundaryCondition3D(; 
            left   = Bx <: Periodic ? PeriodicBC() : ( (isnothing(left)   || bc_type(left)   <:PeriodicBCType ) ? error("Missing/incompatible boundary condition for left side.")   : left   ),
            right  = Bx <: Periodic ? PeriodicBC() : ( (isnothing(right)  || bc_type(right)  <:PeriodicBCType ) ? error("Missing/incompatible boundary condition for right side.")  : right  ),
            front  = By <: Periodic ? PeriodicBC() : ( (isnothing(front)  || bc_type(front)  <:PeriodicBCType ) ? error("Missing/incompatible boundary condition for front side.")  : front  ),
            back   = By <: Periodic ? PeriodicBC() : ( (isnothing(back)   || bc_type(back)   <:PeriodicBCType ) ? error("Missing/incompatible boundary condition for back side.")   : back   ),
            bottom = Bz <: Periodic ? PeriodicBC() : ( (isnothing(bottom) || bc_type(bottom) <:PeriodicBCType ) ? error("Missing/incompatible boundary condition for bottom side.") : bottom ),
            top    = Bz <: Periodic ? PeriodicBC() : ( (isnothing(top)    || bc_type(top)    <:PeriodicBCType ) ? error("Missing/incompatible boundary condition for top side.")    : top    )
    )
end

"""
    display_bc(bc::BoundaryCondition2D)

Displays the boundary condition for each side of a 2D rectangular domain
"""
function display_bc(bc::BoundaryCondition2D)
    println("Bottom BC: ", bc_type(bc.bottom))
    println("Top BC:    ", bc_type(bc.top))
    println("Left BC:   ", bc_type(bc.left))
    println("Right BC:  ", bc_type(bc.right))
end

# boundary condition list for all variables
"""
    AbstractBoundaryConditions

Abstract type for collection of boundary conditions for all variables on a domain
"""
abstract type AbstractBoundaryConditions end

"""
    BoundaryConditions2D(keys, conditions)

Collection of boundary conditions for all variables on a 2D rectangular domain

e.g. `conditions = (; u = u_bc, w = w_bc)`,
where `u_bc` and `w_bc` are of type `BoundaryCondition2D`
"""
struct BoundaryConditions2D <: AbstractBoundaryConditions
    keys::Tuple
    conditions::NamedTuple
end

"""
    BoundaryConditions3D(keys, conditions)

Collection of boundary conditions for all variables on a 3D rectangular domain

e.g. `conditions = (; u = u_bc, w = w_bc)`,
where `u_bc`` and `w_bc` are of type `BoundaryCondition3D`
"""
struct BoundaryConditions3D <: AbstractBoundaryConditions
    keys::Tuple
    conditions::NamedTuple
end

"""
    BoundaryConditions2D(domain, bc_list::NamedTuple; display::Bool = false)
    BoundaryConditions2D(bc_list::NamedTuple; display::Bool = false)

Constructs collection of boundary conditions for all variables on a 2D rectangular domain from
a NamedTuple of boundary conditions on each side for each variable.

`domain` is required if incomplete boundary conditions are provided, to set default BCs based on domain topology.

If `display` is true, the boundary conditions are printed to the console.
"""
function BoundaryConditions2D(domain, bc_list::NamedTuple; display::Bool = false)
    conditions = (; )
    for (key, val) in pairs(bc_list)
        if isa(val, BoundaryCondition2D)
            val_bc = val
        elseif ~isa(val, BoundaryCondition2D) && isa(val, NamedTuple)
            val_bc = BoundaryCondition2D(domain, val)
        else
            error("Boundary condition for variable $key must be of type BoundaryCondition2D or NamedTuple.")
        end
        conditions = (; conditions..., (key => val_bc))
    end
    output = BoundaryConditions2D(keys(conditions), conditions)
    check_bcs(domain, output)
    if display
        display_bcs(output)
    end
    return output
end

function BoundaryConditions2D(bc_list::NamedTuple; display::Bool = false)
    conditions = (; )
    for (key, val) in pairs(bc_list)
        if isa(val, BoundaryCondition2D)
            val_bc = val
        elseif ~isa(val, BoundaryCondition2D) && isa(val, NamedTuple)
            val_bc = BoundaryCondition2D(val)
        else
            error("Boundary condition for variable $key must be of type BoundaryCondition2D or NamedTuple.")
        end
        conditions = (; conditions..., (key => val_bc))
    end
    output = BoundaryConditions2D(keys(conditions), conditions)
    check_bcs(output)
    if display
        display_bcs(output)
    end
    return output
end

"""
    display_bcs(bcs::BoundaryConditions2D)
    
Displays the boundary conditions for all variables on a 2D rectangular domain
"""
function display_bcs(bcs::BoundaryConditions2D)
    for key in bcs.keys
        println("$key:")
        display_bc(bcs.conditions[key])
    end
end

"""
    check_bcs(domain, BC::AbstractBoundaryConditions)
    check_bcs(BC::AbstractBoundaryConditions)

Checks that the boundary conditions are compatible with the domain topology.
Raises an error if incompatible boundary conditions are found.
"""
check_bcs(::Box2D{M, N, Hx, Hz, Bx, Bz}, BC::AbstractBoundaryConditions) where {M, N, Hx, Hz, Bx, Bz} = 
    check_bcs(BC; Bx, Bz)
check_bcs(::Box3D{M, N, P, Hx, Hy, Hz, Bx, By, Bz}, BC::AbstractBoundaryConditions) where {M, N, P, Hx, Hy, Hz, Bx, By, Bz} = 
    check_bcs(BC; Bx, By, Bz)
function check_bcs(BC::AbstractBoundaryConditions; Bx = Nothing, By = Nothing, Bz = Nothing)
    (; conditions) = BC
    # loop over all variables
    for (var, bc) in pairs(conditions)
        # check that periodic BCs are only used where domain is periodic
        for side in intersect((:left, :right, :top, :bottom, :front, :back), fieldnames(typeof(bc)))
            if side in (:left, :right)
                bc_side = getproperty(bc, side)
                bc_type_side = bc_type(bc_side)
                if ((bc_type_side === PeriodicBCType) && Bx <: Bounded) || ((bc_type_side !== PeriodicBCType) && Bx <: Periodic)
                    error("The $side boundary condition for variable $var is of type $bc_type_side but domain has $Bx topology in x-direction.")
                end
            elseif side in (:top, :bottom)
                bc_side = getproperty(bc, side)
                bc_type_side = bc_type(bc_side)
                if ((bc_type_side === PeriodicBCType) && Bz <: Bounded) || ((bc_type_side !== PeriodicBCType) && Bz <: Periodic)
                    error("The $side boundary condition for variable $var is of type $bc_type_side but domain has $Bz topology in z-direction.")
                end
            elseif side in (:front, :back)
                bc_side = getproperty(bc, side)
                bc_type_side = bc_type(bc_side)
                if ((bc_type_side === PeriodicBCType) && By <: Bounded) || ((bc_type_side !== PeriodicBCType) && By <: Periodic)
                    error("The $side boundary condition for variable $var is of type $bc_type_side but domain has $By topology in y-direction.")
                end
            end
        end
        if var == :w
            # check that w has Dirichlet BCs on bounded vertical sides
            for side in (:top, :bottom)
                bc_side = getproperty(bc, side)
                bc_type_side = bc_type(bc_side)
                if bc_type_side !== DirichletBCType && Bz <: Bounded
                    error("The $side boundary condition for variable $var must be DirichletBCType on bounded vertical sides.")
                end
            end
        end
        if var == :u
            # check that u has Dirichlet BCs on bounded horizontal sides
            for side in (:left, :right)
                bc_side = getproperty(bc, side)
                bc_type_side = bc_type(bc_side)
                if bc_type_side !== DirichletBCType && Bx <: Bounded
                    error("The $side boundary condition for variable $var must be DirichletBCType on bounded horizontal sides.")
                end
            end
        end        
    end
end
# Unimplemented code:
# # discrete boundary condition
# struct DiscreteBoundaryConditionSide{BCType} <: AbstractBoundaryConditionSide{BCType}
#     bc_vec::Vector
# end

# # constant boundary condition
# struct ConstantBoundaryConditionSide{BCType} <: AbstractBoundaryConditionSide{BCType}
#     bc_value::Number
# end

# whole domain boundary condition with discrete or periodic BCs only


# discretize_bcs(bc::BoundaryCondition2D, domain::Box2D) = 
#     DiscreteBoundaryCondition2D(
#         discretize_bc_side(bc.bottom, xpoints),
#         discretize_bc_side(bc.top, xpoints),
#         discretize_bc_side(bc.left, zpoints),
#         discretize_bc_side(bc.right, zpoints),
#     )

# # evaluate functional BCs at given points to create discrete BCs
# discretize_bc_side(bc_side::FunctionalBoundaryConditionSide{BCType}, points::Vector) where BCType = 
#     DiscreteBoundaryConditionSide{BCType}(collect(bc_side.bc_func.(points)))
# # evalute constant BCs at given points to create discrete BCs
# # discretize_bc_side(bc_side::ConstantBoundaryConditionSide{BCType}, points::Vector) where BCType = 
#     # DiscreteBoundaryConditionSide{BCType}(fill(bc_side.bc_value, length(points)))
# # pass through discrete BCs
# discretize_bc_side(bc_side::DiscreteBoundaryConditionSide{BCType}, points::Vector) where BCType = bc_side
# # pass through periodic BCs
# discretize_bc_side(bc_side::PeriodicBoundaryConditionSide{BCType}, points::Vector) where BCType = bc_side
