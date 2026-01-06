# type definitions for advection scheme choices

"""
    AdvectionScheme
Abstract type for advection schemes.
"""
abstract type AdvectionScheme end

"""
    AdvEnergyCons
Energy-conserving advection scheme, see Sadourny (1975)
"""
struct AdvEnergyCons <: AdvectionScheme end

"""
    AdvEnstrophyCons
Enstrophy-conserving advection scheme, see Sadourny (1975)
"""
struct AdvEnstrophyCons <: AdvectionScheme end
