# boundary topology types

""" 
    BoundaryTopology

Abstract type for boundary topology of a domain boundary.
"""
abstract type BoundaryTopology end

""" 
    Periodic <: BoundaryTopology
    Bounded  <: BoundaryTopology
Types representing periodic and bounded boundary topologies.
"""
struct Periodic <: BoundaryTopology end
struct Bounded  <: BoundaryTopology end
