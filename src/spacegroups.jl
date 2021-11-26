"""
    SetSpaceGroup{N}

Space group with information about the origin and orientation of the main rotational axis.
"""
struct SetSpaceGroup{N}
    # The actual space group
    spgrp::Crystalline.SpaceGroup{N}
    # Normally this is the point of inversion
    # But this varies for non-centrosymmetric space groups
    origin::SVector{N,Float64}
    # Number of the axis
    index::SVector{N,Int}
end