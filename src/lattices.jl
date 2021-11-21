#=
"""
    AbstractLattice

Supertype for all types that describe crystal lattices.
"""
abstract type AbstractLattice <: Any
end

"""
    Lattice{M,N}

The lattice for an M-dimensional structure with N repeating dimensions.
"""
struct Lattice{M,N}
    # N static vectors of dimension N
    vecs::NTuple{N,NTuple{M,Float64}}
    # Inner constructor
    function Lattice(vecs::NTuple{N,NTuple{M,Float64}})
        @assert M >= N "Number of repeating dimensions ($N) exceeds number of dimensions ($M)"
        # Verify that the vectors are linearly independent
        @assert _is_linearly_independent(vecs) "Vectors are not linearly independent"
        new(vecs)
    end
end

"""
    Lattice{M,N}

The lattice for an M-dimensional structure with N repeating dimensions.
"""
const CrystalLattice{N} = Lattice{N,N}

# Make Lattice{M,N} behave more like a Tuple
Base.getindex(l::Lattice{M,N}, i::Integer) where {M,N} = getindex(l.vecs, i)
=#