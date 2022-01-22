"""
    AbstractLattice

Contains information about a crystal lattice.
"""
abstract type AbstractLattice
end

"""
    RealLattice{N}

Describes an N-dimensional real space crystal lattice, with an associated primitive and 
conventional lattice. These lattices may be the same depending on the lattice type.
"""
struct RealLattice{N} <: AbstractLattice
    # Primitive and conventional lattice vectors are stored together
    prim::SVector{N,SVector{N,Float64}}
    conv::SVector{N,SVector{N,Float64}}
end

# TODO: Implement lattice checking and conversion
# Some of this functionality might already be present in Crystalline.jl
# However, we need to be able to recognize primitive lattices with their own unique conditions

"""
    nuniqlen(basis::AbstractVector{<:AbstractVector{<:Real}}) -> Int

Returns the number of vectors that are the same length.
"""
function nuniqlen(basis::AbstractVector{<:AbstractVector{<:Real}})
    # unique!() does better than unique() here (per @btime)
    return length(unique!(norm.(basis)))
end

# Same thing for matrices
nuniqlen(basis::AbstractMatrix{<:Real}) = nuniqlen(_tovectors(basis))

"""
    northog(basis::AbstractVector{<:AbstractVector{<:Real}}) -> Int

Returns the number of pairs of orthogonal vectors.
"""
function northog(basis::AbstractVector{<:AbstractVector{<:Real}})
    # FIXME: This doesn't work for arbitrary dimensions
    # It needs to use the upper portion of the dot product matrix
    return count(iszero, dot.(basis, basis[circshift(1:length(basis), 1)]))
end

northog(basis::AbstractMatrix{<:Real}) = northog(_tovectors(basis))

# TODO: Write a number of tests for checking if lattices meet certain criteria

#=
Criteria for primitive unit cells in 3 dimensions:
    aP: 6 parameters (a, b, c, α, β, γ)
    mP: 4 parameters (a, b, c, β), α = γ = 90°
    mS: 4 parameters (a, c, α, γ), a = b, α = β
    oP: 3 parameters (a, b, c), α = β = γ = 90°
    oS: 3 parameters (a, c, γ), a = b, α = β = 90°
    oI: 3 parameters (a, α, γ), a = b = c, |cos(α)| + |cos(β)| + |cos(γ)| = 1
    oF: 3 parameters (a, b, c), sin(α)/a = sin(β)/b = sin(γ)/c, α + β + γ = 180°
        3 parameters (a, c, γ), [forgot the rest of this alternate test]
    tP: 2 parameters (a, c), a = b, α = β = γ = 90°
    tI: 2 parameters (a, α), a = b = c, |cos(α)| + |cos(β)| + |cos(γ)| = 1
    hP: 2 parameters (a, c), a = b, α = β = 90°, γ = 120°
    hR: 2 parameters (a, α), a = b = c, α = β = γ
    cP: 1 parameter (a), a = b = c, α = β = γ = 90°
    cI: 1 parameter (a), a = b = c, α = β = γ = arccos(-1/3)°
    cF: 1 parameter (a), a = b = c, α = β = γ = 60°

Rules that need to be checked:
    Number of unique lengths
    Number of orthogonal vectors
    Number of unique angles
    Cosine rule: |cos(α)| + |cos(β)| + |cos(γ)| = 1
    Sine rule: sin(α)/a = sin(β)/b = sin(γ)/c, α + β + γ = 180°
    
    The last three are probably going to be trickier to implement, or at the very least there
    might be some unintuitive linear algebra tricks to make those checks happen faster.
=#



# TODO: what are the crtieria for a "good" lattice basis?
# Potentially worth looking at the Niggli reduction algorithm to do this

# TODO: Tools for working with reciprocal space lattices