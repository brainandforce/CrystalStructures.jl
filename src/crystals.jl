"""
    AbstractCrystal

Supertype for all types that describe crystal structures.
"""
abstract type AbstractCrystal{N} where N <: Any
end

"""
    CrystalStructure{N}

A crystal structure in `N` dimensions. Contains information about the lattice vectors, the space
group (and origin setting), the list of atoms and their positions, and whether the lattice vectors
are primitive and/or conventional.
"""
struct CrystalStructure{N} <: AbstractCrystal{N}
    # Lattice vectors
    basis::Crystalline.DirectBasis{N}
    # Space group - set to 0 if there is no information provided
    spgrp::Int
    # Origin of the space group within the cell
    origin::SVector{N, Float64}
    # Names of atoms
    atoms::Vector{String}
    # Generating set of atoms given the space group and lattice vectors
    sites::Vector{SVector{N, Float64}}
    # Whether the lattice is the conventional lattice
    conv::Bool
    # Whether the lattice is primitive
    # Both can be true if the space group has primitive centering
    # Both can be false if a non-standard setting is used (e.g. hexagonal cell for cubic crystal)
    prim::Bool

    # Inner constructor
    function CrystalStructure{N}(;
        basis::AbstractVector{<:AbstractVector{<:Real}} = _nobasis{N}(),
        spgrp::Integer = 0,
        origin::AbstractVector{<:AbstractVector{<:Real}} = zeros(SVector{3, Float64})
        atoms::AbstractVector{<:AbstractString}
        sites::AbstractVector{<:AbstractVector{<:Real}}
        conv::AbstractBool = true,
        prim::AbstractBool = false,
    )
        # Verify the dimensionality of the data
        basis_vl = [length(v) for v in basis]
        sites_vl = [length(v) for v in sites]
        @assert _allsame(basis_vl) "inconsistent lattice vector lengths"
        @assert _allsame(atomset_vl) "inconsistent site position vector lengths"
        @assert N == first(basis_vl) "expected $N-dimensional basis vectors, " *
            "got $(first(basis_vl))-dimensional vectors"
        @assert N == first(sites_vl) "expected $N-dimensional site position vectors, " *
            "got $(first(sites_vl))-dimensional vectors"
        @assert N == length(basis) "expected $N basis vectors, got $(length(basis))"
        # Make sure that site and origin vector components are between 0 and 1
        origin = origin - floor.(origin)
        sites = [v .- floor.(v) for v in sites]
        # Make sure the space group is valid
        spgrp in NO_SGS[N] || throw(ErrorException("$spgrp is not a valid space group number"))
        return new(basis, spgrp, origin, atoms, sites, conv, prim)
    end
end

"""
    (xtal::CrystalStructure{N};
        basis::Union{AbstractVector{<:AbstractVector{<:Real}}, Nothing} = nothing,
        spgrp::Union{Integer, Nothing} = nothing,
        atoms::Union{AbstractVector{<:AbstractString}, Nothing} = nothing,
        sites::Union{AbstractVector{<:AbstractVector{<:Real}}, Nothing} = nothing,
        conv::Union{AbstractBool, Nothing} = nothing,
        prim::Union{AbstractBool, Nothing} = nothing,
    ) 

Updates a `CrystalStructure` with new parameters. By default, unset parameters are `nothing`; any
set parameters will be used to construct the new `CrystalStructure`.
"""
function CrystalStructure{N}(xtal::CrystalStructure{N};
    basis::Union{AbstractVector{<:AbstractVector{<:Real}}, Nothing} = nothing,
    spgrp::Union{Integer, Nothing} = nothing,
    origin::Union{AbstractVector{<:AbstractVector{<:Real}}, Nothing} = nothing,
    atoms::Union{AbstractVector{<:AbstractString}, Nothing} = nothing,
    sites::Union{AbstractVector{<:AbstractVector{<:Real}}, Nothing} = nothing,
    conv::Union{AbstractBool, Nothing} = nothing,
    prim::Union{AbstractBool, Nothing} = nothing,
)
    # TODO: There's probably a nicer way to do this, isn't there.
    basis === nothing && (basis = xtal.basis)
    spgrp === nothing && (spgrp = xtal.spgrp)
    origin === nothing && (origin = xtal.spgrp)
    atoms === nothing && (atoms = xtal.spgrp)
    sites === nothing && (sites = xtal.sites)
    conv === nothing && (conv = xtal.conv)
    prim === nothing && (spgrp = xtal.prim)
    return CrystalStructure{N}(
        basis=basis,
        spgrp=spgrp,
        origin=origin,
        atoms=atoms,
        sites=sites,
        conv=conv,
        prim=prim,
    )
end

"""
    _enumerate_unique_elements(itr)

Returns a `Vector{Tuple{<:Any, Int}}` that associates every element in `itr` with a numerical index
which counts the number of times that element has occurred counting from the first element.

This is used to count numbers of atoms in a crystal.

# Examples
```jldoctest
julia> _enumerate_unique_elements([2, 2, 3, 2, 3, 3, 1, 1, 2])
9-element Vector{Int}:
 1
 2
 1
 3
 2
 3
 1
 2
 4
```
"""
function _enumerate_unique_elements(v::AbstractVector)
    items = unique(v)
    count = zeros(size(items))
    index = Vector{Int}(undef, length(v))
    for (m, el) in enumerate(v)
        n = findfirst(isequal(el), items)
        index[m] = (count[n] += 1)
    end
    return collect(zip(v, index))
end

function _printsites(io::IO, xtal::CrystalStructure{N})
    output_array = 
    [
        atom
    ]
    println(io, prefix)
end

function Base.display(xtal::CrystalStructure{N}) where N
    println(string(N) * "-dimensional crystal structure")
    if xtal.basis == _nobasis{N}()
        println("Aperiodic (no basis vectors specified)")
    end
end

"""
    CrystalStructureWithData{N,T}

A combination of a `CrystalStructure{N}` with `N`-dimensional volumetric datasets. The datasets are
stored in a `Dict{String, Array{N,T}}`.

When a `CrystalStructureWithData{N,T}` is generated with a single dataset, the dataset may be
accessed using the empty string `""`.
"""
struct CrystalStructureWithData{N,T} <: AbstractCrystal{N}
    xtal::CrystalStructure{N}
    data::Dict{String, Array{N,T}}
end

"""
    
"""
function CrystalStructureWithData{N,T}(xtal::CrystalStructure{N}, data::Array{N,T})
    return CrystalStructureWithData{N,T}(xtal, Dict{String, Array}("" => data))
end