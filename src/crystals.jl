"""
    AbstractCrystal{N}

Supertype for all types that describe crystal structures.

All subtypes of `AbstractCrystalData{N,T}` consist of a basis, space group with origin, atomic
sites, and information about the choice of cell vectors.
"""
abstract type AbstractCrystal{N}
end

"""
    CrystalStructure{N}

A crystal structure in `N` dimensions. Contains information about the lattice vectors, the space
group (and origin setting), the list of atoms and their positions, and whether the lattice vectors
are primitive and/or conventional.
"""
struct CrystalStructure{N} <: AbstractCrystal{N}
    # Lattice vectors
    basis::SVector{N, SVector{N, Float64}}
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
    prim::Bool
    # Both can be true if the space group has primitive centering
    # Both can be false if a non-standard setting is used (e.g. hexagonal cell for cubic crystal)

    # Inner constructor
    function CrystalStructure{N}(;
        basis::AbstractVector{<:AbstractVector{<:Real}} = _basis{N}(),
        spgrp::Integer = 0,
        origin::AbstractVector{<:Real} = zeros(SVector{3, Float64}),
        atoms::AbstractVector{<:AbstractString} = Vector{String}[],
        sites::AbstractVector{<:AbstractVector{<:Real}} = Vector{Float64}[],
        conv::Bool = true,
        prim::Bool = false,
    ) where N
        # Verify the dimensionality of the data
        basis_vl = [length(v) for v in basis]
        sites_vl = [length(v) for v in sites]
        @assert _allsame(basis_vl) "inconsistent lattice vector lengths"
        @assert _allsame(sites_vl) "inconsistent site position vector lengths"
        @assert N == first(basis_vl) "expected $N-dimensional basis vectors, " *
            "got $(first(basis_vl))-dimensional vectors"
        @assert N == first(sites_vl) "expected $N-dimensional site position vectors, " *
            "got $(first(sites_vl))-dimensional vectors"
        @assert N == length(basis) "expected $N basis vectors, got $(length(basis))"
        # Make sure that site and origin vector components are between 0 and 1
        origin = origin - floor.(origin)
        sites = [v .- floor.(v) for v in sites]
        # Make sure the space group is valid
        spgrp <= NO_SGS[N] || throw(ErrorException("$spgrp is not a valid space group number"))
        return new(basis, spgrp, origin, atoms, sites, conv, prim)
    end
end

"""
    (xtal::CrystalStructure{N};
        basis::Union{AbstractVector{<:AbstractVector{<:Real}}, Nothing} = nothing,
        spgrp::Union{Integer, Nothing} = nothing,
        atoms::Union{AbstractVector{<:AbstractString}, Nothing} = nothing,
        sites::Union{AbstractVector{<:AbstractVector{<:Real}}, Nothing} = nothing,
        conv::Union{Bool, Nothing} = nothing,
        prim::Union{Bool, Nothing} = nothing,
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
    conv::Union{Bool, Nothing} = nothing,
    prim::Union{Bool, Nothing} = nothing,
) where N
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

function Base.display(xtal::CrystalStructure{N}) where N
    # New method for printing floats
    # Avoiding the @printf/@sprintf macros because dynamic formatting is difficult
    tostr(x::Number) = lpad(@sprintf("%f", x), 12)
    println(string(N) * "-dimensional crystal structure")
    # Print the space group
    if xtal.spgrp == 0
        println("No space group specified")
    else
        println("Space group " * string(xtal.spgrp))
    end
    println()
    # Print the basis vectors
    if xtal.basis == _basis{N}()
        println("Aperiodic (no basis vectors specified)")
    else
        println("Basis vectors:")
        for n in 1:N
            # '`' is 0x60, and all letters come after this in Unicode
            # e.g. 'a' is 0x61, 'b' is 0x62, 'c' is 0x63...
            println(('`' + n) * ':' * join(tostr.(xtal.basis[n])))
        end
    end
    # Print the atomic sites:
    println("\nAtomic sites:")
    for n in 1:length(xtal.atoms)
        println(rpad(xtal.atoms[n], 4), join(tostr.(xtal.sites[n])))
    end
end

"""
    formula(xtal::CrystalStructure{N})

Returns the formula string for a crystal.
"""
function formula(xtal::CrystalStructure{N}) where N
    
end

#=
function _printsites(io::IO, xtal::CrystalStructure{N})
    output_array = 
    [
        atom
    ]
    println(io, prefix)
end
=#

"""
    AbstractCrystalData{N,T}

Supertype for all types that describe crystal structures with datasets of type `T`.

All subtypes of `AbstractCrystalData{N,T}` must consist of a field containing crystal data and a
field which contains datasets - usually, stored in a dictionary.
"""
abstract type AbstractCrystalData{N,T} <: AbstractCrystal{N}
end

"""
    CrystalStructureWithData{N,T}

A combination of a `CrystalStructure{N}` with `N`-dimensional volumetric datasets. The datasets are
stored in a `Dict{String, DataGrid{N,T}}`.

When a `CrystalStructureWithData{N,T}` is generated with a single dataset, the dataset may be
accessed using the empty string `""`.
"""
struct CrystalStructureWithData{N,T} <: AbstractCrystalData{N,T}
    xtal::CrystalStructure{N}
    data::Dict{String, DataGrid{N,T}}
end

function Base.display(xtaldata::CrystalStructureWithData{N,T}) where {N,T}
    display(xtaldata.xtal)
    println()
    # List the datasets
    if isempty(xtaldata.data)
        println("No datasets specified")
    else
        str = "Datasets (of type "
        l = lastindex(str) - 1
        println(str * string(T) * "):")
        for (k,v) in xtaldata.data
            sz = size(v.data)
            println(lpad(k, l) * ":   " * join(string.(sz), '×'))
            # TODO: provide info for dataset vectors and origin
        end
    end
end
