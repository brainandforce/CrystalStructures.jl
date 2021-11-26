"""
    AbstractCrystal{N}

Supertype for all types that describe crystal structures.
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
    println(string(N) * "-dimensional crystal structure\n")
    # Print the basis vectors
    if xtal.basis == _basis{N}()
        println("Aperiodic (no basis vectors specified)")
    else
        println("Basis vectors:")
        # Print the space group
        if xtal.spgrp == 0
            print("No space group")
        end
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
    DataGrid{N,T}

Used to store `N`-dimensional volumetric data of type `T` for a crystal, giving both the basis 
vectors defining the datagrid and the data itself. This has the advantage of being able to store
data that is defined relative to the primitive cell vectors (common in the case of computational 
data) or the conventional cell vectors.
"""
struct DataGrid{N,T}
    # TODO: Does this need to handle aperiodic structures?
    # Data basis vectors - may or may not match crystal
    basis::SVector{N,SVector{N,Float64}}
    # Origin, defined in reduced coordinates
    origin::SVector{N,Float64}
    data::Array{T,N}
    function DataGrid{N,T}(
        basis::AbstractVector{<:AbstractVector{<:Number}},
        origin::AbstractVector{<:Number},
        data::Array{T,N}
    ) where {N<:Any, T<:Any}
        new{N,T}(basis, origin, data)
    end
end

function DataGrid{N,T}(
    basis::AbstractVector{<:AbstractVector{<:Number}},
    data::Array{T,N}
) where {N<:Any, T<:Any}
    return DataGrid{N,T}(basis, [0, 0, 0], data)
end

"""
    CellMismatch([msg])

The two unit cells defined by the different objects are inequivalent.
"""
struct CellMismatch <: Exception
end

"""
    _check(g1::DataGrid{N,T}, g2::DataGrid{N,T})

Checks that the datagrids are compatible with binary mathematical operations, like addition and
multiplication. Throws a `CellMismatch` if the unit cells are incompatible for addition, and throws
a `DimensionMismatch` if the datagrids are of different dimensionality.

When interpolation functions are included to operate on datagrids (for resizing) these checks may
be altered or removed as needed.
"""
function _check(g1::DataGrid{N,T1}, g2::DataGrid{N,T2}) where {N,T1,T2}
    # TODO: can we do some sort of automatic conversion of unit cells?
    g1.basis == g2.basis || throw(CellMismatch("Unit cells are not the same."))
    # TODO: try adjusting the origin automatically
    g1.origin == g2.origin || throw(CellMismatch("Unit cell origins are not the same."))
    # TODO: try using interpolation to add grids of inequivalent size
    size(g1.data) == size(g2.data) || throw(DimensionMismatch("Data grids differ in dimension."))
    return nothing
end

# Adds two datagrids together
function Base.:+(g1::DataGrid{N,T}, g2::DataGrid{N,T}) where {N,T<:Number}
    _check(g1, g2)
    newdata = g1.data .+ g2.data
    return DataGrid{N,T}(g1.basis, g1.origin, newdata)
end

# Element-wise multiplication of two datagrids
function Base.:*(g1::DataGrid{N,T}, g2::DataGrid{N,T}) where {N,T<:Number}
    _check(g1, g2)
    newdata = g1.data .* g2.data
    return DataGrid{N,T}(g1.basis, g1.origin, newdata)
end

# Scalar multiplication of a datagrid
function Base.:*(s::Number, g::DataGrid{N,T}) where {N,T<:Number}
    return DataGrid{N,T}(g.basis, g.origin, s * g.data)
end

# Same but with reversed order of arguments
function Base.:*(g::DataGrid{N,T}, s::Number) where {N,T<:Number}
    return Base.*(s, g)
end

# Sum up all of the elements of the datagrid
function sum(g::DataGrid{N,T}) where {N,T<:Number}
    return sum(g.data)
end

"""
    CrystalStructureWithData{N,T}

A combination of a `CrystalStructure{N}` with `N`-dimensional volumetric datasets. The datasets are
stored in a `Dict{String, DataGrid{N,T}}`.

When a `CrystalStructureWithData{N,T}` is generated with a single dataset, the dataset may be
accessed using the empty string `""`.
"""
struct CrystalStructureWithData{N,T} <: AbstractCrystal{N}
    xtal::CrystalStructure{N}
    data::Dict{String, DataGrid{N,T}}
end