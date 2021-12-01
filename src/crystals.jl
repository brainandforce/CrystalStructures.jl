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
        @assert length(sites) == length(atoms) "number of atoms and sites do not match"
        # Make sure that site and origin vector components are between 0 and 1
        # This may not be needed for CrystalStructure{N}
        # but would be required for CrystalGenerator{N}
        #=
        origin = origin - floor.(origin)
        sites = [v .- floor.(v) for v in sites]
        =#
        # Make sure the space group is valid
        @assert 0 <= spgrp <= NO_SGS[N] "$spgrp is not a valid space group number in $N dimensions"
        return new(basis, spgrp, origin, atoms, sites, conv, prim)
    end
end

# TODO: Equality checking for CrystalStructure{N} is going to be quite complex!
# Two structures may be equal even if their basis vectors and atomic positions are not the same
# because they may have different settings or different choice of basis vectors.
# This is a preliminary attempt but the details are going to take a while.
function Base.:(==)(xtal1::CrystalStructure{N}, xtal2::CrystalStructure{N}) where N
    # Number of atoms should be equal
    size(xtal1.atoms) == size(xtal2.atoms) || return false
    # Types of atoms should be checked next
    if xtal1.basis != xtal2.basis
        # TODO: check that the basis vectors can be converted
        return false
    end
    # If all checks pass, return true
    return true
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

printbasis(io::IO, xtal::CrystalStructure{N}) where N = printbasis(io, xtal.basis)

"""
    printsites(io::IO, xtal::CrystalStructure{N}) where Names

Prints the atomic sites for a crystal structure in the following format:

```
Atomic sites:
Ti      0.333333    0.666667    0.250000
Ti      0.666667    0.333333    0.750000
```
"""
function printsites(io::IO, xtal::CrystalStructure{N}) where N
    tostr(x::Number) = lpad(@sprintf("%f", x), 12)
    println(io, "Atomic sites:")
    for n in 1:length(xtal.atoms)
        println(io, rpad(xtal.atoms[n], 4), join(tostr.(xtal.sites[n])))
    end
end

# Display the data in a CrystalStructure{N}
function Base.display(xtal::CrystalStructure{N}) where N
    println(string(typeof(xtal)) * ":")
    # Print the space group
    if xtal.spgrp == 0
        println("No space group specified")
    else
        println("Space group " * string(xtal.spgrp))
    end
    # Print the basis vectors
    println()
    printbasis(stdout, xtal)
    # Print the atomic sites:
    println()
    printsites(stdout, xtal)
end

"""
    formula(xtal::CrystalStructure{N})

Returns the formula string for a crystal.
"""
function formula(xtal::CrystalStructure{N}) where N
    
end

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

"""
    printdatasets(io, d::AbstractDict{<:Any, DataGrid{N,T}}) where {N,T}

Prints brief info about the datasets contained by a dictionary of datagrids.
"""
function printdatasets(io, d::AbstractDict{<:Any, DataGrid{N,T}}) where {N,T}
    if isempty(d)
        println(io, "No datasets specified")
    else
        str = "Datasets (of type "
        l = lastindex(str) - 1
        println(io, str * string(T) * "):")
        for (k,v) in d
            sz = size(v.data)
            println(io, lpad(k, l) * ":   " * join(string.(sz), '×'))
            # TODO: provide info for dataset vectors and origin
        end
    end
    return nothing
end

function printbasis(io::IO, xtaldata::CrystalStructureWithData{N,T}) where {N,T}
    printbasis(io, xtaldata.xtal.basis)
end

function printsites(io::IO, xtaldata::CrystalStructureWithData{N,T}) where {N,T}
    printsites(io, xtaldata.xtal)
end

function printdatasets(io::IO, xtaldata::CrystalStructureWithData{N,T}) where {N,T}
    printdatasets(io, xtaldata.data)
end

# Display the data in a CrystalStructureWithData{N,T}
function Base.display(xtaldata::CrystalStructureWithData{N,T}) where {N,T}
    println(string(typeof(xtaldata)) * ":")
    # Print the space group
    if xtaldata.xtal.spgrp == 0
        println("No space group specified")
    else
        println("Space group " * string(xtaldata.xtal.spgrp))
    end
    # Print the basis vectors:
    println()
    printbasis(stdout, xtaldata)
    # Print the atomic sites:
    println()
    printsites(stdout, xtaldata)
    # Print the datasets:
    println()
    printdatasets(stdout, xtaldata)
end
