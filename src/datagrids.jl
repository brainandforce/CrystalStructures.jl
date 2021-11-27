"""
    CellMismatch([msg])

The two unit cells defined by the different objects are inequivalent.
"""
struct CellMismatch <: Exception
end

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

function Base.display(g::DataGrid{N,T}) where {N,T}
    # Used to format the number
    tostr(x::Number) = lpad(@sprintf("%f", x), 12)
    # Header
    println(join(string.(size(g.data)), '×') * " " * string(typeof(g)) * ":\n")
    # Basis vectors
    if g.basis == _basis{N}()
        println("Aperiodic (no basis vectors specified)")
    else
        println("Basis vectors:")
        for n in 1:N
            # '`' is 0x60, and all letters come after this in Unicode
            # e.g. 'a' is 0x61, 'b' is 0x62, 'c' is 0x63...
            println(('`' + n) * ':' * join(tostr.(g.basis[n])))
        end
    end
    println("\nOrigin:\n  " * join(tostr.(g.origin)))
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

# Adds two scalar datagrids together
function Base.:+(g1::DataGrid{N,T}, g2::DataGrid{N,T}) where {N,T<:Number}
    _check(g1, g2)
    newdata = g1.data .+ g2.data
    return DataGrid{N,T}(g1.basis, g1.origin, newdata)
end

# Element-wise scalar multiplication of two datagrids
function Base.:*(g1::DataGrid{N,T}, g2::DataGrid{N,T}) where {N,T<:Number}
    _check(g1, g2)
    newdata = g1.data .* g2.data
    return DataGrid{N,T}(g1.basis, g1.origin, newdata)
end

# Scalar multiplication of a whole datagrid
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
    _make_general(a::AbstractArray{T,N})

Converts a periodic grid to a general grid. This is accomplished by repeating data at the initial
boundary at the final boundary.

This is necessary to correctly format data in some file formats, like XCrysDen XSF files.
"""
function _make_general(a::AbstractArray{T,N}) where {T,N}
    # Reference sizes
    sz = size(a)
    # Loop through every coordinate of newsz
    return [a[c .% sz .+ 1...] for c in Iterators.product((0:b for b in sz)...)]
end

"""
    _make_general(a::AbstractArray{T,N})

Converts a general grid to a periodic grid. Some file formats will specify their datagrids in a
general grid, repeating at the boundaries.

This function checks that the grid can actually converted to a general grid. It is faster to make
a general grid from `A_general` with the following code (works for any dimensionality):

```julia
A_periodic = A_general[(1:n-1 for n in size(a))...]
```
"""
function _make_periodic(a::AbstractArray{T,N}) where {T,N}
    # Check that the data matches at the boundaries
    # Get the grid size
    sz = size(a)
    # TODO: implement boundary checking
    # This might require some weird tuple usage...
    return a[(1:n-1 for n in sz)...]
end