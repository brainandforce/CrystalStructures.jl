const _basis{N} = SVector{N,SVector{N,Float64}}

"""
    _basis{N}()

Returns a set of `N` `N`-dimensional zero vectors, which can serve to indicate that the basis
vectors of a dataset arecons missing, or that the structure is aperiodic.
"""
function _basis{N}() where N
    return zeros(SVector{N,SVector{N,Float64}})
end

"""
    _allsame(itr)

Returns `true` if all the elements of an iterator are identical.
"""
_allsame(itr) = all(x -> x == itr[1], itr)

"""
    _is_linearly_independent(vecs::AbstractMatrix{<:Number}) -> Bool

Returns `true` if a matrix consists of linearly independent vectors, `false` otherwise.

This test is performed by checking that the rank of the matrix is equal to the smallest dimension
of the matrix.
"""
function _is_linearly_independent(vecs::AbstractMatrix{<:Number})
    return LinearAlgebra.rank(vecs) == minimum(size(vecs))
end

"""
    _is_linearly_independent(vecs::NTuple{N,NTuple{M,<:Number}}) -> Bool
    _is_linearly_independent(vecs::NTuple{M,<:Number}...) -> Bool
    _is_linearly_independent(vecs::AbstractVector{<:AbstractVector{<:Number}}) -> Bool

Returns `true` if a set of vectors is linearly independent, false otherwise. 
"""
function _is_linearly_independent(vecs::NTuple{N,NTuple{M,<:Number}}) where {M,N}
    # if M is greater than N, they cannot be linearly independent
    M < N && return false
    mat = [vecs[n][m] for m = 1:M, n = 1:N]
    return _is_linearly_independent(mat)
end

_is_linearly_independent(vecs::NTuple{M,<:Number}...) where M = _is_linearly_independent(vecs)

function _is_linearly_independent(vecs::AbstractVector{<:AbstractVector{<:Number}})
    # Get number of vectors
    N = length(vecs)
    # Get length of each vector in the vector
    M = let l = [length(v) for v in vecs]
        _allsame(l) ? l[1] : throw(ErrorException("Vectors have inconsistent dimensions"))
    end
    # if M is greater than N, they cannot be linearly independent
    M < N && return false
    # Convert to a matrix for rank checking
    mat = [vecs[m][n] for n = 1:N, m = 1:M]
    return _is_linearly_independent(mat)
end

"""
    _reducecoords(B::AbstractMatrix{<:Number}, v::AbstractVector{<:Number})

Converts a vector `v` from Cartesian coordinates to reduced coordinates relative to the basis
vectors `B` of a unit cell, expressed as a collection of vectors.
"""
function _reducecoords(b::AbstractVector{AbstractVector{<:Number}}, v::AbstractVector{<:Number})
    B = hcat(b...)
    return inv(B) * v
end

"""
    printbasis(io::IO, basis::AbstractVector{<:AbstractVector{<:Real}})

Prints the basis vectors in the following format:

```
Basis vectors:
a:    2.847216    0.000000    0.000000
b:   -1.423608    2.465762   -0.000000
c:    0.000000    0.000000    4.538054
```
"""
function printbasis(io::IO, basis::AbstractVector{<:AbstractVector{<:Real}})
    tostr(x::Number) = lpad(@sprintf("%f", x), 12)
    # Check if basis vectors are empty or zero
    if isempty(basis) || any(all.(iszero, v for v in basis))
        println(io, "Aperiodic (no basis vectors specified)")
    else
        println(io, "Basis vectors:")
        for (n,v) in enumerate(basis)
            # '`' is 0x60, and all letters come after this in Unicode
            # e.g. 'a' is 0x61, 'b' is 0x62, 'c' is 0x63...
            println(io, ('`' + n) * ':' * join(tostr.(v)))
        end
    end
    return nothing
end