"""
    _nobasis{N}()

Returns a set of `N` `N`-dimensional zero vectors, which can serve to indicate that the basis
vectors of a dataset are missing, or that the structure is aperiodic.
"""
# TODO: should probably use zeros(SVector{3, Float64}) somewhere here
_nobasis{N}() = SVector{N,SVector{N,Float64}}(SVector{N,Float64}(0 for n=1:N) for n = 1:N)

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