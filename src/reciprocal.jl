# reciprocal.jl: data structures needed to stor reciprocal space data
# TODO: I'm sure there's a Julia package that will generate paths in k-space for band structures
# perhaps we should try to integrate that functionality?

"""
    AbstractKPoints

Supertype for data structures that store k-point information associated with a structure.
"""
abstract type AbstractKPoints
end

"""
    KPointGrid{N} <: AbstractKPoints

Stores information about a k-point grid.

The storage of a k-point grid matches the convention used in abinit's `kptrlatt` variable. With 
this convention, the k-point grid can be interpreted as defining a supercell in real space. Each 
vector comprising the grid is given in terms of the real space primitive basis vectors. The shift
of the k-point grid off the reciprocal space origin (Γ) is also given.
"""
struct KPointGrid{N} <: AbstractKPoints
    # k-point grid
    grid::SVector{SVector{N,Int}}
    # shift of the k-point grid off Γ
    shift::SVector{N, Float64}

    function KPointGrid{N}(
        grid::AbstractVector{<:AbstractVector{<:Real}},
        shift::AbstractVector{<:Real}
    )
        # Verify that the grid vectors are linearly independent
        # TODO: this might need to be done independently for exterior constructors
        @assert _is_linearly_independent(grid...) "k-point grid vectors are not linearly independent"
        # Bring the shift within the first Brillouin zone
        shift = shift - floor.(shift)
        return new(grid, shift)
    end
end

"""
    KPointList{N}

Stores information about individual k-points used in a calculation. This may be useful for
calculations of band structures.
"""
struct KPointList{N} <: AbstractKPoints
    list::Vector{SVector{N,Float64}}
end


"""
    nkpt(k::KPointList{N}) -> Int

Returns the number of k-points in a `KPointList`.
"""
nkpt(k::KPointList{N}) = size(k.list)


function nkpt(k::KPointGrid{N})

end