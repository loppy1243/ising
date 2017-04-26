# The module name is unfortunately not SpinGrid since currently in Julia you cannot have something
# in a module with the same name as the module because the module has itself in scope.
"""
All spin related stuff.
"""
module SpinMod
import Loppy.Indexing.modindex, Loppy.Indexing.modgetindex

export SPIN_UP, SPIN_DOWN
const SPIN_UP = true
const SPIN_DOWN = false

export Spin
"""
Represents spin-up or spin-down. Canonically, these are true and false, respectively.
"""
typealias Spin Bool

export SpinGrid
"""
A grid of spin states.
"""
immutable SpinGrid <: AbstractArray{Spin, 2}
    array::BitArray{2}

    """
        SpinGrid(iter)

    Construct a `SpinGrid` from an iterator.
    """
    SpinGrid(iter) = new(BitArray([x for x in iter]))
    """
        SpinGrid(::AbstractArray{Bool, 2})

    Construct a `SpinGrid` directly from an array.
    """
    SpinGrid(a::AbstractArray{Bool, 2}) = new(a)
end
"""
    SpinGrid(s, isize, jsize)

Construct an `isize` by `jsize` `SpinGrid` with every site having spin `s`.
"""
(::Type{SpinGrid})(s::Spin, isize::Int, jsize::Int) =
    SpinGrid(spinup(s) ? trues(isize, jsize) : falses(isize, jsize))

export RandomSpinGrid
"""
    RandomSpinGrid(isize, jsize)

Construct a `SpinGrid` of size `(isize, jsize)` where each site is initialized with a random spin.
"""
RandomSpinGrid(isize::Int, jsize::Int) = SpinGrid(rand(Bool) for i = 1:isize, j = 1:jsize)
"""
    RandomSpinGrid(s, p_spin, isize, jsize)

Construct a `SpinGrid` of size `(isize, jsize)` where each site has a probability `p_spin` of having
spin `s`.
"""
RandomSpinGrid(s::Spin, p_spin::Real, isize::Int, jsize::Int) =
    SpinGrid(rand(Bool) < p_spinup ? s : flipspin(s) for i = 1:isize, j = 1:jsize)

export MaybeRandomSpinGrid
"""
    MaybeRandomspinGrid(s, p_spin::Nullable(), isize, jsize)

Construct a `SpinGrid` of size `(isize, jsize)` where each site has spin `s` if `p_spin` is null,
otherwise initialize each site to `s` with a probability of `p_spin`.
"""
function MaybeRandomSpinGrid{T<:Real}(s::Spin, p_spin::Nullable{T}, isize::Int, jsize::Int)
    if isnull(p_spin)
        SpinGrid(s, isize, jsize)
    else
        RandomSpinGrid(s, get(p_spin), isize, jsize)
    end

end

export spinup, spindown
"""
    spinup(s)

Tests whether `s` is spin-up.
"""
spinup(s::Spin) = s
"""
    spindown(s)

Tests whether 's` is spin-down.
"""
spindown(s::Spin) = !s

export flipspin, flipspin!
"""
    flipspin(s)

Return the spin that is opposite of `s`.
"""
flipspin(s::Spin) = !s
"""
    flipspin(sg::SpinGrid, i[, j])

Change sg so that the specified spin site now has the opposite spin.

The spin site to consider is specified by `i` and `j`. If `j` is not present, then `i` may be either
a 2-tuple or a linear index into `sg` so that the spin site under consideration is either `sg[i...]`
or `sg[i]`. If `j` is present, then the spin site under consideration is `sg[i, j]`.
"""
flipspin!(sg::SpinGrid, i::Int, j::Int) = (sg[i, j] = flipspin(sg[i, j]))
flipspin!(sg::SpinGrid, ixs::NTuple{2, Int}) = flipspin!(sg, ixs...)
flipspin!(sg::SpinGrid, ixl::Int) = flipspin!(sg, ind2sub(sg, ixl)...)

## Various array-related methods.
Base.copy(sg::SpinGrid) = SpinGrid(copy(sg.array))
# Note that on these next two we make sg into a torus!
Base.getindex(sg::SpinGrid, i::Int, j::Int) = modgetindex(sg.array, i, j)
Base.setindex!(sg::SpinGrid, val::Spin, i::Int, j::Int) =
    sg.array[modindex(size(sg), i, j)...] = val
Base.size(sg::SpinGrid) = size(sg.array)
Base.length(sg::SpinGrid) = length(sg.array)
Base.similar{N, S}(sg::SpinGrid, ::Type{S}, dims::NTuple{N, Int}) = similar(sg.array, S, dims)

# For some reason having this uncommented causes it to be called and makes the code slower...?
#function Base.convert(::Type{SpinGrid}, a::AbstractArray{Bool, 2})
#    display(stacktrace())
#    SpinGrid(convert(BitArray{2}, a))
#end

## Delegate printing methods to the inner array.
Base.show(sg::SpinGrid) = show(sg.array)
Base.display(sg::SpinGrid) = display(sg.array)

## Write out sg as a raw RGB image. SPIN_UP => red, SPIN_DOWN => Blue.
function Base.write(stream::IO, sg::SpinGrid)
    sum(sg) do s
        write(stream, spinup(s) ? UInt8[0xff, 0x00, 0x00] : UInt8[0x00, 0x00, 0xff])
    end
end

## Bring in Neighborhood and make it a submodule. It has to be at the end here because it uses
## things from SpinMod.
include("Neighborhood.jl")
import .Neighborhood.neighborhood
export neighborhood

end # SpinMod
