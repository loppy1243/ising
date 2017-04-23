module SpinMod
import Loppy.Indexing.modindex, Loppy.Indexing.modgetindex

export SPIN_UP, SPIN_DOWN
const SPIN_UP = true
const SPIN_DOWN = false

export Spin
typealias Spin Bool

export SpinGrid
immutable SpinGrid <: AbstractArray{Spin, 2}
    array::BitArray{2}

    SpinGrid(iter) = new(BitArray([x for x in iter]))
    SpinGrid(a::AbstractArray{Bool, 2}) = new(a)
end
(::Type{SpinGrid})(s::Spin, isize::Int, jsize::Int) =
    SpinGrid(spinup(s) ? trues(isize, jsize) : falses(isize, jsize))

export RandomSpinGrid
RandomSpinGrid(isize::Int, jsize::Int) = SpinGrid(rand(Bool) for i = 1:isize, j = 1:jsize)
RandomSpinGrid(s::Spin, p_spin::Number, isize::Int, jsize::Int) =
    SpinGrid(rand(Bool) < p_spinup ? s : flipspin(s) for i = 1:isize, j = 1:jsize)

export MaybeRandomSpinGrid
function MaybeRandomSpinGrid(s::Spin, p_spin::Nullable{Float64}, isize::Int, jsize::Int)
    if isnull(p_spin)
        SpinGrid(s, isize, jsize)
    else
        RandomSpinGrid(s, get(p_spin), isize, jsize)
    end

end

export spinup, spindown
spinup(s::Spin) = s
spindown(s::Spin) = !s

export flipspin, flipspin!
flipspin(s::Spin) = !s
flipspin!(sg::SpinGrid, i::Int, j::Int) = (sg[i, j] = flipspin(sg[i, j]))
flipspin!(sg::SpinGrid, ixl::Int) = flipspin!(sg, ind2sub(sg, ixl)...)

Base.copy(sg::SpinGrid) = SpinGrid(copy(sg.array))
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

Base.show(sg::SpinGrid) = show(sg.array)
Base.display(sg::SpinGrid) = display(sg.array)

function Base.write(stream::IO, sg::SpinGrid)
    sum(sg) do s
        write(stream, spinup(s) ? UInt8[0xff, 0x00, 0x00, 0xff] : UInt8[0x00, 0x00, 0xff, 0xff])
    end
end

include("Neighborhood.jl")
import .Neighborhood.neighborhood
export neighborhood

end # SpinMod
