module Ising
using Unitful
export Spin, SpinGrid, RandomSpinGrid, spinup, spindown, SPIN_UP, SPIN_DOWN, flipspin, logprob,
       spinflipaff

# Inverse temperature of the system.
const INV_TEMP = 0.5u"1/eV"
# Interaction strength of the system (assuming interactions are homogenous and isotropic).
const INT_STR = 1e0u"eV"
const SPIN_UP = true
const SPIN_DOWN = false
const NEIGH_SIZE = 3

typealias Spin Bool

spinup(s::Spin) = s
spindown(s::Spin) = !s

immutable SpinGrid <: AbstractArray{Spin, 2}
    array::BitArray{2}

    SpinGrid(iter) = new(BitArray([x for x in iter]))
    SpinGrid(a::AbstractArray{Bool, 2}) = new(a)
end

function (::Type{SpinGrid})(s::Spin, isize::Int, jsize::Int)
    SpinGrid(spinup(s) ? trues(isize, jsize) : falses(isize, jsize))
end

# For some reason having this uncommented causes it to be called and makes the code slower...?
#function Base.convert(::Type{SpinGrid}, a::AbstractArray{Bool, 2})
#    display(stacktrace())
#    SpinGrid(convert(BitArray{2}, a))
#end
Base.show(sg::SpinGrid) = show(sg.array)
Base.display(sg::SpinGrid) = display(sg.array)

function RandomSpinGrid(isize::Int, jsize::Int)
    SpinGrid(rand(Bool) for i = 1:isize, j = 1:jsize)
end

function RandomSpinGrid(p_spinup::Number, isize::Int, jsize::Int)
    SpinGrid(rand(Bool) < p_spinup ? SPIN_UP : SPIN_DOWN for i = 1:isize, j = 1:jsize)
end

Base.copy(sg::SpinGrid) = SpinGrid(copy(sg.array))
Base.getindex(sg::SpinGrid, i::Int, j::Int) = modgetindex(sg.array, i, j)

function Base.setindex!(sg::SpinGrid, val::Spin, i::Int, j::Int)
    sg.array[modindex(size(sg), i, j)...] = val
end

flipspin(s::Spin) = !s
flipspin(sg::SpinGrid, i::Int, j::Int) = (sg[i, j] = flipspin(sg[i, j]))

Base.size(sg::SpinGrid) = size(sg.array)
Base.length(sg::SpinGrid) = length(sg.array)
Base.similar{N, S}(sg::SpinGrid, ::Type{S}, dims::NTuple{N, Int}) = similar(sg.array, S, dims)

function Base.write(stream::IO, sg::SpinGrid)
    sum(sg) do s
        write(stream, spinup(s) ? UInt8[0xff, 0x00, 0x00, 0xff] : UInt8[0x00, 0x00, 0xff, 0xff])
    end
end

@generated function modindex{N}(dims::NTuple{N, Int}, ixs::Vararg{Int, N})
    js(i) = :(mod(ixs[$i], dims[$i]))

    quote
        tuple($((:($(js(i)) == 0 ? dims[$i] : $(js(i))) for i = 1:N)...))
    end
end

function modindex{T, N}(a::AbstractArray{T, N}, ixs::Vararg{Int, N})
    modindex(size(a), ixs...)
end

@generated function modgetindex{T, N}(a::AbstractArray{T, N}, ixs::Vararg{Int, N})
    js(i) = :(mod(ixs[$i], size(a, $i)))

    quote
        a[$((:($(js(i)) == 0 ? size(a, $i) : $(js(i))) for i = 1:N)...)]
    end
end

#modindex{N}(x, ixs::NTuple{N, Int}) = modindex(x, ixs...)
#modgetindex{N}(x, ixs::NTuple{N, Int}) = modgetindex(x, ixs...)

# Don't know is this is needed...
#immutable NeighborhoodIter
#end

immutable NeighborhoodIndexIter
    grid::SpinGrid
    current::NTuple{2, Int}
    start::NTuple{2, Int}
    ending::NTuple{2, Int}
end

function neighborhood(sg::SpinGrid, i::Int, j::Int, n::Int)
    NeighborhoodIndexIter(sg, (i, j), (i-n, j-n), (i+n, j+n))
end

function Base.start(iter::NeighborhoodIndexIter)
    iter.start
end

function Base.next(iter::NeighborhoodIndexIter, state::NTuple{2, Int})
    next_state = if state[1] < iter.ending[1]
        (state[1] + 1, state[2])
    else
        (iter.start[1], state[2] + 1)
    end

    if state == iter.current
        next(iter, next_state)
    else
        (modindex(iter.grid, state...), next_state)
    end
end

function Base.done(iter::NeighborhoodIndexIter, state::NTuple{2, Int})
    state[2] > iter.ending[2]
end

function hamil(sg::SpinGrid)
    visited = falses(size(sg))

    sum = 0.0
    for i in 1:endof(sg)
        if visited[i]
            continue
        end

        for jxs in neighborhood(sg, ind2sub(sg, i)..., NEIGH_SIZE)
            if visited[jxs...]
                continue
            end

            sum += sg[i] $ sg[jxs...]
        end
        visited[i] = true
    end

    -INT_STR*sum
end

function localhamil(sg::SpinGrid, i::Int, j::Int)
    -INT_STR * sum(neighborhood(sg, i, j, NEIGH_SIZE)) do ixs
        sg[i, j] $ sg[ixs...]
    end
end

logprob(sg::SpinGrid) = -INV_TEMP*hamil(sg)

function spinflipaff(sg::SpinGrid, i::Int, j::Int)
    h1 = localhamil(sg, i, j)
    flipspin(sg, i, j)
    h2 = localhamil(sg, i, j)
    flipspin(sg, i, j)

    -INV_TEMP*(h2 - h1)
end

end # Ising
