module Neighborhood
import Loppy.Indexing.modindex
using SpinMod

immutable IndexIter
    gridsize::NTuple{2, Int}
    current::NTuple{2, Int}
    start::NTuple{2, Int}
    ending::NTuple{2, Int}
    n::Int
end

export neighborhood
neighborhood(sg::SpinGrid, n::Int, i::Int, j::Int) =
    IndexIter(size(sg), (i, j), (i-n, j-n), (i+n, j+n), n)
neighborhood(sg::SpinGrid, n::Int, ixs::NTuple{2, Int}) = neighborhood(sg, n, ixs...)
neighborhood(sg::SpinGrid, n::Int, ixl::Int) = neighborhood(sg, n, ind2sub(sg, ixl)...)

Base.start(iter::IndexIter) = iter.start
Base.done(iter::IndexIter, state::NTuple{2, Int}) =
    state[2] == iter.ending[2] && state[1] > iter.current[1]

function Base.next(iter::IndexIter, state::NTuple{2, Int})
    while true
        dist = hypot(state[1] - iter.current[1], state[2] - iter.current[2])
        if state != iter.current && dist <= iter.n
            return (modindex(iter.gridsize, state), (state[1] + 1, state[2]))
        elseif state[1] > iter.current[1]
            state = (iter.start[1], state[2] + 1)
        elseif state[2] > iter.ending[2]
            throw("IMPOSSIBLE")
        else
            state = (state[1] + 1, state[2])
        end
    end
end

end # Neighborhood
