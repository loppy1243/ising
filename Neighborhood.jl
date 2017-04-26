"""
Submodule of SpinMod. Exports a function that returns an iterator over the neighborhood of a spin
site in a SpinGrid.
"""
module Neighborhood
import Loppy.Indexing.modindex
using SpinMod

"""
Iterator over a neighborhood that gives an index at each step.
"""
immutable IndexIter
    # Size of the spin grid.
    gridsize::NTuple{2, Int}
    # The neighborhood size.
    n::Int
    # The site which whose neighborhood we're considering.
    current::NTuple{2, Int}
    # The starting corner of the rectangle enclosing the neighborhood.
    start::NTuple{2, Int}
    # The ending corner of the rectangle enclosing the neighborhood.
    ending::NTuple{2, Int}
end

export neighborhood
"""
    neighborhood(sg::SpinGrid, neigh_size, i[, j])

Return an iterator over the neighborhood of size `neigh_size` of the specified spin site in `sg`

The spin site to consider is specified by `i` and `j`. If `j` is not present, then `i` may be either
a 2-tuple or a linear index into `sg` so that the spin site under consideration is either `sg[i...]`
or `sg[i]`. If `j` is present, then the spin site under consideration is `sg[i, j]`.
"""
neighborhood(sg::SpinGrid, n::Int, i::Int, j::Int) =
    IndexIter(size(sg), n, (i, j), (i-n, j-n), (i+n, j+n))
neighborhood(sg::SpinGrid, n::Int, ixs::NTuple{2, Int}) = neighborhood(sg, n, ixs...)
neighborhood(sg::SpinGrid, n::Int, ixl::Int) = neighborhood(sg, n, ind2sub(sg, ixl)...)

## The necessary methods to override to implement the iterator interface.

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
