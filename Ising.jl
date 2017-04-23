module Ising
using SpinMod

# Not (currently) used.
#function hamil(sg::SpinGrid, J::Function)
#    visited = falses(size(sg))
#
#    sum = 0.0
#    for i in 1:endof(sg)
#        if visited[i]
#            continue
#        end
#
#        for jxs in neighborhood(sg, NEIGH_SIZE, ind2sub(sg, i)...)
#            if visited[jxs...]
#                continue
#            end
#
#            sum += J(sg, ind2sub(sg, i)..., jxs...) * (sg[i] $ sg[jxs...])
#        end
#        visited[i] = true
#    end
#
#    # Please don't forget the minus sign here.
#    -sum
#end

localhamil(sg::SpinGrid, neigh_size::Int, J::Array{Float64, 2}, ixs::NTuple{2, Int}) =
    localhamil(sg, neigh_size, J, ixs...)
localhamil(sg::SpinGrid, neigh_size::Int, J::Array{Float64, 2}, ixl::Int) =
    localhamil(sg, neigh_size, J, ind2sub(sg, ixl)...)
function localhamil(sg::SpinGrid, neigh_size::Int, J::Array{Float64, 2}, i::Int, j::Int)
    # Please don't forget the minus sign here.
    -(sum(neighborhood(sg, neigh_size, i, j)) do ixs
        (J[ixs...] * (sg[i, j] $ sg[ixs...]))
    end)
end

export hamildiff
hamildiff(sg::SpinGrid, neigh_size::Int, J::Function, ixl::Int) =
    hamildiff(sg, neigh_size, J, ind2sub(sg, ixl)...)

let Jinit = true, J = Array{Float64, 2}()
global hamildiff
function hamildiff(sg::SpinGrid, neigh_size::Int, Jmap::Function, i::Int, j::Int)
    if Jinit
        J = Array{Float64, 2}(size(sg))
    end

    for ixl = 1:endof(J)
        J[ixl] = Jmap(sg, i, j, ind2sub(sg, ixl)...)
    end

    h1 = localhamil(sg, neigh_size, J, i, j)
    flipspin!(sg, i, j)
    h2 = localhamil(sg, neigh_size, J, i, j)
    flipspin!(sg, i, j)

    h2 - h1
end
end

end # Ising
