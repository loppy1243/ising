"""
The 
"""
module Ising
using SpinMod

# Not (currently) used.
#"""
#    hamil(sg::SpinGrid, J::Function)
#
#Computes the Hamiltonian of `sg` with interaction strengths determined by `J`.
#
#`J` should take `sg` and two spin sites and return a `Real` with the following signature:
#    J(::SpinGrid, ::Int, ::Int, ::Int, ::Int)::Real
#"""
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

"""
    localhamil(sg::SpinGrird, neigh_size, J::Function, i[, j])

Computes the contribution of a neighborhood of size `neigh_size` to the Hamiltonian of `sg`.

`J` should take `sg` and two spin sites and return a `Real` with the following signature:
    J(::SpinGrid, ::Int, ::Int, ::Int, ::Int)::Real
The spin site to consider is specified by `i` and `j`. If `j` is not present, then `i` may be either
a 2-tuple or a linear index into `sg` so that the spin site under consideration is either `sg[i...]`
or `sg[i]`. If `j` is present, then the spin site under consideration is `sg[i, j]`.
"""
localhamil(sg::SpinGrid, neigh_size::Int, J::Function, ixs::NTuple{2, Int}) =
    localhamil(sg, neigh_size, J, ixs...)
localhamil(sg::SpinGrid, neigh_size::Int, J::Function, ixl::Int) =
    localhamil(sg, neigh_size, J, ind2sub(sg, ixl)...)
function localhamil(sg::SpinGrid, neigh_size::Int, J::Function, i::Int, j::Int)
    # Please don't forget the minus out front.
    #
    # Each site in sg is a Bool, but we want to treat them as +-1 and take the product of two
    # of them; this is exactly xor-ing the two sites. The result then gets converted to an
    # appropriate Integer to be multiplied by J().
    -(sum(neighborhood(sg, neigh_size, i, j)) do ixs
        (J(sg, i, j, ixs...) * xor(sg[i, j], sg[ixs...]))
    end)
end

export hamildiff
"""
    hamildiff(sg::SpinGrid, neigh_size, J::Function, i[, j])

Computes the difference induced in the Hamiltonian of `sg` if one spin site is flipped.

`J` should take `sg` and two spin sites and return a `Real` with the following signature:
    J(::SpinGrid, ::Int, ::Int, ::Int, ::Int)::Real
This difference only depends on the neighborhood, of size `neigh_size`, of the site being flipped.
The spin site to consider is specified by `i` and `j`. If `j` is not present, then `i` may be either
a 2-tuple or a linear index into `sg` so that the spin site under consideration is either `sg[i...]`
or `sg[i]`. If `j` is present, then the spin site under consideration is `sg[i, j]`.
"""
hamildiff(sg::SpinGrid, neigh_size::Int, J::Function, ixl::Int) =
    hamildiff(sg, neigh_size, J, ind2sub(sg, ixl)...)
function hamildiff(sg::SpinGrid, neigh_size::Int, J::Function, i::Int, j::Int)
    h1 = localhamil(sg, neigh_size, J, i, j)
    flipspin!(sg, i, j)
    h2 = localhamil(sg, neigh_size, J, i, j)
    flipspin!(sg, i, j)

    h2 - h1
end

end # Ising
