include("includes.jl")

"""
A Monte Carlo simulation following the Ising model without magnetic field.
"""
module IsingMain
import SpinMod
using SimulationMod, SpinMod

function main()
    print("Running simulation... ")
    # Choose the simulation to run here.
    sim = constantβ(20, 20, 1e15, 15000)
    println("Done.")

    print("Writing video as mp4...")
    SimulationMod.write_mp4(sim, 400, "sim.mp4")
    println("Done.")
end

"""
    constantβ(isize, jsize, β, trans; initspin=SPIN_DOWN, J=1.0, p_spin=Nullable(), neigh_size=6)

Construct and run a simulation over `trans` transitions with an `isize` by `jsize` grid, constant
``β`` and ``J``, and a neighborhood size of `neigh_size`.

`initspin` specifies the initial spin state of the system. If `p_spin` is not `Nullable()`, then the
simulation grid is initialized randomly with a given site having probability `p_spin` of being in
state `initspin`.
"""
function constantβ(isize::Int, jsize::Int, β::Real, trans::Int;
                   initspin::Spin=SPIN_DOWN, J::Real=1.0, p_spin=Nullable{Float64}(),
                   neigh_size::Int=6)
    init_sg = MaybeRandomSpinGrid(initspin, p_spin, isize, jsize)

    @simulate! Simulation(init_sg, neigh_size, trans, J, β)
end

"""
    rangedβ(isize, jsize, βr::Range; initspin=SPIN_DOWN, J=1.0, p_spin=Nullable(), neigh_size=6)

Similar to `constantβ`, except for each transition of the system ``β`` steps through `βr`.
"""
function rangedβ{T<:Real}(isize::Int, jsize::Int, βr::Range{T};
                            initspin::Spin=SPIN_DOWN, J::Real=1.0, p_spin=Nullable{Float64}(),
                            neigh_size::Int=6)
    init_sg = MaybeRandomSpinGrid(initspin, p_spin, isize, jsize)
    sim = Simulation(init_sg, neigh_size, length(βr), J, βr.start)

    @simulate! sim begin s, dh
        sim.βmap .+= βr.step
    end
end

"""
    inverseJ(isize, jsize, J_pow, β, trans; initspin=SPIN_DOWN, p_spin=Nullable(), neigh_size=6)

Similar to `constantβ`, except J now depends on the `J_pow`-th power of the distance between two
interacting sites.

For one site with indices `(i, j)` and another with indices `(k, l)`, the interaction
strength between them is ``J = ((i-k)^2 + (j-l)^2)^{\frac{\text{J\_pow}}{2}}``.
"""
function inverseJ(isize::Int, jsize::Int, J_pow::Int, β::Real, trans::Int;
                  initspin::Spin=SPIN_DOWN, p_spin=Nullable{Float64}(), neigh_size::Int=6)
    init_sg = MaybeRandomSpinGrid(initspin, p_spin, isize, jsize)

    function Jmap(sg::SpinGrid, i::Int, j::Int, k::Int, l::Int)
        hypot(i - k, j - l)^J_pow
    end

    @simulate! Simulation(init_sg, neigh_size, trans, Jmap, β)
end


"""
    βflow(isize, jsize, β0, trans; initspin=SPIN_DOWN, J=1.0, p_spin=Nullable(), neigh_size=6)

Similar to `constantβ`, except every spin site has its own ``β`` which models a flow of energy.

Each spin site starts with ``β = \tt{β0}``. If a site transitions, then it evenly divides the
resulting change in the Hamiltonian amongst its neighbors and adds that to each of their ``β``.
"""
function βflow(isize::Int, jsize::Int, β0::Real, trans::Int;
               initspin::Spin=SPIN_DOWN, J::Real=1.0, p_spin=Nullable{Float64}(),
               neigh_size::Int=6)
    init_sg = MaybeRandomSpinGrid(initspin, p_spin, isize, jsize)
    sim = Simulation(init_sg, neigh_size, trans, J, β0)

    @simulate! sim begin s, dh
        for ixs in SpinMod.neighborhood(sim.grid, sim.neigh_size, s)
            sim.βmap[ixs...] = 1 / (1/sim.βmap[ixs...] - dh/trunc(pi*neigh_size^2))
        end
    end
end

end # IsingMain

# Run the program if we're not at the REPL.
if !isinteractive(); IsingMain.main() end
