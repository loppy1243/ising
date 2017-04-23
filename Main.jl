include("includes.jl")

module IsingMain
import SpinMod
using SimulationMod, SpinMod

function main()
    print("Running simulation... ")
    sim = inverseJ(100, 100, -2, 1e-10, 150000)
    println("Done.")

    print("Writing raw video to file... ")
    write("sim.raw", sim)
    println("Done.")

    println("Rendering video...")
    SimulationMod.render(sim, 100, "sim.raw", "sim.mp4")
    println("Done.")
end

function constantβ(isize::Int, jsize::Int, β::Number, trans::Int;
                   initspin::Spin=SPIN_DOWN, J::Number=1.0, p_spin=Nullable{Float64}(),
                   neigh_size::Int=6)
    init_sg = MaybeRandomSpinGrid(initspin, p_spin, isize, jsize)

    @simulate! Simulation(init_sg, neigh_size, trans, J, β)
end

function rangedβ{T<:Number}(isize::Int, jsize::Int, βr::Range{T};
                            initspin::Spin=SPIN_DOWN, J::Number=1.0, p_spin=Nullable{Float64}(),
                            neigh_size::Int=6)
    init_sg = MaybeRandomSpinGrid(initspin, p_spin, isize, jsize)
    sim = Simulation(init_sg, neigh_size, length(βr), J, βr.start)

    @simulate! sim begin s, dh
        sim.βmap .+= βr.step
    end
end

function inverseJ(isize::Int, jsize::Int, J_pow::Int, β::Number, trans::Int;
                  initspin::Spin=SPIN_DOWN, p_spin=Nullable{Float64}(), neigh_size::Int=6)
    init_sg = MaybeRandomSpinGrid(initspin, p_spin, isize, jsize)

    function Jmap(sg::SpinGrid, i::Int, j::Int, k::Int, l::Int)
        hypot(i - k, j - l)^J_pow
    end

    @simulate! Simulation(init_sg, neigh_size, trans, Jmap, β)
end

function βflow(isize::Int, jsize::Int, J::Float64, β0::Number, trans::Int;
               initspin::Spin=SPIN_DOWN, p_spin=Nullable{Float64}(), neigh_size::Int=6)
    init_sg = MaybeRandomSpinGrid(initspin, p_spin, isize, jsize)
    sim = Simulation(init_sg, neigh_size, trans, J, β0)

    @simulate! sim begin s, dh
        for ixs in SpinMod.neighborhood(sim.grid, sim.neigh_size, s)
            sim.βmap[ixs...] = 1 / (1/sim.βmap[ixs...] - dh/trunc(pi*neigh_size^2))
        end
    end
end

end # IsingMain

if !isinteractive(); IsingMain.main() end
