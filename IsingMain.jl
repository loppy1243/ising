# ffmpeg rendering command:
# ffmpeg -f rawvideo # Input codec.
#        -pix_fmt rgba
#        -s ${IMG_WIDTH}x${IMG_HEIGHT} # Frame size.
#        -framerate ${1/DURATION} # This tells ffmpeg what the simulation speed is.
#        -i sim.raw # Raw video input file.
#        -frames ${DURATION/TIME_STEP + 1} # +1 because of the initial state of the world.
#        -f h264 # Output codec.
#        -r 24 # Output framerate.
#        -y sim.mp4 # Output file.

module IsingMain
include("Ising.jl")
using .Ising

immutable Simulation
    init::SpinGrid
    grid::SpinGrid
    trans::Array{Nullable{Int}, 1}
    # Maps a two spin sites to their interaction strength.
    Jmap::Function
    # Maps a spin site to its β.
    βmap::Array{Float64, 2}

    Simulation(init::SpinGrid, trans::Int, Jmap::Array{Array{Float64, 2}, 2}, βmap::Array{Float64, 1})
        trans = Array{Nullable{Int}, 1}(trans+1)
        trans[1] = init
        new(init, copy(init), Jmap, βmap)
    end
end

function (::Type{Simulation})(init::SpinGrid, trans::Int, J::Number, β::Number)
    Simulation(init, trans, fill(fill(J, size(init)), size(init)), fill(β, size(init)))
end

function (::Type{Simulation})(init::SpinGrid, trans::Int, Jmap::Array{Array{Float64, 2}, 2}, β::Number)
    Simulation(init, trans, Jmap, fill(β, size(init)))
end

function (::Type{Simulation})(init::SpinGrid, trans::Int, J::Number, βmap::Array{Float64, 2})
    Simulation(init, trans, fill(fill(J, size(init)), size(init)), βmap)
end

function main()
    print("Running simulation... ")
    sim = rangedβ(100, 100, linspace(0.01, 1, 50000))
    println("Done.")

    print("Writing raw video to file... ")
    write("sim.raw", sim)
    println("Done.")

    println("Rendering video...")
    render(sim, 100, "sim.raw", "sim.mp4")
    println("Done.")
end

function simulate!(sim::Simulation)
    simulate!(x -> nothing, sim)
end

function simulate!(update::Function, sim::Simulation)
    for i = 1:endof(sim.trans)-1
        update(i)
        s = rand(1:endof(sim.grid))
        aff = spinflipaff(sg, β, J, s)

        sim.trans[i+1] = if rand() < exp(aff)
            transpin!(sim.grid, s)
            Nullable(pos)
        else
            Nullable{NTuple{2, Int}}()
        end
    end

    sim
end

function constantβ(isize::Int, jsize::Int, β::Number, trans::Int;
                   initspin::Spin=SPIN_DOWN, J::Number=1.0, p_spin=Nullable{Float64}())
    init_sg = MaybeRandomSpinGrid(initspin, p_spin, isize, jsize)
    sim = Simulation(init_sg, trans, J, β)

    simulate!(sim, trans)
end

function rangedβ(isize::Int, jsize::Int, βr::Range{Number};
                 initspin::Spin=SPIN_DOWN, J::Number=1.0, p_spin=Nullable{Float64}())
    init_sg = MaybeRandomSpinGrid(initspin, p_spin, isize, jsize)
    sim = Simulation(init_sg, length(βr), J, βr.start)

    simulate!(sim) do _
        sim.βmap .+= βr.step
    end
end

function inverseJ(isize::Int, jsize::Int, β::Number, J_pow::Int;
                  initspin::Spin=SPIN_DOWN, p_spin=Nullable{Float64}())
    init_sg = MaybeRandomSpinGrid(initspin, p_spin, isize, jsize)
end

function Base.write(stream::IO, sim::Simulation)
    sim.grid .= sim.init

    sum(sim.trans) do f
        if !isnull(f)
            transpin!(sg, get(f)...)
        end
    
        write(stream, sg)
    end
end

# tps -- Transitions per second.
function render(sim::Simulation, tps::Int, infile::AbstractString, outfile::AbstractString)
    ffmpeg = `ffmpeg -loglevel warning -f rawvideo -pixel_format rgba
                     -video_size $(size(sim.init, 2))x$(size(sim.init, 1)) -framerate $(tps)
                     -i $(infile) -frames $(length(sim.trans)) -f h264 -r 24 -y $(outfile)`
    run(ffmpeg)
end

end # IsingMain

if !isinteractive()
    IsingMain.main()
end
