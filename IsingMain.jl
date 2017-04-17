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
    flips::Array{Nullable{NTuple{2, Int}}, 1}
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

function constantβ(isize::Int, jsize::Int, β::Number, trans::Int;
                   initspin::Spin=SPIN_DOWN, J::Number=1.0)
    init_sg = SpinGrid(initspin, isize, jsize)
    sg = copy(init_sg)
    states = Array{Nullable{NTuple{2, Int}}, 1}(trans+1)
    states[1] = Nullable{NTuple{2, Int}}()

    for i = 1:trans
        pos = ind2sub(sg, rand(1:endof(sg)))
        aff = spinflipaff(sg, pos...)

        states[i+1] = if aff > 0 || rand() < exp(aff)
            flipspin(sg, pos...)
            Nullable(pos)
        else
            Nullable{NTuple{2, Int}}()
        end
    end

    Simulation(init_sg, states)
end

function rangedβ(isize::Int, jsize::Int, iter; initspin::Spin=SPIN_DOWN, J::Number=1.0)
    init_sg = SpinGrid(initspin, isize, jsize)
    sg = copy(init_sg)
    states = Array{Nullable{NTuple{2, Int}}, 1}(length(iter) + 1)
    states[1] = Nullable{NTuple{2, Int}}()

    for (i, β) in enumerate(iter)
        pos = ind2sub(sg, rand(1:endof(sg)))
        aff = spinflipaff(sg, β, J, pos...)

        states[i+1] = if aff > 0 || rand() < exp(aff)
            flipspin(sg, pos...)
            Nullable(pos)
        else
            Nullable{NTuple{2, Int}}()
        end
    end

    Simulation(init_sg, states)
end

function Base.write(stream::IO, sim::Simulation)
    sg = copy(sim.init)

    sum(sim.flips) do f
        if !isnull(f)
            flipspin(sg, get(f)...)
        end
    
        write(stream, sg)
    end
end

# tps -- Transitions per second.
function render(sim::Simulation, tps::Int, infile::AbstractString, outfile::AbstractString)
    ffmpeg = `ffmpeg -loglevel warning -f rawvideo -pixel_format rgba
                     -video_size $(size(sim.init, 2))x$(size(sim.init, 1)) -framerate $(tps)
                     -i $(infile) -frames $(length(sim.flips)) -f h264 -r 24 -y $(outfile)`
    run(ffmpeg)
end

end # IsingMain

if !isinteractive()
    IsingMain.main()
end
