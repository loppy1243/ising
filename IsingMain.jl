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

    function Simulation(init::SpinGrid, trans::Int, Jmap::Function, βmap::Array{Float64, 2})
        trans = Array{Nullable{Int}, 1}(trans+1)
        trans[1] = Nullable{Int}()
        new(init, copy(init), trans, Jmap, βmap)
    end
end

function (::Type{Simulation})(init::SpinGrid, trans::Int, J::Number, β::Number)
    Simulation(init, trans, (sg, i, j, k, l) -> J, fill(β, size(init)))
end

function (::Type{Simulation})(init::SpinGrid, trans::Int, Jmap::Function, β::Number)
    Simulation(init, trans, Jmap, fill(β, size(init)))
end

function (::Type{Simulation})(init::SpinGrid, trans::Int, J::Number, βmap::Array{Float64, 2})
    Simulation(init, trans, (sg, i, j, k, l) -> J, βmap)
end

function main()
    print("Running simulation... ")
    sim = inverseJ(20, 20, -2, 1e-10, 150000)
    println("Done.")

    #print("Writing raw video to file... ")
    #write("sim.raw", sim)
    #println("Done.")

    #println("Rendering video...")
    #render(sim, 100, "sim.raw", "sim.mp4")
    #println("Done.")
end

macro simulate!(sim)
    quote
        simm = ($sim)::Simulation
        for i = 2:endof(simm.trans)
            s = rand(1:endof(simm.grid))
            dh = hamildiff(simm.grid, simm.Jmap, s)

            simm.trans[i] = if rand() < exp(-simm.βmap[s]*dh)
                flipspin!(simm.grid, s)
                Nullable{Int}(s)
            else
                Nullable{Int}()
            end
        end

        simm
    end
end

macro simulate!(sim, body)
    @assert body.head == :block            &&
            body.args[2].head == :tuple    &&
            length(body.args[2].args) == 2
              
    s_sym = body.args[2].args[1]
    dh_sym = body.args[2].args[2]
    body = body.args[3:end]

    quote
        simm = ($sim)::Simulation
        for i = 2:endof(simm.trans)
            $s_sym = rand(1:endof(simm.grid))
            $dh_sym = hamildiff(simm.grid, simm.Jmap, s)

            simm.trans[i] = if rand() < exp(-simm.βmap[$s_sym]*$dh_sym)
                flipspin!(simm.grid, $s_sym)
                Nullable{Int}($s_sym)
            else
                Nullable{Int}()
            end

            $(body...)
        end

        simm
    end
end

function constantβ(isize::Int, jsize::Int, β::Number, trans::Int;
                   initspin::Spin=SPIN_DOWN, J::Number=1.0, p_spin=Nullable{Float64}())
    init_sg = MaybeRandomSpinGrid(initspin, p_spin, isize, jsize)

    @simulate! Simulation(init_sg, trans, J, β)
end

function rangedβ{T<:Number}(isize::Int, jsize::Int, βr::Range{T};
                            initspin::Spin=SPIN_DOWN, J::Number=1.0, p_spin=Nullable{Float64}())
    init_sg = MaybeRandomSpinGrid(initspin, p_spin, isize, jsize)
    sim = Simulation(init_sg, length(βr), J, βr.start)

    @simulate! sim begin s, dh
        sim.βmap .+= βr.step
    end
end

function inverseJ(isize::Int, jsize::Int, J_pow::Int, β::Number, trans::Int;
                  initspin::Spin=SPIN_DOWN, p_spin=Nullable{Float64}())
    init_sg = MaybeRandomSpinGrid(initspin, p_spin, isize, jsize)

    function Jmap(sg::SpinGrid, i::Int, j::Int, k::Int, l::Int)
        hypot(i - k, j - l)^J_pow
    end

    @simulate! Simulation(init_sg, trans, Jmap, β)
end

function Base.write(stream::IO, sim::Simulation)
    sim.grid .= sim.init

    sum(sim.trans) do s
        if !isnull(s)
            flipspin!(sim.grid, get(s))
        end
    
        write(stream, sim.grid)
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
