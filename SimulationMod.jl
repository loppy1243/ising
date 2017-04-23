module SimulationMod
using SpinMod, Ising

export Simulation
immutable Simulation
    init::SpinGrid
    grid::SpinGrid
    neigh_size::Int
    trans::Array{Nullable{Int}, 1}
    # Maps a two spin sites to their interaction strength.
    Jmap::Function
    # Maps a spin site to its β.
    βmap::Array{Float64, 2}

    function Simulation(init::SpinGrid, neigh_size::Int, trans::Int, Jmap::Function,
                        βmap::Array{Float64, 2})
        trans = Array{Nullable{Int}, 1}(trans+1)
        trans[1] = Nullable{Int}()
        new(init, copy(init), neigh_size, trans, Jmap, βmap)
    end
end
(::Type{Simulation})(init::SpinGrid, neigh_size::Int, trans::Int, J::Number, β::Number) =
    Simulation(init, neigh_size, trans, (sg, i, j, k, l) -> J, fill(β, size(init)))
(::Type{Simulation})(init::SpinGrid, neigh_size::Int, trans::Int, Jmap::Function, β::Number) =
    Simulation(init, neigh_size, trans, Jmap, fill(β, size(init)))
(::Type{Simulation})(init::SpinGrid, neigh_size::Int, trans::Int, J::Number,
                     βmap::Array{Float64, 2}) =
    Simulation(init, neigh_size, trans, (sg, i, j, k, l) -> J, βmap)

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
    # ffmpeg rendering command:
    # ffmpeg -f rawvideo # Input codec.
    #        -pix_fmt rgba
    #        -video_size ${IMG_WIDTH}x${IMG_HEIGHT} # Frame size.
    #        -framerate ${1/DURATION} # This tells ffmpeg what the simulation speed is.
    #        -i sim.raw # Raw video input file.
    #        -frames ${DURATION/TIME_STEP + 1} # +1 because of the initial state of the world.
    #        -f h264 # Output codec.
    #        -r 24 # Output framerate.
    #        -y sim.mp4 # Output file (-y is overwrite).
    ffmpeg = `ffmpeg -loglevel warning -f rawvideo -pixel_format rgba
                     -video_size $(size(sim.init, 2))x$(size(sim.init, 1)) -framerate $(tps)
                     -i $(infile) -frames $(length(sim.trans)) -f h264 -r 24 -y $(outfile)`
    run(ffmpeg)
end

export @simulate!
macro simulate!(sim)
    quote
        simm = $(esc(sim))::Simulation
        for i = 2:endof(simm.trans)
            s = rand(1:endof(simm.grid))
            dh = hamildiff(simm.grid, simm.neigh_size, simm.Jmap, s)

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
        simm = $(esc(sim))::Simulation
        for i = 2:endof(simm.trans)
            $s_sym = rand(1:endof(simm.grid))
            $dh_sym = hamildiff(simm.grid, simm.neigh_size, simm.Jmap, s)

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

end # SimulationMod
