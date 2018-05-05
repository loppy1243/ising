"""
Contains the machinery of the simulation. Same naming caveat as SpinMod.
"""
module SimulationMod
using SpinMod, Ising

export Simulation
"""
Holds all the information necessary to run construct a simulation.
"""
immutable Simulation
    # The initial state.
    init::SpinGrid
    # The intermediary state. This is volatile.
    grid::SpinGrid
    # The size (radius) of the neighborhood this simulation uses.
    neigh_size::Int
    # An array of the (linear) indices of the spin sites that get flipped in chronological order.
    # Always starts with Nullable() to represent the initial state.
    trans::Array{Nullable{Int}, 1}
    # Maps two spin sites to their interaction strength. Should have the signature
    #     Jmap(::SpinGrid, ::Int, ::Int, ::Int, ::Int)
    Jmap::Function
    # Maps a spin site to its inverse temperature.
    βmap::Array{Float64, 2}

    """
        Simulation(init, neigh_size, trans, Jmap, βmap)

    Construct a Simulation with initial state `init`, neighborhood size `neigh_size`, total
    transitions `trans`, and (variable) interaction strength and inverse temperature specified by
    `Jmap` and `βmap`.
    """
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

## Write a Simulation by writing each intermediate state in turn.
function Base.write(stream::IO, sim::Simulation)
    sim.grid .= sim.init

    # sum() so that we ultimately get the total number of bytes written.
    sum(sim.trans) do s
        if !isnull(s)
            flipspin!(sim.grid, get(s))
        end
    
        write(stream, sim.grid)
    end
end

"""
    write_mp4(sim::Simulation, tps::Int, outfile)

Write out H264 video of the simulation to `outfile` as an MP4 file.

This is accomplished by writing out a raw RGB video and then using `ffmpeg` to convert. `tps`
is the number of real-time transitions per second, i.e., the simulation speed.
"""
function render(sim::Simulation, tps::Int, outfile::AbstractString)
    # ffmpeg command:
    # ffmpeg -f rawvideo # Input codec.
    #        -pix_fmt rgb24 # Pixel format.
    #        -video_size ${IMG_WIDTH}x${IMG_HEIGHT} # Frame size.
    #        -framerate ${1/DURATION} # This tells ffmpeg what the simulation speed is.
    #        -i sim.raw # Raw video input file.
    #        -frames ${TRANS} # +1 because of the initial state of the world.
    #        -f h264 # Output codec.
    #        -r 24 # Output framerate.
    #        -y sim.mp4 # Output file (-y is overwrite).

    mktemp() do rawfile, rawfile_stream
        ffmpeg = `ffmpeg -loglevel warning -f rawvideo -pixel_format rgb24
                         -video_size $(size(sim.init, 2))x$(size(sim.init, 1))
                         -framerate $(tps) -i $(rawfile) -frames $(length(sim.trans)) -f h264
                         -r 24 -y $(outfile)`
        write(rawfile_stream, sim)
        run(ffmpeg)
    end
end

export @simulate!
"""
    @simulate! sim

Run the simulation specified by `sim`.

Returns the simulation for convenience. Cobbles `sim.grid`.
"""
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

"""
    @simulate! sim begin s, dh
        [block...]
    end

Run the simulation specified by `sim`, running `block` after each iteration.

`s` is the linear index of the spin site that was just checked and `dh` is the change in the
Hamiltonian induced.
"""
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
