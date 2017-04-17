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

const OUT_FILE = "sim"
const TRANS = 50000
const TRANS_PER_SEC = 100 

function main()
    init_sg = SpinGrid(SPIN_DOWN, 10, 10)
    sg = copy(init_sg)
    states = Array{Nullable{NTuple{2, Int}}, 1}(TRANS+1)
    states[1] = Nullable{NTuple{2, Int}}()

    print("Running simulation... ")
    for i = 1:TRANS
        println("Iteration: ", i)
        pos = ind2sub(sg, rand(1:endof(sg)))
        aff = spinflipaff(sg, pos...)

        states[i+1] = if aff > 0 || rand() < exp(aff)
            flipspin(sg, pos...)
            Nullable(pos)
        else
            Nullable{NTuple{2, Int}}()
        end
    end
    println("Done.")

    ffmpeg = `ffmpeg -loglevel warning -f rawvideo -pix_fmt rgba
                     -s $(size(init_sg, 2))x$(size(init_sg, 1)) -framerate $(TRANS_PER_SEC)
                     -i $(OUT_FILE).mp4 -frames $(TRANS+1) -f h264 -r 24 -y $(OUT_FILE).mp4`

    print("Writing raw video to file... ")
    open("$(OUT_FILE).raw", "w") do out
        for s in states
            if !isnull(s)
                flipspin(init_sg, get(s)...)
            end

            write(out, init_sg)
        end
    end
    println("Done.")

    println("Converting raw video to h264...")
    run(ffmpeg)
    println("Done.")
end

end # IsingMain

if !isinteractive()
    IsingMain.main()
end
