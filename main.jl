# Load your local modules
include("src/parameters.jl")
include("src/constants.jl")
include("src/target.jl")
include("src/swerling.jl")
include("src/radarprocessing.jl")
include("src/Radar.jl")
include("src/RadarUI.jl")
# include("src/utils.jl")

# using Radar
# using RadarProcessingChain
# using RadarUI
# using Targets
using REPL.TerminalMenus
using GLMakie

println("\e[38;5;83mWelcome to my radar simulator, Radar.jl!")
mode = while(true)
    options = ["Spawn", "Game", "Quit"]
    selection = request("Choose a Radar mode:", RadioMenu(options))
    chosen_mode = options[selection]
    
    if options[selection] == "Quit"
        println("Exiting the radar simulator.\e[37m")
        exit(0)
    elseif options[selection] == "Spawn"
        println("Spawn mode activated!")
    else
        println("Game mode activated!")
    end
    
    return chosen_mode
end


liveTargets = Vector{Targets.Target}(
                    [Targets.CartesianTarget(randn() * 100e3, randn() * 100e3, randn() * 5000, randn() * 5000; swerlingModel = :Swerling0, rcs = 1.0, status = :unknown) for _ in 1:20]  # Randomly generated targets
                ) # Live targets to be processed

GC.enable(true)
# Initialize system
# radar = Radar.initRadar()
radar = Radar.pulseRadar()

# targets = generateTargets()
fig, fig2, azimuth, beamPoints, targetDots, echo, compression, doppler = RadarUI.setupRadarDisplay(radar, liveTargets; initialAzimuth = radar.azimuth)

@async begin
    try 
        RadarUI.updateDisplay(azimuth, beamPoints, targetDots, echo, compression, doppler, radar, liveTargets)
    catch e
        println("Error in updateDisplay: $e")
    end
end
# Show the UI

# Wow, this works!
display(GLMakie.Screen(), fig)
display(GLMakie.Screen(), fig2)


# Start live radar
Radar.runRadar(radar, liveTargets)

@async RadarUI.updateTargets(liveTargets)

println("Press Enter to close the radar display...")

readline()
println("Radar display closed.\e[37m")
