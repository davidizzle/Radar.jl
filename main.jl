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
mode = ""
while(true)
    global mode
    options = ["Visualize & Spawn", "Game", "Quit"]
    selection = request("Choose a Radar mode:", RadioMenu(options))
    chosen_mode = options[selection]
    for _ = 1:4
        print("\e[A\e[2K")  # Clear entire current line
    end 
    if options[selection] == "Quit"
        println("Exiting the radar simulator.\e[37m")
        exit(0)
    elseif options[selection] == "Visualize & Spawn"
        println("V&S mode activated!")
        println("Left-click to spawn targets, Right-Click to toggle omniscient mode.")
    end
    
    mode = options[selection]
    break
end

if mode == "Game"
    while(true)
        global mode
        options = ["Incoming!", "Spy Hunt", "TBC...", "Quit"]
        selection = request("Choose a game mode:", RadioMenu(options))
        chosen_mode = options[selection]
        for _ = 1:5
            print("\e[A\e[2K")  # Clear entire current line
        end 

        if options[selection] == "Quit"
            println("Exiting the radar simulator.\e[37m")
            exit(0)
        elseif options[selection] == "Incoming!"
            println("Incoming! challenge activated!")
            println("Left-click to shoot down enemies, Right-Click to toggle clockwise/counter-clockwise.")
        elseif options[selection] == "Spy Hunt"
            println("Spy Hunt challenge activated!")
        else
            print("Work in progress... ")
            continue
        end
        mode = options[selection]
        break
    end
end

if mode == "Visualize & Spawn"
    liveTargets = Vector{Targets.Target}(
                    [Targets.CartesianTarget(randn() * 100e3, randn() * 100e3, randn() * 5000, randn() * 5000; swerlingModel = :Swerling0, rcs = 1.0, status = :unknown) for _ in 1:20]  # Randomly generated targets
                ) # Live targets to be processed
elseif mode == "Incoming!"
    liveTargets = Vector{Targets.Target}(
                    [Targets.toCartesian(Targets.PolarTarget(rand() * 100e3 + 100e3, rand() * 360, -1000.0, 0.0; swerlingModel = :Swerling0, rcs = 1.0, status = :foe)) for _ in 1:20]  # Randomly generated targets
                ) # Live targets to be processed
elseif mode == "Spy Hunt"
    println("not yet supported...")
    quit(0)
end

GC.enable(true)
# Initialize system
# radar = Radar.initRadar()
radar = Radar.pulseRadar()

# targets = generateTargets()
PPI, plots, azimuth, beamPoints, targetDots, echo, compression, doppler, omniscience = RadarUI.setupRadarDisplay(radar, liveTargets, mode; initialAzimuth = radar.azimuth)

@async begin
    try 
        RadarUI.updateDisplay(azimuth, beamPoints, targetDots, echo, compression, doppler, radar, liveTargets, omniscience)
    catch e
        println("Error in updateDisplay: $e")
    end
end
# Show the UI

# Wow, this works!
display(GLMakie.Screen(), PPI)
display(GLMakie.Screen(), plots)


# Start live radar
Radar.runRadar(radar, liveTargets)

@async RadarUI.updateTargets(liveTargets, mode)

println("Press Enter to quit...")

readline()
println("Radar.jl quit.\e[37m")