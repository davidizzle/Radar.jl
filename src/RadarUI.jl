module RadarUI
    using GLMakie
    using Plots
    using GeometryBasics
    using Observables
    using Dates
    using FileIO
    using ..Radar
    using ..Targets
    using ..Parameters
    using ..Constants
    # using .Radar

    mutable struct TargetDetection
        target::Targets.Target
        position::Point2f
        detection_time::Float64 # Using `time()` for this
        outerColor::RGBAf
        innerColor::RGBAf
        initial_alpha::Float64  # Store the initial transparency
        lifetime::Float64       # How long the dot should be visible (in seconds)
    end

    function TargetDetection(target::Targets.Target, detection_time::Float64, outerColor::RGBAf, innerColor::RGBAf, initial_alpha::Float64, lifetime::Float64)
        tempTarget = Targets.toCartesian(target)
        position = Point2f(tempTarget.x/1e5, tempTarget.y/1e5)
        return TargetDetection(target, position, detection_time, outerColor, innerColor, initial_alpha, lifetime)
    end

    const TRAIL_LENGTH = 5
    beamHistory = Observable(Vector{Vector{Point2f}}())

    function polarToCartesian(θ_deg, r)
        θ_rad = deg2rad(θ_deg)
        # println(θ_rad)
        return Point2f(r * cos(θ_rad), r * sin(θ_rad))
    end
    
    function updateBeamTrail!(currentAzimuth::Float64, maxRange::Float64)
        newBeam = [Point2f(0, 0), polarToCartesian(currentAzimuth, maxRange)]

        if length(beamHistory[]) >= TRAIL_LENGTH
            popfirst!(beamHistory[])
        end

        push!(beamHistory[], newBeam)
        notify(beamHistory)
    end

    function setupRadarDisplay(radar::Radar.pulseRadar, targets, mode::String; initialAzimuth=0.0)
        set_theme!(theme_black())
        fig = Figure(size=(1920, 1080), backgroundcolor=:black)
        maxRange = Parameters.radarCoverageRange / 1e5  # Convert to kilometers
        # Echo Section
        ax_echo = Axis(fig[1, 1],
        xlabel="Range (Km)", ylabel="Amplitude",
        title="Raw Echo", limits=(0, radar.datacube.rangeBins, 1e-2, 1e2), yscale=log10
        )
        

        # echo_data = Observable(radar.datacube.data[:, 1, 1]) 
        echo_data = Observable(radar.echo) 
        # Create a lifted observable of absolute values
        echo = lift(x -> abs.(x), echo_data)
        lines!(ax_echo, 1:length(echo[]), echo, color=:blue)
        
        # Pulse Compression Section
        ax_pc = Axis(fig[1, 2],
        xlabel="Range (Km)", ylabel="Amplitude",
        title="Compressed Echo", limits=(0, radar.datacube.rangeBins, 1e-2, 1e2), yscale=log10
        )
        pc_data = Observable(radar.compressedEcho) # Example: 500 points
        pc = lift(x -> abs.(x), pc_data)
        lines!(ax_pc, 1:length(pc[]), pc, color=:blue)

        # Doppler Processing Section
        ax_dopp = Axis(fig[2, 1], xlabel="Range (Km)", ylabel="Radial Velocity",
        title="Doppler Processed")
        ax_dopp_3d = Axis3(fig[2, 2], xlabel="Range (Km)", ylabel="Radial Velocity",
        title="Doppler Processed (3D)", limits=(0, radar.datacube.rangeBins, 0, radar.datacube.nPulses, 1e-2, 1e2))
        # dopp = Observable(abs.(radar.datacube.data[:, :, 1])) 
        shift = round(Int, size(radar.dopplerEcho, 2) / 2)
        dopp = Observable(radar.dopplerEcho) 
        # abs_dopp = lift(x -> circshift!(x, abs.(x), (0, shift)), dopp)
        abs_dopp = lift(x -> circshift(abs.(x), (0, shift)), dopp)
        # abs_dopp = lift(x -> 10*log10.(abs.(x)), dopp)
        # @show typeof(radar.datacube.data[:, :, 1])   # should be Matrix{ComplexF64}
        # @show size(radar.datacube.data[:, :, 1])  
        λ = Constants.lightSpeed / Parameters.radarCarrierFrequency
        # 1
        PRF = Parameters.radarCoverageRange * ( 1 + radar.dutyCycle ) * 2 
        GLMakie.heatmap!(ax_dopp, abs_dopp, colormap=:viridis, colorrange=(0, 50), colorscale=log10)

        GLMakie.surface!(ax_dopp_3d, abs_dopp, colormap=:viridis, colorrange=(0, 50), colorscale=log10)
 

        # Radar PPI Section

        # Vector of observable dots
        fig2 = Figure(size=(1080, 1080), backgroundcolor=:black)
        
        ax = Axis(fig2[1,1],
            backgroundcolor=:black,
            xgridcolor=:gray,
            ygridcolor=:gray,
            xlabel="X", ylabel="Y",
            limits = (-maxRange, maxRange, -maxRange, maxRange),
            )
        interactions(ax)
        deregister_interaction!(ax, :rectanglezoom)
        deregister_interaction!(ax, :dragpan)
        deregister_interaction!(ax, :scrollzoom)
        deregister_interaction!(ax, :limitreset)

        color_map = Dict(
            :friend  => RGBAf(0.2, 0.4, 0.8, 0.0),
            :foe     => RGBAf(0.8, 0.2, 0.2, 0.0),
            :unknown => RGBAf(1.0, 1.0, 0.0, 0.0)
        )


        detectionDots = Observable(TargetDetection[]) 
        detectionDots[] = [TargetDetection(t, time(), RGBAf(0.85, 0.88, 0.92, 1.0), color_map[t.status], 1.0, 5.0) for t in targets]
        positions = @lift [td.position for td in $detectionDots]
        inner_colors = @lift [td.innerColor for td in $detectionDots]
        outer_colors = @lift [td.outerColor for td in $detectionDots]
        # colors = @lift begin
        #     now = time()
        #     [RGBAf(td.outerColor.r, td.outerColor.g, td.outerColor.b,
        #         max(0, td.initial_alpha * (1.0 - (now - td.detection_time)/td.lifetime)))
        #     for td in $detectionDots]
        # end
        GLMakie.scatter!(ax, positions, color=inner_colors, strokecolor=outer_colors, markersize=10, strokewidth=2)

        scene = ax.scene

        omniscience = Observable(false)
        
        if (mode == "Visualize & Spawn")
            on(scene.events.mousebutton) do event
                if event.button == Mouse.left && event.action == Mouse.press
                    pos = GLMakie.mouseposition(scene)
                    if !isnothing(pos)
                        x, y = Tuple(pos)
                        # Scale to match target coordinates
                        # @show pos
                        x *= 1e5
                        y *= 1e5
                        new_target = Targets.CartesianTarget(x, y, randn() * 500, randn() * 500)
                        push!(targets, new_target)
                    end
                end
            end
            register_interaction!(ax, :rightclick) do event::MouseEvent, axis
                if event.type == MouseEventTypes.rightclick

                    omniscience[] = !(omniscience[])

                end
            end
        elseif mode == "Incoming!"
            register_interaction!(ax, :leftclick) do event::MouseEvent, axis
                if event.type == MouseEventTypes.leftclick
                    # println(event.data)
                    filter!(t -> abs(Targets.toCartesian(t).x/1e5 - event.data[1]) >= 3e-2 || abs(Targets.toCartesian(t).y/1e5 - event.data[2]) >= 3e-2, targets)
                end
            end
            register_interaction!(ax, :rightclick) do event::MouseEvent, axis
                if event.type == MouseEventTypes.rightclick
                    radar.clockwise = !(radar.clockwise)
                end
            end
        end
        # register_interaction!(ax, :my_interaction) do event::MouseEvent, GLMakie.Axis
        #     if event.type == MouseEventTypes.leftclick
        #         println("Left click at: ", event.data)
        #         x, y = Tuple(pos)
        #         x *= 1e5
        #         y *= 1e5
        #         new_target = Targets.CartesianTarget(x, y, randn() * 5000, randn() * 5000)
        #         push!(targets, new_target)
        #     end
        # end
        # Vector of detected dots
        # tempPos = [Point2f(Targets.toCartesian(t).x/1e5, Targets.toCartesian(t).y/1e5)  for t in targets]

        # dot_position_lift = @lift [td.position for td in $detectionDots]
        # dot_color_lift = @lift begin
        #     current_time = time()
        #     colors = Vector{RGBAf}()
        #     for td in $detectionDots
        #         time_since_detection = current_time - td.detection_time
        #         remaining_lifetime_ratio = 1.0 - (time_since_detection / td.lifetime)
        #         alpha = max(0.0, td.initial_alpha * remaining_lifetime_ratio) # Fade out
        #         push!(colors, RGBAf(1.0, 1.0, 0.0, alpha)) # Yellow dot, fading alpha
        #     end
        #     colors
        # end
        
        # GLMakie.scatter!(ax, tempPos, 
        # # color=dot_color_lift,
        # # strokecolor=dot_color_lift,
        #  markersize=15, strokewidth=1, strokecolor=:white)


        # Draw concentric circles
        step = maxRange / 5
        for r in step:step:2*maxRange
            lines!(ax, [polarToCartesian(θ, r) for θ in 0:1:    360], color=:gray, linewidth=0.5)
        end

        # General background section
        # Draw radial spokes
        for θ in 0:30:330
            lines!(ax, [Point2f(0, 0), polarToCartesian(θ, maxRange)], color=:gray, linewidth=0.5)
        end

        # Background
        image_idx = rand(1:7)
        my_image_data = FileIO.load("assets/bkgrnd$(image_idx).png") 
        my_image_data = rotr90(my_image_data)
        image!(ax, -maxRange*1.1..maxRange*1.1, -maxRange*1.1..maxRange*1.1, my_image_data, alpha=0.6)

        for i in 1:TRAIL_LENGTH
            alpha = 0.5 - (i - 1) / TRAIL_LENGTH / 2
            lineplot = lift(beamHistory) do beams
                i ≤ length(beams) ? beams[end - i + 1] : [Point2f(0, 0), Point2f(0, 0)]
            end
            lines!(ax, lineplot, color=(:green, alpha), linewidth=8)
        end

        azimuthDeg = Observable(initialAzimuth)
        beam_start = Point2f(0, 0)
        beam_end = polarToCartesian(initialAzimuth, maxRange)
        beamPoints = Observable([beam_start, beam_end])

        lines!(ax, beamPoints, color=:green, linewidth=5)

        return fig, fig2, azimuthDeg, beamPoints, detectionDots, echo_data, pc_data, dopp, omniscience
    end

    function updateDisplay(az::Observable, beamPoints::Observable, detectionDots::Observable, echo::Observable, compression::Observable, doppler::Observable, radar::Radar.pulseRadar, targets::Vector{Targets.Target}, omniscience::Observable{Bool})

        range = Float64(sqrt(beamPoints[][2][1]^2 + beamPoints[][2][2]^2))
        while true
            # echo[] = radar.datacube.data[:, 1, 1] # Update echo data
            echo[] = radar.echo # Update echo data
            compression[] = radar.compressedEcho # Update compressed echo data
            # pc_data[] = radar.datacube.data[:, 1, 1] # Update
            az[] = radar.azimuth
            notify(az)
            
            # Update beam endpoint
            newEnd = polarToCartesian(az[], range)
            # beamPoints[][2] = newEnd
            # println("Updating display 3...")
            beamPoints[] = [Point2f(0,0), newEnd]
            # beamPoints[][2] =  newEnd
            notify(beamPoints)
            updateBeamTrail!(az[], range)

            doppler[] .= radar.dopplerEcho # Update Doppler data
            notify(doppler)
            # doppler[] .= abs.(radar.datacube.data[:, :, 1]) # Update Doppler data

            # Easiest thing to do is for targets to all check if they have been spotted
            # if so update detection time, and compute alpha of inner from here. 
            # Alpha of outer should be configurable from menu (God mode)
            # detectionDots[] = [TargetDetection(t, time(), RGBAf(1.0, 1.0, 1.0, 1.0), RGBAf(1.0, 0.0, 0.0, 0.0), 1.0, 5.0) for t in targets]
            # detectionDots[] = [TargetDetection(t, time(), RGBAf(1.0, 1.0, 1.0, 1.0), RGBAf(1.0, 0.0, 0.0, 0.0), 1.0, 5.0) for t in targets]

            
            syncDetections!(detectionDots, targets, omniscience)
            for dot in detectionDots[]
                dummy = Targets.toPolar(dot.target)

                # Could have been busted
                if abs(dummy.azimuth - radar.azimuth) < 0.4
                    rangeBin = round(Int, dummy.distance * 2 / Constants.lightSpeed * Parameters.radarSamplingFrequency)
                    dot.detection_time = time()
                    dot.innerColor = RGBAf(dot.innerColor.r, dot.innerColor.g, dot.innerColor.b, 1.0)
                    if rangeBin < radar.datacube.rangeBins && (any(radar.detections[rangeBin-5:rangeBin+5]))
                        # Busted!
                    end
                end
            end

            notify(detectionDots)
            sleep(1e-2)
        end
    end

    function syncDetections!(detectionDots::Observable{Vector{TargetDetection}}, targets::Vector{Targets.Target}, omniscience::Observable{Bool})
        existing = detectionDots[]
        new_dots = Vector{TargetDetection}()
        outer_alpha = omniscience[] ? 1.0 : 0.0

        for tgt in targets
            # Try to find a matching detection by some identity (use objectid or a unique ID field)
            found = findfirst(dot -> dot.target === tgt, existing)

            if isnothing(found)
                # New target -> create a new detection
                if tgt.status == :friend
                    innerC = RGBAf(0.2, 0.4, 0.8, outer_alpha)
                elseif tgt.status == :foe
                    innerC = RGBAf(0.8, 0.2, 0.2, outer_alpha)
                else
                    innerC = RGBAf(1.0, 1.0, 0.0, outer_alpha)
                end

                cart = Targets.toCartesian(tgt)
                new_dot = TargetDetection(
                    tgt,
                    Point2f(cart.x / 1e5, cart.y / 1e5),
                    time(),
                    RGBAf(0.85, 0.88, 0.92, 1.0),
                    innerC,
                    1.0,
                    5.0,
                )
                push!(new_dots, new_dot)
            else
                # Existing target -> update position, keep color/alpha/time
                dot = existing[found]
                cart = Targets.toCartesian(tgt)
                time_elapsed_from_detection = time() - dot.detection_time
                new_alpha = max(0.0, 1.0 - time_elapsed_from_detection / 1.0)
                dot.innerColor = RGBAf(dot.innerColor.r, dot.innerColor.g, dot.innerColor.b, new_alpha)
                dot.outerColor = RGBAf(dot.outerColor.r, dot.outerColor.g, dot.outerColor.b, outer_alpha)
                dot.position = Point2f(cart.x / 1e5, cart.y / 1e5)
                push!(new_dots, dot)
            end
        end

        # Replace with the updated detection list
        detectionDots[] = new_dots
    end

    function updateTargets(targets::Vector{Targets.Target}, mode::String)
        try
            prevTime = time()
            while true
                for t in targets
                        Targets.updatePosition(t)

                        if mode == "Incoming!"
                            if Targets.toPolar(t).distance <= 0.05
                                println("You got blown up!")
                                sleep(1)
                                quit(0)
                            end
                        end
                end
                # Remove targets that are out of range
                filter!(t -> Targets.toPolar(t).distance <= Parameters.radarCoverageRange || Targets.toPolar(t).radialSpeed > 0, targets)

                if mode == "Incoming!"
                    # Generate a foe 50% prob every second
                    if time() > prevTime + 1 && rand() < 0.5
                        push!(targets, Targets.toCartesian(Targets.PolarTarget(rand() * 100e3 + 100e3, rand() * 360, -1000.0, 0.0; swerlingModel = :Swerling0, rcs = 1.0, status = :foe)))
                        prevTime = time()
                    end
                elseif mode == "Spy Hunt"

                end

                sleep(1e-2)  # Sleep to avoid busy-waiting
            end
        catch e
            @error "Error in updateTargets: $e"
        end
    end

    export setupRadarDisplay, runRadar
end