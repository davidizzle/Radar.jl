module RadarUI
    using GLMakie
    using Plots
    using GeometryBasics
    using Observables
    using Dates
    using ..Radar
    using ..Targets
    # using .Radar

    struct TargetDetection
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

    const TRAIL_LENGTH = 30
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

    function setupRadarDisplay(radar::Radar.pulseRadar, targets; maxRange=5.0, initialAzimuth=0.0)
        fig = Figure(size=(1920, 1080), backgroundcolor=:black)

        # Echo Section
        ax_echo = Axis(fig[1, 1],
        xlabel="Range (Km)", ylabel="Amplitude",
        title="Echo", limits=(0, radar.datacube.rangeBins, 1e-2, 1e2), yscale=log10
        )

        # echo_data = Observable(radar.datacube.data[:, 1, 1]) 
        echo_data = Observable(radar.echo) 
        # Create a lifted observable of absolute values
        echo = lift(x -> abs.(x), echo_data)
        lines!(ax_echo, 1:length(echo[]), echo, color=:blue)
        
        # Pulse Compression Section
        ax_pc = Axis(fig[1, 2],
        xlabel="Range (Km)", ylabel="Amplitude",
        title="Echo", limits=(0, radar.datacube.rangeBins, 1e-2, 1e2), yscale=log10
        )
        pc_data = Observable(radar.compressedEcho) # Example: 500 points
        pc = lift(x -> abs.(x), pc_data)
        lines!(ax_pc, 1:length(pc[]), pc, color=:blue)

        # Doppler Processing Section
        ax_dopp = Axis(fig[1, 3], xlabel="Range (Km)", ylabel="Velocity",
        title="Echo")
        # dopp = Observable(abs.(radar.datacube.data[:, :, 1])) 
        dopp = Observable(radar.dopplerEcho) 
        abs_dopp = lift(x -> abs.(x), dopp)
        # @show typeof(radar.datacube.data[:, :, 1])   # should be Matrix{ComplexF64}
        # @show size(radar.datacube.data[:, :, 1])  
        GLMakie.heatmap!(ax_dopp, abs_dopp, colormap=:viridis, colorrange=(0, 50))


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
        

        detectionDots = Observable(TargetDetection[]) 
        detectionDots[] = [TargetDetection(t, time(), RGBAf(1.0, 0.0, 0.0, 1.0), RGBAf(1.0, 1.0, 0.0, 1.0), 1.0, 5.0) for t in targets]
        positions = @lift [td.position for td in $detectionDots]
        colors = @lift begin
            now = time()
            [RGBAf(td.outerColor.r, td.outerColor.g, td.outerColor.b,
                max(0, td.initial_alpha * (1.0 - (now - td.detection_time)/td.lifetime)))
            for td in $detectionDots]
        end
        GLMakie.scatter!(ax, positions, color=colors, markersize=10, strokewidth=1)

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
        #  markersize=15, strokewidth=1, strokecolor=:white)


        # Draw concentric circles
        step = maxRange / 5
        for r in step:step:2*maxRange
            lines!(ax, [polarToCartesian(θ, r) for θ in 0:1:360], color=:gray, linewidth=0.5)
        end

        # General background section
        # Draw radial spokes
        for θ in 0:30:330
            lines!(ax, [Point2f(0, 0), polarToCartesian(θ, maxRange)], color=:gray, linewidth=0.5)
        end

        for i in 1:TRAIL_LENGTH
            alpha = 0.5 - (i - 1) / TRAIL_LENGTH / 2
            lineplot = lift(beamHistory) do beams
                i ≤ length(beams) ? beams[end - i + 1] : [Point2f(0, 0), Point2f(0, 0)]
            end
            lines!(ax, lineplot, color=(:green, alpha), linewidth=2)
        end

        azimuthDeg = Observable(initialAzimuth)
        beam_start = Point2f(0, 0)
        beam_end = polarToCartesian(initialAzimuth, maxRange)
        beamPoints = Observable([beam_start, beam_end])

        lines!(ax, beamPoints, color=:green, linewidth=2)

        return fig, fig2, azimuthDeg, beamPoints, detectionDots, echo_data, pc_data, dopp
    end

    function updateDisplay(az::Observable, beamPoints::Observable, detectionDots::Observable, echo::Observable, compression::Observable, doppler::Observable, radar::Radar.pulseRadar, targets::Vector{Targets.Target})

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
            # updateBeamTrail!(az[], range)

            doppler[] .= radar.dopplerEcho # Update Doppler data
            notify(doppler)
            # doppler[] .= abs.(radar.datacube.data[:, :, 1]) # Update Doppler data

            detectionDots[] = [TargetDetection(t, time(), RGBAf(1.0, 0.0, 0.0, 1.0), RGBAf(1.0, 1.0, 0.0, 1.0), 1.0, 5.0) for t in targets]
            
            sleep(1e-2)
        end
    end

    function updateTargets(targets::Vector{Targets.Target})
        try
            while true
                for t in targets
                        Targets.updatePosition(t)
                end
                sleep(1e-2)  # Sleep to avoid busy-waiting
            end
        catch e
            @error "Error in updateTargetDots: $e"
        end
    end

    export setupRadarDisplay, runRadar
end