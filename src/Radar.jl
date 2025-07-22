
module Radar
    using Distributions
    using Random
    using Observables
    using LinearAlgebra
    using DSP
    using ..Parameters
    using ..Constants
    using ..Targets

    mutable struct DataCube <: AbstractArray{ComplexF64, 3}
        data::Array{ComplexF64, 3}
        rangeBins::Int
        nPulses::Int
        nChannels::Int
    end

    function DataCube(rangeBins::Int, nPulses::Int, nChannels::Int)
        # Initialize all to zeros
        @assert rangeBins > 0 "Number of range bins must be positive"
        @assert nPulses > 0 "Number of pulses must be positive"
        @assert nChannels > 0 "Number of channels must be positive"
        return DataCube(zeros(rangeBins, nPulses, nChannels), rangeBins, nPulses, nChannels)
    end

    Base.size(dc::DataCube) = size(dc.data) # Delegates to the internal array's size

    Base.getindex(dc::DataCube, i::Int, j::Int, k::Int) = getindex(dc.data, i, j, k)
    Base.getindex(dc::DataCube, I::Vararg{Int,N}) where N = getindex(dc.data, I...) # For arbitrary number of indices

    Base.setindex!(dc::DataCube, v, i::Int, j::Int, k::Int) = setindex!(dc.data, v, i, j, k)
    Base.setindex!(dc::DataCube, v, I::Vararg{Int,N}) where N = setindex!(dc.data, v, I...) # For arbitrary number of indices

    # Implement methods for StridedArray compatibility needed for FFTW In-Place
    Base.parent(dc::DataCube) = dc.data
    Base.parentindices(dc::DataCube) = map(Base.Slice, axes(dc.data)) # Or (Base.Slice(1:size(A,1)), Base.Slice(1:size(A,2)), Base.Slice(1:size(A,3)))
    Base.strides(dc::DataCube) = strides(dc.data)

    clear(dc::DataCube) = fill!(dc.data, ComplexF64(0.0))

    mutable struct pulseRadar
        datacube::DataCube
        carrierFrequency::Float64
        bandwidth::Float64
        samplingFrequency::Float64
        chirp::Vector{ComplexF64}
        matchedFilter::Vector{ComplexF64}
        dutyCycle::Float64
        # azimuth::Observable{Float64}
        azimuth::Float64
        elevation::Float64
        azimuthStep::Float64
        angularSpeed::Float64
        clockwise::Bool
        # For charting output
        echo::Vector{ComplexF64}
        compressedEcho::Vector{ComplexF64}
        dopplerEcho::AbstractMatrix{ComplexF64}
        detections::Vector{Bool}
        detectionVals::Vector{ComplexF64}
    end

    """ A simple pulse radar system that generates a chirp signal with linear frequency modulation, 
        sampling frequency of 5 MHz, and a carrier frequency of 3 GHz.
        Every range bin is then 60 m apart, and assuming a Radar coverage range of 240 km,
        the number of range bins by default is 4000.

    """
    function pulseRadar()

        coverageRange = Parameters.radarCoverageRange
        carrierFreq = Parameters.radarCarrierFrequency
        BW = Parameters.radarChirpBW
        samplingFreq = Parameters.radarSamplingFrequency
        dutyCycle = Parameters.radarDutyCycle / 100.0  # Convert percentage to fraction
        clockwise = Parameters.clockwise

        rangeBins = Int(floor(coverageRange / (Constants.lightSpeed / (2 * samplingFreq))))  # Number of range bins based on coverage range and sampling frequency
        dc = DataCube(rangeBins, Parameters.nPulsesPerBurst, Parameters.nChannels)
        chirpWindowLength = floor(Int, dutyCycle * rangeBins)
        chirpWindow = range(0, chirpWindowLength / samplingFreq, length = chirpWindowLength)

        # chirp = 3 .* exp.(im * (2 * π * carrierFreq * chirpWindow .+ π * BW / chirpWindowLength .* chirpWindow .^ 2))
        chirp = 3 * exp.(im * π * BW / chirpWindow[end] * ( chirpWindow .- chirpWindow[end] / 2).^2)
        matchedFilter = conj(reverse(chirp)) ./ norm(chirp, 2)

        # Here is where we might apply SLL
        matchedFilter .*= DSP.kaiser(length(matchedFilter), Parameters.kaiserSLL)
        matchedFilter = vcat(zeros(ComplexF64, rangeBins - length(matchedFilter)), matchedFilter)

        pulseRadar(dc, carrierFreq, BW, samplingFreq, chirp, matchedFilter, dutyCycle, Parameters.radarInitialAzimuth,
         Parameters.radarInitialElevation, Parameters.radarAzimuthStep, Parameters.radarAngularVelocity, clockwise,
         zeros(ComplexF64, rangeBins), zeros(ComplexF64, rangeBins), zeros(ComplexF64, rangeBins, Parameters.nPulsesPerBurst * 3),
         falses(rangeBins), zeros(ComplexF64, rangeBins))
    end


    ### RADAR PROCESSING CHAIN
    module RadarProcessingChain
        using FFTW
        using Distributions
        using Plots
        using Random
        using DSP
        using ..Radar
        using ...Targets
        using ...SwerlingModels
        using ...Parameters
        using ...Constants

        function calculateTargetEcho(target::Targets.Target, rangeBins::Int, az::Float64, fs::Float64, fc::Float64, chirpedSignal::Vector{ComplexF64})       
            
            # println("AZ = $az, Target AZ = $(target.azimuth), Δ = $(abs(az - target.azimuth))")
            if abs(az - target.azimuth) >= 0.5
                return nothing  # No echo if target is outside the azimuth coverage
            end
            # if abs(az - target.azimuth) < 0.5
            #     println("Target! At azimuth: $az") # No echo if target is outside the azimuth coverage
            # end

            
            # Calculate the time delay based on the target's radial speed and distance
            delayTime = 2 * target.distance / Constants.lightSpeed
            delayRangeBins = round(Int, delayTime * fs)
            
            # Doppler frequency shift based on the target's radial speed
            t = (0:length(chirpedSignal)-1) / fs
            dopplerShift = 2 * target.radialSpeed / Constants.lightSpeed * fc
            dopplerPhase = exp.(im .* 2 .* π .* dopplerShift .* t)
            
            # Calculate attenuation
            Pt = Parameters.radarPower
            Pn = Constants.boltzmann * Parameters.radarTemperature * Parameters.radarChirpBW  # Noise power
            G = 10^(Parameters.radarAntennaGain / 10) # Given in dB
            λ = Constants.lightSpeed / fc

            attenuation = (Pt * G^2 * λ^2 * target.rcs) / ((4π)^3 * (target.distance^4 + 1e-6) * Pn)
            # @show attenuation
            # attenuation = 1e10 * SwerlingModels.swerlingAmplitude(target.rcs, target.swerlingModel) / (target.distance^2 + 1e-6)
            
            echo = chirpedSignal .* attenuation .* dopplerPhase
            paddedEcho = zeros(ComplexF64, rangeBins)
            paddedEcho[1:length(echo)] .= echo  
            # plt = Plots.plot(abs.(circshift(paddedEcho, delayRangeBins)), title="Echo", xlabel="Bin", ylabel="|Amplitude|")
            # display(plt)
            # Generate the echo signal with Doppler shift and noise
            return circshift(paddedEcho, delayRangeBins)
        end

        function calculateClutterEcho(rangeBins::Int, p::Int; strength=0.5, dopplerSpread=0.01)
            taper = exp.(-range(0, rangeBins-1) ./ (rangeBins * 4))

            clutter = strength * randn(ComplexF64, rangeBins) .* taper
            phaseVel = 2π * dopplerSpread .* rand(rangeBins)
            clutter .*= exp.(im .* phaseVel .* p)

            return clutter
        end

        # const calculateClutterEcho = let
        #     clutter_buf = nothing
        #     taper = nothing
        #     phase_buf = nothing
        #     prev_N = 0  # track the size used to know when to reallocate

        #     function (rangeBins::Int, p::Int; strength=0.5, dopplerSpread=0.01)
        #         if clutter_buf === nothing || rangeBins != prev_N
        #             # (Re)initialize on first use or when size changes
        #             clutter_buf = Vector{ComplexF64}(undef, rangeBins)
        #             taper = exp.(-range(0, rangeBins - 1) ./ (rangeBins * 4))
        #             phase_buf = Vector{Float64}(undef, rangeBins)
        #             prev_N = rangeBins
        #         end

        #         randn!(clutter_buf)
        #         @. clutter_buf = strength * clutter_buf * taper
        #         rand!(phase_buf)
        #         @. phase_buf = 2π * dopplerSpread * phase_buf
        #         @. clutter_buf *= cis(phase_buf * p)

        #         return clutter_buf
        #     end
        # end

        function calculateTerrainMask(mask::Vector{ComplexF64}, az::Float64, samples::Int, fs::Float64)

            for element in Parameters.terrainOccludedVisibility
                if !(az in element[1])
                    continue
                end
                
                # If we are in the azimuth range of the terrain occlusion
                kmRange = element[2]
                visibility = element[3]
                
                binsize = Constants.lightSpeed / (2 * fs)  # Range bin size in meters
                
                # Convert km range to range bins
                startRange = first(kmRange) * 1000
                endRange   = last(kmRange) * 1000

                startBin = max(0, round(Int, startRange / binsize))
                endBin = min(samples, round(Int, endRange / binsize))

                mask[startBin:endBin] .= visibility  # Set visibility in the range bins
            end

            return mask
        end

        """    Function to process the radar echo signal using pulse compression.
        It applies a matched filter to the echo signal to enhance the target detection.
        Parameters:
        - r: pulseRadar object containing radar parameters and chirp signal.
        - echo: Vector of ComplexF64 representing the received echo signal.
        Returns:
        - Compressed echo signal as a Vector of ComplexF64.
        """
        function pulseCompression(echo::Vector{ComplexF64}, matchedFilter::Vector{ComplexF64})
            # Pad matched filter
            # @show size(echo), size(paddedFilter) 
            return ifft(fft(echo) .* fft(matchedFilter))
            
            # # Apply matched filter to the echo signalq
            # res = conv(echo, matchedFilter)

            # # Normalize the compressed echo
            # res ./= maximum(abs.(res))

            # # Center-trim
            # N = length(echo)
            # halfLength = length(matchedFilter) ÷ 2
            # res[halfLength:halfLength + N - 1]
        end

        function dopplerProcessing!(compressedEchoes::AbstractMatrix{ComplexF64}; window::Vector{Float64} = nothing)
            if window !== nothing
                compressedEchoes .*= window'
            end
            
            # FFT along the columns
            fft!(compressedEchoes, 2)  # FFT to process the echo in frequency domain
            compressedEchoes .= conj.(compressedEchoes) ./ sqrt(Parameters.nPulsesPerBurst)  # Renormalize after gain

            # Out-of-place circshift along rows (1st dim)
            # shift_amount = size(compressedEchoes, 1) ÷ 2
            # compressedEchoes .= circshift(compressedEchoes, (shift_amount, 0))
        end

        """
        Performs Cell-Averaging CFAR on a 2D matrix (e.g., range-Doppler map).

        # Arguments:
        - data: 2D matrix (ComplexF64 or Real)
        - guard: number of guard cells (per side)
        - ref: number of reference cells (per side)
        - K: scaling factor (e.g., 1.4 for 13 reference cells)

        # Returns:
        - detections: BitMatrix of size(data), where true = detection
        """
        function caCFAR(data::AbstractMatrix; guard::Int=1, ref::Int=50, K::Float64=2.4)
            rows, cols = size(data)
            detections = falses(rows, cols)

            # Work on magnitude
            mag = abs.(data)

            averageThreshold = K * mean(mag)

            detections .= mag .> averageThreshold
            
            # totalWindow = guard + ref
            # for r in (1 + totalWindow):(rows - totalWindow)
            #     for c in 1:cols
            #         # Define reference window
            #         rStart, rEnd = r - totalWindow, r + totalWindow

            #         # Exclude guard + CUT window
            #         rGuardStart, rGuardEnd = r - guard, r + guard

            #         # Get reference cells
            #         window = mag[rStart:rEnd, c]
            #         window[rGuardStart - rStart + 1 : rGuardEnd - rStart + 1] .= NaN

            #         # Compute threshold
            #         noise = mean(skipmissing(window))
            #         threshold = K * noise
                    

            #         if mag[r, c] > threshold
            #             detections[r, c] = true
            #         end
            #     end
            # end

            return detections
        end

        function processingRoutine(r::Radar.pulseRadar, targets::Vector{Targets.Target}, lastUpdateTime::Float64)

            # r.datacube.data = zeros(ComplexF64, r.datacube.rangeBins, r.datacube.nPulses, r.datacube.nPulses) # Clear the data cube before processing
            Radar.clear(r.datacube)
            for pulse in 1:r.datacube.nPulses

                # Refresh azimuth
                # Time elapsed update
                currentTime = time()
                timeElapsed = currentTime - lastUpdateTime
                lastUpdateTime = currentTime
                
                # Azimuth spanned update
                deltaAz = timeElapsed * r.angularSpeed
                deltaAz = r.clockwise ? deltaAz : -deltaAz
                r.azimuth = mod(r.azimuth + deltaAz, 360)
                
                # Simulate listening for echoes  
                # 2 ms == ~300km coverage range
                # Radars can actually be smart and listen for residual echoes in following pulses.
                # However, for now, this is beyond the scope of this simple tool. 2 ms is a dense enough
                # grid to scan for targets, i.e. ~ 0.06 deg at 12 seconds per revolution
                # TODO: Listen for residual echoes, cut listening time by Npulses
                sleep(Parameters.radarCoverageRange / (2 * Constants.lightSpeed))
                
                for channel in 1:r.datacube.nChannels

                    r.datacube[:, pulse, channel] .= craftComplexEcho!(r, r.datacube[:, pulse, channel], targets, pulse)
                    r.echo .= r.datacube[:, pulse, channel]  # Update echo for charting
                    
                    # FAST-TIME: Apply pulse compression
                    r.datacube[:, pulse, channel] .= pulseCompression(r.datacube[:, pulse, channel], r.matchedFilter)
                    r.compressedEcho .=  r.datacube[:, pulse, channel]  # Update compressed echo for charting
                    
                end
            end
            
            # r.echo .= r.datacube[:, 1, 1]  # Update echo for charting
            # # FAST-TIME: Apply pulse compression
            # compressedEcho = pulseCompression!(r.datacube[], r.matchedFilter)

            # SLOW-TIME: Apply Doppler processing
            r.dopplerEcho .= hcat(r.datacube.data[:, :, 1], zeros(r.datacube.rangeBins, r.datacube.nPulses * 2))  # Matrix of Doppler data
            dopplerProcessing!(r.dopplerEcho; window=DSP.kaiser(size(r.dopplerEcho, 2), Parameters.kaiserSLL))
            r.datacube[:, :, 1] .= dopplerProcessing!(r.datacube[:, :, 1]; window=DSP.kaiser(r.datacube.nPulses, Parameters.kaiserSLL))
            # if maximum(abs.(r.datacube[:, :, 1])) > 30
            #     println("Max Doppler data: $(maximum(abs.(r.datacube[:, :, 1])))")
            # end
            r.detectionVals .= vec(maximum(abs.(r.datacube[:, :, 1]), dims=2))
            r.detections .= vec(any(caCFAR(r.datacube[:, :, 1]), dims=2))

            # TODO: Implement monopulse with added channels

            return lastUpdateTime
        end

        function craftComplexEcho!(r::Radar.pulseRadar, out::Vector{ComplexF64}, targets::Vector{Targets.Target}, pulseNum::Int)
            
            # Preallocate auxiliary vectors
            mask = ones(ComplexF64, r.datacube.rangeBins)

            for target in targets
                # Calculate echo for each target
                targetEcho = calculateTargetEcho(Targets.toPolar(target), r.datacube.rangeBins, r.azimuth, r.samplingFrequency, r.carrierFrequency, r.chirp) 
                if targetEcho !== nothing
                    out .+= targetEcho  # Add the echo to the output
                end
            end
            
            # Add clutter 
            out .+= calculateClutterEcho(r.datacube.rangeBins, pulseNum)
            # Add terrain mask
            # plt = Plots.plot(abs.(out), title="Echo2", xlabel="Bin", ylabel="|Amplitude|")
            # display(plt)  # Or `gui()` to force window display
            out .*= calculateTerrainMask(mask, r.azimuth, r.datacube.rangeBins, r.samplingFrequency)

            # Add noise to the echo
            # out .+= randn(ComplexF64, r.datacube.rangeBins) * 0.1  # Add some noise
            noisyBackground = rand(Normal(0, 1), r.datacube.rangeBins) # Simulated noise in the background
            out .+= noisyBackground  # Add noise to the echo
        end
    end

       
    function runRadar(
        radar::Radar.pulseRadar,
        targets::Vector{Targets.Target}
    )
        @async begin
            try
                lastUpdateTime = time()
                while true

                    # Update detection routine
                    lastUpdateTime = RadarProcessingChain.processingRoutine(radar, targets, lastUpdateTime)

                    # This is too frequent, we ideally want to only reliqnquish control when we are not updating the UI
                    # yield() 
                    # Workaround: sleep a fine-tuned amount
                    sleep(1e-2)
                end
            catch e
                @error "Radar processing failed!" exception=(e, catch_backtrace())
            end
        end
    end

end 

