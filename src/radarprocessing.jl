

module RadarProcessingChain
    using FFTW
    include("radar.jl")
    using .Radar

    function calculateTargetEcho(r::pulseRadar, target::Target)

        noisyBackground = rand(Normal(0, 1), r.rangeBins) # Simulated noise in the background
        
        if abs(r.azimuth - target.azimuth) >= 0.5
            return noisyBackground  # No echo if target is outside the azimuth coverage
        end

        # Calculate the time delay based on the target's radial speed and distance
        delayTime = 2 * target.distance / Constants.lightSpeed
        delayRangeBins = round(Int, delayTime * r.samplingFrequency)

        # Doppler frequency shift based on the target's radial speed
        dopplerShift = 2 * target.radialSpeed / Constants.lightSpeed * r.carrierFrequency

        # Generate the echo signal with Doppler shift and noise
        echoSignal = r.chirp .* exp(im * 2 * π * dopplerShift * (0:length(r.chirp)-1) / r.samplingFrequency)
        echoSignal = circshift(echoSignal, delayRangeBins) + noisyBackground
    end

    """    Function to process the radar echo signal using pulse compression.
    It applies a matched filter to the echo signal to enhance the target detection.
    Parameters:
    - r: pulseRadar object containing radar parameters and chirp signal.
    - echo: Vector of ComplexF64 representing the received echo signal.
    Returns:
    - Compressed echo signal as a Vector of ComplexF64.
    """
    function pulseCompression(r::pulseRadar, echo::Vector{ComplexF64})
        # Apply matched filter to the echo signal
        compressedEcho = conv(echo, r.matchedFilter, mode = :same)

        # Normalize the compressed echo
        normalizedEcho = compressedEcho ./ maximum(abs.(compressedEcho))

        return normalizedEcho
    end

    function dopplerProcessing!(r::pulseRadar, compressedEchoes::Matrix{ComplexF64}, target::Target)
        # FFT along the columns
        compressedEchoes = fft(compressedEchoes, 2)  # FFT to process the echo in frequency domain
        compressedEchoes = adjoint(compressedEchoes) / sqrt(size(compressedEchoes)[2])  # Renormalize after gain
    end

    function processingRoutine(r::pulseRadar, targets::Vector{Target})

        for target in targets
            calculateTargetEcho(r, target)

            # Simulate the echo signal for the target
            echo = simulateEcho(r, target)

        end

        # Apply pulse compression
        compressedEcho = pulseCompression(r, echo)

        # Apply Doppler processing
        dopplerProcessedEcho = dopplerProcessing(r, compressedEcho, target)

        return dopplerProcessedEcho
    end

end