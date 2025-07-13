
module Radar
    using Distributions
    using Random
    include("constants.jl")
    using .Constants

    struct DataCube <: AbstractArray{ComplexF64, 3}
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

    struct pulseRadar
        datacube::DataCube
        carrierFrequency::Float64
        bandwidth::Float64
        samplingFrequency::Float64
        chirp::Vector{ComplexF64}
        matchedFilter::Vector{ComplexF64}
        dutyCycle::Float64
        azimuth::Float64
        elevation::Float64
        azimuthStep::Float64
        radialSpeed::Float64
    end

    """ A simple pulse radar system that generates a chirp signal with linear frequency modulation, 
        sampling frequency of 5 MHz, and a carrier frequency of 3 GHz.
        Every range bin is then 60 m apart, and assuming a Radar coverage range of 240 km,
        the number of range bins by default is 4000.

    """
    function pulseRadar(rangeBins::Int = 4000, carrierFreq::Float64 = 3e9, BW::Float64 = 2e6, samplingFreq::Float64 = 5e6, dutyCycle::Float64 = 0.1)
        
        dc = DataCube(rangeBins, Constants.nPulsesPerBurst, Constants.nChannels)

        chirpWindowLength = dutyCycle * rangeBins
        chirpWindow = range(0, chirpWindowLength / samplingFreq, length = chirpWindowLength)

        chirp = 3 .* exp.(im * (2 * π * carrierFreq * chirpWindow .+ π * BW / chirpWindowLength .* chirpWindow .^ 2))
        matchedFilter = conj(reverse(chirp)) ./ norm(chirp, 2)

        pulseRadar(dc, carrierFreq, BW, samplingFreq, chirp, matchedFilter, dutyCycle, Constants.radarInitialAzimuth,
         Constants.radarInitialElevation, Constants.radarAzimuthStep, Constants.radarAngularVelocity)
    end

end 

