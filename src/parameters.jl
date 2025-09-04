module Parameters

    const radarCoverageRange = 200e3  # Coverage range in meters
    const radarPower = 1e6 # Radar power in Watts, 1 MW
    const radarTemperature = 290  # Noise temperature in Kelvin
    const radarAntennaGain = 36 # Antenna gain in dB
    const radarAngularVelocity = 30  # deg per s
    const radarAzimuthStep = 0.5  # degrees
    const radarInitialElevation = 0  # degrees
    const radarInitialAzimuth = 0  # degrees
    const nPulsesPerBurst = 10  # Number of pulses per burst
    const nChannels = 1  # Number of channels in the radar system
    const radarCarrierFrequency = 9e9  # Carrier frequency in Hz
    const radarSamplingFrequency = 5e6  # Sampling frequency in Hz
    const radarChirpBW = 2e6  # Chirp bandwidth in Hz
    const radarDutyCycle = 10  # Duty cycle in percentage
    const clockwise = true
    # Terrain occluded visibility: (azStart:azEnd [deg], rangeStart:rangeEnd [km], visibility)
    const terrainOccludedVisibility = [(31:67, 60:90, 0.5), (142:158, 100:200, 0), (253:289, 80:100, 0.2), (290:326, 200:400, 0.6)]  # Example visibility ranges for terrain occlusion

    # SLL 
    const kaiserSLL = 5.47 # 40 dB
    # const kaiserSLL = 6.76 # 50 dB

end