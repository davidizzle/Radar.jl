

module Targets

    struct Target
        distance::Float64  # Distance to the target in meters
        radialSpeed::Float64  # Radial speed of the target in m/s
        azimuth::Float64  # Azimuth angle of the target in radians
        swerlingModel::Symbol  # Swerling model type
        rcs::Float64  # Radar Cross Section (RCS) in square meters
    end

    function Target(distance::Float64, radialSpeed::Float64, azimuth::Float64, swerlingModel::Symbol = :Swerling0, rcs::Float64 = 1.0)
        @assert distance >= 0 "Distance must be non-negative"
        # @assert radialSpeed >= 0 "Radial speed must be non-negative"
        @assert azimuth >= 0 && azimuth < 2π "Azimuth must be in the range [0, 2π)"
        @assert swerlingModel in [:Swerling0, :Swerling1, :Swerling2, :Swerling3, :Swerling4] "Unsupported Swerling model: $swerlingModel"
        @assert rcs >= 0 "RCS must be non-negative"

        return Target(distance, radialSpeed, azimuth, swerlingModel, rcs)
    end

end