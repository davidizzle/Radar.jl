

module Targets

    include("swerling.jl")
    using .SwerlingModels

    abstract type Target end

    mutable struct PolarTarget <: Target
        distance::Float64  # Distance to the target in meters
        azimuth::Float64  # Azimuth angle of the target in radians
        radialSpeed::Float64  # Radial speed of the target in m/s, +Approaching, -Receding
        tangentialSpeed::Float64  # Tangential speed of the target in m/s, +CW, -CCW
        swerlingModel::Symbol  # Swerling model type
        rcs::Float64  # Radar Cross Section (RCS) in square meters
        lastUpdateTime::Float64  # Last update time for the target
        status::Symbol  # Status of the target, e.g., :friend, :foe, :unknown
    end
    
    function PolarTarget(distance::Float64, azimuth::Float64, radialSpeed::Float64, tangentialSpeed::Float64; swerlingModel::Symbol = :Swerling0, rcs::Float64 = 1.0, status = :unknown)
        @assert distance >= 0 "Distance must be non-negative"
        # @assert radialSpeed >= 0 "Radial speed must be non-negative"
        # @assert azimuth >= 0 && azimuth < 360 "Azimuth must be in the range [0, 360)"
        @assert swerlingModel in [:Swerling0, :Swerling1, :Swerling2, :Swerling3, :Swerling4] "Unsupported Swerling model: $swerlingModel"
        @assert rcs >= 0 "RCS must be non-negative"
        
        azimuth = mod(azimuth, 360)
        return PolarTarget(distance, azimuth, radialSpeed, tangentialSpeed, swerlingModel, rcs, time(), status)
    end
    
    mutable struct CartesianTarget <: Target
        x::Float64  # X coordinate of the target in meters
        y::Float64  # Y coordinate of the target in meters
        xVel::Float64  # Velocity along X axis
        yVel::Float64  # Velocity along Y axis
        swerlingModel::Symbol  # Swerling model type
        rcs::Float64  # Radar Cross Section (RCS) in square meters
        lastUpdateTime::Float64  # Last update time for the target
        status::Symbol  # Status of the target, e.g., :friend, :foe, :unknown
    end
    
    function CartesianTarget(x::Float64, y::Float64, xVel::Float64, yVel::Float64; swerlingModel::Symbol = :Swerling0, rcs::Float64 = 1.0, status = :unknown)
        @assert swerlingModel in [:Swerling0, :Swerling1, :Swerling2, :Swerling3, :Swerling4] "Unsupported Swerling model: $swerlingModel"
        @assert rcs >= 0 "RCS must be non-negative"

        return CartesianTarget(x, y, xVel, yVel, swerlingModel, rcs, time(), status)
    end

    function updatePosition(target::PolarTarget)
        # Update the target's position based on its radial speed and the current time
        currentTime = time()
        timeElapsed = currentTime - target.lastUpdateTime
        target.lastUpdateTime = currentTime

        azimuth_rad = target.azimuth * π / 180  # Convert azimuth to radians
        x = target.distance * cos(azimuth_rad)
        y = target.distance * sin(azimuth_rad)
        vx = target.radialSpeed * cos(azimuth_rad) - target.tangentialSpeed * sin(azimuth_rad)
        vy = target.radialSpeed * sin(azimuth_rad) + target.tangentialSpeed * cos(azimuth_rad)

        traversedOrigin = false
        x_old = x
        y_old = y
        x += vx * timeElapsed
        y += vy * timeElapsed
        
        traversedOrigin = (x * x_old <= 0) && (y * y_old <= 0) 

        target.distance = sqrt(x^2 + y^2)
        azimuth_rad = atan(y, x) 
        target.radialSpeed = vx * cos(azimuth_rad) + vy * sin(azimuth_rad)
        target.tangentialSpeed = -vx * sin(azimuth_rad) + vy * cos(azimuth_rad)
        target.azimuth = mod(azimuth_rad * 180 / π, 360)  # Convert back to degrees

        @show target.azimuth

        # Ensure distance does not go negative
        if traversedOrigin 
            target.radialSpeed = -target.radialSpeed  # Reverse direction if distance goes negative
        end
    end

    function updatePosition(target::CartesianTarget)
        # Update the target's position based on its velocities
        currentTime = time()
        timeElapsed = currentTime - target.lastUpdateTime
        target.lastUpdateTime = currentTime
        
        target.x += target.xVel * timeElapsed
        target.y += target.yVel * timeElapsed
    end

    function switchTarget(target::PolarTarget)
        target.azimuth *= 2π / 360  # Convert azimuth to radians
        x = target.distance * cos(target.azimuth)
        y = target.distance * sin(target.azimuth)
        vx = target.radialSpeed * cos(target.azimuth) - target.tangentialSpeed * sin(target.azimuth)
        vy = target.radialSpeed * sin(target.azimuth) + target.tangentialSpeed * cos(target.azimuth)
        return CartesianTarget(x, y, vx, vy, target.swerlingModel, target.rcs, target.lastUpdateTime, target.status)
    end

    function switchTarget(target::CartesianTarget)
        distance = sqrt(target.x^2 + target.y^2)
        azimuth = atan(target.y, target.x) 
        radialSpeed = target.xVel * cos(azimuth) + target.yVel * sin(azimuth)
        tangentialSpeed = -target.xVel * sin(azimuth) + target.yVel * cos(azimuth)
        azimuth = mod(azimuth * 180 / π, 360)
        return PolarTarget(distance, azimuth, radialSpeed, tangentialSpeed, target.swerlingModel, target.rcs, target.lastUpdateTime, target.status)
    end

    function toPolar(target::CartesianTarget)
        return switchTarget(target)
    end
    function toPolar(target::PolarTarget)
        return target
    end
    
    function toCartesian(target::CartesianTarget)
        return target
    end
    function toCartesian(target::PolarTarget)
        return switchTarget(target)
    end

end