using Distributions

module SwerlingModels

    function swerlingAmplitude(rcs::Float64, model::Symbol, n::Int = 1)
        """
        Calculate the Swerling amplitude based on the RCS, model type, and number of pulses.
        
        Parameters:
        - rcs: Radar Cross Section (RCS) in square meters.
        - model: Swerling model type (e.g., :Swerling1, :Swerling2, etc.).
        - n: Number of pulses (default is 1).
        
        Returns:
        - Amplitude as a Float64.
        """

        @assert rcs >= 0 "RCS must be non-negative"
        @assert model in [:Swerling0, :Swerling1, :Swerling2] "Unsupported Swerling model: $model"

        if model == :Swerling0
            # Constant amplitude, no fluctuation
            return sqrt(rcs)
        elseif model == :Swerling1 || model == :Swerling2
            # Fluctuates per pulse or scan, exponential distribution i.e. (χ², 2 DoF)
            return sqrt(rcs * rand(Exponential(1.0)))
        elseif model == :Swerling2 || model == :Swerling4
            # Fluctuates per pulse or scan, χ² with 4 DoF
            return sqrt(rcs * (rand(Exponential(1.0)) + rand(Exponential(1.0))))
        else
            error("Unsupported Swerling model: $model")
        end
    end

end