mutable struct SDASettings{T<:Real}
    maxit::Int64
    tol::T
    verbose::Bool
    freq::Int64
    stepsize::T
    mu_init::T
    mu_threshold::T
    mu_freq::Int64
    function SDASettings{T}(; kwargs...) where T
        new(
            get(kwargs, :maxit, 1e+6),
            get(kwargs, :tol, 1e-6),
            get(kwargs, :verbose, true),
            get(kwargs, :freq, 1000),
            get(kwargs, :stepsize, 1e-4),
            get(kwargs, :mu_init, 0.0),
            get(kwargs, :mu_threshold, 1e-8),
            get(kwargs, :mu_freq, 1000),
        )
    end
    function SDASettings(; kwargs...)
        new{Float64}(kwargs...)
    end
end

# default setting
default_settings = SDASettings{Float64}()