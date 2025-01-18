mutable struct SDASettings{T<:Real}
    maxit::Int64
    tol::T
    verbose::Bool
    freq::Int64
    stepsize::T
    function SDASettings{T}(;maxit, tol, verbose, freq, stepsize)  where {T<:Real}
        new(maxit, tol, verbose, freq, stepsize)
    end
    function SDASettings(;maxit, tol, verbose, freq, stepsize)
        new{Float64}(maxit, tol, verbose, freq, stepsize)
    end
end

# default setting
default_settings = SDASettings{Float64}(maxit=100000, tol=1e-4, verbose=false, freq=100, stepsize=1e-3)