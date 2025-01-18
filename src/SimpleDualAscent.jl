module SimpleDualAscent

using SparseArrays
using LinearAlgebra

M, V = AbstractMatrix, AbstractVector
function solve_sda(A::M{T}, b::V{T}, c::V{T}, l::V{T}, u::V{T}, settings) where {T<:Real}
    started = time()

    y = zeros(size(A, 1))  
    
    z = c - A' * y

    x = zeros(size(A, 2))
    x[z .< 0] .= u[z .< 0]
    x[z .> 0] .= l[z .> 0]
    x[z .== 0] .= (u[z .== 0] + l[z .== 0]) / 2
    r = (b - A * x)

    α = settings.stepsize
    τ = settings.tol
    Nmax = settings.maxit
    i = 1
    while LinearAlgebra.norm(r) > τ && i < Nmax
        i += 1
        
        y = y + α .* r
        
        if settings.verbose && i % settings.freq == 0
            @info "Iteration $i: residual=$(LinearAlgebra.norm(r))"
        end

        z = c - A' * y

        x[z .< 0] .= u[z .< 0]
        x[z .> 0] .= l[z .> 0]
        x[z .== 0] .= (u[z .== 0] + l[z .== 0]) / 2

        r = (b - A * x)
    end

    if settings.verbose
        @info "Finished $i: residual=$(LinearAlgebra.norm(r))"
    end

    zₗ = max.(z, 0)
    zᵤ = max.(-z, 0)
    
    finished = time()
    return x, y, zₗ, zᵤ, finished - started, r
end

include("settings.jl")
include("MOI_wrapper.jl")
include("utils.jl")

end
