module SimpleDualAscent

using SparseArrays
using LinearAlgebra

function solve_sda(A, b, c, l, u, settings)
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
            @info "Iteration $i: $(LinearAlgebra.norm(r)) $(b - A * x) $z, $x, $y"
        end

        z = c - A' * y

        x[z .< 0] .= u[z .< 0]
        x[z .> 0] .= l[z .> 0]
        x[z .== 0] .= (u[z .== 0] + l[z .== 0]) / 2

        r = (b - A * x)
    end

    @info "Finished $i: $(LinearAlgebra.norm(r)) $(b - A * x) $z, $x, $y"

    zₗ = max.(z, 0)
    zᵤ = max.(-z, 0)
    
    finished = time()
    return x, y, zₗ, zᵤ, finished - started
end

include("settings.jl")
include("MOI_wrapper.jl")
include("utils.jl")

end
