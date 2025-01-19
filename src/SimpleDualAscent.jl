module SimpleDualAscent

using SparseArrays
using LinearAlgebra

include("SDA_utils.jl")
include("SDA_settings.jl")
function solve_sda(A::M{T}, b::V{T}, c::V{T}, l::V{T}, u::V{T}, settings) where {T<:Real}
    started = time()

    x = zeros(size(A, 2))
    xb = zeros(size(A, 2))
    y = zeros(size(A, 1))
    z = zeros(size(A, 2))
    r = zeros(size(A, 1))
    rb = zeros(size(A, 1))
    μ = settings.mu_init

    z = recover_z!(z, c, A, y)
    xb = recover_x_barrier!(xb, z, l, u, μ)
    rb = residual!(rb, b, A, xb)

    α = settings.stepsize
    τ = settings.tol
    Nmax = settings.maxit

    i = 1
    while LinearAlgebra.norm(rb) > τ && i < Nmax
        y = update_y!(y, α, rb)
        iteration_log(µ, i, r, rb, settings)
        z = recover_z!(z, c, A, y)
        xb = recover_x_barrier!(xb, z, l, u, μ)
        rb = residual!(rb, b, A, xb)
        
        x = recover_x!(x, z, l, u)
        r = residual!(r, b, A, x)
        μ = update_μ(μ, r, i, settings)
        i += 1
    end
    zₗ, zᵤ = extract_zlzu(z)
    finished = time()

    if settings.verbose
        @info "Finished $i: r=$(LinearAlgebra.norm(r)) rb=$(LinearAlgebra.norm(rb)) µ=$μ"
    end

    return x, y, zₗ, zᵤ, finished - started, r
end

include("MOI_wrapper.jl")
include("MOI_utils.jl")

end
