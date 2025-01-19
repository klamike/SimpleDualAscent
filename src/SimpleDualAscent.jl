module SimpleDualAscent

using SparseArrays
using LinearAlgebra

include("SDA_utils.jl")
include("SDA_settings.jl")
function solve_sda(A::M{T}, b::V{T}, c::V{T}, l::V{T}, u::V{T}, settings) where {T<:Real}
    started = time()

    x = zeros(size(A, 2))
    y = zeros(size(A, 1))
    z = zeros(size(A, 2))
    r = zeros(size(A, 1))
    μ = settings.mu_init

    recover_z!(z, c, A, y)
    recover_x_barrier!(x, z, l, u, μ)
    residual!(r, b, A, x)

    α = settings.stepsize
    τ = settings.tol
    Nmax = settings.maxit

    i = 1
    while LinearAlgebra.norm(r) > τ && i < Nmax
        update_y!(y, α, r)
        iteration_log(µ, i, r, settings)
        recover_z!(z, c, A, y)
        recover_x_barrier!(x, z, l, u, μ)
        residual!(r, b, A, x)
        update_μ!(μ, r, i, settings)
        i += 1
    end
    recover_x!(x, z, l, u)
    zₗ, zᵤ = extract_zlzu(z)
    finished = time()

    if settings.verbose
        @info "Finished $i: r=$(LinearAlgebra.norm(r)) µ=$μ"
    end

    return x, y, zₗ, zᵤ, finished - started, r
end

include("MOI_wrapper.jl")
include("MOI_utils.jl")

end
