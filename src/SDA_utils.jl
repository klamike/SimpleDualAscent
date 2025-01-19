M, V = AbstractMatrix, AbstractVector

function recover_z!(z, c, A, y)
    z .= c - A' * y
end

function recover_x!(x, z, l, u)
    x[z.<0] .= u[z.<0]
    x[z.>0] .= l[z.>0]
    x[z.==0] .= (u[z.==0] + l[z.==0]) / 2

    return x
end

function residual!(r, b, A, x)
    r .= b - A * x
end

function update_y!(y, α, r)
    y .= y + α * r
end

function iteration_log(µ, i, r, settings)
    if settings.verbose && i % settings.freq == 0
        @info "Iteration $i: r=$(LinearAlgebra.norm(r)) µ=$µ"
    end
end

function extract_zlzu(z)
    max.(z, 0), max.(-z, 0)
end

function recover_x_barrier!(x, z, l, u, µ)
    if μ > 0
        z₌₀ = z .== 0
        nz = .!z₌₀
        v = µ ./ z[nz]
        x[z₌₀] .= (u[z₌₀] + l[z₌₀]) / 2
        x[nz] .= (u[nz] + l[nz]) / 2 + v - sign.(z[nz]) .* hypot.(v, (u[nz] - l[nz]) / 2)

        return x
    else
        return recover_x!(x, z, l, u)
    end
end

function update_μ!(μ, r, i, settings)
    # TODO: monitor r to adjust μ
    if i % settings.mu_freq == 0
        μ /= 2
    end
    if μ < settings.mu_threshold
        μ = 0
    end
    return μ
end
