function bound_slacks!(A, l, u)
    no_upper_bound = findall(isinf, u)
    no_lower_bound = findall(isinf, l)

    Ai, Aj, Av = findnz(A)
    Am, An = size(A)

    for j in no_upper_bound
        i = Ai[findall(Aj .== j)]
        if isempty(i)
            error("Unbounded variable $j does not appear in A")
        end
        @assert length(i) == 1 "Unbounded variable $j appears multiple times in A"
        i = first(i)

        pos_x, neg_x = [], []
        for k in 1:An
            if k == j
                continue
            elseif A[i, k] > 0
                push!(pos_x, k)
            elseif A[i, k] < 0
                push!(neg_x, k)
            end
        end

        u[j] = (sum(A[i, neg_x] .* l[neg_x]) + sum(A[i, pos_x] .* u[pos_x]))

        @info "Deduced upper bound $(u[j]) for variable $j, so the bounds are $(l[j]) ≤ x[$j] ≤ $(u[j])"
    end

    for j in no_lower_bound
        i = Ai[findall(Aj .== j)]
        if isempty(i)
            error("Unbounded variable $j does not appear in A")
        end
        @assert length(i) == 1 "Unbounded variable $j appears multiple times in A"
        i = first(i)

        pos_x, neg_x = [], []
        for k in 1:An
            if k == j
                continue
            elseif A[i, k] > 0
                push!(pos_x, k)
            elseif A[i, k] < 0
                push!(neg_x, k)
            end
        end

        l[j] = (sum(A[i, pos_x] .* l[pos_x]) + sum(A[i, neg_x] .* u[neg_x]))

        @info "Deduced lower bound $(l[j]) for variable $j, so the bounds are $(l[j]) ≤ x[$j] ≤ $(u[j])"
    end

    @assert all(isfinite, l) "Could not deduce bounds for all variables."
    @assert all(isfinite, u) "Could not deduce bounds for all variables."
    @assert all(l .<= u) "Lower bounds are greater than upper bounds."
end


function cache_to_data(cache::OptimizerCache{T}) where {T<:Real}
    A = convert(
        SparseArrays.SparseMatrixCSC{T,Int},
        cache.constraints.coefficients,
    )
    b = cache.constraints.constants
    b = -b # because @odow set Ax+b ∈ {0}
    c = zeros(T, size(A, 2))
    offset = cache.objective.scalar_affine.constant
    for term in cache.objective.scalar_affine.terms
        c[term.variable.value] += term.coefficient
    end
    if cache.objective.sense == MOI.MAX_SENSE
        c *= -1
    end
    l = cache.variables.lower
    u = cache.variables.upper

    bound_slacks!(A, l, u)
    return A, b, c, l, u
end


function populate_dest!(dest::Optimizer{T}, src, index_map, x, y, zₗ, zᵤ, solve_time, ts, ps) where {T<:Real}
    for i in MOI.get(src, MOI.ListOfVariableIndices())
        dest.x_primal[i] = x[index_map[i].value]
    end
    for i in MOI.get(src, MOI.ListOfConstraintIndices{MOI.VectorAffineFunction{T},MOI.Zeros}())
        dest.y_dual[i] = y[index_map[i].value]
    end
    for i in MOI.get(src, MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.GreaterThan{T}}())
        dest.zl_dual[i] = zₗ[index_map[i].value]
    end
    for i in MOI.get(src, MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.LessThan{T}}())
        dest.zu_dual[i] = zᵤ[index_map[i].value]
    end
    dest.termination_status = ts
    dest.primal_status = ps
    dest.dual_status = MOI.FEASIBLE_POINT
    dest.solve_time = solve_time
end