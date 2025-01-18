## MOI interface for SimpleDualAscent
import MathOptInterface as MOI

MOI.Utilities.@product_of_sets(RHS, MOI.Zeros)


const OptimizerCache{T <: Real} = MOI.Utilities.GenericModel{
    T,
    MOI.Utilities.ObjectiveContainer{T},
    MOI.Utilities.VariablesContainer{T},
    MOI.Utilities.MatrixOfConstraints{T,
        MOI.Utilities.MutableSparseMatrixCSC{T, Int, MOI.Utilities.OneBasedIndexing},
        Vector{T}, RHS{T},
    },
}


mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    x_primal::Dict{MOI.VariableIndex,T}
    y_dual::Dict{MOI.ConstraintIndex,T}
    zl_dual::Dict{MOI.ConstraintIndex,T}
    zu_dual::Dict{MOI.ConstraintIndex,T}
    termination_status::MOI.TerminationStatusCode
    primal_status::MOI.ResultStatusCode
    dual_status::MOI.ResultStatusCode
    settings::SDASettings
    solve_time::Float64

    function Optimizer{T}() where {T <: Real}
        return new(
            Dict{MOI.VariableIndex,T}(), Dict{MOI.ConstraintIndex,T}(), Dict{MOI.ConstraintIndex,T}(), Dict{MOI.ConstraintIndex,T}(),
            MOI.OPTIMIZE_NOT_CALLED, MOI.UNKNOWN_RESULT_STATUS, MOI.UNKNOWN_RESULT_STATUS, deepcopy(SimpleDualAscent.default_settings), 0.0,
        )
    end

    function Optimizer()
        return Optimizer{Float64}()
    end
end


MOI.supports_constraint(::Optimizer{T}, ::Type{MOI.VectorAffineFunction{T}}, ::Type{MOI.Zeros}) where {T <: Real} = true
MOI.supports_constraint(::Optimizer{T}, ::Type{MOI.VariableIndex}, ::Type{MOI.GreaterThan{T}}) where {T <: Real} = true
MOI.supports_constraint(::Optimizer{T}, ::Type{MOI.VariableIndex}, ::Type{MOI.LessThan{T}}) where {T <: Real} = true
MOI.supports(::Optimizer{T}, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}})  where {T <: Real} = true


function MOI.optimize!(dest::Optimizer{T}, src::MOI.ModelLike) where {T<:Real}
    cache = OptimizerCache{T}()
    index_map = MOI.copy_to(cache, src)
    
    A, b, c, l, u = cache_to_data(cache)

    x, y, zₗ, zᵤ, st, r = solve_sda(A, b, c, l, u, dest.settings)

    rnorm = LinearAlgebra.norm(r)
    ts, ps = if rnorm ≈ 0
        MOI.OPTIMAL, MOI.FEASIBLE_POINT
    elseif rnorm < dest.settings.tol
        MOI.ALMOST_OPTIMAL, MOI.NEARLY_FEASIBLE_POINT
    elseif rnorm > dest.settings.tol
        MOI.ITERATION_LIMIT, MOI.INFEASIBLE_POINT
    end

    populate_dest!(dest, src, index_map, x, y, zₗ, zᵤ, st, ts, ps)

    return index_map, false
end

MOI.set(model::Optimizer, attr::MOI.RawOptimizerAttribute, value) = setfield!(model.settings, Symbol(attr.name), value)
MOI.set(model::Optimizer, ::MOI.Silent, value) = setfield!(model.settings, :verbose, !value)

MOI.get(::Optimizer{T}, ::MOI.SolverName) where {T <: Real} = T !== Float64 ? "SimpleDualAscent{$(T)}" : "SimpleDualAscent"
MOI.get(model::Optimizer, attr::MOI.RawOptimizerAttribute) = getfield(model.settings, Symbol(attr.name))
MOI.get(model::Optimizer, ::MOI.Silent) = !model.settings.verbose
MOI.get(model::Optimizer, ::MOI.SolveTimeSec) = model.solve_time
MOI.get(model::Optimizer, ::MOI.VariablePrimal, x::MOI.VariableIndex) = model.x_primal[x]
MOI.get(model::Optimizer{T}, ::MOI.ConstraintDual, x::MOI.ConstraintIndex{MOI.VectorAffineFunction{T},MOI.Zeros}) where {T <: Real} = model.y_dual[x]
MOI.get(model::Optimizer{T}, ::MOI.ConstraintDual, x::MOI.ConstraintIndex{MOI.VariableIndex,MOI.GreaterThan{T}})  where {T <: Real} = model.zl_dual[x]
MOI.get(model::Optimizer{T}, ::MOI.ConstraintDual, x::MOI.ConstraintIndex{MOI.VariableIndex,MOI.LessThan{T}})  where {T <: Real} = model.zu_dual[x]
MOI.get(model::Optimizer, ::MOI.ResultCount) = model.termination_status == MOI.OPTIMAL ? 1 : 0
MOI.get(model::Optimizer, ::MOI.RawStatusString) = "$(model.termination_status)"
MOI.get(model::Optimizer, ::MOI.TerminationStatus) = model.termination_status
MOI.get(model::Optimizer, ::MOI.PrimalStatus) = model.primal_status
MOI.get(model::Optimizer, ::MOI.DualStatus) = model.dual_status

MOI.is_empty(model::Optimizer) = (
    isempty(model.x_primal) && isempty(model.y_dual) && isempty(model.zl_dual) && isempty(model.zu_dual) &&
    model.termination_status == MOI.OPTIMIZE_NOT_CALLED && model.solve_time == 0.0 &&
    model.primal_status == MOI.UNKNOWN_RESULT_STATUS && model.dual_status == MOI.UNKNOWN_RESULT_STATUS
)

function MOI.empty!(model::Optimizer)
    empty!(model.x_primal); empty!(model.y_dual); empty!(model.zl_dual); empty!(model.zu_dual)
    model.termination_status = MOI.OPTIMIZE_NOT_CALLED; model.solve_time = 0.0
    model.primal_status = MOI.UNKNOWN_RESULT_STATUS; model.dual_status = MOI.UNKNOWN_RESULT_STATUS
    return # NOTE: settings are kept!
end