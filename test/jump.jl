using JuMP


function simple_model()
    model = Model(SimpleDualAscent.Optimizer)
    @variable(model, x ∈ MOI.Interval(0, 1))
    @variable(model, y ∈ MOI.Interval(0, 1))
    @constraint(model, eq, x + y == 1)
    @objective(model, Min, x + 2y)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x) == 1.0
    @test value(y) == 0.0
    @test dual(eq) ≈ 1.0

    model = Model(SimpleDualAscent.Optimizer)
    @variable(model, 0 ≤ x ≤ 0)
    @variable(model, 0 ≤ y ≤ 1)
    @constraint(model, eq, x + y == 1)
    @objective(model, Min, x + 2y)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x) == 0.0
    @test value(y) == 1.0
    @test isapprox(dual(eq), 2.0, atol=2e-3)

    model = Model(SimpleDualAscent.Optimizer)
    # set_silent(model)
    @variable(model, [0,0][i] ≤ x[i in 1:2] ≤ [1,1][i])
    @constraint(model, eq, x[1] + x[2] == 1)
    @constraint(model, eq2, 2x[1] - x[2] == 2)
    @objective(model, Min, x[1] + 2x[2])
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value.(x) == [1.0, 0.0]
    @test dual(eq) ≈ 0.2
    @test dual(eq2) ≈ 0.4
end

@testset "JuMP" begin
    simple_model()
end
