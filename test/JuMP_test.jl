import Pkg
Pkg.activate(".")
push!(LOAD_PATH, "/Users/guilhermebodin/LPIP.jl/src")
using LPIP, Test

using MathOptInterface, JuMP

const MOI = MathOptInterface

A = [
    2. 1 
    1 2
]
b = [4.;4]
c = [4; 3]
d = 3.

model = Model(with_optimizer(LPIP.Optimizer))
@variable(model, x[1:2] >= 0)
@constraint(model, con, A*x .<= b)
@objective(model, Max, c'x + d)
@time optimize!(model)
JuMP.value.(x)
JuMP.objective_value(model)