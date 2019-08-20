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
c = [-4.; -3]
d = 3.

lpip_pb = LPIPLinearProblem{Float64}(A, b, c, d, 0)
params = LPIP.Params(;rho = 0.01, verbose = true)
@time lp = LPIP.interior_points(lpip_pb, params)

model = Model(with_optimizer(LPIP.Optimizer, rho = 0.01, verbose = true))
@variable(model, x[1:2])
@constraint(model, A * x - b in MOI.Nonpositives(2))
@objective(model, Min, c'x + d)
optimize!(model)
JuMP.value.(x)
JuMP.objective_value(model)

backend(model)