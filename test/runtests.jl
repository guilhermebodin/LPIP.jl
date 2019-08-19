import Pkg
Pkg.activate(".")
push!(LOAD_PATH, "/Users/guilhermebodin/LPIP.jl/src")
using LPIP, Test


A = [
    2. 1 
    1 2
]
b = [4.;4]
c = [-4.; -1]
d = 3.

lpip_pb = LPIPLinearProblem{Float64}(A, b, c, d, 0)
params = LPIP.Params(;rho = 0.01, verbose = true)

@time lp = LPIP.interior_points(lpip_pb, params)

using JuMP, GLPK

model = Model(with_optimizer(GLPK.Optimizer))
@variable(model, x[1:2] >= 0)
@constraint(model, A*x .== b)
@objective(model, Min, c'x + d)
@time optimize!(model)
JuMP.termination_status(model)
JuMP.value.(x)
JuMP.objective_value(model)