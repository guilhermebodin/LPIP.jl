import Pkg
Pkg.activate(".")
push!(LOAD_PATH, "/home/guilhermebodin/Documents/Github/LPIP.jl/src")
using LPIP, Test, GLPK

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

A = rand(10, 10)
b = rand(10)
c = rand(10)

model = Model(with_optimizer(LPIP.Optimizer, rho = 0.01, verbose = true))
# model = Model(with_optimizer(GLPK.Optimizer))
@variable(model, x[1:10])
@constraint(model, A * x - b in MOI.Zeros(10))
@objective(model, Min, c'x + 3)
optimize!(model)
JuMP.value.(x)
JuMP.objective_value(model)
JuMP.termination_status(model)