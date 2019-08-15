push!(LOAD_PATH, "/Users/guilhermebodin/LPIP.jl/src")
using LPIP


# A = [
#     2 1 
#     1 2
# ]
# b = [4;4]
# c = [-4; -3]
# d = 3

n = 100
A = rand(n, n)
b = rand(n)
c = -rand(n)
d = 0

linear_problem = RawLinearProblem{Float64}(A, b, c, d)
lpip_pb = LPIP.LPIPLinearProblem{Float64}(linear_problem)
params = LPIP.Params()

@time lp = LPIP.interior_points(lpip_pb, params)

using JuMP, Clp

model = Model(with_optimizer(Clp.Optimizer))
@variable(model, x[1:n])
@constraint(model, A*x .<= b)
@objective(model, Min, c'x)
@time optimize!(model)