push!(LOAD_PATH, "/Users/guilhermebodin/LPIP.jl/src")
using LPIP

A = rand(100, 100)
b = rand(100)
c = -rand(100)

linear_problem = RawLinearProblem{Float64}(A, b, c)
params = LPIP.Params()

lp = LPIP.interior_points(linear_problem, params)