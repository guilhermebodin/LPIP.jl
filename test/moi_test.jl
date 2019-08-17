import Pkg
Pkg.activate(".")
push!(LOAD_PATH, "/Users/guilhermebodin/LPIP.jl/src")
using LPIP, Test

using MathOptInterface

const MOI = MathOptInterface
const CI = MOI.ConstraintIndex
const VI = MOI.VariableIndex

const MOIU = MOI.Utilities

MOIU.@model LPIPModelData () () (MOI.Zeros, MOI.Nonnegatives) () (MOI.SingleVariable,) (MOI.ScalarAffineFunction,) () (MOI.VectorAffineFunction,)

const optimizer = MOIU.CachingOptimizer(LPIPModelData{Float64}(), LPIP.Optimizer())

MOI.empty!(optimizer)
@test MOI.is_empty(optimizer)
#= 
    min -3x - 2y - 4z
s.t.    
    x +  y +  z == 3
    y +  z == 2
    x>=0 y>=0 z>=0
=#

v = MOI.add_variables(optimizer, 3)

vov = MOI.VectorOfVariables(v)

c = MOI.add_constraint(optimizer, MOI.VectorAffineFunction(MOI.VectorAffineTerm.([1,1,1,2,2], MOI.ScalarAffineTerm.(1.0, [v;v[2];v[3]])), [-3.0,-2.0]), MOI.Zeros(2))

MOI.set(optimizer, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.([-3.0, -2.0, -4.0], v), 0.0))
MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)

MOI.optimize!(optimizer)

obj = MOI.get(optimizer, MOI.ObjectiveValue())

@test obj ≈ -9.33333 atol = 1e-2

Xr = MOI.get(optimizer, MOI.VariablePrimal(), X)

@test Xr ≈ [1.3333, 1.3333] atol = 1e-2