using MathOptInterface

const MOI = MathOptInterface
const CI = MOI.ConstraintIndex
const VI = MOI.VariableIndex

const MOIU = MOI.Utilities

const SF = Union{MOI.SingleVariable, MOI.ScalarAffineFunction{Float64}, MOI.VectorOfVariables, MOI.VectorAffineFunction{Float64}}
const SS = Union{MOI.EqualTo{Float64}, MOI.GreaterThan{Float64}, MOI.LessThan{Float64}, MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives}

mutable struct MOISolution
    status::Int
    primal::Vector{Float64} # primal of variables
    dual::Vector{Float64} # dual of constraints
    objval::Float64
    dual_objval::Float64
    gap::Float64
    time::Float64
end
MOISolution() = MOISolution(0, Float64[], Float64[], Float64[], NaN, NaN, NaN, NaN, NaN, NaN, 0)

mutable struct Optimizer <: MOI.AbstractOptimizer
    lpip_pb::LPIPLinearProblem{Float64}
    maxsense::Bool
    sol::MOISolution

    params::Params
    function Optimizer(args)
        new()
    end
end
function Optimizer(;args...)
    return Optimizer(args)
end

function MOI.is_empty(optimizer::Optimizer)
    !optimizer.maxsense && optimizer.data === nothing
end
function MOI.empty!(optimizer::Optimizer)
    optimizer.maxsense = false
    optimizer.data = nothing # It should already be nothing except if an error is thrown inside copy_to
    optimizer.sol.ret_val = 0
end

function MOI.supports(::Optimizer,
                      ::Union{MOI.ObjectiveSense,
                              MOI.ObjectiveFunction{MOI.SingleVariable},
                              MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}})
    return true
end

function MOI.supports_constraint(optimizer::Optimizer, F::Type{<:SF}, S::Type{<:SS})
    return true
end

function MOI.supports_constraint(::Optimizer, ::Type{MOI.AbstractFunction}, ::Type{MOI.AbstractSet})
    return false 
end

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; copy_names::Bool = true)
    return MOIU.allocate_load(dest, src, copy_names)
end

# MOI getters 
MOI.get(::Optimizer, ::MOI.SolverName) = "LPIP"
MOI.get(optimizer::Optimizer, ::MOI.SolveTime) = optimizer.lpip_pb.time
MOI.get(optimizer::Optimizer, ::MOI.PrimalStatus) = optimizer.sol.status == 1 ? MOI.FEASIBLE_POINT : MOI.INFEASIBLE_POINT
MOI.get(optimizer::Optimizer, ::MOI.DualStatus) = optimizer.sol.status == 1 ? MOI.FEASIBLE_POINT : MOI.INFEASIBLE_POINT

function MOI.get(optimizer::Optimizer, ::MOI.TerminationStatus)
    s = optimizer.sol.status
    @assert 0 <= s <= 4
    if s == 0
        return MOI.OPTIMIZE_NOT_CALLED
    elseif s == 1
        return MOI.OPTIMAL
    elseif s == 2
        return MOI.INFEASIBLE_OR_UNBOUNDED
    elseif s == 3
        return MOI.TIME_LIMIT
    elseif s == 4
        return MOI.ITERATION_LIMIT
    end
end
function MOI.get(optimizer::Optimizer, ::MOI.ResultCount)
    if MOI.get(optimizer, MOI.TerminationStatus()) == MOI.INFEASIBLE_OR_UNBOUNDED
        return 0
    else
        return 1
    end
end

# MOI setters

# MOI optimize
function MOI.optimize!(optimizer::Optimizer)
    
end