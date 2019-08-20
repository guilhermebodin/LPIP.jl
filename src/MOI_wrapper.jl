using MathOptInterface

const MOI = MathOptInterface
const CI = MOI.ConstraintIndex
const VI = MOI.VariableIndex

const MOIU = MOI.Utilities

const SF = MOI.VectorAffineFunction{Float64}
const SS = Union{MOI.Zeros, MOI.Nonnegatives}

mutable struct MOISolution
    status::Int
    primal::Vector{Float64} # primal of variables
    slack::Vector{Float64}
    dual::Vector{Float64} # dual of constraints
    obj_val::Float64
    dual_objval::Float64
    solve_time::Float64
    iter::Int
end
MOISolution() = MOISolution(0, Float64[], Float64[], Float64[], NaN, NaN, NaN, 0)

mutable struct ModelData
    m::Int # Numbero of constraints
    n::Int # Number of variables
    I::Vector{Int} # List of rows
    J::Vector{Int} # List of cols
    V::Vector{Float64} # List of coefficients
    b::Vector{Float64} # constants
    obj_constant::Float64 # The objective is min c'x + obj_constant
    c::Vector{Float64} # Objective coefficients
end

mutable struct ConeData
    f::Int # length of the zero cone (equality constraints)
    l::Int # length of the nonnegatives cone (<= constraints)
    function ConeData()
        return new(0, 0)
    end
end

mutable struct Optimizer <: MOI.AbstractOptimizer
    cone::ConeData
    data::Union{Nothing, ModelData} # only non-Void between MOI.copy_to and MOI.optimize!
    maxsense::Bool
    sol::MOISolution
    params::Params
    function Optimizer(kwargs...)
        return new(ConeData(), nothing, false, MOISolution(), Params(kwargs...))
    end
end
function Optimizer(;args...)
    return Optimizer(args)
end


function MOI.is_empty(optimizer::Optimizer)
    return !optimizer.maxsense && optimizer.data === nothing
end
function MOI.empty!(optimizer::Optimizer)
    optimizer.maxsense = false
    optimizer.data = nothing # It should already be nothing except if an error is thrown inside copy_to
    optimizer.sol.status = 0
    return
end

MOIU.supports_allocate_load(model::Optimizer, copy_names::Bool) = true

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

function MOI.supports(::Optimizer, ::MOI.VariablePrimalStart,
                      ::Type{MOI.VariableIndex})
    return false
end

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; kws...)
    return MOIU.automatic_copy_to(dest, src; kws...)
end

# Computes cone dimensions
constroffset(cone::ConeData, ci::CI{<:MOI.AbstractFunction, MOI.Zeros}) = ci.value
function _allocate_constraint(cone::ConeData, f::MOI.AbstractFunction, s::MOI.Zeros)
    ci = cone.f
    cone.f += MOI.dimension(s)
    return ci
end
constroffset(cone::ConeData, ci::CI{<:MOI.AbstractFunction, MOI.Nonnegatives}) = cone.f + ci.value
function _allocate_constraint(cone::ConeData, f::MOI.AbstractFunction, s::MOI.Nonnegatives)
    ci = cone.l
    cone.l += MOI.dimension(s)
    return ci
end
constroffset(optimizer::Optimizer, ci::CI) = constroffset(optimizer.cone, ci::CI)
function MOIU.allocate_constraint(optimizer::Optimizer, f::F, s::S) where {F <: MOI.AbstractFunction, S <: MOI.AbstractSet}
    return CI{F, S}(_allocate_constraint(optimizer.cone, f, s))
end
constrrows(s::MOI.AbstractVectorSet) = 1:MOI.dimension(s)

output_index(t::MOI.VectorAffineTerm) = t.output_index
variable_index_value(t::MOI.ScalarAffineTerm) = t.variable_index.value
variable_index_value(t::MOI.VectorAffineTerm) = variable_index_value(t.scalar_term)
coefficient(t::MOI.ScalarAffineTerm) = t.coefficient
coefficient(t::MOI.VectorAffineTerm) = coefficient(t.scalar_term)

function MOIU.load_constraint(optimizer::Optimizer, ci::CI, f::MOI.VectorAffineFunction, s::MOI.AbstractVectorSet)
    A = sparse(output_index.(f.terms), variable_index_value.(f.terms), coefficient.(f.terms))
    # sparse combines duplicates with + but does not remove zeros created so we call dropzeros!
    dropzeros!(A)
    I, J, V = findnz(A)

    offset = constroffset(optimizer, ci)
    rows = constrrows(s)
    i = offset .+ rows

    b = f.constants
    optimizer.data.b[i] = b
    append!(optimizer.data.I, offset .+ I)
    append!(optimizer.data.J, J)
    append!(optimizer.data.V, -V)
    return
end

function MOIU.allocate_variables(optimizer::Optimizer, nvars::Integer)
    optimizer.cone = ConeData()
    return VI.(1:nvars)
end

# allocate block
function MOIU.allocate(::Optimizer, ::MOI.VariablePrimalStart,
                       ::MOI.VariableIndex, ::Union{Nothing, Float64})
end
function MOIU.allocate(::Optimizer, ::MOI.ConstraintPrimalStart,
                       ::MOI.ConstraintIndex,
                       ::Union{Nothing, AbstractVector{Float64}})
end
function MOIU.allocate(::Optimizer, ::MOI.ConstraintDualStart,
                       ::MOI.ConstraintIndex,
                       ::Union{Nothing, AbstractVector{Float64}})
end
function MOIU.allocate(optimizer::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    optimizer.maxsense = sense == MOI.MAX_SENSE
end
function MOIU.allocate(::Optimizer, ::MOI.ObjectiveFunction,
                       ::MOI.Union{MOI.SingleVariable,
                                   MOI.ScalarAffineFunction{Float64}})
end

function MOIU.load_variables(optimizer::Optimizer, nvars::Integer)
    cone = optimizer.cone
    m = cone.f + cone.l
    I = Int[]
    J = Int[]
    V = Float64[]
    b = zeros(m)
    c = zeros(nvars)
    optimizer.data = ModelData(m, nvars, I, J, V, b, 0., c)
    return 
end

function MOIU.load(::Optimizer, ::MOI.ObjectiveSense, ::MOI.OptimizationSense)
end
function MOIU.load(optimizer::Optimizer, ::MOI.ObjectiveFunction,
                   f::MOI.SingleVariable)
    MOIU.load(optimizer,
              MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
              MOI.ScalarAffineFunction{Float64}(f))
end
function MOIU.load(optimizer::Optimizer, ::MOI.ObjectiveFunction,
                   f::MOI.ScalarAffineFunction)
    c0 = Vector(sparsevec(variable_index_value.(f.terms), coefficient.(f.terms),
                          optimizer.data.n))
    optimizer.data.obj_constant = f.constant
    optimizer.data.c = optimizer.maxsense ? -c0 : c0
    return nothing
end

# MOI getters 
MOI.get(::Optimizer, ::MOI.SolverName) = "LPIP"
MOI.get(optimizer::Optimizer, ::MOI.SolveTime) = optimizer.sol.solve_time
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

function MOI.get(optimizer::Optimizer, ::MOI.VariablePrimal, vi::VI)
    return optimizer.sol.primal[vi.value]
end
function MOI.get(optimizer::Optimizer, a::MOI.VariablePrimal, vi::Vector{VI})
    return MOI.get.(optimizer, a, vi)
end
MOI.get(optimizer::Optimizer, ::MOI.ObjectiveValue) = optimizer.sol.obj_val
function MOI.get(optimizer::Optimizer, ::MOI.ConstraintPrimal,
                 ci::CI{<:MOI.AbstractFunction, S}) where S <: MOI.AbstractSet
    offset = constroffset(optimizer, ci)
    rows = constrrows(optimizer, ci)
    primal = optimizer.sol.slack[offset .+ rows]
    return primal
end
function MOI.get(optimizer::Optimizer, ::MOI.ConstraintDual,
                 ci::CI{<:MOI.AbstractFunction, S}) where S <: MOI.AbstractSet
    offset = constroffset(optimizer, ci)
    rows = constrrows(optimizer, ci)
    dual = optimizer.sol.dual[offset .+ rows]
    return dual
end

# MOI setters

# MOI optimize
function MOI.optimize!(optimizer::Optimizer)
    cone = optimizer.cone
    A = sparse(optimizer.data.I, optimizer.data.J, optimizer.data.V)
    b = optimizer.data.b
    obj_constant = optimizer.data.obj_constant
    c = optimizer.data.c
    # optimizer.data = nothing # Allows GC to free optimizer.data before A is loaded to LPIP
    
    lpip_pb = LPIP.LPIPLinearProblem{Float64}(A, b, c, obj_constant, cone.l)

    t0 = time()
    sol = LPIP.interior_points(lpip_pb, optimizer.params)
    solve_time = time() - t0

    # Query solutions
    status = sol.status
    primal = sol.variables.x
    slack = sol.variables.x[end - cone.l:end]
    dual = sol.variables.p
    obj_val = sol.obj_val
    dual_obj_val = sol.dual_obj_val
    iter = sol.iter
    optimizer.sol = MOISolution(status, primal, slack, dual, obj_val, dual_obj_val, solve_time, iter)
    return
end