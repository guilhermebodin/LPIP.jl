export RawLinearProblem

struct RawLinearProblem{T}
    A::AbstractMatrix{T}
    b::AbstractVector{T}
    c::AbstractVector{T}
    d::T
end

mutable struct LPIPVariables{T} 
    x::Vector{T}
    s::Vector{T}
    p::Vector{T}

    function LPIPVariables{T}(m::Int, n::Int) where T
        return new{T}(
                    ones(T, n),
                    ones(T, n),
                    ones(T, m)
                )
    end
end

struct LPIPLinearProblem{T}
    A::AbstractMatrix{T}
    b::AbstractVector{T}
    c::AbstractVector{T}
    d::T
    variables::LPIPVariables{T}
    m::Int # Number of equality constraints
    n::Int # Number of variables

    function LPIPLinearProblem{T}(linear_problem::RawLinearProblem{T}) where T
        m, n = size(linear_problem.A)
        return new{T}(
                    linear_problem.A,
                    linear_problem.b,
                    linear_problem.c,
                    linear_problem.d,
                    LPIPVariables{T}(m, n),
                    m,
                    n,
               )
    end
end

mutable struct Result{T}
    variables::LPIPVariables{T}
    obj_val::T
    dual_obj_val::T
    status::Int
    iter::Int
    time::Float64

    function Result{T}(vars::LPIPVariables{T}, obj_val::T, dual_obj_val::T, 
                       status::Int, iter::Int, time::Float64) where T
        return new(vars, obj_val, dual_obj_val, status, iter, time)
    end
end

function Result(lpip_pb::LPIPLinearProblem{T}, status::Int, iter::Int, time::Float64) where T
    obj_val = dot(lpip_pb.c, lpip_pb.variables.x) + lpip_pb.d
    dual_obj_val = dot(lpip_pb.b, lpip_pb.variables.p) + lpip_pb.d
    return Result{T}(lpip_pb.variables, obj_val, dual_obj_val, status, iter, time)
end

mutable struct Params
    verbose::Bool
    max_iter::Int
    rho::Float64
    alpha::Float64
    tol::Float64
    time_limit::Float64

    function Params(;verbose::Bool = false, max_iter::Int = 100, rho::Float64 = 0.1,
                    alpha::Float64 = 0.9, tol::Float64 = 1e-6, time_limit::Float64 = 1e10)
        (rho <= 0 || rho > 1) && error("rho must be between in (0, 1]")
        (alpha <= 0 || alpha >= 1) && error("alpha must be in (0, 1)")
        return new(verbose, max_iter, rho, alpha, tol, time_limit)
    end
end