export LPIPLinearProblem

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
    m::Int # Number of constraints
    n::Int # Number of variables
    f::Int # Number of equality constraints
    l::Int # Number of nonnegative constraints

    function LPIPLinearProblem{T}(A::AbstractMatrix{T}, b::AbstractVector{T}, c::AbstractVector{T},
                                  d::T, l::Int) where T

        # The number of equality constraints must be 
        # number of constraints - number of nonnegative
        f = l > size(A, 1) ? error("The number of nonnegative constraints must be smaller than the number of rows") : size(A, 1) - l
        extended_A = extend_A(A, f, l)
        extended_c = extend_c(c, l)
        m, n = size(extended_A)
        return new{T}(
                    extended_A, b,
                    extended_c, d,
                    LPIPVariables{T}(m, n),
                    m, n, f, l
               )
    end
end

function extend_A(A::AbstractMatrix{T}, f::Int, l::Int) where T
    if l > 0 && f != 0 
        return [A [zeros(f); Matrix{T}(I, l, l)]]
    elseif f == 0
        return [A Matrix{T}(I, l, l)]
    end
    return A
end

function extend_c(c::AbstractVector{T}, l::Int) where T
    if l > 0
        return [c; zeros(T, l)]
    end
    return c
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