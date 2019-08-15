export RawLinearProblem

struct RawLinearProblem{T}
    A::AbstractMatrix{T}
    b::AbstractVector{T}
    c::AbstractVector{T}
end

mutable struct LPIPVariables{T} 
    x::Vector{T}
    s::Vector{T}
    p::Vector{T}

    function LPIPVariables{T}(m::Int, n::Int) where T
        return new{T}(
                    ones(T, m+n),
                    ones(T, m+n),
                    ones(T, m)
                )
    end
end

struct LPIPLinearProblem{T}
    A::AbstractMatrix{T}
    b::AbstractVector{T}
    c::AbstractVector{T}
    variables::LPIPVariables{T}
    m::Int # Number of constraints
    n::Int # Number of variables
    time::Float64

    function LPIPLinearProblem{T}(linear_problem::RawLinearProblem{T}) where T
        m, n = size(linear_problem.A)
        return new{T}(
                    extended_A(linear_problem.A, m),
                    linear_problem.b,
                    extended_c(linear_problem.c, m),
                    LPIPVariables{T}(m, n),
                    m,
                    n,
                    0.0
               )
    end
end

mutable struct Result{T}
    variables::LPIPVariables{T}
    status::Int
    iter::Int
    time::Float64

    function Result{T}(vars::LPIPVariables{T}, status::Int, iter::Int, time::Float64) where T
        return new(vars, status, iter, time)
    end
end

function extended_A(A::AbstractMatrix{T}, m::Int) where T
    return [A Matrix{T}(I, m, m)]
end

function extended_c(c::AbstractVector{T}, m::Int) where T
    return [c; zeros(T, m)]
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