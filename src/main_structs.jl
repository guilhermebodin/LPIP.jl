struct RawLinearProblem{T}
    A::AbstractMatrix{T}
    b::AbstractVector{T}
    c::AbstractVector{T}
end

mutable struct LPIPVariables{T}
    x::Vector{T}
    s::Vector{T}
    p::Vector{T}
end

struct LPIPLinearProblem{T}
    A::AbstractMatrix{T}
    b::AbstractVector{T}
    c::AbstractVector{T}
    variables::LPIPVariables{T}
    status::Int
    m::Int # Number of constraints
    n::Int # Number of variables
end

function build_lpip_problem(linear_problem::RawLinearProblem{T}) where T

end

mutable struct Params
    verbose::Bool
    max_iter::Int
    rho::Float64
    alpha::Float64
    tol::Float64
end