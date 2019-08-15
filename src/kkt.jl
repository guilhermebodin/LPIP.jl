# file to solve the kkt system
mutable struct NewtonDirections{T}
    dx::AbstractVector{T}
    ds::AbstractVector{T}
    dp::AbstractVector{T}

    function NewtonDirections{T}(d::AbstractVector{T}, m::Int, n::Int) where T
        return new(
            zeros(T, m+n), # x
            zeros(T, m+n), # s
            zeros(T, m)    # p
        )
    end
end

mutable struct NewtonSystem{T}
    J::AbstractMatrix{T}
    F::AbstractVector{T}
    d::NewtonDirections{T}

    function NewtonSystem{T}(lpip_pb::LPIPLinearProblem{T}, params::Params) where T
        m = lpip_pb.m
        n = lpip_pb.n
        J = zeros(T, 2(m + n) + m, 2(m+n) + m) 
        F = zeros(T, 2(m + n) + m) 
        return new(J, F, NewtonDirections{T}(zeros(T, 2(m + n) + m), m, n))
    end
end

function build_J_F!(newton_system::NewtonSystem{T}, lpip_pb::LPIPLinearProblem{T}, params::Params) where T
    # rename m, n
    m = lpip_pb.m
    n = lpip_pb.n
    
    # Build diagonal matrices
    X = Diagonal(lpip_pb.variables.x)
    S = Diagonal(lpip_pb.variables.s)
    # Auxiliary e
    e = ones(T, m + n) # m + n
    mi = params.rho * dot(lpip_pb.variables.x, lpip_pb.variables.s) / (m + n)

    # Build J
    build_J!(newton_system, lpip_pb.A, X, S, m, n)

    # Build F
    build_F!(newton_system, lpip_pb.A, lpip_pb.b, 
             lpip_pb.c, X, S, lpip_pb.variables.x, lpip_pb.variables.s,
             lpip_pb.variables.p, e, mi, m, n)

    return 
end

function build_J!(newton_system::NewtonSystem{T}, A::AbstractMatrix{T}, X::AbstractMatrix{T}, 
                 S::AbstractMatrix{T}, m::Int, n::Int) where T

    newton_system.J = [A zeros(m, m) zeros(m, m + n);
                       zeros(m + n, m + n) A' Matrix(I, m + n, m + n);
                       S zeros(m + n, m) X]
    return
end

function build_F!(newton_system::NewtonSystem{T}, A::AbstractMatrix{T}, b::AbstractVector{T},
                  c::AbstractVector{T}, X::AbstractMatrix{T}, S::AbstractMatrix{T},
                  x::AbstractVector{T}, s::AbstractVector{T}, p::AbstractVector{T},
                  e::AbstractVector{T}, mi::T, m::Int, n::Int) where T

    newton_system.F = [A * x - b;
                       A' * p + s - c;
                       X * S * e - mi * e]
    return
end

function solve_kkt(newton_system::NewtonSystem{T}, lpip_pb::LPIPLinearProblem{T}, params::Params) where T
    build_J_F!(newton_system, lpip_pb, params)
    d = newton_system.J \ - newton_system.F # solve the system

    # Organize the solution
    m = lpip_pb.m
    n = lpip_pb.n

    newton_system.d.dx = d[1:m + n]
    newton_system.d.dp = d[(m + n + 1):(m + n + m)]
    newton_system.d.ds = d[(m + n + m + 1):(2 * (m + n) + m)]
    return 
end