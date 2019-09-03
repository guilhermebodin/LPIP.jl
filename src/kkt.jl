# file to solve the kkt system
mutable struct NewtonDirections{T}
    dx::AbstractVector{T}
    ds::AbstractVector{T}
    dp::AbstractVector{T}

    function NewtonDirections{T}(d::AbstractVector{T}, m::Int, n::Int) where T
        return new(
            zeros(T, n), # x
            zeros(T, n), # s
            zeros(T, m)  # p
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
        J = zeros(T, 2n + m, 2n + m)
        build_first_J!(J, lpip_pb.A, lpip_pb.variables.x, lpip_pb.variables.s, m, n)
        F = zeros(T, 2n + m) 
        return new(J, F, NewtonDirections{T}(zeros(T, 2n + m), m, n))
    end
end

function build_J_F!(newton_system::NewtonSystem{T}, lpip_pb::LPIPLinearProblem{T}, params::Params) where T
    # rename m, n
    m = lpip_pb.m
    n = lpip_pb.n
    
    mi = params.rho * dot(lpip_pb.variables.x, lpip_pb.variables.s) / n

    # Build J
    build_J!(newton_system, lpip_pb.A, lpip_pb.variables.x, lpip_pb.variables.s, m, n)

    # Build F
    build_F!(newton_system, lpip_pb.A, lpip_pb.b, 
             lpip_pb.c, lpip_pb.variables.x, lpip_pb.variables.s,
             lpip_pb.variables.p, mi, m, n)

    return 
end

function build_first_J!(J::AbstractMatrix{T}, A::AbstractMatrix{T}, x::AbstractVector{T}, 
                        s::AbstractVector{T}, m::Int, n::Int) where T
    
    for j in 1:n
        J[m + n + j, j] = s[j]
        J[m + j, m + n + j] = one(T)
        J[m + n + j, m + n + j] = x[j]
        for i in 1:m
            J[i, j] = A[i, j]
            J[m + j, n + i] = A[i, j]
        end
    end
    return 
end

function build_J!(newton_system::NewtonSystem{T}, A::AbstractMatrix{T}, x::AbstractVector{T}, 
                  s::AbstractVector{T}, m::Int, n::Int) where T

    for j in 1:n
        newton_system.J[m + n + j, j] = s[j]
        newton_system.J[m + n + j, m + n + j] = x[j]
    end
    return
end

function build_F!(newton_system::NewtonSystem{T}, A::AbstractMatrix{T}, b::AbstractVector{T},
                  c::AbstractVector{T}, x::AbstractVector{T}, s::AbstractVector{T},
                  p::AbstractVector{T}, mi::T, m::Int, n::Int) where T

    dx_aux = A * x - b
    dp_aux = A' * p + s - c
    for i in 1:n
        newton_system.F[n + m + i] = x[i]*s[i] - mi
        newton_system.F[m + i] = dp_aux[i]
    end
    for i in 1:m
        newton_system.F[i] = dx_aux[i]
    end
    return
end

function solve_kkt(newton_system::NewtonSystem{T}, lpip_pb::LPIPLinearProblem{T}, params::Params) where T
    build_J_F!(newton_system, lpip_pb, params)
    d = newton_system.J \ - newton_system.F # solve the system

    # Organize the solution
    m = lpip_pb.m
    n = lpip_pb.n

    newton_system.d.dx = d[1:n]
    newton_system.d.dp = d[(n + 1):(n + m)]
    newton_system.d.ds = d[(n + m + 1):(n + m + n)]
    return 
end