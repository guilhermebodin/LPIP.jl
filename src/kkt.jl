# file to solve the kkt system
mutable struct NewtonSystem{T}
    J::AbstractMatrix{T}
    F::AbstractVector{T}
    d_result::AbstractVector{T}
    dx_aux::AbstractVector{T}
    dp_aux::AbstractVector{T}

    function NewtonSystem{T}(lpip_pb::LPIPLinearProblem{T}, params::Params) where T
        J = zeros(T, 2lpip_pb.n + lpip_pb.m, 2lpip_pb.n + lpip_pb.m)
        build_first_J!(J, lpip_pb.A, lpip_pb.variables.x, lpip_pb.variables.s, lpip_pb.m, lpip_pb.n)
        return new(J, 
                   zeros(T, 2lpip_pb.n + lpip_pb.m) , 
                   zeros(T, 2lpip_pb.n + lpip_pb.m) , 
                   zeros(T, lpip_pb.m), 
                   zeros(T, lpip_pb.n))
    end
end

function build_J_F!(newton_system::NewtonSystem{T}, lpip_pb::LPIPLinearProblem{T}, params::Params) where T
    mi = params.rho * dot(lpip_pb.variables.x, lpip_pb.variables.s) / lpip_pb.n

    # Build J
    build_J!(newton_system.J, lpip_pb.A, lpip_pb.variables.x, lpip_pb.variables.s, lpip_pb.m, lpip_pb.n)

    # Build F
    build_F!(newton_system.F, newton_system.dx_aux, newton_system.dp_aux,
             lpip_pb.A, lpip_pb.b, lpip_pb.c, 
             lpip_pb.variables.x, lpip_pb.variables.s,
             lpip_pb.variables.p, mi, lpip_pb.m, lpip_pb.n)

    return 
end

function build_first_J!(J::AbstractMatrix{T}, A::AbstractMatrix{T}, x::AbstractVector{T}, 
                        s::AbstractVector{T}, m::Int, n::Int) where T
    
    @inbounds for j in 1:n
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

function build_J!(J::AbstractMatrix{T}, A::AbstractMatrix{T}, x::AbstractVector{T}, 
                  s::AbstractVector{T}, m::Int, n::Int) where T

    @inbounds for j in 1:n
        J[m + n + j, j] = s[j]
        J[m + n + j, m + n + j] = x[j]
    end
    return
end

function build_F!(F::AbstractVector{T}, dx_aux::AbstractVector{T}, dp_aux::AbstractVector{T},
                  A::AbstractMatrix{T}, b::AbstractVector{T},
                  c::AbstractVector{T}, x::AbstractVector{T}, s::AbstractVector{T},
                  p::AbstractVector{T}, mi::T, m::Int, n::Int) where T

    LinearAlgebra.BLAS.gemv!('N', 1.0, A, x, 0.0, dx_aux) # newton_system.d.dx_aux = Ax
    LinearAlgebra.BLAS.gemv!('T', 1.0, A, p, 0.0, dp_aux) # newton_system.d.dp_aux = A'p
    @. F[(n + m + 1):(n + m + n)] = x*s - mi
    @. F[(m + 1):(m + n)] = dp_aux + s - c
    @. F[1:m] = dx_aux - b
    return
end

function solve_kkt(newton_system::NewtonSystem{T}, lpip_pb::LPIPLinearProblem{T}, params::Params) where T
    build_J_F!(newton_system, lpip_pb, params)
    newton_system.d_result = newton_system.J \ - newton_system.F # solve the system
    return 
end