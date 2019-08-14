# file to solve the kkt system
mutable struct NewtonDirections{T}
    dx::AbstractVector{T}
    ds::AbstractVector{T}
    dp::AbstractVector{T}

    function NewtonDirections{T}(d::AbstractVector{T}, lpip_pb::LPIPLinearProblem{T}) where T

        dx = [d[i] for i=1:m+n] # m + n
        dp = [d[i] for i=(m+n)+1:(m+n)+m] # m
        ds = [d[i] for i=(m+n+m)+1:(m+n+m)+m+n] # m + n
    end
end

mutable struct NewtonSystem{T}
    J::AbstractMatrix{T}
    F::AbstractVector{T}
    d::NewtonDirections{T}

    function NewtonSystem{T}(lpip_pb::LPIPLinearProblem{T}) where T
        
    end
end

function fill_newton_system!(newton_system::NewtonSystem{T}, lpip_pb::LPIPLinearProblem{T}) where T
    # # Build newton system (J and F)
    # μ = ρ * s' * x / (m + n)
    # X = Diagonal(x)
    # S = Diagonal(s)
    # # Solve the linear system
    # J = [A zeros(m, m) zeros(m, m + n);
    #     zeros(m + n, m + n) A' Matrix(I, m + n, m + n);
    #     S zeros(m + n, m) X]

    # F = [A*x-b;
    #      A'*p+s-c;
    #      X*S*e-μ*e]
    return 
end

function solve_kkt(newton_system::NewtonSystem{T}) where T
    d = newton_system.J \ - newton_system.F


    return 
end

       