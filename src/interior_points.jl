"""
 0 - NotSolved
 1 - Optimal
 2 - Infeasible or unbounded
 3 - Time limit
 4 - Iteration limit
"""
function interior_points(A::AbstractMatrix{T}, b::AbstractVector{T}, c::AbstractVector{T},
                         d::T, l::Int, params::Params) where T
    lpip_pb = LPIPLinearProblem{T}(A, b, c, d, l)
    return interior_points(lpip_pb, params)
end

function interior_points(lpip_pb::LPIPLinearProblem{T}, params::Params) where T
    t0 = time() # Start the timer
    status = 0 # Not solved
    newton_system = NewtonSystem{T}(lpip_pb, params)
    # Interior points iteration
    @inbounds for i in 1:params.max_iter
        # Optimality test
        check_optimality(lpip_pb, params.tol) == 1 && return Result(lpip_pb, 1, i, time() - t0)
        # Solve the system and fill newton directions
        solve_kkt(newton_system, lpip_pb, params)
        # Update x, s, p
        update_lpip_pb_vars(newton_system, lpip_pb, params)
        # Check if problem is ubounded or infeasible
        check_infeasible_or_unbounded(lpip_pb) == 2 && return Result(lpip_pb, 2, i, time() - t0)
        # Check if time limit was reached
        check_time_limit(params.time_limit, t0) == 3 && return Result(lpip_pb, 3, i, time() - t0)
    end
    # Iteration limit
    return Result(lpip_pb, 4, params.max_iter, time() - t0)
end

function check_optimality(lpip_pb::LPIPLinearProblem{T}, tol::Float64) where T
    if abs(dot(lpip_pb.variables.x, lpip_pb.variables.s)) < tol
        return 1 # Optimal
    end
    return 0 # Not solved
end

function check_time_limit(time_limit::Float64, t0::Float64)
    if time() - t0 > time_limit
        return 3 # Time limit
    end
    return 0 # Not solved
end

function check_infeasible_or_unbounded(lpip_pb::LPIPLinearProblem{T}) where T
    
end

function update_lpip_pb_vars(newton_system::NewtonSystem{T}, lpip_pb::LPIPLinearProblem{T}, params::Params) where T
    # Find step lenghts
    beta_primal = 1
    beta_dual = 1

    for i in 1:lpip_pb.n
        if newton_system.d.dx[i] <= 0
            if beta_primal > - lpip_pb.variables.x[i] / newton_system.d.dx[i]
                beta_primal = - lpip_pb.variables.x[i] / newton_system.d.dx[i]
            end
        end
        if newton_system.d.ds[i] <= 0
            if beta_dual > - lpip_pb.variables.s[i] / newton_system.d.ds[i]
                beta_dual = - lpip_pb.variables.s[i] / newton_system.d.ds[i]
            end
        end
    end

    # Update variables
    @. lpip_pb.variables.x += params.alpha * beta_primal * newton_system.d.dx
    @. lpip_pb.variables.p += params.alpha * beta_dual * newton_system.d.dp
    @. lpip_pb.variables.s += params.alpha * beta_dual * newton_system.d.ds
    return
end