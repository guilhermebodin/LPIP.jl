"""
-1 - NotStarted
 0 - NotSolved
 1 - Optimal
 2 - Ilimited
 3 - Unbounded
 4 - Time limit
 5 - Iteration limit
"""
function interior_points(linear_problem::RawLinearProblem{T}, params::Params) where T
    
    status = -1 # Not solved
    lpip_pb = LPIP.LPIPLinearProblem{T}(linear_problem)
    newton_system = NewtonSystem{T}(lpip_pb, params)
    t0 = time()
    # Interior points iteration
    @inbounds for i in 1:params.max_iter
        # Optimality test
        if check_optimality(lpip_pb, params.tol) == 1
            return Result{T}(lpip_pb.variables, status, i, time() - t0)
        end
        # Solve the system and fill newton directions
        solve_kkt(newton_system, lpip_pb, params)
        # Update x, s, p
        update_lpip_pb_vars(newton_system, lpip_pb, params)
        # Check if problem is ubounded or infeasible
        check_ilimited_or_unbounded() # Decide how to do this

        if check_time_limit(params.time_limit, t0) == 4
            return Result{T}(lpip_pb.variables, status, i, time() - t0)
        end
    end
    # Iteration limit
    return Result{T}(lpip_pb.variables, 5, params.max_iter, time() - t0)
end

function check_optimality(lpip_pb::LPIPLinearProblem{T}, tol::Float64) where T
    if abs(dot(lpip_pb.variables.x, lpip_pb.variables.s)) < tol
        return 1 # Optimal
    end
    return 0 # Not solved
end

function check_time_limit(time_limit::Float64, t0::Float64)
    if time() - t0 > time_limit
        return 4
    end
    return 0
end

function check_ilimited_or_unbounded()
    
end

function update_lpip_pb_vars(newton_system::NewtonSystem{T}, lpip_pb::LPIPLinearProblem{T}, params::Params) where T
    ## Find step lenghts
    beta_primal = 1
    beta_dual = 1

    for i in 1:(lpip_pb.m + lpip_pb.n)
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

    # Multiply betas to alpha
    beta_primal = params.alpha * beta_primal
    beta_dual = params.alpha * beta_dual
    # Update variables
    @. lpip_pb.variables.x = lpip_pb.variables.x + beta_primal * newton_system.d.dx
    @. lpip_pb.variables.p = lpip_pb.variables.p + beta_dual * newton_system.d.dp
    @. lpip_pb.variables.s = lpip_pb.variables.s + beta_dual * newton_system.d.ds
    return
end