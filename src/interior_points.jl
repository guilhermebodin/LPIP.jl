"""
 0 - NotSolved
 1 - Optimal
 2 - Infeasible
 3 - Unbounded
 4 - Time limit
 5 - Iteration limit
"""
const LPIP_NOT_SOLVED = 0
const LPIP_OPTIMAL = 1
const LPIP_INFEASIBLE = 2
const LPIP_UNBOUNDED = 3
const LPIP_TIME_LIMIT = 4
const LPIP_ITERATION_LIMIT = 5

function interior_points(A::AbstractMatrix{T}, b::AbstractVector{T}, c::AbstractVector{T},
                         d::T, l::Int, params::Params) where T
    lpip_pb = LPIPLinearProblem{T}(A, b, c, d, l)
    return interior_points(lpip_pb, params)
end

function interior_points(lpip_pb::LPIPLinearProblem{T}, params::Params) where T
    print_header(params.verbose)
    t0 = time() # Start the timer
    check_unbounded(lpip_pb) == LPIP_UNBOUNDED && return Result(lpip_pb, LPIP_UNBOUNDED, 0, time() - t0)
    check_infeasible(lpip_pb) == LPIP_INFEASIBLE && return Result(lpip_pb, LPIP_INFEASIBLE, 0, time() - t0)
    newton_system = NewtonSystem{T}(lpip_pb, params)
    # Interior points iteration
    @inbounds for i in 1:params.max_iter
        # Optimality test
        check_optimality(lpip_pb, params.tol) == LPIP_OPTIMAL && return Result(lpip_pb, LPIP_OPTIMAL, i, time() - t0)
        # Solve the system and fill newton directions
        solve_kkt(newton_system, lpip_pb, params)
        # Update x, s, p
        update_lpip_pb_vars(newton_system, lpip_pb, params)
        # Check if time limit was reached
        check_time_limit(params.time_limit, t0) == LPIP_TIME_LIMIT && return Result(lpip_pb, LPIP_TIME_LIMIT, i, time() - t0)
    end
    # Iteration limit
    return Result(lpip_pb, LPIP_ITERATION_LIMIT, params.max_iter, time() - t0)
end

function check_optimality(lpip_pb::LPIPLinearProblem{T}, tol::Float64) where T
    if abs(dot(lpip_pb.variables.x, lpip_pb.variables.s)) < tol
        return LPIP_OPTIMAL
    end
    return LPIP_NOT_SOLVED
end

function check_time_limit(time_limit::Float64, t0::Float64)
    if time() - t0 > time_limit
        return LPIP_TIME_LIMIT # Time limit
    end
    return LPIP_NOT_SOLVED
end

function check_infeasible(lpip_pb::LPIPLinearProblem{T}) where T
    # Detecting Infeasibility in Infeasible-Interior-Point Methods for Optimization
    # solve A'y + s = 0
    # if s is all positive
    # if b'y is positive
    # infeasible
    # We might want to set s = ones beforehand
end

function check_unbounded(lpip_pb::LPIPLinearProblem{T}) where T
    # Detecting Infeasibility in Infeasible-Interior-Point Methods for Optimization
    # solve Ax = 0
    # if x is all positive
    # if c'x is negative
    # unbounded
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