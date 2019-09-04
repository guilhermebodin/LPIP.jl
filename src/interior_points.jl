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
const LPIP_DUAL_INFEASIBLE = 3
const LPIP_TIME_LIMIT = 4
const LPIP_ITERATION_LIMIT = 5

function interior_points(A::AbstractMatrix{T}, b::AbstractVector{T}, c::AbstractVector{T},
                         d::T, l::Int, params::Params) where T
    lpip_pb = LPIPLinearProblem{T}(A, b, c, d, l)
    return interior_points(lpip_pb, params)
end

function interior_points(lpip_pb::LPIPLinearProblem{T}, params::Params) where T
    print_header(params, lpip_pb)
    t0 = time() # Start the timer
    newton_system = NewtonSystem{T}(lpip_pb, params)
    # Interior points iteration
    for i in 1:params.max_iter
        # Optimality test
        check_optimality(lpip_pb, params) == LPIP_OPTIMAL && return Result(lpip_pb, LPIP_OPTIMAL, i, time() - t0)
        check_unbounded(lpip_pb, params) == LPIP_DUAL_INFEASIBLE && return Result(lpip_pb, LPIP_DUAL_INFEASIBLE, 0, time() - t0)
        check_infeasible(lpip_pb, params) == LPIP_INFEASIBLE && return Result(lpip_pb, LPIP_INFEASIBLE, 0, time() - t0)
        # Solve the system and fill newton directions
        solve_kkt(newton_system, lpip_pb, params)
        # Update x, s, p
        update_lpip_pb_vars!(lpip_pb.variables.x, lpip_pb.variables.s, lpip_pb.variables.p,
                             newton_system.d_result, params.alpha, lpip_pb.n, lpip_pb.m)
        # Check if time limit was reached
        check_time_limit(params.time_limit, t0) == LPIP_TIME_LIMIT && return Result(lpip_pb, LPIP_TIME_LIMIT, i, time() - t0)
        # Print result of the iteration
        print_evolution(params, lpip_pb, i, t0)
    end
    # Iteration limit
    return Result(lpip_pb, LPIP_ITERATION_LIMIT, params.max_iter, time() - t0)
end

function check_optimality(lpip_pb::LPIPLinearProblem{T}, params::Params) where T
    if abs(dot(lpip_pb.variables.x, lpip_pb.variables.s)) < params.tol
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

function check_infeasible(lpip_pb::LPIPLinearProblem{T}, params::Params) where T
    if abs(dot(lpip_pb.b, lpip_pb.variables.p)) >= params.infeasible_unbounded_tol
        return LPIP_INFEASIBLE
    end
    return LPIP_NOT_SOLVED
end

function check_unbounded(lpip_pb::LPIPLinearProblem{T}, params::Params) where T
    if abs(dot(lpip_pb.c, lpip_pb.variables.x)) >= params.infeasible_unbounded_tol
        return LPIP_DUAL_INFEASIBLE
    end
    return LPIP_NOT_SOLVED
end

function update_lpip_pb_vars!(x::AbstractVector{T}, s::AbstractVector{T}, p::AbstractVector{T},
                              d::AbstractVector{T}, alpha::Float64, n::Int, m::Int) where T
    # Find step lenghts
    beta_primal = 1
    beta_dual = 1

    @inbounds for i in 1:n
        if d[i] <= 0
            if beta_primal > - x[i] / d[i]
                beta_primal = - x[i] / d[i]
            end
        end
        if d[n + m + i] <= 0
            if beta_dual > - s[i] / d[n + m + i]
                beta_dual = - s[i] / d[n + m + i]
            end
        end
    end

    alpha_beta_primal = alpha * beta_primal
    alpha_beta_dual = alpha * beta_dual

    @inbounds  for i in 1:n
        x[i] += alpha_beta_primal*d[i]
        s[i] += alpha_beta_dual*d[n + m + i]
    end
    @inbounds for i in 1:m
        p[i] += alpha_beta_dual*d[n + i]
    end
    return
end