function interior_points(linear_problem::LinearProblem{T}, params::Params) where T
    
    lpip_pb = build_lpip_problem(linear_problem)
    newton_system = NewtonSystem{T}(lpip_pb)

    # Interior points iteration
    for i in 1 : params.max_iter
        # Optimality test
        check_optimality(lpip_pb, params.tol)
        if lpip_pb.status == 1 # Optimal solution
            break
        end
        
        # Build kkt matrixes
        fill_newton_system!(newton_system, lpip_pb)
        # Solve the system and fill newton directions
        solve_kkt(newton_system, lpip_pb)

        ## Find step lenghts
        minx = 10e10
        mins = 10e10

        # Find the minimum - x / dx
        for i = 1 : m + n
            if dx[i] <= 0
                aux = - x[i] / dx[i]
                if minx > aux
                    minx = aux
                end
            end
        end

        minx = α * minx

        # Find the minimum - s / ds
        for i = 1 : m + n
            if ds[i] <= 0
                aux = - s[i] / ds[i]
                if mins > aux
                    mins = aux
                end
            end
        end

        mins = α * mins

        # beta primal
        bp = min(1, minx)

        # beta dual
        bd = min(1, mins)

        ## Solution update
        x_new = x + bp * dx
        p_new = p + bd * dp
        s_new = s + bd * ds

        # Update new values
        x = x_new
        p = p_new
        s = s_new

        # Check if the problem is unlimited
    end
    return
end

function check_optimality(lpip_pb::LPIPLinearProblem{T}, tol::Float64)
    if dot(lpip_pb.varaibles.x, lpip_pb.varaibles.s) < tol
        lpip_pb.status = 1
    end
    return
end