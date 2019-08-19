function print_header(params::Params, lpip_pb::LPIPLinearProblem{T}) where T
    if params.verbose 
        println("                                                   ")
        println(" LPIP.jl v0.1.0  (c) Guilherme Bodin,  2019        ")
        println("                                                   ")
        println("Problem with: ", lpip_pb.n - lpip_pb.l, " variables,")
        println("              ", lpip_pb.f, " equality constraints,")
        println("              ", lpip_pb.l, " nonnegative constraints.")
        println("                                                   ")
        println("Parameters:  rho = ", params.rho, " alpha = ", params.alpha, " tol = ", params.tol)
        println("             max_iter = ", params.max_iter, " time limit (s) = ", params.time_limit)
        println("                                                   ")
        println("||     obj    |    gap    | iter |   time (s) ||")
    end
end

function print_evolution(params::Params, lpip_pb::LPIPLinearProblem{T}, iter::Int, t0::Float64) where T
    if params.verbose
        t1 = time() - t0
        obj = dot(lpip_pb.c, lpip_pb.variables.x) + lpip_pb.d
        gap = abs(dot(lpip_pb.variables.x, lpip_pb.variables.s))

        str_obj = " " * @sprintf("%.6f", obj) * "  |"
        str_gap   = " " * @sprintf("%.6f", gap) * "  |"
        str_iter = @sprintf("%d", iter) * " |"
        str_time   = "   " * @sprintf("%.4f", t1) * "   ||"

        str = "||"
        str *= " "^max(0, 11 - length(str_obj)) * str_obj
        str *= " "^max(0, 11 - length(str_gap)) * str_gap
        str *= " "^max(0, 7 - length(str_iter)) * str_iter
        str *= " "^max(0, 12 - length(str_time)) * str_time

        println(str)
    end
end