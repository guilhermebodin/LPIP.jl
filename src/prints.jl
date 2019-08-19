function print_header(verbose::Bool)
    if verbose 
        println("===================================================")
        println(" LPIP.jl v0.1.0  (c) Guilherme Bodin,  2019        ")
        # Problem with n variables x equality constraints and k non-negative constraints
        println("||    obj    |    gap    |   iter    |   time    ||")
    end
end