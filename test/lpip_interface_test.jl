@testset "LPIP Interface" begin
    @testset "Production problem" begin
        A = [
            2. 1 
            1 2
        ]
        b = [4.;4]
        c = [-4.; -3]
        d = 3.

        lpip_pb = LPIPLinearProblem{Float64}(A, b, c, d, 0)
        params = LPIP.Params()
        lp = LPIP.interior_points(lpip_pb, params)

        @test lp.obj_val ≈ -6.3333 atol = 1e-4 rtol = 1e-4
        @test lp.dual_obj_val ≈ -6.3333 atol = 1e-4 rtol = 1e-4
        @test lp.variables.x ≈ [1.3333; 1.3333] atol = 1e-4 rtol = 1e-4
        @test lp.variables.p ≈ [-1.6666; -0.6666] atol = 1e-4 rtol = 1e-4
        @test lp.variables.s ≈ [0.0; 0.0] atol = 1e-4 rtol = 1e-4
        @test lp.status == 1
    end
end