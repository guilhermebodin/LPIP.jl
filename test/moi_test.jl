using MathOptInterface

const MOI = MathOptInterface
const MOIU = MOI.Utilities
const MOIT = MOI.Test
const MOIB = MOI.Bridges

const CI = MOI.ConstraintIndex
const VI = MOI.VariableIndex


MOIU.@model LPIPModelData () () (MOI.Zeros, MOI.Nonnegatives) () () (MOI.ScalarAffineFunction,) () (MOI.VectorAffineFunction,)

universal_fallback = MOIU.UniversalFallback(LPIPModelData{Float64}())
const optimizer = MOIU.CachingOptimizer(universal_fallback, LPIP.Optimizer())
# const optimizer = MOIU.CachingOptimizer(LPIPModelData{Float64}(), LPIP.Optimizer())
const config = MOIT.TestConfig(atol=1e-5, rtol=1e-5, infeas_certificates = false)

@testset "MOI tests" begin
    @testset "SolverName" begin
        @test MOI.get(optimizer, MOI.SolverName()) == "LPIP"
    end
    @testset "Unit" begin
        MOIT.unittest(MOIB.full_bridge_optimizer(optimizer, Float64), config,[
            # Quadratic functions are not supported
            "solve_qcp_edge_cases", "solve_qp_edge_cases",
            # Integer and ZeroOne sets are not supported
            "solve_integer_edge_cases", "solve_objbound_edge_cases"
            ]
        )
    end
    @testset "MOI Continuous Linear" begin
        MOIT.contlineartest(MOIB.full_bridge_optimizer(optimizer, Float64), config, [
            # infeasible/unbounded
            # "linear8a", "linear8b", "linear8c", "linear12",
            # linear10 is poorly conditioned
            "linear10",
            # linear9 is requires precision
            "linear9",
            # primalstart not accepted
            "partial_start"
            ]
        )
        # MOIT.linear9test(MOIB.SplitInterval{Float64}(optimizer_lin_hd), config)
    end
end

