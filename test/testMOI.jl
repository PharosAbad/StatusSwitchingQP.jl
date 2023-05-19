
"run test for MOI"
module TestMOI
#=
much thanks to [Oscar Dowson](https://github.com/odow), who provide this code, and help me walk through
=#

using Test

import MathOptInterface as MOI
import StatusSwitchingQP

function runtests()
    for name in names(@__MODULE__; all=true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_runtests()
    config = MOI.Test.Config(
        atol=1e-6,
        rtol=1e-6,
        optimal_status=MOI.OPTIMAL,
        exclude=Any[
            MOI.VariableBasisStatus,
            MOI.ConstraintBasisStatus,
            MOI.ObjectiveBound,
            MOI.ConstraintDual,
            MOI.DualObjectiveValue,
            #MOI.DualStatus,
        ],
    )
    model = MOI.instantiate(
        StatusSwitchingQP.Optimizer;
        with_bridge_type=Float64,
        with_cache_type=Float64
    )
    MOI.Test.runtests(
        model,
        config;
        exclude=String[
            "test_model_LowerBoundAlreadySet",
            "test_model_UpperBoundAlreadySet",
            #= "test_model_supports_constraint_ScalarAffineFunction_EqualTo",
            "test_constraint_ScalarAffineFunction_EqualTo",
            "test_modification_const_scalar_objective",
            "test_modification_const_vectoraffine_zeros",
            "test_objective_ObjectiveFunction",
            "test_objective_qp_ObjectiveFunction_edge_cases",
            "test_quadratic_nonhomogeneous",
            "test_solve_TerminationStatus_DUAL_INFEASIBLE",
            "test_solve_result_index",  =#
        ]
    )
    return
end

function test_SolverName()
    @test MOI.get(StatusSwitchingQP.Optimizer(), MOI.SolverName()) == "StatusSwitchingQP"
    return
end

end # module TestMOI

# This line at the end of the file runs all the tests!
TestMOI.runtests()

