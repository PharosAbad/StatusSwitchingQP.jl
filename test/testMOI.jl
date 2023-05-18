#=
# ============================ /test/MOI_wrapper.jl ============================
module TestMOI

import StatusSwitchingQP
using Test

import MathOptInterface as MOI

#=
const _EXPLICIT_METHOD_FAILURES = [
    "test_objective_qp_ObjectiveFunction_edge_cases",
    "test_objective_qp_ObjectiveFunction_zero_ofdiag",
    "test_quadratic_duplicate_terms",
    "test_quadratic_integration",
    "test_quadratic_nonhomogeneous",
    "test_linear_Semicontinuous_integration",
    "test_linear_Semiinteger_integration",
    "test_attribute_RawStatusString",
    "test_attribute_SolveTimeSec",
    "test_constraint_ScalarAffineFunction_EqualTo",
    "test_constraint_ScalarAffineFunction_GreaterThan",
    "test_constraint_ScalarAffineFunction_Interval",
    "test_constraint_ScalarAffineFunction_LessThan",
    "test_constraint_ScalarAffineFunction_duplicate",
    "test_constraint_VectorAffineFunction_duplicate",
    "test_linear_FEASIBILITY_SENSE",
    "test_linear_VariablePrimalStart_partial",
    "test_model_supports_constraint_ScalarAffineFunction_EqualTo",
    "test_modification_affine_deletion_edge_cases",
    "test_modification_coef_scalar_objective",
    "test_modification_coef_scalaraffine_lessthan",
    "test_modification_const_scalar_objective",
    "test_modification_const_vectoraffine_nonpos",
    "test_modification_delete_variable_with_single_variable_obj",
    "test_modification_delete_variables_in_a_batch",
    "test_modification_func_scalaraffine_lessthan",
    "test_modification_func_vectoraffine_nonneg",
    "test_modification_multirow_vectoraffine_nonpos",
    "test_modification_set_scalaraffine_lessthan",
    "test_modification_set_singlevariable_lessthan",
    "test_modification_transform_singlevariable_lessthan",
    "test_objective_FEASIBILITY_SENSE_clears_objective",
    "test_objective_ObjectiveFunction_VariableIndex",
    "test_objective_ObjectiveFunction_blank",
    "test_objective_ObjectiveFunction_constant",
    "test_objective_ObjectiveFunction_duplicate_terms",
    "test_objective_ObjectiveSense_FEASIBILITY_SENSE",
    "test_solve_TerminationStatus_DUAL_INFEASIBLE",
    "test_solve_result_index",
    "test_variable_solve_with_lowerbound",
    "test_variable_solve_with_upperbound",
] =#

const _EXPLICIT_METHOD_FAILURES = [
#=    "test_attribute_RawStatusString",
    "test_attribute_SolveTimeSec",
    "test_constraint_ScalarAffineFunction",
    "test_constraint_VectorAffineFunction_duplicate",
    "test_linear_FEASIBILITY_SENSE",
    "test_linear_VariablePrimalStart_partial",
    "test_linear_complex",
    #"test_linear_open_intervals",
    "test_linear_variable", =#
    "test_model_LowerBoundAlreadySet",
    "test_model_UpperBoundAlreadySet",
    "test_model_supports_constraint_ScalarAffineFunction_EqualTo",
#=    "test_modification_",
    "test_objective_FEASIBILITY_SENSE_clears_objective",
    "test_objective_ObjectiveFunction",
    "test_objective_ObjectiveSense_FEASIBILITY_SENSE",
    "test_objective_qp_",
    "test_quadratic_duplicate_terms",
    "test_quadratic_integration",
    "test_quadratic_nonhomogeneous",
    "test_solve_VariableIndex",
    "test_solve_TerminationStatus_DUAL_INFEASIBLE",
    "test_solve_result_index",
    "test_unbounded_M",
    "test_variable_solve_with_",    =#
]


const OPTIMIZER = MOI.instantiate(
    MOI.OptimizerWithAttributes(StatusSwitchingQP.Optimizer, MOI.Silent() => true),
)

const BRIDGED = MOI.instantiate(
    MOI.OptimizerWithAttributes(StatusSwitchingQP.Optimizer, MOI.Silent() => true),
    with_bridge_type=Float64,
)

# See the docstring of MOI.Test.Config for other arguments.
const CONFIG = MOI.Test.Config(
    # Modify tolerances as necessary.
    atol=1e-6,
    rtol=1e-6,
    # Use MOI.LOCALLY_SOLVED for local solvers.
    optimal_status=MOI.OPTIMAL,
    # Pass attributes or MOI functions to `exclude` to skip tests that
    # rely on this functionality.
    exclude=Any[MOI.VariableName, MOI.delete,
        MOI.ConstraintName,
        MOI.VariableBasisStatus,
        MOI.ConstraintBasisStatus,
        MOI.ObjectiveBound],
)

"""
    runtests()

This function runs all functions in the this Module starting with `test_`.
"""
function runtests()
    for name in names(@__MODULE__; all=true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
end

"""
    test_runtests()

This function runs all the tests in MathOptInterface.Test.

Pass arguments to `exclude` to skip tests for functionality that is not
implemented or that your solver doesn't support.
"""
function test_runtests()
    MOI.Test.runtests(
        BRIDGED,
        CONFIG,
        exclude=[
            "test_attribute_NumberOfThreads",
            "test_quadratic_",
            _EXPLICIT_METHOD_FAILURES...,
        ],
        # This argument is useful to prevent tests from failing on future
        # releases of MOI that add new tests. Don't let this number get too far
        # behind the current MOI release though. You should periodically check
        # for new tests to fix bugs and implement new features.
        exclude_tests_after=v"0.10.5",
    )
    return
end

"""
    test_SolverName()

You can also write new tests for solver-specific functionality. Write each new
test as a function with a name beginning with `test_`.
"""
function test_SolverName()
    @test MOI.get(StatusSwitchingQP.Optimizer(), MOI.SolverName()) == "StatusSwitchingQP"
    return
end

end # module TestMOI

# This line at the end of the file runs all the tests!
TestMOI.runtests()

=#




module TestMOI

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
            "test_model_supports_constraint_ScalarAffineFunction_EqualTo",
            "test_constraint_ScalarAffineFunction_EqualTo",
            "test_modification_const_scalar_objective",
            "test_modification_const_vectoraffine_zeros",
            "test_objective_ObjectiveFunction",
            "test_objective_qp_ObjectiveFunction_edge_cases",
            "test_quadratic_nonhomogeneous",
            "test_solve_TerminationStatus_DUAL_INFEASIBLE",
            "test_solve_result_index",
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

