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
        ],
    )
    model = MOI.instantiate(
        StatusSwitchingQP.Optimizer;
        with_bridge_type = Float64,
        with_cache_type = Float64,
    )
    MOI.Test.runtests(model, config)
    return
end

function test_SolverName()
    @test MOI.get(StatusSwitchingQP.Optimizer(), MOI.SolverName()) == "StatusSwitchingQP"
    return
end

end # module TestMOI

# This line at the end of the file runs all the tests!
TestMOI.runtests()
