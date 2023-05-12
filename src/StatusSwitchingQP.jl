"Status Switching Method for Quadratic Programming and Linear Programming"
module StatusSwitchingQP
using LinearAlgebra
using PrecompileTools
#PrecompileTools.verbose[] = true


include("./types.jl")
export Status, IN, DN, UP, OE, EO, Event
export LP, QP, Settings

include("./utils.jl")
export getRowsGJ, getRowsGJr

include("./Simplex.jl")
using .Simplex
export SimplexLP

include("./SSQP.jl")
using .SSQP
export solveQP

include("./SSLP.jl")
using .SSLP
export solveLP



# #=
@setup_workload begin
    # Putting some things in `@setup_workload` instead of `@compile_workload` can reduce the size of the
    # precompile file and potentially make loading faster
    V = [1/100 1/80 1/100
        1/80 1/16 1/40
        1/100 1/40 1/25]
    #E = [109 / 100; 23 / 20; 119 / 100]

    up = [0.7; +Inf; 0.7]     #Inf means no bounded
    #P0 = QP(V; u=up)
    #P = QP(QP(V; u=up), E)


    c = [-3.0, -2]
    #d = [0.0, 0]
    #u = [Inf, Inf]
    G = [-1.0 3; 1 -5]
    g = [12.0; 5]
    A = zeros(0, 2)
    b = zeros(0)


    @compile_workload begin
        # all calls in this block will be precompiled, regardless of whether
        # they belong to your package or not (on Julia 1.8 and higher)
        #P = QP(P0, E, 0.0)
        P = QP(V; u=up)
        z, S, iter = solveQP(P)

        Q = LP(c, A, b; G=G, g=g)
        #Q = LP(c, A, b; d=d, u=u, G=G, g=g)
        #sol = solveLP(Q)
        res = SimplexLP(Q)
    end
end
# =#

end
