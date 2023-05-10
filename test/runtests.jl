using StatusSwitchingQP
using Test

@testset "StatusSwitchingQP.jl" begin

    #SSLP, unbounded
    c = [-3.0, -2]
    d = [0.0, 0]
    u = [Inf, Inf]
    G = [-1.0 3; 1 -5]
    g = [12.0; 5]
    A = zeros(0, 2)
    b = zeros(0)

    Q = LP(c, A, b; d=d, u=u, G=G, g=g)
    x, S, status = solveLP(Q)
    @test status == 3   #infinitely many solutions



    #SSQP
    V = [1/100 1/80 1/100
        1/80 1/16 1/40
        1/100 1/40 1/25]
    E = [109 / 100; 23 / 20; 119 / 100]

    up = [0.7; +Inf; 0.7]     #Inf means no bounded
    P = QP(QP(V; u=up), E)
    z, Sp, iter = solveQP(P)
    @test Status[UP, IN, IN] == Sp

end
