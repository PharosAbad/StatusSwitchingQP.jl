"Status Switching Method for QP"
module SSQP

using LinearAlgebra
using StatusSwitchingQP: Status, IN, DN, UP, OE, EO, Event, LP, QP, getRowsGJ, getRowsGJr, Settings
export solveQP
using StatusSwitchingQP.Simplex: cDantzigLP, maxImprvLP, stpEdgeLP #, Simplex


@inline function polishSz!(S, z, d, u, G, g, N, J, tol)
    #for final result
    @inbounds for k in 1:N
        if S[k] == DN
            z[k] = d[k]
        elseif S[k] == UP
            z[k] = u[k]
        else    #IN
            if abs(z[k] - d[k]) < tol
                z[k] = d[k]
                S[k] = DN
            elseif abs(z[k] - u[k]) < tol
                z[k] = u[k]
                S[k] = UP
            end
        end
    end

    @inbounds for j in 1:J
        S[N+j] = abs(g[j] - z' * G[j, :]) < tol ? EO : OE
    end

end


@inline function freeK!(S, z, V, q, N, tol)  #for K=0
    #modify: S
    p = V * z + q
    S0 = copy(S)

    t = true   #hit optimal
    @inbounds for k in 1:N
        if (p[k] >= -tol && S[k] == UP) || (p[k] <= tol && S[k] == DN)
            #if (p[k] > tol && S[k] == UP) || (p[k] < -tol && S[k] == DN)
            S[k] = IN
            t = false
        end
    end
    if t
        return 1
    else
        ip = findall(S .== IN)
        if length(ip) > 0 && norm(p[ip], Inf) <= tol  #all movable are optimal
            S[ip] = S0[ip]  #restore the status
            return 1
        end
        return -1
    end

end

@inline function aStep!(p, z::Vector{T}, S, F, Og, alpha, G, g, d, u, fu, fd, N, J, tol) where {T}
    #compute step
    Lo = Vector{Event{T}}(undef, 0)
    ik = findall(F)
    @inbounds for k in eachindex(alpha)
        j = ik[k]
        t = p[k]
        h = z[j]
        dL = (d[j] - h) / t
        uL = (u[j] - h) / t
        if t > tol && fu[j]
            push!(Lo, Event{T}(IN, UP, j, uL))
        elseif t < -tol && fd[j]
            push!(Lo, Event{T}(IN, DN, j, dL))
        end
    end

    @inbounds if J > 0
        zo = g[Og] - G[Og, :] * z
        po = G[Og, F] * p
        ik = findall(Og)
        for k in eachindex(zo)
            j = ik[k]
            t = po[k]
            if t > tol
                push!(Lo, Event{T}(OE, EO, j, zo[k] / t))
            end
        end
    end

    L1::T = 1.0
    nL = length(Lo)
    if nL > 0
        sort!(Lo, by=x -> x.L)
        L1 = Lo[1].L
    end

    @inbounds if L1 < 1.0
        #adding blocking constraints
        z[F] .+= L1 * p
        #multi blocking
        for i in 1:lastindex(Lo)
            Lt = Lo[i]
            if Lt.L - L1 > tol
                break
            end
            k = Lt.id
            To = Lt.To
            if To == EO
                k += N
            end
            S[k] = To
            if k <= N
                z[k] = To == DN ? d[k] : u[k]
            end
        end
        #= ik = findall(F)
        for k in ik
            if abs(z[k] - d[k]) < tol
                z[k] = d[k]
                S[k] = DN
            elseif abs(z[k] - u[k]) < tol
                z[k] = u[k]
                S[k] = UP
            end
        end =#
        return -1
    else
        # if step size L1 == 1.0, maybe some z_i hit bounds, leave it for KKTchk
        z[F] = alpha
        return 1
    end

end

@inline function KKTchk!(S, F, B, Eg, gamma, alphaL, AE, GE, idAE, ra, N, M, tolG::T) where {T}
    ib = findall(B)
    Li = Vector{Event{T}}(undef, 0)
    @inbounds for k in eachindex(gamma)
        j = ib[k]
        t = gamma[k]
        if S[j] == UP && t > tolG
            push!(Li, Event{T}(UP, IN, j, -t))
        elseif S[j] == DN && t < -tolG
            push!(Li, Event{T}(DN, IN, j, t))
        end
    end

    JE = size(GE, 1)
    @inbounds if JE > 0
        iE = zeros(Int, JE)
        iR = findall(ra .> M)
        iE[idAE[ra[iR]]] = iR
        Lda = zeros(T, JE)
        for j in 1:JE
            k = iE[j]
            if k == 0
                x = AE' \ GE[j, F]
                Lda[j] = alphaL' * x
            else
                Lda[j] = alphaL[k]
            end
        end

        ib = findall(Eg)
        for k in 1:JE
            t = Lda[k]
            if t < -tolG
                push!(Li, Event{T}(EO, OE, ib[k], t))
            end
        end
    end

    nL = length(Li)
    @inbounds if nL > 0   #entering one only, the tighest one
        sort!(Li, by=x -> x.L)
        Lt = Li[1]
        k = Lt.id
        To = Lt.To
        if To == OE
            k += N
        end
        S[k] = To
        return -1
    else
        return 1
    end
end

"""

        solveQP(Q::QP{T}; settings, settingsLP) where T
        solveQP(V::Matrix{T}, q::Vector{T}, A, b, G, g, d, u; settings, settingsLP) where T
        solveQP(Q::QP{T}, S::Vector{Status}, x0; settings) where T

for quadratic programming problems: the initial feasible point (S, x0) if given,  can be obtained from SimplexLP

```math
    min   (1/2)z′Vz+q′z
    s.t.   Az=b ∈ R^{M}
           Gz≤g ∈ R^{J}
           d≤z≤u ∈ R^{N}
```

Outputs

    z               : solution,  N x 1 vector
    S               : Vector{Status}, (N+J)x1
    status          : > 0 if successful (=iter_count); = 0 if infeasibility detected; < 0 fail (=-1 if numerical error, =-maxIter if not converged)

See also [`QP`](@ref), [`SimplexLP`](@ref), [`Settings`](@ref), [`initQP`](@ref), [`initSSQP`](@ref)
"""
function solveQP(V::Matrix{T}, q::Vector{T}, A::Matrix{T}, b::Vector{T}, G::Matrix{T}, g::Vector{T},
    d::Vector{T}, u::Vector{T}; settings=Settings{T}(), settingsLP=settings) where {T}

    #N::Int32 = length(q)
    N = length(q)
    M = length(b)
    J = length(g)
    Q = QP(V, A, G, q, b, g, d, u, N, M, J)
    solveQP(Q; settings=settings, settingsLP=settingsLP)
end

function solveQP(Q::QP{T}; settings=Settings{T}(), settingsLP=settings) where {T}
    #x0, S, status = initSSQP(Q, settingsLP)
    if Q.mc <= 0
        return zeros(T, Q.N), fill(DN, Q.N), -1
    end
    x0, S, status = initQP(Q, settingsLP)
    if status <= 0  #infeasible or numerical error
        return x0, S, status
    end
    solveQP(Q, S, x0; settings=settings)
end


function solveQP(Q::QP{T}, S, x0; settings=Settings{T}()) where {T}
    #(; V, A, G, q, b, g, d, u, N, M, J) = Q
    V = Q.V
    A = Q.A
    G = Q.G
    q = Q.q
    b = Q.b
    g = Q.g
    d = Q.d
    u = Q.u
    N = Q.N
    M = Q.M
    J = Q.J

    #(; maxIter, tol, tolN, tolG) = settings
    #(; maxIter, tol, tolG, pivot) = settings
    #(; maxIter, tol, tolG) = settings
    maxIter = settings.maxIter
    tol = settings.tol
    tolG = settings.tolG

    #refineRows = pivot == :column ? getRowsGJ : getRowsGJr
    #refineRows = getRowsGJr

    fu = u .< Inf   #finite upper bound
    fd = d .> -Inf   #finite lower bound

    Sz = @view S[1:N]
    Se = @view S[(N.+(1:J))]
    z = copy(x0)

    iter = 0
    @inbounds while true
        #while true
        iter += 1
        if iter > maxIter
            return z, S, -iter
        end

        F = (Sz .== IN)
        K = sum(F)
        if K == 0
            status = freeK!(S, z, V, q, N, tol)
            if status > 0
                return z, S, iter
            else
                continue
            end
        end

        B = .!F
        Eg = (Se .== EO)
        Og = (Se .== OE)
        GE = @view G[Eg, :]
        AE = vcat(A[:, F], GE[:, F])
        idAE = vcat(axes(A, 1), axes(GE, 1)) # id of each constraint
        zB = z[B]
        AB = vcat(A[:, B], GE[:, B])
        bE = vcat(b, g[Eg]) - AB * zB

        #= ra = getRows(AE, tol)
        W = length(ra)
        if W < length(bE)
            rb = getRows([AE bE], tol)
            if W != length(rb)
                return z, S, 0    #infeasible
            end
            AE = AE[ra, :]
            bE = bE[ra]
            AB = AB[ra, :]
        end =#

        #ra, la = refineRows([AE bE], tol)
        ra, la = getRowsGJr([AE bE], tol)
        W = length(ra)
        if W < length(bE)
            if W != la
                return z, S, -1
            end
            AE = AE[ra, :]
            bE = bE[ra]
            AB = AB[ra, :]
        end


        iV = inv(cholesky(V[F, F]))
        VBF = V[B, F]
        c = VBF' * zB + q[F]
        mT = iV * AE'   #T=V_{I}⁻¹A_{I}′
        C = AE * mT
        C = (C + C') / 2
        C = inv(cholesky(C))
        TC = mT * C
        VQ = iV - mT * TC'    #Q=I-A_{I}′CT′   V_{I}⁻¹Q=V_{I}⁻¹-TCT′
        alpha = TC * bE - VQ * c    #α=TCb_{E}-V_{I}⁻¹Qc
        p = alpha - z[F]

        #direction p ≠ 0
        if norm(p, Inf) > tolG  #= IMPORTANT: in theory, norm(p, Inf)=0  == norm(p)=0, but in numerical computation, NO!  p = (e/2)*1,  norm(p) = (e/2)*sqrt(N) > e if N>4 =#
            #if norm(p) > tolN
            status = aStep!(p, z, S, F, Og, alpha, G, g, d, u, fu, fd, N, J, tol)
            if status < 0
                continue
            end
            #else
            #z[F] = alpha
        end

        #= if iter >= 61
            display(iter)
        end =#

        #direction p = 0
        #α_{λ}=-C(T′c+b_{E})
        alphaL = -(TC' * c + C * bE)
        gamma = VBF * alpha + V[B, B] * zB + q[B] + AB' * alphaL
        #z[F] = alpha
        status = KKTchk!(S, F, B, Eg, gamma, alphaL, AE, GE, idAE, ra, N, M, tolG)
        if status > 0
            #= ik = findall(F)
            for k in ik #check fake IN
                if abs(z[k] - d[k]) < tol
                    z[k] = d[k]
                    S[k] = DN
                elseif abs(z[k] - u[k]) < tol
                    z[k] = u[k]
                    S[k] = UP
                end
            end =#
            polishSz!(S, z, d, u, G, g, N, J, tol)


            #should we compute the final analytical z[I]?
            #z1 = copy(z)
            #z1[F] = alphaCal(F, z, Se, V, A, G, q, b, g, tol)
            #z[F] = (z[F]+z1[F])/2

            return z, S, iter
        end
    end
end



"""
        x0, S, status = initSSQP(Q::QP{T}, settingsLP)

performing Phase-I Simplex, do not handle free variables such that -∞ < x < +∞. and d should be finite. OK for EfficientFrontier
"""
function initSSQP(Q::QP{T}, settingsLP) where {T}
    #(; A, G, b, g, d, u, N, M, J) = Q
    A = Q.A
    G = Q.G
    b = Q.b
    g = Q.g
    d = Q.d
    u = Q.u
    N = Q.N
    M = Q.M
    J = Q.J

    #(; tol, rule) = settingsLP
    tol = settingsLP.tol
    rule = settingsLP.rule

    solveLP = cDantzigLP
    if rule == :maxImprovement
        solveLP = maxImprvLP
    elseif rule == :stpEdgeLP
        solveLP = stpEdgeLP
    end

    #An initial feasible point by performing Phase-I Simplex on the polyhedron
    Ms = M + J  #convert Gz<=g to equality
    Ns = N + J
    N1 = Ms + Ns
    S = fill(DN, N1)
    B = collect(Ns .+ (1:Ms))
    S[B] .= IN

    As = [A zeros(T, M, J)
        G Matrix{T}(I, J, J)]
    invB = Matrix{T}(I, Ms, Ms)
    bs = [b; g]
    ds = [d; zeros(T, J)]
    us = [u; fill(Inf, J)]

    q = As * ds
    for j in 1:Ms
        invB[j, j] = bs[j] >= q[j] ? one(T) : -one(T)
    end
    #q = abs.(As * ds - bs)
    q = abs.(q - bs)
    c1 = [zeros(T, Ns); fill(one(T), Ms)]   #灯塔的　模型　是　min
    A1 = [As invB]
    b1 = bs
    d1 = [ds; zeros(T, Ms)]
    u1 = [us; fill(Inf, Ms)]

    _iH, x, _invB = solveLP(c1, A1, b1, d1, u1, B, S; invB=invB, q=q, tol=tol)

    x0 = x[1:N]
    S = S[1:N+J]
    for k in N+1:N+J
        S[k] = S[k] == IN ? OE : EO
    end

    f = sum(x[Ns+1:end])
    if f > tol
        #error("feasible region is empty")
        return x0, S, 0
    end


    return x0, S, 1

end

"""
        x0, S, status = initQP(Q::QP{T}, settingsLP)

performing Phase-I Simplex on the polyhedron {Az=b, Gz≤g, d≤z≤u} to find an initial feasible point
allowing free variables such that -∞ < z < +∞
"""
function initQP(Q::QP{T}, settingsLP) where {T}
    #An initial feasible point by performing Phase-I Simplex on the polyhedron
    #(; A, G, b, g, d, u, N, M, J) = Q
    A = Q.A
    G = Q.G
    b = Q.b
    g = Q.g
    d = Q.d
    u = Q.u
    N = Q.N
    M = Q.M
    J = Q.J
    #(; tol, rule) = settingsLP
    tol = settingsLP.tol
    rule = settingsLP.rule

    solveLP = cDantzigLP
    if rule == :maxImprovement
        solveLP = maxImprvLP
    elseif rule == :stpEdgeLP
        solveLP = stpEdgeLP
    end

    #convert free variable: -∞ < x < +∞
    fu = u .== Inf   #no upper bound
    fd = d .== -Inf   #no lower bound
    fv = fu .& fd  #free variable
    iv = findall(fv)
    n = length(iv)
    id = findall(fd .& (.!fv))   # (-∞, u]

    #add slack variables for Gz<=g , and 2nd part of free variables
    M0 = M + J
    N0 = N + J + n
    A0 = [A zeros(T, M, J) -A[:, iv]
        G Matrix{T}(I, J, J) -G[:, iv]]
    b0 = [b; g]
    d0 = [d; zeros(T, J + n)]
    u0 = [u; fill(Inf, J + n)]


    #(-∞, +∞)  -> two copy of [0, +∞)
    d0[iv] .= 0     #the 1st part of free variable
    #u0[iv] .= Inf

    #(-∞, u]  -> [-u, +∞)
    d0[id] .= -u0[id]
    u0[id] .= Inf
    A0[:, id] .= -A0[:, id]

    N1 = M0 + N0
    S = fill(DN, N1)
    B = collect(N0 .+ (1:M0))
    S[B] .= IN

    invB = Matrix{T}(I, M0, M0)
    q = A0 * d0
    for j in 1:M0
        invB[j, j] = b0[j] >= q[j] ? one(T) : -one(T)
    end
    q = abs.(q - b0)
    c1 = [zeros(T, N0); fill(one(T), M0)]   #灯塔的　模型　是　min
    A1 = [A0 invB]
    b1 = b0
    d1 = [d0; zeros(T, M0)]
    u1 = [u0; fill(Inf, M0)]


    #f, x, _q, _invB, _iH = EfficientFrontier.Simplex.cDantzigLP(c1, A1, b1, d1, u1, B, S; invB=invB, q=q, tol=tol)
    _iH, x0, _invB = solveLP(c1, A1, b1, d1, u1, B, S; invB=invB, q=q, tol=tol)
    x = x0[1:N]
    S = S[1:N+J]
    f = sum(x0[N0+1:end])
    if f > tol
        #error("feasible region is empty")
        return x, S, 0   #0 infeasible
    end


    for k in N+1:N+J    #inequalities
        S[k] = S[k] == IN ? OE : EO
    end

    if n > 0    #free variable
        x[iv] .-= x0[N+J+1:N+J+n]
        S[iv] .= IN
    end

    m = length(id)
    if m > 0   # flip u d
        x[id] = -x[id]
        for k in 1:m
            #S[k] = S[k] == DN ? UP : DN
            if S[k] == DN
                S[k] == UP
            end
        end
    end
    return x, S, 1
end

end

