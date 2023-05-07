"Status Switching Method for LP"
module SSLP

# ------------ 2023-04-07 19:52:11
#考虑允许自由变量 (原始态, 无需转化成 两正变量的差)

using LinearAlgebra
using StatusSwitchingQP:  Status, IN, DN, UP, OE, EO, Event

export auxLP, solveLP


"""

        Settings(P::Problem; kwargs...)        The default Settings to given Problem
        Settings(; kwargs...)       The default Settings is set by Float64 type
        Settings{T<:AbstractFloat}(; kwargs...)

kwargs are from the fields of Settings{T<:AbstractFloat} for Float64 and BigFloat

            tol::T          #2^-26 ≈ 1.5e-8  general scalar
            pivot::Symbol    #pivot for purging redundant rows {:column, :row}

"""
struct Settings{T<:AbstractFloat}
    maxIter::Int    #7777
    tol::T          #2^-26
    tolN::T         #2^-26
    tolG::T         #2^-27 for Greeks (beta and gamma)
    pivot::Symbol    #column pivoting
end

Settings(; kwargs...) = Settings{Float64}(; kwargs...)

function Settings{Float64}(; maxIter=7777,
    tol=2^-26,
    tolN=2^-26,
    tolG=2^-27,
    pivot=:column)
    Settings{Float64}(maxIter, tol, tolN, tolG, pivot)
end

function Settings{BigFloat}(; maxIter=7777,
    tol=BigFloat(2)^-76,
    tolN=BigFloat(2)^-76,
    tolG=BigFloat(2)^-77,
    pivot=:column)
    Settings{BigFloat}(maxIter, tol, tolN, tolG, pivot)
end


function getRowsGJ(X::Matrix{T}, tol=eps(norm(X, Inf))) where {T}
    #Gauss-Jordan elimination, code form rref_with_pivots!
    A = copy(X)
    nr, nc = size(A)
    rows = Vector{Int64}()
    r0 = collect(1:nr)  #original row numb
    nc1 = nc - 1
    l1 = 0  #without the last col
    i = j = 1
    while i <= nr && j <= nc
        (m, mi) = findmax(abs.(A[r0[i:nr], j]))
        mi = mi + i - 1
        if m <= tol
            if j == nc1
                l1 = length(rows)
            end
            j += 1
        else
            append!(rows, r0[mi])
            r0[mi], r0[i] = r0[i], r0[mi]   #keep tracting, not move the data
            n = r0[i]
            d = A[n, j]
            for k = j:nc
                A[n, k] /= d
            end
            for k in 1:nr
                if k != n
                    d = A[k, j]
                    for l = j:nc
                        A[k, l] -= d * A[n, l]
                    end
                end
            end
            if j == nc1
                l1 = length(rows)
            end
            i += 1
            j += 1
        end
    end
    return rows, l1

end

function getRowsGJr(X::Matrix{T}, tol=eps(norm(X, Inf))) where {T}      # row poviting
    #Gauss-Jordan elimination, code form rref_with_pivots!
    A = copy(X)
    nr, nc = size(A)
    rows = Vector{Int64}()
    c0 = collect(1:nc)  #original col numb
    nc1 = nc - 1
    l1 = 0  #without the last col
    i = j = 1
    while i <= nr && j <= nc
        (m, mj) = findmax(abs.(A[i, c0[j:nc]]))
        mj = mj + j - 1
        if m <= tol
            i += 1
        else
            append!(rows, i)
            c0[mj], c0[j] = c0[j], c0[mj]   #keep tracting, not move the data
            n = c0[j]
            d = A[i, n]
            for k = c0[j:nc]
                A[i, k] /= d
            end
            for k in 1:nr
                if k != i
                    d = A[k, n]
                    for l in c0[j:nc]
                        A[k, l] -= d * A[i, l]
                    end
                end
            end
            l1 = j
            i += 1
            j += 1
        end
    end
    return rows, l1

end

function freeK!(S, U, D, c, N, tol)  #for K=0
    #modify: S
    t = true   #hit optimal
    for k in 1:N
        if (c[k] >= -tol && U[k]) || (c[k] <= tol && D[k])
        #if (c[k] > tol && U[k]) || (c[k] < -tol && D[k])
            S[k] = IN
            t = false
        end
    end
    return t ? 1 : -1

end

function freeK1!(S, U, D, c, N, tol)  #for K=0
    #modify: S
    S0 = copy(S)
    t = true   #hit optimal
    for k in 1:N
        if (c[k] >= -tol && U[k]) || (c[k] <= tol && D[k])
            #if (c[k] > tol && U[k]) || (c[k] < -tol && D[k])
            S[k] = IN
            t = false
        end
    end
    #return t ? 1 : -1
    if t
        return 1
    else
        ip = findall(S .== IN)
        if length(ip) > 0 && norm(c[ip], Inf) <= tol  #all movable are optimal
            S[ip] = S0[ip]  #restore the status
            return 1    #to do: chance to be infinitely many solutions, or let it handle by KKTchk?
        end
        return -1
    end

end

function aStep!(S, x, p, F, Og, c::Vector{T}, G, g, d, u, tol) where {T}  #for K>0 (allow W = 0)
    #function aStep!(S, x, p, F, Og, c::Vector{T}, G, g, d, u, fu, fd, tol) where {T}  #for K>0 (allow W = 0)
    #modify: S, x
    N = length(c)
    J = length(g)
    Lo = Vector{Event{T}}(undef, 0)
    s::T = 0.0

    ik = findall(F)
    for k in eachindex(p)   #p: Kx1
        t = p[k]
        j = ik[k]
        if t < -tol #&& fd[j]
            push!(Lo, Event{T}(IN, DN, j, (d[j] - x[j]) / t))
            #= if t < -tol
                if fd[j]
                    push!(Lo, Event{T}(IN, DN, j, (d[j] - x[j]) / t))
                else
                    push!(Lo, Event{T}(IN, IN, j, -Inf))
                end =#
        elseif t > tol #&& fu[j]
            push!(Lo, Event{T}(IN, UP, j, (u[j] - x[j]) / t))
            #= elseif t > tol
                if fu[j]
                    push!(Lo, Event{T}(IN, UP, j, (u[j] - x[j]) / t))
                else
                    push!(Lo, Event{T}(IN, IN, j, Inf))
                end =#
        end
    end

    if J > 0
        zo = g[Og] - G[Og, :] * x
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

    nL = length(Lo)
    if nL > 0
        sort!(Lo, by=x -> x.L)
        #= Lt = Lo[1]
        s = Lt.L
        Ls = [Lt]   =#
        s = Lo[1].L
    else
        return 1    #W=0? or Inf many sol ? or return -1 ?
        #return -1
    end

    if isinf(s)
        return 3
    end

    x[F] .+= s * p

    for i in 1:lastindex(Lo)
        Lt = Lo[i]
        #=k = Lt.id
        To = Lt.To
        if Lt.L - s < tol
            if To == EO
                k += N
            end
            S[k] = To
            if k <= N
                x[k] = To == DN ? d[k] : u[k]
            end
        else
            break
        end =#
        if Lt.L - s > tol
            break
        end
        k = Lt.id
        To = Lt.To
        if To == EO
            k += N
        end
        S[k] = To
        if k <= N
            x[k] = To == DN ? d[k] : u[k]
        end
    end

    return -1
end

function KKTchk!(S, F, B, U, hW, AB, Eg, c::Vector{T}, AE, GE, idAE, ra, M, JE, tolG) where {T}

    N = length(c)
    Li = Vector{Event{T}}(undef, 0)
    ib = findall(B)
    nb = length(ib)
    if nb > 0
        #if length(ib) > 0
        hB = c[B] + AB' * hW
        for k in eachindex(hB)
            j = ib[k]
            t = hB[k]
            if U[j]
                if t > tolG
                    push!(Li, Event{T}(UP, IN, j, -t))
                end
            else
                if t < -tolG
                    push!(Li, Event{T}(DN, IN, j, t))
                end
            end
        end
    end

    if JE > 0
        iE = zeros(Int, JE)
        iR = findall(ra .> M)
        iE[idAE[ra[iR]]] = iR
        Lda = zeros(T, JE)
        for j in 1:JE
            k = iE[j]
            if k == 0
                x = AE' \ GE[j, F]
                Lda[j] = hW' * x
            else
                Lda[j] = hW[k]
            end
        end

        ie = findall(Eg)
        for k in 1:JE
            t = Lda[k]
            if t < -tolG
                push!(Li, Event{T}(EO, OE, ie[k], t))
            end
        end
    end

    nL = length(Li)
    if nL > 0
        #sort!(Li, by=x -> x.L, rev=true)
        sort!(Li, by=x -> x.L)  #argmin
        Lt = Li[1]
        k = Lt.id
        To = Lt.To
        if To == OE
            k += N
        end
        S[k] = To
    else    #to do: indicate infinitely many solutions
        #= if any(abs.(hB) .< tol)      #when all IN, no hB
            return 2    #infinitely many solutions
        end =#
        if nb > 0 && any(abs.(hB) .< tolG)
            return 2
        end
        return 1    #hit optimal
    end
    return -1
end



function solveLP(Q::LP{T}; nS=Settings{T}()) where {T}
    (; c, A, b, G, g, d, u, M, J) = Q

    if J == 0 && M == 0 #a box
        #display("W = 0, K = 0, J = 0, a box")
        return boxLP(c, d, u, N; tol=nS.tol)
    end

    #S, x0 = initLP(Q, nS)
    S, x0 = auxLP(Q, nS)
    x, f, status = solveLP(c, A, b, G, g, d, u, S, x0; nS=nS)
    return S, x, f, status  #status:  1 unique; 0 infeasible; 2 infinitely many sol; 3 unbounded ; -1 in process
end

function solveLP(Q::LP{T}, S, x0; nS=Settings{T}()) where {T}
    (; c, A, b, G, g, d, u) = Q

    x, f, status = solveLP(c, A, b, G, g, d, u, S, x0; nS=nS)
    return S, x, f, status  #status:  1 unique; 0 infeasible; 2 infinitely many sol; 3 unbounded ; -1 in process
end

function solveLP(c::Vector{T}, A, b, G, g, d, u, S, x0; nS=Settings{T}()) where {T}
    #S is updated

    (; maxIter, tol, tolG, pivot) = nS
    N = length(c)

    refineRows = pivot == :column ? getRowsGJ : getRowsGJr

    M = length(b)
    J = length(g)
    if J == 0 && M == 0 #a box
        Sb, x, f, status = boxLP(c, d, u, N; tol=tol)
        S .= Sb
        return x, f, status
    end

    #fu = u .< Inf   #finite upper bound
    #fd = d .> -Inf   #finite lower bound

    Sx = @view S[1:N]
    Se = @view S[N+1:N+J]
    x = copy(x0)

    #s::T = 1.0  #stepsize
    iter = 0
    iCal = 0
    #@inbounds     while true
    while true
        iter += 1
        if iter > maxIter
            f = c' * x
            status = -iter
            return x, f, status
        end

        F = (Sx .== IN) #free variable are always IN
        K = sum(F)


        #= if mod(iter, 10) == 0
            display(iter)
            #display( (iter, -c'*x) )
        end =#

        #= if iter >= 12
            display(("check ", iter))
        end =#

        U = (Sx .== UP)
        D = (Sx .== DN)

        if K == 0
            status = freeK!(S, U, D, c, N, tol)
            if status >= 0
                f = c' * x
                return x, f, status
            else
                continue
            end
        end

        # B = U .|| D     #some free variable are not bounded
        B = .!F

        Eg = (Se .== EO)
        Og = (Se .== OE)
        GE = @view G[Eg, :]
        AE = vcat(A[:, F], GE[:, F])
        idAE = vcat(axes(A, 1), axes(GE, 1)) # id of each constraint
        JE = size(GE, 1)

        x[D] = d[D]
        x[U] = u[U]
        #xRefine!(S, x, d, u, N, tol)

        zB = x[B]
        AB = vcat(A[:, B], GE[:, B])
        bE = vcat(b, g[Eg]) - AB * zB

        #ra1 = getRows([AE bE], tol)

        # #=
        ra, la = refineRows([AE bE], tol)
        W = length(ra)

        #display((length(ra1), length(ra), K, W, iter))

        if W < length(bE)
            if W != la
                error("infeasible")
            end
            AE = AE[ra, :]
            bE = bE[ra]
            AB = AB[ra, :]
        end # =#


        # #=
        if K > 1 && K == W && iter - iCal > 77
            #x[F] = inv(AE) * bE
            #x[F] = AE\bE
            x[F] = inv(lu(AE)) * bE
            #display(("--- refine ", iter))
            iCal = iter
        end # =#

        if W == 0   # W = 0, K > 0
            #display("W = 0, K > 0, J > 0")
            p = -c[F]
        else    # W > 0, K > 0
            #working set, no chance for redundant?      for the feasible set is convex
            #有可能，但那只通过 单个顶点， 与对面的平行 2023-04-10 14:34:57

            cF = @view c[F]
            #=
            C = AE * AE'
            invC = inv(cholesky(C))
            #        invC = pinv(C)
            #p = AE' * (AE' \ cF) -cF #much slow
            p = AE' * invC * AE * cF - cF
            =#

            # #=
            hW = AE' \ cF   #
            p = AE' * hW - cF
            # =#

            #=
            #infeasible: grow7.mps, seba.mps
            C = AE * AE'
            C = (C + C') / 2
            invC = inv(cholesky(C))
            hW = invC * AE * cF
            p = AE' * hW - cF
            =#
        end

        #direction p ≠ 0
        if norm(p, Inf) > tolG
            #if norm(p, Inf) > tolG * ( abs( log2(cond(invC))) +1)  #when M>100, compute invC is a heavy burden and not stable  2023-04-29 22:32:18
            #status = aStep!(S, x, p, F, Og, c, G, g, d, u, fu, fd, tol)
            status = aStep!(S, x, p, F, Og, c, G, g, d, u, tol)
            #display(("--- ", iter, norm(A*x-b, Inf)))
            if status >= 0
                f = c' * x
                return x, f, status
            else
                continue
            end
        end


        #direction p = 0
        #hW = -(invC * AE * cF)
        #status = KKTchk!(S, F, B, U, hW, AB, Eg, c, AE, GE, idAE, ra, M, JE, tolG)
        status = KKTchk!(S, F, B, U, -hW, AB, Eg, c, AE, GE, idAE, ra, M, JE, tolG)
        if status >= 0
            f = c' * x
            #display((W, K))
            if W < K
                status = 2  #infinitely many solutions
            end
            return x, f, status
        end

    end
end


function initLP(Q::LP{T}, nS) where {T}
    #An initial feasible point by performing Phase-I Simplex on the polyhedron
    #(; c, A, b, G, g, d, u, N, M, J) = Q
    (; A, b, G, g, d, u, N, M, J) = Q
    #(; maxIter, tol, tolN) = nS
    tol = nS.tol

    #convert free variable: -∞ < x < +∞
    fu = u .== Inf   #no upper bound
    fd = d .== -Inf   #no lower bound
    fv = fu .&& fd  #free variable
    iv = findall(fv)
    n = length(iv)
    id = findall(fd .&& .!fv)   # (-∞, u]
    #display(id)

    #add slack variables for Gz<=g , and 2nd part of free variables
    M0 = M + J
    N0 = N + J + n
    A0 = [A zeros(T, M, J) -A[:, iv]
        G Matrix{T}(I, J, J) -G[:, iv]]
    b0 = [b; g]
    d0 = [d; zeros(T, J + n)]
    u0 = [u; fill(Inf, J + n)]

    #= c0 = [c; zeros(T, J); -c[iv]]   #not needed
    c0[id] .= -c0[id]   =#
    #display(c0')

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


    #to do: for SSLP, some q may not be IN

    G1 = zeros(T, 0, N1)
    g1 = zeros(T, 0)
    x1 = copy(d1)
    x1[B] = q
    x0, f, status = solveLP(c1, A1, b1, G1, g1, d1, u1, S, x1; nS=nS)

    #S0 = copy(S)
    #display((x0, S0, f, status))


    if abs(f) > tol
        #display(f)
        error("feasible region is empty")
    end

    #variables: restore free, or flip back
    x = x0[1:N]
    S = S[1:N+J]

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
    return S, x
end

#simple bound only
function boxLP(Q::LP{T}; nS=Settings{T}()) where {T}
    (; c, d, u, N, M, J) = Q
    tol = nS.tol

    M + J == 0 || error("Not a box LP (only simple bound)")

    return boxLP(c, d, u, N; tol=tol)
end

function boxLP(c, d, u, N; tol=2^-26)
    #fu = u .< Inf   #finite upper bound
    #fd = d .> -Inf   #finite lower bound
    fu = u .== Inf   #no upper bound
    fd = d .== -Inf   #no lower bound
    status = -1  # 1: unique; 0 infeasible; 2 infinitely many sol; 3 unbounded; -1 computing
    x = copy(d)
    S = fill(DN, N)

    for k in 1:N
        if abs(c[k]) <= tol  #c[k] = 0
            status = 2
            #break
            continue    #to find a solution
        end

        if c[k] > tol   # c > 0
            if fd[k]
                status = 3
                break
            end
        else    # c < 0
            x[k] = u[k]
            if fu[k]
                status = 3
                break
            end
            S[k] = UP
        end
    end
    f = c' * x
    return S, x, f, status
end


function auxLP(Q::LP{T}, nS) where {T}
    #An initial feasible point by performing Phase-I Simplex on the polyhedron
    (; A, b, G, g, d, u, N, M, J) = Q
    tol = nS.tol

    fu = u .== Inf   #no upper bound
    fd = d .== -Inf   #no lower bound
    fv = fu .&& fd  #free variable
    iv = findall(fv)
    id = findall(fd .&& .!fv)   # (-∞, u]

    #add slack variables for Gz<=g, and artificial variables for Ax=b, and s for aux
    M0 = M + J
    N0 = N + M0 + 1
    A0 = [A zeros(T, M, J) Matrix{T}(I, M, M) -ones(T, M)
        G Matrix{T}(I, J, J) zeros(T, J, M) -ones(T, J)]
    b0 = [b; g]
    d0 = [d; zeros(T, J + M + 1)]
    u0 = [u; fill(Inf, J + M + 1)]

    #compute x0
    x = copy(d0)
    x[iv] .= 0
    x[id] .= u0[id]
    Aa = @view A0[:, 1:N]
    xa = @view x[1:N]
    s, is = findmax(Aa * xa - b0)
    if s < 0
        s = 0.0
        is = 0
    end

    x[end] = s
    x[N+1:N+J] .= (g - G * xa) .+ s
    x[N+J+1:N+J+M] .= (b - A * xa) .+ s

    #fill the status
    S = fill(DN, N0)    #all are equalities
    S[iv] .= IN
    S[id] .= UP

    S[N+1:N+J+M+1] .= IN    #for moving on
    if is == 0
        S[end] = DN
    else
        is = is > M ? is - M : is + J
        S[N+is] = DN
    end

    c0 = zeros(T, N0)
    c0[N+J+1:end] .= 1.0
    G0 = zeros(T, 0, N0)
    g0 = zeros(T, 0)
    x0, f, status = solveLP(c0, A0, b0, G0, g0, d0, u0, S, x; nS=nS)

    if abs(f) > tol #* 3000
        error("feasible region is empty")
    end

    x = x0[1:N]
    S = S[1:N+J]

    for k in N+1:N+J    #inequalities
        S[k] = S[k] == IN ? OE : EO
    end

    return S, x

end



end



