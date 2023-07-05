"Simplex algorithm"
module Simplex
using LinearAlgebra
using StatusSwitchingQP: Status, IN, DN, UP, OE, EO, Event, LP, Settings, getRowsGJr
export SimplexLP, cDantzigLP, maxImprvLP, stpEdgeLP

#=
REMARK: writing the next basis as a product of the current basis times an easily invertible matrix can be extended over several iterations. We do not adopt this for accuracy
    If this method is adopt, how often should one recompute an inverse of the current basis?
=#



"""
        stpEdgeLP(c, A, b, d, u, B, S; invB, q, tol=2^-26)

using  steepest edge pivot rule to solve LP (switch to Bland's rule if iters > N)

```math
        min   f=c′x
        s.t.  Ax=b
              d≤x≤u
```
the native model requires that `d` is finite
B    : index set of basic variable, always sorted
S    : Vector{Status}, Nx1
invB : inverse of the basic matrix
q    : x[B]

Outputs             B and S in caller is changed
    status          : 1 unique; 0 infeasible; 2 infinitely many sol; 3 unbounded
    x               : solution,  N x 1 vector
    invB            : inverse of the basic matrix


"""
function stpEdgeLP(c::Vector{T}, A, b, d, u, B, S; invB, q, tol=2^-26) where {T}
    #B is always sorted. B and S in caller is changed, compute invB and q=x[B] each step, switch Dantzig to Bland rule when iter > N

    N = length(c)
    M = length(b)
    F = trues(N)
    F[B] .= false
    gt = zeros(T, M)    #step size
    ip = zeros(Int, M)    #tracking the rows of p
    Sb = fill(DN, M)    #State of leaving to be

    ud = u - d
    du = -ud
    fu = u .< Inf   #finite upper bound

    x = copy(d)
    @inbounds begin
        iu = findall(S .== UP)
        x[iu] = u[iu]

        Y = invB * A[:, F]
        h = c[F] - Y' * c[B]
        id = S[F] .== DN
        h[id] .= -h[id]
        ih = findall(h .> tol)
        hp = h[ih]
        iH = findall(F)[ih]
    end
    nH = length(iH)
    Edge = true
    #Edge = false
    loop = 0
    @inbounds while nH > 0
        loop += 1
        if loop > N
            Edge = false    #switch to Bland rule
        end
        #display((Edge, B))
        if Edge             #steepest edge
            #=
            se = zeros(T, nH)
            for n in eachindex(hp)
                y = Y[:, ih[n]]
                #se[n] = hp[n]^2 / (1 + y' * y)
                se[n] = hp[n]^2 / (1 + sum(y .^ 2))
            end
            =#

            # #=
            #Yh = Y[:, ih] .^2
            #y = sum(Yh, dims=1)
            #y = vec(sum(Y[:, ih] .^ 2, dims=1))
            y = vec(sum(Y[:, ih] .^ 2, dims=1)) .+ 1
            #se = hp .^2
            #se ./= (1 .+ vec(y))
            #se = hp .^ 2 ./ (1 .+ y)
            se = hp .^ 2 ./ y
            # =#

            k0 = argmax(se)
        else
            k0 = 1
            #k0 = argmax(hp)
        end
        k = iH[k0]

        #display((Edge, B, hp, se, k, iH[argmax(hp)]))
        #display((Edge, B, hp, k, iH[argmax(hp)]))

        p = invB * A[:, k]
        kd = S[k] == DN
        m = 0
        if kd
            for j in 1:M
                i = B[j]
                if p[j] > tol
                    m += 1
                    gt[m] = (q[j] - d[i]) / p[j]
                    ip[m] = j
                    Sb[m] = DN
                elseif p[j] < -tol #&& fu[i]
                    m += 1
                    gt[m] = (q[j] - u[i]) / p[j]
                    ip[m] = j
                    Sb[m] = UP
                end
            end

            if m == 0   # p=0 => A[:,k]=0
                if fu[k]    #DN -> UP
                    l = -1
                else    # unbounded
                    x[B] = q
                    return 3, x, invB
                end
            else
                (gl, l) = findmin(gt[1:m])  #gl>0
                if fu[k]
                    if gl >= ud[k]   #DN -> UP
                        l = -1
                    else
                        Sl = Sb[l]
                        l = ip[l]
                    end
                else
                    if isinf(gl) #unbounded
                        x[B] = q
                        return 3, x, invB
                    end
                    Sl = Sb[l]
                    l = ip[l]
                end
            end

        else    #UP
            for j in 1:M
                i = B[j]
                if p[j] > tol #&& fu[i]
                    m += 1
                    gt[m] = (q[j] - u[i]) / p[j]
                    ip[m] = j
                    Sb[m] = UP
                elseif p[j] < -tol
                    m += 1
                    gt[m] = (q[j] - d[i]) / p[j]
                    ip[m] = j
                    Sb[m] = DN
                end
            end

            if m == 0   # p=0 => A[:,k]=0
                l = -2  #UP -> DN
            else
                (gl, l) = findmax(gt[1:m])  #gl<0
                if gl <= du[k]  #du[k] is finite, for S[k]==UP and d is finite
                    l = -2  #UP -> DN
                else
                    Sl = Sb[l]
                    l = ip[l]
                end
            end
        end


        if l == -1  #flip the state
            S[k] = UP
            x[k] = u[k]
            #q -= g0 * p
        elseif l == -2  #flip the state
            S[k] = DN
            x[k] = d[k]
            #q -= g0 * p
        elseif l > 0
            m = l
            l = B[l]       #leaving index
            F[k] = false
            F[l] = true
            B[m] = k

            #B = sort(B)
            sort!(B)
            #invB = inv(A[:, B])
            invB = inv(lu(A[:, B]))

            S[k] = IN
            S[l] = Sl
            x[l] = Sl == DN ? d[l] : u[l]
            Y = invB * A[:, F]
            #q = invB * b - Y * x[F]
        end

        q = invB * b - Y * x[F]
        h = c[F] - Y' * c[B]
        id = S[F] .== DN
        h[id] .= -h[id]
        #ih = h .> tol
        ih = findall(h .> tol)
        hp = h[ih]
        iH = findall(F)[ih]
        nH = length(iH)
    end

    x[B] = q
    ms = any(abs.(h) .< tol)
    status = ms ? 2 : 1
    #display((status, length(h), h))
    return status, x, invB
end

#=
#much slow
# (switch to Bland's rule if zero progress temporarily)
function stpEdgeLP(c::Vector{T}, A, b, d, u, B, S; invB, q, tol=2^-26) where {T}
    #B is always sorted. B and S in caller is changed, compute invB and q=x[B] each step, switch Dantzig to Bland rule when iter > N
    #=
    not able to prove that steepest edge is anti-cycling. But all cycling example are unbounded, that is not real cycling, as unbounded is detected always.
    Hence, all possible cycling are at the same point (step size is zero)
    =#


    N = length(c)
    M = length(b)
    F = trues(N)
    F[B] .= false
    gt = zeros(T, M)    #step size
    ip = zeros(Int, M)    #tracking the rows of p
    Sb = fill(DN, M)    #State of leaving to be

    ud = u - d
    du = -ud
    fu = u .< Inf   #finite upper bound

    x = copy(d)
    @inbounds begin
        iu = findall(S .== UP)
        x[iu] = u[iu]

        Y = invB * A[:, F]
        h = c[F] - Y' * c[B]
        id = S[F] .== DN
        h[id] .= -h[id]
        ih = findall(h .> tol)
        hp = h[ih]
        iH = findall(F)[ih]
    end
    nH = length(iH)
    Edge = true
    #Edge = false
    #loop = 0
    gl = zero(T)
    @inbounds while nH > 0
        #=
        loop += 1
        if loop > N
            Edge = false    #switch to Bland rule
        end
        =#
        if Edge             #steepest edge
            y = vec(sum(Y[:, ih] .^ 2, dims=1)) .+ 1
            se = hp .^ 2 ./ y
            k0 = argmax(se)
        else
            k0 = 1
            #k0 = argmax(hp)
        end
        k = iH[k0]

        #display((Edge, B, hp, se, k, iH[argmax(hp)]))
        #display((Edge, B, hp, k, iH[argmax(hp)]))

        p = invB * A[:, k]
        kd = S[k] == DN
        m = 0
        if kd
            for j in 1:M
                i = B[j]
                if p[j] > tol
                    m += 1
                    gt[m] = (q[j] - d[i]) / p[j]
                    ip[m] = j
                    Sb[m] = DN
                elseif p[j] < -tol #&& fu[i]
                    m += 1
                    gt[m] = (q[j] - u[i]) / p[j]
                    ip[m] = j
                    Sb[m] = UP
                end
            end

            if m == 0   # p=0 => A[:,k]=0
                if fu[k]    #DN -> UP
                    l = -1
                else    # unbounded
                    x[B] = q
                    return 3, x, invB
                end
            else
                (gl, l) = findmin(gt[1:m])  #gl>0
                if fu[k]
                    if gl >= ud[k]   #DN -> UP
                        l = -1
                    else
                        Sl = Sb[l]
                        l = ip[l]
                    end
                else
                    if isinf(gl) #unbounded
                        x[B] = q
                        return 3, x, invB
                    end
                    Sl = Sb[l]
                    l = ip[l]
                end
            end

        else    #UP
            for j in 1:M
                i = B[j]
                if p[j] > tol #&& fu[i]
                    m += 1
                    gt[m] = (q[j] - u[i]) / p[j]
                    ip[m] = j
                    Sb[m] = UP
                elseif p[j] < -tol
                    m += 1
                    gt[m] = (q[j] - d[i]) / p[j]
                    ip[m] = j
                    Sb[m] = DN
                end
            end

            if m == 0   # p=0 => A[:,k]=0
                l = -2  #UP -> DN
            else
                (gl, l) = findmax(gt[1:m])  #gl<0
                if gl <= du[k]  #du[k] is finite, for S[k]==UP and d is finite
                    l = -2  #UP -> DN
                else
                    Sl = Sb[l]
                    l = ip[l]
                end
            end
        end


        if l == -1  #flip the state
            S[k] = UP
            x[k] = u[k]
            #q -= g0 * p
        elseif l == -2  #flip the state
            S[k] = DN
            x[k] = d[k]
            #q -= g0 * p
        elseif l > 0

            if Edge && abs(gl) < tol    #temporarily switch to Bland rule if zero stepsize
                Edge = false
                continue
            end
            Edge = true

            m = l
            l = B[l]       #leaving index
            F[k] = false
            F[l] = true
            B[m] = k

            #B = sort(B)
            sort!(B)
            #invB = inv(A[:, B])
            invB = inv(lu(A[:, B]))

            S[k] = IN
            S[l] = Sl
            x[l] = Sl == DN ? d[l] : u[l]
            Y = invB * A[:, F]
            #q = invB * b - Y * x[F]
        end

        q = invB * b - Y * x[F]
        h = c[F] - Y' * c[B]
        id = S[F] .== DN
        h[id] .= -h[id]
        #ih = h .> tol
        ih = findall(h .> tol)
        hp = h[ih]
        iH = findall(F)[ih]
        nH = length(iH)
    end

    x[B] = q
    ms = any(abs.(h) .< tol)
    status = ms ? 2 : 1
    return status, x, invB
end
=#




"""
        cDantzigLP(c, A, b, d, u, B, S; invB, q, tol=2^-26)

using  combined Dantzig's pivot rule to solve LP (combine the largest-distance rule with Dantzig's pivot rule, and switch to Bland's rule if iters > N)

```math
        min   f=c′x
        s.t.  Ax=b
              d≤x≤u
```
the native model requires that `d` is finite
B    : index set of basic variable, always sorted
S    : Vector{Status}, Nx1
invB : inverse of the basic matrix
q    : x[B]

Outputs             B and S in caller is changed
    status          : 1 unique; 0 infeasible; 2 infinitely many sol; 3 unbounded
    x               : solution,  N x 1 vector
    invB            : inverse of the basic matrix


"""
function cDantzigLP(c::Vector{T}, A, b, d, u, B, S; invB, q, tol=2^-26) where {T}
    #B is always sorted. B and S in caller is changed, compute invB and q=x[B] each step, switch Dantzig to Bland rule when iter > N

    N = length(c)
    M = length(b)
    F = trues(N)
    F[B] .= false
    gt = zeros(T, M)    #step size
    ip = zeros(Int, M)    #tracking the rows of p
    Sb = fill(DN, M)    #State of leaving to be

    ud = u - d
    du = -ud
    fu = u .< Inf   #finite upper bound

    cA = zeros(T, N)    #norm of each column of A
    x = copy(d)
    @inbounds begin
        for k in 1:N
            cA[k] = norm(A[:, k])
        end
        #ldr = all(cA .> tol)    #some cols maybe zero. EfficientFrontier do not have this prob

        #x = zeros(T, N)
        #x[S.==DN] .= d[S.==DN]
        #x[S.==UP] .= u[S.==UP]
        iu = findall(S .== UP)
        x[iu] = u[iu]
        #x[B] = q

        Y = invB * A[:, F]
        h = c[F] - Y' * c[B]
        ih = S[F] .== DN
        h[ih] .= -h[ih]
        ih = h .> tol
        hp = h[ih]
        iH = findall(F)[ih]
    end
    nH = length(iH)
    Bland = false
    loop = 0
    @inbounds while nH > 0
        loop += 1
        if loop > N
            Bland = true    #Dantzig rule switch to Bland rule
        end
        #display((Bland, B))

        #Dantzig rule or Bland rule
        #k0 = Bland ? 1 : argmax(hp)
        k0 = Bland ? 1 : argmax(hp ./ cA[iH]) #Largest-Distance Rule
        k = iH[k0]
        p = invB * A[:, k]
        kd = S[k] == DN
        m = 0
        if kd
            for j in 1:M
                i = B[j]
                if p[j] > tol
                    m += 1
                    gt[m] = (q[j] - d[i]) / p[j]
                    ip[m] = j
                    Sb[m] = DN
                elseif p[j] < -tol #&& fu[i]
                    m += 1
                    gt[m] = (q[j] - u[i]) / p[j]
                    ip[m] = j
                    Sb[m] = UP
                end
            end

            if m == 0   # p=0 => A[:,k]=0
                if fu[k]    #DN -> UP
                    l = -1
                else    # unbounded
                    x[B] = q
                    return 3, x, invB
                end
            else
                (gl, l) = findmin(gt[1:m])  #gl>0
                if fu[k]
                    if gl >= ud[k]   #DN -> UP
                        l = -1
                    else
                        Sl = Sb[l]
                        l = ip[l]
                    end
                else
                    if isinf(gl) #unbounded
                        x[B] = q
                        return 3, x, invB
                    end
                    Sl = Sb[l]
                    l = ip[l]
                end
            end

        else    #UP
            for j in 1:M
                i = B[j]
                if p[j] > tol #&& fu[i]
                    m += 1
                    gt[m] = (q[j] - u[i]) / p[j]
                    ip[m] = j
                    Sb[m] = UP
                elseif p[j] < -tol
                    m += 1
                    gt[m] = (q[j] - d[i]) / p[j]
                    ip[m] = j
                    Sb[m] = DN
                end
            end

            if m == 0   # p=0 => A[:,k]=0
                l = -2  #UP -> DN
            else
                (gl, l) = findmax(gt[1:m])  #gl<0
                if gl <= du[k]  #du[k] is finite, for S[k]==UP and d is finite
                    l = -2  #UP -> DN
                else
                    Sl = Sb[l]
                    l = ip[l]
                end
            end
        end


        if l == -1  #flip the state
            S[k] = UP
            x[k] = u[k]
            #q -= g0 * p
        elseif l == -2  #flip the state
            S[k] = DN
            x[k] = d[k]
            #q -= g0 * p
        elseif l > 0
            m = l
            l = B[l]       #leaving index
            F[k] = false
            F[l] = true
            B[m] = k

            #B = sort(B)
            sort!(B)
            #invB = inv(A[:, B])
            invB = inv(lu(A[:, B]))

            S[k] = IN
            S[l] = Sl
            x[l] = Sl == DN ? d[l] : u[l]
            Y = invB * A[:, F]
            #q = invB * b - Y * x[F]
        end

        q = invB * b - Y * x[F]
        h = c[F] - Y' * c[B]
        ih = S[F] .== DN
        h[ih] .= -h[ih]
        ih = h .> tol
        hp = h[ih]
        iH = findall(F)[ih]
        nH = length(iH)
    end


    x[B] = q

    ms = any(abs.(h) .< tol)
    status = ms ? 2 : 1
    return status, x, invB
end


"""
        maxImprvLP(c, A, b, d, u, B, S; invB, q, tol=2^-26)

using  max improvement pivot rule to solve LP (no cycle since it becomes Bland's rule if the improvement is 0)

```math
        min   f=c′x
        s.t.  Ax=b
              d≤x≤u
```
the native model requires that `d` is finite
B    : index set of basic variable, always sorted
S    : Vector{Status}, Nx1
invB : inverse of the basic matrix
q    : x[B]

Outputs             B and S in caller is changed
    status          : 1 unique; 0 infeasible; 2 infinitely many sol; 3 unbounded
    x               : solution,  N x 1 vector
    invB            : inverse of the basic matrix


"""
function maxImprvLP(c::Vector{T}, A, b, d, u, B, S; invB, q, tol=2^-26) where {T}
    #greatest improvement, B is always sorted. B and S in caller is change, compute invB and q=x[B] each step
    N = length(c)
    M = length(b)
    F = trues(N)
    F[B] .= false
    gt = zeros(T, M)    #theta
    ip = zeros(Int, M)    #tracking the rows of p
    Sb = fill(DN, M)    #State of leaving to be

    ud = u - d
    du = -ud
    fu = u .< Inf   #finite upper bound

    #x = zeros(T, N)
    #x[S.==UP] .= u[S.==UP]
    #x[S.==DN] .= d[S.==DN]
    x = copy(d)
    @inbounds begin
        iu = findall(S .== UP)
        x[iu] = u[iu]
        #x[B] = q

        Y = invB * A[:, F]
        h = c[F] - Y' * c[B]
        vS = S[F]   #State of leaving to be, for each candidate k
        g = zeros(T, N - M)    #theta for each candidate k
        ig = zeros(Int, N - M)    #min subscrip for theta for each candidate k
        #vD = falses(N - M)  #DN or not, for each candidate k


        ih = S[F] .== DN
        h[ih] .= -h[ih]
        ih = h .> tol
        iH = findall(F)[ih]
    end
    nH = length(iH)
    @inbounds while nH > 0
        P = @view Y[:, ih]
        for n in 1:nH
            k = iH[n]
            p = P[:, n]
            kd = S[k] == DN
            #vD[n] = kd
            m = 0
            if kd
                for j in 1:M
                    i = B[j]
                    if p[j] > tol
                        m += 1
                        gt[m] = (q[j] - d[i]) / p[j]
                        ip[m] = j
                        Sb[m] = DN
                    elseif p[j] < -tol #&& u[i] < Inf
                        m += 1
                        gt[m] = (q[j] - u[i]) / p[j]
                        ip[m] = j
                        Sb[m] = UP
                    end
                end

                g0 = ud[k]
                if m == 0
                    if fu[k]    #DN -> UP
                        g[n] = g0
                        l = 1
                        ip[1] = -1
                    else    # unbounded
                        x[B] = q
                        return 3, x, invB
                    end
                else
                    (g[n], l) = findmin(gt[1:m])

                    if fu[k]
                        if g[n] >= g0 #DN -> UP
                            g[n] = g0
                            ip[l] = -1
                        end
                    else
                        if isinf(g[n])  #unbounded
                            x[B] = q
                            return 3, x, invB
                        end
                    end
                end

            else    #UP
                for j in 1:M
                    i = B[j]
                    if p[j] > tol && u[i] < Inf
                        m += 1
                        gt[m] = (q[j] - u[i]) / p[j]
                        ip[m] = j
                        Sb[m] = UP
                    elseif p[j] < -tol
                        m += 1
                        gt[m] = (q[j] - d[i]) / p[j]
                        ip[m] = j
                        Sb[m] = DN
                    end
                end
                g0 = du[k]  #finite
                if m == 0
                    g[n] = g0
                    l = 1
                    ip[1] = -2  #UP -> DN
                else
                    (g[n], l) = findmax(gt[1:m])
                    if g[n] <= g0
                        g[n] = g0
                        ip[l] = -2
                    end
                end
            end
            ig[n] = ip[l]
            vS[n] = Sb[l]
        end
        k = getfield(findmax(abs.(h[ih] .* g[1:nH])), 2)
        l = ig[k]
        #= if l < 0 && isinf(u[k])
            error("infeasible or unbounded")
        end =#

        #kd = vD[k]
        p = P[:, k]
        #gl = g[k]
        Sl = vS[k]
        k = iH[k]   #entering index
        if l == -1  #flip the state
            S[k] = UP
            x[k] = u[k]
            #q -= gl * p
        elseif l == -2  #flip the state
            S[k] = DN
            x[k] = d[k]
            #q -= gl * p
        elseif l > 0
            m = l
            l = B[l]       #leaving index
            F[k] = false
            F[l] = true
            B[m] = k

            #B = sort(B)
            sort!(B)
            #invB = inv(A[:, B])
            invB = inv(lu(A[:, B]))

            S[k] = IN
            S[l] = Sl
            x[l] = Sl == DN ? d[l] : u[l]
            Y = invB * A[:, F]
            #q = invB * b - Y * x[F]
        end

        q = invB * b - Y * x[F]
        h = c[F] - Y' * c[B]
        ih = S[F] .== DN
        h[ih] .= -h[ih]
        ih = h .> tol
        iH = findall(F)[ih]
        nH = length(iH)
    end


    x[B] = q

    ms = any(abs.(h) .< tol)
    status = ms ? 2 : 1
    return status, x, invB

end



"""

        SimplexLP(P::LP{T}; settings=Settings{T}(), min=true, Phase1=false)
        SimplexLP(c, A, b, d, u; settings=Settings{T}(), min=true, Phase1=false)

solve LP by simplex method. If `min=false`, we maximize the objective function

Outputs
    x               : solution,  N x 1 vector
    S               : Vector{Status}, (N+J)x1
    status          : 1 unique; 0 infeasible; 2 infinitely many sol; 3 unbounded

See also [`Status`](@ref), [`LP`](@ref), [`Settings`](@ref), [`cDantzigLP`](@ref), [`maxImprvLP`](@ref)
"""
function SimplexLP(P::LP{T}; settings=Settings{T}(), min=true, Phase1=false) where {T}
    #t0 = time()
    #An initial feasible point by performing Phase-I Simplex on the polyhedron
    #(; c, A, b, G, g, d, u, N, M, J, mc) = P
    c = P.c
    A = P.A
    b = P.b
    G = P.G
    g = P.g
    d = P.d
    u = P.u
    N = P.N
    M = P.M
    J = P.J
    mc = P.mc

    if mc <= 0
        return zeros(T, N), fill(DN, N), -1
    end
    #(; tol, rule) = settings
    tol = settings.tol
    rule = settings.rule

    solveLP = cDantzigLP
    if rule == :maxImprovement
        solveLP = maxImprvLP
    elseif rule == :stpEdgeLP
        solveLP = stpEdgeLP
    end

    #convert free variable: -∞ < x < +∞
    fu = (u .== Inf)   #no upper bound
    fd = (d .== -Inf)   #no lower bound
    fv = (fu .& fd)  #free variable
    iv = findall(fv)
    n = length(iv)
    #id = findall(fd .&& (.!fv))   # (-∞, u]
    id = findall(fd .& (.!fv))   # (-∞, u]  #ver 1.6

    #add slack variables for Gz<=g , and 2nd part of free variables
    nj = N + J
    M0 = M + J
    N0 = nj + n
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

    #purge redundancy
    m0 = rank(A0)
    if m0 < M0
        ra, la = getRowsGJr([A0 b0], tol)
        if length(ra) != la
            return zeros(T, N), fill(DN, N), 0   #0 infeasible
        end
        if m0 != la
            return zeros(T, N), fill(DN, N), -1   # -1 numerical errors
        end
        M0 = m0
        A0 = A0[ra, :]
        b0 = b0[ra]
    end


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

    iH, x, invB = solveLP(c1, A1, b1, d1, u1, B, S; invB=invB, q=q, tol=tol)
    f = sum(x[N0+1:end])
    if abs(f) > tol
        #error("feasible region is empty")
        return x[1:N], S[1:nj], 0   #0 infeasible
    end

    #t1 = time()
    #display(t1 - t0)

    #=
    if Phase1
        if n > 0    #free variable
            x[iv] .-= x0[nj+1:nj+n]
            S[iv] .= IN
        end
        m = length(id)
        if m > 0   # flip u d
            x[id] = -x[id]
            for k in 1:m
                if S[k] == DN
                    S[k] == UP
                end
            end
        end
        return x[1:N], S[1:nj], iH
    end
    =#


    #Phase II    --- --- phase 2 --- ---
    if Phase1
        x = x[1:N]
    else
        q = x[B]
        c0 = [c; zeros(T, J + n)]
        c0[id] .= -c0[id]
        c0[nj+1:end] .= -c0[iv]

        #Artificial Variable in Basis, driving out AV
        ib = findall(B .<= N0)
        m = length(ib)
        if m < M0
            iB = B[ib]
            F = trues(N0)
            F[iB] .= false
            iF = findall(F)
            ic = [iB; iF]
            #ra, la = getRowsGJr(copy(A0[:, ic]'), tol)
            ra, la = getRowsGJr(A0[:, ic]', tol)
            B = sort(ic[ra])
            iA = setdiff(B, iB)
            invB = inv(lu(A0[:, B]))
            S[iA] .= IN
            q = x[B]
        end

        #display(time() - t1)

        S = S[1:N0]
        if !min
            c0 = -c0
        end
        iH, x0, invB = solveLP(c0, A0, b0, d0, u0, B, S; invB=invB, q=q, tol=tol)
        x = x0[1:N]
    end

    S = S[1:nj]
    #variables: restore free, or flip back
    for k in N+1:nj    #inequalities
        S[k] = S[k] == IN ? OE : EO
    end

    if n > 0    #free variable
        x[iv] .-= x0[nj+1:nj+n]
        S[iv] .= IN

        #always infinitely many solutions when free variables
        for k in 1:M0
            if B[k] > nj
                B[k] = iv[B[k]-nj]
            end
        end

        F = trues(nj)
        F[B] .= false
        invB = inv(lu(A0[:, B]))
        Y = invB * A0[:, findall(F)]
        c0 = c0[1:nj]
        h = c0[F] - Y' * c0[B]
        ih = abs.(h) .< tol
        iH = sum(ih) > 0 ? 2 : 1
    end

    m = length(id)
    if m > 0   # flip u d
        x[id] = -x[id]
        for k in 1:m
            if S[k] == DN
                S[k] == UP
            end
        end
    end
    return x, S, iH
end

function SimplexLP(c::Vector{T}, A, b, d, u; settings=Settings{T}(), min=true, Phase1=false) where {T}
    tol = settings.tol
    rule = settings.rule
    M, N = size(A)

    solveLP = cDantzigLP
    if rule == :maxImprovement
        solveLP = maxImprvLP
    elseif rule == :stpEdgeLP
        solveLP = stpEdgeLP
    end

    #convert free variable: -∞ < x < +∞
    fu = (u .== Inf)   #no upper bound
    fd = (d .== -Inf)   #no lower bound
    fv = (fu .& fd)  #free variable
    iv = findall(fv)
    n = length(iv)
    id = findall(fd .& (.!fv))   # (-∞, u]  #ver 1.6

    #add slack variables for Gz<=g , and 2nd part of free variables

    M0 = M
    N0 = N + n
    A0 = [A -A[:, iv]]
    b0 = b
    d0 = [d; zeros(T, n)]
    u0 = [u; fill(Inf, n)]

    #(-∞, +∞)  -> two copy of [0, +∞)
    d0[iv] .= 0     #the 1st part of free variable
    #u0[iv] .= Inf

    #(-∞, u]  -> [-u, +∞)
    d0[id] .= -u0[id]
    u0[id] .= Inf
    A0[:, id] .= -A0[:, id]

    #purge redundancy
    m0 = rank(A0)
    if m0 < M0
        ra, la = getRowsGJr([A0 b0], tol)
        if length(ra) != la
            return zeros(T, N), fill(DN, N), 0   #0 infeasible
        end
        if m0 != la
            return zeros(T, N), fill(DN, N), -1   # -1 numerical errors
        end
        M0 = m0
        A0 = A0[ra, :]
        b0 = b0[ra]
    end


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

    iH, x, invB = solveLP(c1, A1, b1, d1, u1, B, S; invB=invB, q=q, tol=tol)
    f = sum(x[N0+1:end])
    if abs(f) > tol
        return x[1:N], S[1:N], 0   #0 infeasible
    end

    #Phase II    --- --- phase 2 --- ---
    if Phase1
        x = x[1:N]
    else
        q = x[B]
        c0 = [c; zeros(T, n)]
        c0[id] .= -c0[id]
        c0[N+1:end] .= -c0[iv]

        #Artificial Variable in Basis, driving out AV
        ib = findall(B .<= N0)
        m = length(ib)
        if m < M0
            iB = B[ib]
            F = trues(N0)
            F[iB] .= false
            iF = findall(F)
            ic = [iB; iF]
            #ra, la = getRowsGJr(copy(A0[:, ic]'), tol)
            ra, la = getRowsGJr(A0[:, ic]', tol)
            B = sort(ic[ra])
            iA = setdiff(B, iB)
            invB = inv(lu(A0[:, B]))
            S[iA] .= IN
            q = x[B]
        end

        #display(time() - t1)

        S = S[1:N0]
        if !min
            c0 = -c0
        end
        iH, x0, invB = solveLP(c0, A0, b0, d0, u0, B, S; invB=invB, q=q, tol=tol)
        x = x0[1:N]
    end

    S = S[1:N]
    #variables: restore free, or flip back

    if n > 0    #free variable
        x[iv] .-= x0[N+1:N+n]
        S[iv] .= IN

        #always infinitely many solutions when free variables
        for k in 1:M0
            if B[k] > N
                B[k] = iv[B[k]-N]
            end
        end

        F = trues(N)
        F[B] .= false
        invB = inv(lu(A0[:, B]))
        Y = invB * A0[:, findall(F)]
        c0 = c0[1:N]
        h = c0[F] - Y' * c0[B]
        ih = abs.(h) .< tol
        iH = sum(ih) > 0 ? 2 : 1
    end

    m = length(id)
    if m > 0   # flip u d
        x[id] = -x[id]
        for k in 1:m
            if S[k] == DN
                S[k] == UP
            end
        end
    end
    return x, S, iH


end



end

