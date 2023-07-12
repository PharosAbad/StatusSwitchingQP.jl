"Criss-Cross Algorithm for LP"
module CCLP

using LinearAlgebra
using StatusSwitchingQP: Status, IN, DN, UP, OE, EO, Event, LP, QP, getColumnsQR, getRowsGJr, Settings, getRowsQR, cAb, boxLP
export solveLP

function initSx(c::Vector{T}, A, b; settings=Settings{T}()) where {T}    # 2023-07-11 18:57:37   外点法
    #d = 0, upper bound is u=+∞
    tol = settings.tol
    N = length(c)
    #M = length(b)

    #=
    p = A' * (A' \ c) - c
    #ip = sortperm(abs.(p))
    #ip = sortperm(p)
    ip = sortperm(p; rev=true)
    =#
    p = c - A' * (A' \ c)
    ip = sortperm(p)


    S = fill(DN, N)

    ra, la = getRowsGJr(A[:, ip]', tol)
    #ra = 1:length(b)
    #ra = getColumnsQR(A[:, ip], tol)

    ib = ip[ra]
    S[ib] .= IN

    #iz = findall(S .!= IN)
    x = lu(A[:, ib]) \ b

    #display(norm(A * x - b, Inf))
    return x, S, ib

end



"""

        solveLP(Q::LP{T}; settings=Settings{T}())
        solveLP(c::Vector{T}, A, b; settings)

using Least-Index Criss-Cross Method (smallest subscript criss-cross algorithm), solve the following LP defined by Q::LP

```math
min   f=c′x
s.t.  Ax=b  ∈R^{M}
      Gx≤g  ∈R^{J}
      d≤x≤u ∈R^{N}
```

Outputs
    x               : solution,  N x 1 vector
    S               : Vector{Status}, (N+J)x1
    status          : 1 unique; 0 infeasible; 2 infinitely many sol; 3 unbounded or dual infeasible; -1 numerical errors; -maxIter, not done

solveLP(c, A, b; settings):  for

```math
min   f=c′x
s.t.  Ax=b  ∈R^{M}
      x≥0 ∈R^{N}
```


See also [`Status`](@ref), [`LP`](@ref), [`Settings`](@ref)

"""
function solveLP(c::Vector{T}, A, b; settings=Settings{T}()) where {T}   #2023-07-11 19:06:38
    tol = settings.tol
    maxIter = settings.maxIter

    x, S, B = initSx(c, A, b; settings=settings)

    N = length(c)
    F = trues(N)
    F[B] .= false

    #display((N, length(B)))

    iter = 0
    #maxIter = 17
    @inbounds while true
        iter += 1
        if iter > maxIter
            return x, S, -iter
        end

        B = findall(.!F)
        invB = inv(lu(A[:, B]))
        #Y = invB * A[:, F]
        x = invB * b

        #h = c[F] - Y' * c[B]
        AF = @view A[:, F]
        #AF = A[:, F]
        h = c[F] - AF' * (invB' * c[B])
        ik = findall(h .< -tol)
        ih = findall(F)
        K = ih[ik]

        D = x .< -tol
        il = findall(D)
        L = B[il]

        nK = length(K)
        nL = length(L)
        q = false
        w = 0
        if nK > 0
            if nL > 0
                if L[1] > K[1]
                    w = ik[1]
                else
                    w = il[1]
                    q = true
                end
            else
                w = ik[1]
            end
        elseif nL > 0
            w = il[1]
            q = true
        end

        if w == 0
            break
        end

        if q   #leaving
            #y = Y[w, :]

            #y = invB[w, :]' * AF
            #J = findall(vec(y) .< -tol)
            y = AF' * invB[w, :]
            J = findall(y .< -tol)

            #println("leaving ", y)
            #println("leaving ", J)

            if length(J) == 0
                return x, S, 0 #1 unique; 0 infeasible; 2 infinitely many sol; 3 unbounded; -1 numerical errors; -maxIter, not done       -2 dual infeasible
            end

            k = ih[J[1]]
            l = B[w]
        else    #entering
            #y = Y[:, w]
            k = ih[w]
            y = invB * A[:, k]
            J = findall(y .> tol)
            #println("entering ", y)
            #println("entering ", J)
            if length(J) == 0
                #return x, S, -2 #1 unique; 0 infeasible; 2 infinitely many sol; 3 unbounded; -1 numerical errors; -maxIter, not done       -2 dual infeasible
                return x, S, 3 #1 unique; 0 infeasible; 2 infinitely many sol; 3 unbounded or dual infeasible; -1 numerical errors; -maxIter, not done
            end
            #k = ih[w]
            l = B[J[1]]
        end

        #display((k, l))

        F[l] = true
        F[k] = false
        S[k] = IN
        S[l] = DN

    end
    return x, S, 1

end

#=
function solveLP(c::Vector{T}, A, b; settings=Settings{T}()) where {T}   #2023-07-12 20:42:53
    tol = settings.tol
    maxIter = settings.maxIter

    x, S, B = initSx(c, A, b; settings=settings)

    N = length(c)
    F = trues(N)
    F[B] .= false

    iter = 0
    Least = false
    @inbounds while true
        iter += 1
        if iter > maxIter
            return x, S, -iter
        end
        if iter > N
            Least = true    #switch to Least-Index rule
        end

        B = findall(.!F)
        invB = inv(lu(A[:, B]))
        x = invB * b
        AF = @view A[:, F]

        h = c[F] - AF' * (invB' * c[B])
        ik = findall(h .< -tol)
        ih = findall(F)
        K = ih[ik]

        D = x .< -tol
        il = findall(D)
        L = B[il]

        nK = length(K)
        nL = length(L)
        q = false
        w = 0
        if nK > 0
            if nL > 0
                if L[1] > K[1]
                    w = ik[1]
                else
                    w = il[1]
                    q = true
                end
            else
                w = ik[1]
            end
        elseif nL > 0
            w = il[1]
            q = true
        end

        if w == 0
            break
        end

        if q   #leaving
            y = AF' * invB[w, :]
            J = findall(y .< -tol)
            nJ = length(J)

            if nJ == 0
                return x, S, 0 #1 unique; 0 infeasible; 2 infinitely many sol; 3 unbounded; -1 numerical errors; -maxIter, not done       -2 dual infeasible
            else
                if nJ == 1
                    j1 = 1
                else
                    #j1 = Least ? 1 : argmin(y[J])
                    j1 = Least ? 1 : argmax(y[J])
                end
            end

            k = ih[J[j1]]
            l = B[w]
        else    #entering
            k = ih[w]
            y = invB * A[:, k]
            J = findall(y .> tol)
            nJ = length(J)
            if nJ == 0
                #return x, S, -2 #1 unique; 0 infeasible; 2 infinitely many sol; 3 unbounded; -1 numerical errors; -maxIter, not done       -2 dual infeasible
                return x, S, 3 #1 unique; 0 infeasible; 2 infinitely many sol; 3 unbounded or dual infeasible; -1 numerical errors; -maxIter, not done
            else
                if nJ == 1
                    j1 = 1
                else
                    #j1 = Least ? 1 : argmax(y[J])
                    j1 = Least ? 1 : argmin(y[J])
                end
            end
            #l = B[J[1]]
            l = B[J[j1]]
        end


        F[l] = true
        F[k] = false
        S[k] = IN
        S[l] = DN
    end
    return x, S, 1

end
=#

# DN for free variable means x=0, not x = -∞
function solveLP(Q::LP{T}; settings=Settings{T}()) where {T}    #2023-07-11 19:09:57
    #(; N, J) = Q
    N = Q.N
    J = Q.J
    nj = N + J
    #nj = N + Q.J

    tol = settings.tol


    if J == 0 && Q.M == 0 #a box
        #return boxLP(Q.c, Q.d, Q.u, N; tol=tol)
        return boxLP(Q; settings=settings)
    end


    status, c, A, b, d0, c0, A0, iv, id, ib = cAb(Q; tol=tol)
    if status <= 0
        return zeros(T, N), fill(DN, N), status
    end

    #t0 = time_ns()
    y, S, status = solveLP(c, A, b; settings=settings)
    #t1 = time_ns()
    #println("optm ", (t1 - t0) / 1e9, " seconds")

    @inbounds begin
        x = zeros(T, length(c))
        B = findall(S .== IN)
        x[B] = y
        x[1:nj] .+= d0

        #display((S, x))
        #return x, S, status
        if status <= 0
            return x[1:N], S[1:nj], status
        end


        #B = findall(S .== IN)
        n = length(iv)
        if n > 0    #free variable
            x[iv] .-= x0[nj+1:nj+n]
            M0 = size(A0, 1)
            for k in 1:M0
                t = B[k]
                if t > nj
                    B[k] = iv[t-nj]
                    S[B[k]] = IN
                end
            end

        end

        m = length(id)
        if m > 0   # flip u d
            x[id] = -x[id]
            for k in 1:m
                i = id[k]
                if S[i] == DN
                    S[i] = UP
                end
            end
        end

        m = length(ib)
        if m > 0   # upper bound
            L0 = nj + n
            for k in 1:m
                if S[k+L0] == DN
                    S[ib[k]] = UP
                    #=
                    t = ib[k]
                    S[t] = UP
                    B = setdiff(B, t)
                    =#
                end
            end
        end


        #detect infinitely many solution
        S = S[1:nj]
        #F = findall(S .!= IN)
        B = findall(S .== IN)
        F = trues(nj)
        F[B] .= false
        invB = inv(lu(A0[:, B]))
        #Y = invB * A0[:, F]
        #h = c0[F] - Y' * c0[B]
        h = c0[F] - A0[:, F]' * (invB' * c0[B])
        ih = abs.(h) .< tol
        status = sum(ih) > 0 ? 2 : 1

        x = x[1:N]

        for k in N+1:nj
            S[k] = S[k] == DN ? EO : OE
        end
    end
    return x, S, status
end

end