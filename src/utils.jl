
@inline function getRowsGJ(X::AbstractMatrix{T}, tol=eps(norm(X, Inf))) where {T}
    #@inline function getRowsGJ(X::Matrix{T}, tol=eps(norm(X, Inf))) where {T}
    #Gauss-Jordan elimination, code form rref_with_pivots!
    A = copy(X)
    nr, nc = size(A)
    rows = Vector{Int64}()
    r0 = collect(1:nr)  #original row numb
    nc1 = nc - 1
    l1 = 0  #without the last col
    i = j = 1
    @inbounds while i <= nr && j <= nc
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

#@inline function getRowsGJr(X::Matrix{T}, tol=eps(norm(X, Inf))) where {T}      # row poviting
#function getRowsGJr(X::Matrix{T}, tol=eps(norm(X, Inf))) where {T}      # row poviting
function getRowsGJr(X::AbstractMatrix{T}, tol=2^-33) where {T}     # row poviting
    #Gauss-Jordan elimination, code form rref_with_pivots!
    A = copy(X)
    nr, nc = size(A)
    rows = Vector{Int64}()
    c0 = collect(1:nc)  #original col numb
    #nc1 = nc - 1
    l1 = 0  #without the last col
    i = j = 1
    @inbounds while i <= nr && j <= nc
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

#see https://discourse.julialang.org/t/matrix-division-vs-compute-by-hand-why-such-big-difference/98288
function getRows(A::AbstractMatrix{T}, tol=eps(norm(X, Inf))) where {T}
    #function getRows(A::Matrix{T}, tol=sqrt(eps(T))) where {T}
    #indicate the non-redundant rows, the begaining few rows can be zeros (redundant)
    M, N = size(A)
    if N == 0
        @warn "zero columns" size(A)
        return collect(axes(A, 1))
    end
    R = falses(M)
    if M == 0
        return findall(R)
    end

    r1 = M + 1
    #find the 1st non-zero row
    @inbounds for r in 1:M
        v = @view A[r, :]
        if norm(v, Inf) <= tol
            continue
        else
            R[r] = true
            r1 = r + 1
            break
        end
    end

    #rows after the 1st non-zero row
    H = @view A[R, :]
    @inbounds for r in r1:M
        #v = @view A[r:r, :]
        v = @view A[r, :]
        if norm(v, Inf) > tol && norm(v - H' * (H' \ v), Inf) > tol
            R[r] = true
            H = @view A[R, :]
        end
    end
    return findall(R)
end

#see https://www.mathworks.com/matlabcentral/answers/574543-algorithm-to-extract-linearly-dependent-columns-in-a-matrix#answer_474601
function getRowsQR(A::AbstractMatrix{T}, tol=2^-33) where {T}
    #Extract linearly independent subset of matrix rows
    ci = Vector{Int64}()
    G = qr(A', ColumnNorm())
    d = abs.(diag(G.R))
    #d1 = d[1]
    #if d1 > tol
    #r = findlast(d .>= tol * d1)
    if d[1] > tol
        r = findlast(d .>= tol)
        ci = G.p[1:r]
    end
    #return ci, G
    return ci
end

function getColumnsQR(A::AbstractMatrix{T}, tol=2^-33) where {T}
    #Extract linearly independent subset of matrix columns
    ci = Vector{Int64}()
    G = qr(A, ColumnNorm())
    d = abs.(diag(G.R))
    if d[1] > tol
        r = findlast(d .>= tol)
        ci = G.p[1:r]
    end
    return ci
end


function cAbdu(Q::LP{T}; tol=2^-26) where {T}
    #convert Gx≤g to equalities by adding slack variables
    #(; c, A, b, G, g, d, u, M, J) = P
    c = Q.c
    A = Q.A
    b = Q.b
    G = Q.G
    g = Q.g
    d = Q.d
    u = Q.u
    #N = Q.N
    M = Q.M
    J = Q.J

    # #=
    #purge redundancy in Ax=b
    ir, L = getRowsGJr([A b], tol)
    if length(ir) != L
        #error("infeasible")
        #return zeros(T, N), fill(DN, N), 0   #0 infeasible
        return 0, c, A, b, d, u
    end
    if L != M
        M = L
        A = A[ir, :]
        b = b[ir]
    end

    #g = g[1:2] # not change Q.g
    #g .= G[:,1]    #will change Q.g
    # =#

    #add slack variables, convert to all equalities
    A0 = [A zeros(T, M, J)
        G Matrix{T}(I, J, J)]
    b0 = [b; g]
    d0 = [d; zeros(T, J)]
    u0 = [u; fill(Inf, J)]
    c0 = [c; zeros(T, J)]
    return 1, c0, A0, b0, d0, u0

end


function cAb(Q::LP{T}; tol=2^-26) where {T}
    #convert Gx≤g to equalities by adding slack variables
    #convert free variables, convert (-∞, u]
    #(; c, A, b, G, g, d, u, M, J) = P
    c = Q.c
    A = Q.A
    b = Q.b
    G = Q.G
    g = Q.g
    d = Q.d
    u = Q.u
    N = Q.N
    M = Q.M
    J = Q.J


    #purge redundancy in Ax=b
    ir, L = getRowsGJr([A b], tol)
    if length(ir) != L
        return 0, c, A, b, d, c, A, Int64[], Int64[], Int64[]
    end
    if L != M
        M = L
        A = A[ir, :]
        b = b[ir]
    end

    #free variable: -∞ < x < +∞
    fu = (u .== Inf)   #no upper bound
    fd = (d .== -Inf)   #no lower bound
    fv = (fu .& fd)  #free variable
    iv = findall(fv)
    n = length(iv)
    fb = .!fu   #finite upper bound
    id = findall(fd .& fb)  # (-∞, u]  #ver 1.6
    m = length(id)
    ib = findall(.!fd .& fb)    # finite d and u
    L = length(ib)


    #add slack variables for Gz<=g
    c0 = [c; zeros(T, J)]
    A0 = [A zeros(T, M, J)
        G Matrix{T}(I, J, J)]
    b0 = [b; g]
    d0 = [d; zeros(T, J)]
    #u0 = [u; fill(Inf, J)]

    #flip   (-∞, u]  -> [-u, +∞)
    if m > 0
        c0[id] .= -c0[id]
        A0[:, id] .= -A0[:, id]
        d0[id] .= -u[id]
        #u0[id] .= Inf
    end

    #free variable
    d0[iv] .= 0     #the 1st part of free variable
    #make sure d is finite       d0 is finite now

    #2nd part of free variables
    c1 = [c0; -c0[iv]]
    A1 = [A0 -A0[:, iv]]
    N1 = N + J + n


    #d1 = [d0; zeros(T, n)]
    #u1 = [u0; fill(Inf, n)]

    #x<=u to equalities
    Au = [A1 zeros(T, M + J, L)
        Matrix{T}(I, N1, N1)[ib, :] Matrix{T}(I, L, L)]
    bu = [b0 - A0 * d0; u[ib] - d[ib]]
    cu = [c1; zeros(T, L)]

    return 1, cu, Au, bu, d0, c0, A0, iv, id, ib    #last three: free, up-down, up-box

end

#simple bound only
function boxLP(Q::LP{T}; settings=Settings{T}()) where {T}
    #(; c, d, u, N, M, J) = Q
    c = Q.c
    d = Q.d
    u = Q.u
    N = Q.N
    M = Q.M
    J = Q.J
    tol = settings.tol

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
    #f = c' * x
    return x, S, status
end
