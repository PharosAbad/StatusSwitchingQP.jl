
@inline function getRowsGJ(X::Matrix{T}, tol=eps(norm(X, Inf))) where {T}
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
function getRowsGJr(X::Matrix{T}, tol=eps(norm(X, Inf))) where {T}      # row poviting
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
function getRows(A::Matrix{T}, tol=sqrt(eps(T))) where {T}
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

