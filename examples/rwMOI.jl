
import MathOptInterface as MOI
import MathOptInterface.Utilities as MOIU

#=
    mps = MOIU.Model{Float64}()
    fn = "./grow7.mps.gz"
    MOI.read_from_file(mps, fn)
    P = MOI2LP(mps)
    display(P)
=#

function getConstraints(P::MOIU.Model{T}, N, tol) where {T}

    Ab = Vector{Vector{T}}(undef, 0)
    Gg = Vector{Vector{T}}(undef, 0)  #<=
    d = fill(-Inf, N)
    u = fill(Inf, N)
    for (F, S) in MOI.get(P, MOI.ListOfConstraintTypesPresent())
        for cref in MOI.get(P, MOI.ListOfConstraintIndices{F,S}())
            f = MOI.get(P, MOI.ConstraintFunction(), cref)
            s = MOI.get(P, MOI.ConstraintSet(), cref)
            if F <: MOI.ScalarAffineFunction
                t = zeros(T, N + 1)
                nt = 0
                for term in f.terms
                    t[term.variable.value] = term.coefficient
                    nt += 1
                end
                if nt == 0
                    row_name = MOI.get(P, MOI.ConstraintName(), cref)
                    @warn "skipping redundant rows: " * row_name
                    continue
                end
                if S <: MOI.EqualTo
                    t[end] = s.value
                    push!(Ab, t)
                elseif S <: MOI.GreaterThan
                    t[end] = s.lower
                    push!(Gg, -t)
                elseif S <: MOI.LessThan
                    t[end] = s.upper
                    push!(Gg, t)
                elseif S <: MOI.Interval
                    row_name = MOI.get(P, MOI.ConstraintName(), cref)
                    @warn "repacking AffineFunction range rows: " * row_name s.lower s.upper
                    if abs(s.lower - s.upper) < tol   # an Equality
                        t[end] = (s.lower + s.upper) / 2
                        push!(Ab, t)
                    else    #dispact to Gg
                        t[end] = s.lower
                        push!(Gg, -t)
                        t[end] = s.upper
                        push!(Gg, t)
                    end
                end
            elseif F <: MOI.VariableIndex
                if S <: MOI.EqualTo
                    d[f.value] = s.value
                    u[f.value] = s.value
                    #display((F,S))
                elseif S <: MOI.GreaterThan
                    d[f.value] = s.lower
                elseif S <: MOI.LessThan
                    u[f.value] = s.upper
                elseif S <: MOI.Interval
                    d[f.value] = s.lower
                    u[f.value] = s.upper
                end
            end
        end
    end
    M = lastindex(Ab)
    J = lastindex(Gg)
    A = Matrix{T}(undef, M, N)
    b = Vector{T}(undef, M)
    G = Matrix{T}(undef, J, N)
    g = Vector{T}(undef, J)

    #display(("sizes ", N, M, J))

    for k in 1:M
        A[k, :] = Ab[k][1:end-1]
        b[k] = Ab[k][end]
    end
    for k in 1:J
        G[k, :] = Gg[k][1:end-1]
        g[k] = Gg[k][end]
    end

    return A, b, G, g, d, u, M, J
end

function getConstraintsI(P::MOIU.Model{T}, N, tol) where {T}
    #convert AffineFunction range rows by adding a new MOI.Interval variable
    Ab = Vector{Vector{T}}(undef, 0)
    Gg = Vector{Vector{T}}(undef, 0)  #<=
    d = fill(-Inf, N)
    u = fill(Inf, N)
    nr = 0  #line number
    ir = Vector{Int}()  #rec the line number

    iv = Vector{Vector{T}}(undef, 0)    #Interval [d, u]
    #vd = Vector{T}()
    #vu = Vector{T}()

    for (F, S) in MOI.get(P, MOI.ListOfConstraintTypesPresent())
        for cref in MOI.get(P, MOI.ListOfConstraintIndices{F,S}())
            f = MOI.get(P, MOI.ConstraintFunction(), cref)
            s = MOI.get(P, MOI.ConstraintSet(), cref)
            if F <: MOI.ScalarAffineFunction
                t = zeros(T, N + 1)
                nt = 0
                for term in f.terms
                    t[term.variable.value] = term.coefficient
                    nt += 1
                end
                if nt == 0
                    row_name = MOI.get(P, MOI.ConstraintName(), cref)
                    @warn "skipping redundant rows: " * row_name
                    continue
                end

                if S <: MOI.EqualTo
                    t[end] = s.value
                    push!(Ab, t)
                    nr += 1
                elseif S <: MOI.GreaterThan
                    t[end] = s.lower
                    push!(Gg, -t)
                elseif S <: MOI.LessThan
                    t[end] = s.upper
                    push!(Gg, t)
                elseif S <: MOI.Interval
                    nr += 1
                    row_name = MOI.get(P, MOI.ConstraintName(), cref)
                    @warn "repacking AffineFunction range rows: " * row_name s.lower s.upper
                    if abs(s.lower - s.upper) < tol   # an Equality
                        t[end] = (s.lower + s.upper) / 2
                        #=     push!(Ab, t)
                        else    #dispact to Gg
                            t[end] = s.lower
                            push!(Gg, -t)
                            t[end] = s.upper
                            push!(Gg, t)
                            =#
                    else    #add new variable
                        push!(ir, nr)
                        push!(iv, [s.lower, s.upper])
                        #push!(vd, s.lower)
                        #push!(vu, s.upper)
                    end
                    push!(Ab, t)
                end
            elseif F <: MOI.VariableIndex
                if S <: MOI.EqualTo
                    d[f.value] = s.value
                    u[f.value] = s.value
                    #display((F,S))
                elseif S <: MOI.GreaterThan
                    d[f.value] = s.lower
                elseif S <: MOI.LessThan
                    u[f.value] = s.upper
                elseif S <: MOI.Interval
                    d[f.value] = s.lower
                    u[f.value] = s.upper
                end
            end
        end
    end
    K = length(ir)
    M = lastindex(Ab)
    J = lastindex(Gg)
    #A = Matrix{T}(undef, M, N)
    A = zeros(T, M, N + K)
    b = Vector{T}(undef, M)
    #G = Matrix{T}(undef, J, N)
    G = zeros(T, J, N + K)
    g = Vector{T}(undef, J)

    #display(("sizes ", N, M, J))

    if K == 0
        for k in 1:M
            A[k, :] = Ab[k][1:end-1]
            b[k] = Ab[k][end]
        end
        for k in 1:J
            G[k, :] = Gg[k][1:end-1]
            g[k] = Gg[k][end]
        end

        return A, b, G, g, d, u, M, J
    end

    for k in 1:M
        A[k, 1:N] = Ab[k][1:end-1]
        b[k] = Ab[k][end]
    end
    for k in 1:J
        G[k, 1:N] = Gg[k][1:end-1]
        g[k] = Gg[k][end]
    end

    #d = [d; vd]
    d = [d; zeros(T, K)]
    #u = [u; vu]
    u = [u; zeros(T, K)]

    for k in eachindex(ir)
        n = N + k
        A[ir[k], n] = -1.0
        d[n] = iv[k][1]
        u[n] = iv[k][2]
    end
    #N += K
    return A, b, G, g, d, u, M, J
end

function MOI2LP(P::MOIU.Model{T}, tol=2^-26) where {T}

    #N = P.constraints.num_variables   #for MOI object only
    N = MOI.get(P, MOI.NumberOfVariables())
    f = MOI.get(P, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}())
    c = zeros(T, N)
    for term in f.terms
        c[term.variable.value] = term.coefficient
    end

    A, b, G, g, d, u, M, J = getConstraintsI(P, N, tol) #AffineFunction range rows are converted by adding a new MOI.Interval variable

    #return c, A, b, G, g, d, u
    #return LP(c, A, b, G, g, d, u, N, M, J)


    N1 = length(d)
    #display((N, N1))
    if N != N1
        c = [c; zeros(T, N1 - N)]   #fix c, if new MOI.Interval variable added
    end

    return LP(c, A, b, G, g, d, u, N1, M, J, 1)

end


function LP2MOI(P::LP{T}) where {T}
    #(; c, A, b, G, g, d, u, N, M, J) = P
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

    H = MOIU.Model{T}()
    x = MOI.add_variables(H, N)

    # min c′z
    MOI.set(H, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(H, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}(), MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(c, x), 0.0))

    for m in 1:M
        MOI.add_constraint(H, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(A[m, :], x), 0.0), MOI.EqualTo(b[m]))
    end

    for j in 1:J
        MOI.add_constraint(H, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(G[j, :], x), 0.0), MOI.LessThan(g[j]))
    end

    for k in 1:N
        MOI.add_constraint(H, x[k], MOI.Interval(d[k], u[k]))
    end

    #print(H)
    return H
end


function MOI2QP(P::MOIU.Model{T}, tol=2^-26) where {T}

    N = MOI.get(P, MOI.NumberOfVariables())
    f = MOI.get(P, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{T}}())
    #display(f)

    # f = (1/2)z′Vz+q′z     the 0.5 factor in front of the Q matrix is a common source of bugs
    V = zeros(T, N, N)
    #=
    for term in f.quadratic_terms
        V[term.variable_1.value, term.variable_2.value] = term.coefficient
    end
    for i in 2:N
        for j in 1:i
            V[i, j] = V[j, i]
        end
    end
    #display(V)
    =#
    for term in f.quadratic_terms
        V[term.variable_1.value, term.variable_2.value] += term.coefficient     #duplicate_terms
    end
    for i in 2:N
        for j in 1:i-1
            V[i, j] += V[j, i]  #Duplicate off-diagonal terms
            V[j, i] = V[i, j]   #make sure symmetric
        end
    end

    q = zeros(T, N)
    for term in f.affine_terms
        #q[term.variable.value] = 2 * term.coefficient
        #q[term.variable.value] = term.coefficient
        q[term.variable.value] += term.coefficient  #duplicate_terms
    end
    #display(q)

    A, b, G, g, d, u, M, J = getConstraints(P, N, tol)  #repacking AffineFunction range rows. No need to touch V and q

    return QP(V, A, G, q, b, g, d, u, N, M, J, 1)
end


function QP2MOI(P::QP{T}) where {T}
    #(; V, A, G, q, b, g, d, u, N, M, J) = P
    V = P.V
    A = P.A
    G = P.G
    q = P.q
    b = P.b
    g = P.g
    d = P.d
    u = P.u
    N = P.N
    M = P.M
    J = P.J

    H = MOIU.Model{T}()
    x = MOI.add_variables(H, N)

    # min (1/2)z′Vz+q′z
    MOI.set(H, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    qt = Vector{MOI.ScalarQuadraticTerm{T}}(undef, 0)
    for i in 1:N
        for j in i:N    #upper triangle
            push!(qt, MOI.ScalarQuadraticTerm(V[i, j], x[i], x[j]))
        end
    end
    #at = MOI.ScalarAffineTerm.(q / 2, x)
    at = MOI.ScalarAffineTerm.(q, x)
    f = MOI.ScalarQuadraticFunction(qt, at, 0.0)
    MOI.set(H, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{T}}(), f)
    #display(f)

    for m in 1:M
        MOI.add_constraint(H, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(A[m, :], x), 0.0), MOI.EqualTo(b[m]))
    end

    for j in 1:J
        MOI.add_constraint(H, MOI.ScalarAffineFunction(MOI.ScalarAffineTerm.(G[j, :], x), 0.0), MOI.LessThan(g[j]))
    end

    for k in 1:N
        MOI.add_constraint(H, x[k], MOI.Interval(d[k], u[k]))
    end

    #print(H)
    return H

end
