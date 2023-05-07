
import MathOptInterface as MOI

#=
mps = MOI.FileFormats.Model(format=MOI.FileFormats.FORMAT_MPS)
fn = "./grow7.mps.gz"
#MOI.read_from_file(H, fn)
H = read_from_file(fn)
MOI.copy_to(mps, H)
c, A, b, G, g, d, u = MOI2LP(mps)
=#
function MOI2LP(P, tol=2^-26)

    T = Float64
    #N = P.constraints.num_variables   #for MOI object only
    N = MOI.get(P, MOI.NumberOfVariables())
    f = MOI.get(P, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}())
    c = zeros(T, N)
    for term in f.terms
        c[term.variable.value] = term.coefficient
    end

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

    display(("sizes ", N, M, J))

    for k in 1:M
        A[k, :] = Ab[k][1:end-1]
        b[k] = Ab[k][end]
    end
    for k in 1:J
        G[k, :] = Gg[k][1:end-1]
        g[k] = Gg[k][end]
    end
    return c, A, b, G, g, d, u
end

