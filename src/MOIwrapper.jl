
using TOML

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    #struct Optimizer <: MOI.AbstractOptimizer
    # Fields go here
    Problem::Union{Nothing,StatusSwitchingQP.QP{T},StatusSwitchingQP.LP{T}}
    Settings::StatusSwitchingQP.Settings{T}
    Results::Union{Nothing,Tuple{Vector{T},Vector{Status},Int}}
    Sense::MOI.OptimizationSense
    Silent::Bool
    f0::T

    function Optimizer{T}(; user_settings...) where {T}
        Problem = nothing
        Settings = StatusSwitchingQP.Settings{T}()
        Results = nothing
        Sense = MOI.MIN_SENSE
        Silent = true
        f0 = T(0)
        opt = new(Problem, Settings, Results, Sense, Silent, f0)
        #= for (key, value) in user_settings
            MOI.set(opt, MOI.RawOptimizerAttribute(string(key)), value)
        end =#
        if isempty(user_settings)  #default setting
            return opt
        end
        opt.Settings = StatusSwitchingQP.Settings{T}(; user_settings...)
        return opt
    end
end

Optimizer(args...; kwargs...) = Optimizer{Float64}(args...; kwargs...)

#=
All Optimizers must implement the following methods:

empty!
is_empty

=#


function MOI.empty!(opt::Optimizer{T}) where {T}
    #flush everything, keeping the currently configured settings
    opt.Problem = nothing
    opt.Settings = opt.Settings #preserve settings / no change
    opt.Results = nothing
    opt.Sense = MOI.MIN_SENSE # model parameter, so needs to be reset
end

MOI.is_empty(opt::Optimizer{T}) where {T} = isnothing(opt.Problem)



# Solver Attributes, get/set

MOI.get(opt::Optimizer, ::MOI.SolverName) = "StatusSwitchingQP"
MOI.get(opt::Optimizer, ::MOI.SolverVersion) = TOML.parsefile(joinpath(pkgdir(StatusSwitchingQP), "Project.toml"))["version"]
MOI.get(opt::Optimizer, ::MOI.RawSolver) = StatusSwitchingQP.solveQP
MOI.get(opt::Optimizer, ::MOI.ResultCount) = Int(!isnothing(opt.Results))
MOI.get(opt::Optimizer, ::MOI.SolveTimeSec) = 0.1


MOI.supports(::Optimizer, ::MOI.Silent) = true
MOI.get(opt::Optimizer, ::MOI.Silent) = opt.Silent
MOI.set(opt::Optimizer, ::MOI.Silent, v::Bool) = (opt.Silent = v)



#=
function Base.show(io::IO, model::Optimizer)
    return print(io, "NewSolver with the pointer $(model.ptr)")
end
=#


#=
For each attribute

get gets the current value of the attribute
set sets a new value of the attribute. Not all attributes can be set. For example, the user can't modify the SolverName.
supports returns a Bool indicating whether the solver supports the attribute.

=#


#needed by JuMP
#=
function MOI.supports(::Optimizer{T},
    ::Union{MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}},MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{T}}}) where {T}
    return true
end
#MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{Quadratic}) = true
function MOI.supports(::Optimizer,
    ::Type{<:Union{MOI.ObjectiveFunction{MOI.ScalarAffineFunction},MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction}}})
    return true
end

function MOI.supports(::Optimizer,
    ::Type{MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction}})
    return true
end
=#

#MOI.supports(::Optimizer, ::Type{MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction}}) = true

MOI.supports(::Optimizer{T}, ::MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{T}}) where {T} = true
MOI.supports(::Optimizer{T}, ::MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}) where {T} = true

#=
function MOI.supports_constraint(::Optimizer{T},
    ::Union{MOI.VariableIndex,MOI.ScalarAffineFunction{T}},
    ::Union{MOI.LessThan{T},MOI.GreaterThan{T},MOI.EqualTo{T},MOI.Interval{T}}) where {T}
    return true
end
=#



function MOI.supports_constraint(
    ::Optimizer{T},
    ::Type{MOI.ScalarAffineFunction{T}},
    #::Type{<:Union{MOI.LessThan{T},MOI.GreaterThan{T},MOI.EqualTo{T},MOI.Interval{T}}}) where {T}
    ::Type{<:Union{MOI.LessThan{T},MOI.GreaterThan{T},MOI.EqualTo{T}}}) where {T}
    return true
end

function MOI.supports_constraint(
    ::Optimizer{T},
    ::Type{MOI.VariableIndex},
    #::Type{<:Union{MOI.LessThan{T},MOI.GreaterThan{T},MOI.EqualTo{T},MOI.Interval{T}}}) where {T}
    ::Type{<:Union{MOI.LessThan{T},MOI.GreaterThan{T},MOI.Interval{T}}}) where {T}
    return true
end




#=
If your solver separates data loading and the actual optimization into separate steps, implement the copy_to interface.

copy_to(::ModelLike, ::ModelLike)
optimize!(::ModelLike)

All Optimizers must implement the following attributes:

DualStatus
PrimalStatus
RawStatusString
ResultCount
TerminationStatus

You should also implement the following attributes:

ObjectiveValue
SolveTimeSec
VariablePrimal
=#


function MOI.copy_to(dest::Optimizer{T}, src::MOI.ModelLike) where {T}

    #idxmap = MOIU.IndexMap(dest, src)
    #variables, map = variablesMap(src)
    dest.Sense = MOI.get(src, MOI.ObjectiveSense())
    map, dest.Problem = MOI2QP(dest, src)
    if norm(dest.Problem.V, Inf) == 0   #LP
        #map, dest.Problem = MOI2LP(dest, src)
        Q = dest.Problem
        #dest.Problem = LP(Q.q, Q.A, Q.b, Q.G, Q.g, Q.d, Q.u, Q.N, Q.M, Q.J)
        dest.Problem = LP(Q.q, Q.A, Q.b; G=Q.G, g=Q.g, d=Q.d, u=Q.u)
    end
    return map
end

#=
function MOIU.IndexMap(dest::Optimizer, src::MOI.ModelLike)

    idxmap = MOIU.IndexMap()

    vis_src = MOI.get(src, MOI.ListOfVariableIndices())
    for i in eachindex(vis_src)
        idxmap[vis_src[i]] = MOI.VariableIndex(i)
    end
    i = 0
    for (F, S) in MOI.get(src, MOI.ListOfConstraintTypesPresent())
        MOI.supports_constraint(dest, F, S) || throw(MOI.UnsupportedConstraint{F,S}())
        cis_src = MOI.get(src, MOI.ListOfConstraintIndices{F,S}())
        for ci in cis_src
            i += 1
            idxmap[ci] = MOI.ConstraintIndex{F,S}(i)
        end
    end

    return idxmap
end
=#

#=
function variablesMap(src::MOI.ModelLike)
    #https://github.com/jump-dev/GLPK.jl/blob/master/src/MOI_wrapper/MOI_copy.jl#L79  _init_index_map
    variables = MOI.get(src, MOI.ListOfVariableIndices())
    map = MOIU.IndexMap()
    k = 0
    for x in variables
        k += 1
        map[x] = MOI.VariableIndex(k)
    end
    return variables, map
end
=#


function MOI.optimize!(opt::Optimizer{T}) where {T}

    P = opt.Problem
    if P.mc == -20
        N = P.N
        if P.M > 0
            x = P.A \ P.b
            opt.Results = (x, fill(DN, N), 1)
        else  #no constraints
            if typeof(P) <: QP
                o = norm(P.V, Inf) == 0 && norm(P.q, Inf) == 0
            else
                o = norm(P.c, Inf) == 0
            end
            opt.Results = (zeros(T, N), fill(DN, N), o ? 1 : 3)
        end
        return nothing
    end

    #opt.Results = solveQP(opt.Problem; settings=opt.Settings, settingsLP=settings) #error
    #opt.Results = solveQP(opt.Problem; settings=opt.Settings, settingsLP=opt.Settings)

    #opt.Results = solveQP(opt.Problem; settings=opt.Settings)
    if typeof(opt.Problem) <: QP
        opt.Results = solveQP(opt.Problem; settings=opt.Settings)
    else
        #min = opt.Sense == MOI.MAX_SENSE ? false : true
        #opt.Results = SimplexLP(opt.Problem; settings=opt.Settings, min=min)    # open-intervals, supports `[l, Inf)`, `(-Inf, u]`,  and `(-Inf, Inf)`
        #opt.Results = solveLP(opt.Problem; settings=opt.Settings)
        opt.Results = SimplexLP(opt.Problem; settings=opt.Settings)
    end

    nothing
end


MOI.supports(::Optimizer, ::MOI.TerminationStatus) = true
function MOI.get(opt::Optimizer, ::MOI.TerminationStatus)
    !isnothing(opt.Results) || return MOI.OPTIMIZE_NOT_CALLED
    # > 0 if successful (=iter_count); = 0 if infeasibility detected; < 0 fail (=-1 if numerical error, =-maxIter if not converged)
    st = opt.Results[3]
    if st == 3
        return MOI.INFEASIBLE_OR_UNBOUNDED
    elseif st == 1 || st == 2
        return MOI.OPTIMAL
    elseif st == 0
        return MOI.INFEASIBLE
    elseif st == -1
        return MOI.NUMERICAL_ERROR
    else
        return MOI.ITERATION_LIMIT
    end
end

#=
MOI.OPTIMAL,
    MOI.ITERATION_LIMIT,
    MOI.TIME_LIMIT,
    MOI.INFEASIBLE,
    MOI.DUAL_INFEASIBLE,
    MOI.ALMOST_OPTIMAL,
    MOI.NUMERICAL_ERROR,
    MOI.NUMERICAL_ERROR
=#

function MOI.get(opt::Optimizer, a::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(opt, a)
    x = opt.Results[1]
    if typeof(opt.Problem) <: QP
        f = x' * (opt.Problem.V * x) / 2 + opt.Problem.q' * x
    else
        f = x' * opt.Problem.c
    end
    return (opt.Sense == MOI.MIN_SENSE ? f : -f) + opt.f0
end

MOI.supports(::Optimizer, ::MOI.DualObjectiveValue) = false



MOI.supports(::Optimizer, ::MOI.VariablePrimal) = true
function MOI.get(opt::Optimizer, a::MOI.VariablePrimal, vi::MOI.VariableIndex)
    MOI.check_result_index_bounds(opt, a)
    return opt.Results[1][vi.value]
end
function MOI.get(opt::Optimizer, a::MOI.VariablePrimal)
    MOI.check_result_index_bounds(opt, a)
    return opt.Results[1]
end

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = false

MOI.supports(::Optimizer, a::MOI.DualStatus) = true
#MOI.get(::Optimizer, ::MOI.DualStatus) = MOI.NO_SOLUTION
#MOI.get(::Optimizer, ::MOI.DualStatus) = MOI.FEASIBLE_POINT
function MOI.get(::Optimizer, attr::MOI.DualStatus)
    return attr.result_index == 1 ? MOI.FEASIBLE_POINT : MOI.NO_SOLUTION
end





MOI.supports(::Optimizer, ::MOI.PrimalStatus) = true
function MOI.get(opt::Optimizer, attr::MOI.PrimalStatus)
    if attr.result_index != 1    #988 passed, 4 failed, 127 errored
        return MOI.NO_SOLUTION
    end
    !isnothing(opt.Results) || return MOI.NO_SOLUTION
    st = opt.Results[3]
    if st == 0
        return MOI.INFEASIBLE_POINT
    else
        return MOI.FEASIBLE_POINT
    end
end

MOI.get(opt::Optimizer, ::MOI.RawStatusString) = string(opt.Results[3])

MOI.supports(::Optimizer, ::MOI.ConstraintDual) = false
#= MOI.supports(::Optimizer, ::MOI.ConstraintDual) = true
function MOI.get(opt::Optimizer, a::MOI.ConstraintDual, ci::MOI.ConstraintIndex)

    return opt.Results[1][1]
end =#



function getConstraints(P, N, T)
    #function getConstraints(P, N, tol, T)
    #function getConstraints(P::MOIU.Model{T}, N, tol) where {T}

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
                    #= elseif S <: MOI.Interval
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
                        end =#
                else
                    throw(MOI.UnsupportedConstraint{F,S}())
                end
            elseif F <: MOI.VariableIndex
                #= if S <: MOI.EqualTo
                    d[f.value] = s.value
                    u[f.value] = s.value
                    #display((F,S))
                else =#
                if S <: MOI.GreaterThan
                    d[f.value] = s.lower
                elseif S <: MOI.LessThan
                    u[f.value] = s.upper
                elseif S <: MOI.Interval
                    d[f.value] = s.lower
                    u[f.value] = s.upper
                else
                    throw(MOI.UnsupportedConstraint{F,S}())
                end
            else
                throw(MOI.UnsupportedConstraint{F,S}())
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

function MOI2LP(dest::Optimizer{T}, MP) where {T}
    #function MOI2LP(P::MOIU.Model{T}) where {T}
    #function MOI2LP(P::MOIU.Model{T}, tol=2^-26) where {T}

    P = MOIU.Model{T}()
    map = MOI.copy_to(P, MP)

    #N = P.constraints.num_variables   #for MOI object only
    N = MOI.get(P, MOI.NumberOfVariables())
    f = MOI.get(P, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{T}}())
    c = zeros(T, N)
    for term in f.terms
        c[term.variable.value] = term.coefficient
    end

    if dest.Sense == MOI.MAX_SENSE
        c = -c
    end

    A, b, G, g, d, u, M, J = getConstraints(P, N, T) #getConstraints(P, N, tol, T)

    #return c, A, b, G, g, d, u
    #return LP(c, A, b, G, g, d, u, N, M, J)
    return map, LP(c, A, b; G=G, g=g, d=d, u=u)
end


function LP2MOI(P::LP{T}) where {T}
    (; c, A, b, G, g, d, u, N, M, J) = P

    #= if VERSION ≥ v"1.7.0"
        (; c, A, b, G, g, d, u, N, M, J) = P
    else
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
    end =#

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


function MOI2QP(dest::Optimizer{T}, MP) where {T}
    #function MOI2QP(P::MOIU.Model{T}, tol=2^-26) where {T}
    #tol = dest.Settings.tol
    P = MOIU.Model{T}()
    map = MOI.copy_to(P, MP)

    N = MOI.get(P, MOI.NumberOfVariables())
    f = MOI.get(P, MOI.ObjectiveFunction{MOI.ScalarQuadraticFunction{T}}())
    #display(f)

    # f = (1/2)z′Vz+q′z     the 0.5 factor in front of the Q matrix is a common source of bugs
    V = zeros(T, N, N)
    for term in f.quadratic_terms
        V[term.variable_1.value, term.variable_2.value] += term.coefficient     #duplicate_terms
    end
    for i in 2:N
        for j in 1:i-1
            V[i, j] += V[j, i]
            V[j, i] = V[i, j]   #make sure symmetric
        end
    end
    #V = (V + V') / 2    #Duplicate off-diagonal terms
    #display(V)

    q = zeros(T, N)
    for term in f.affine_terms
        #q[term.variable.value] = 2 * term.coefficient
        q[term.variable.value] += term.coefficient  ##duplicate_terms
    end
    #display(q)

    dest.f0 = f.constant

    if dest.Sense == MOI.MAX_SENSE
        V = -V
        q = -q
    end

    A, b, G, g, d, u, M, J = getConstraints(P, N, T) #getConstraints(P, N, tol, T)

    #return QP(V, A, G, q, b, g, d, u, N, M, J)
    return map, QP(V; A=A, G=G, q=q, b=b, g=g, d=d, u=u)
end


function QP2MOI(P::QP{T}) where {T}
    (; V, A, G, q, b, g, d, u, N, M, J) = P
    #= if VERSION ≥ v"1.7.0"
        (; V, A, G, q, b, g, d, u, N, M, J) = P
    else
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
    end =#

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