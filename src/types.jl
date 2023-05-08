


"""

        @enum Status

Status: assets go IN/OUT(DN or UP), or inequalities go binded (EO, as equality) or not (OE, original ineq), with fields:

            IN  #within the lower and upper bound
            DN  #down, lower bound
            UP  #upper bound
            OE  #original <= not active
            EO  #edge, <= as =

"""
@enum Status begin
    IN
    DN  #down, lower bound
    UP  #upper bound
    OE  #original <=, not active
    EO  #edge, <= as =, active
end



"""

        struct Event{T<:AbstractFloat}

Events that assets go IN/OUT(DN or UP), or inequalities go binded (EO, as equality) or not (OE, original ineq), with fields:

            From::Status
            To::Status
            id::Int     #asset ID
            L::T        #L

"""
struct Event{T<:AbstractFloat}
    From::Status
    To::Status
    id::Int
    L::T
end



"""
        LP(c::Vector{T}, A::Matrix{T}, b::Vector{T}) where {T}

The SSLP takes the following form
min   v=c′x
s.t.  Ax=b  ∈R^{M}
      Gx≤g  ∈R^{J}
      d≤x≤u ∈R^{N}

if free, d=-Inf, u=Inf
"""
struct LP{T<:AbstractFloat}    #standard LP, or structure of LP
    c::Vector{T}
    A::Matrix{T}
    b::Vector{T}
    G::Matrix{T}
    g::Vector{T}
    d::Vector{T}
    u::Vector{T}
    N::Int
    M::Int
    J::Int
end


"""
default SSLP taks the following form
min   v=c′x
s.t.  Ax=b  ∈R^{M}
      x≥0   ∈R^{N}
"""
function LP(c::Vector{T}, A::Matrix{T}, b::Vector{T}; N=length(c),
    u=fill(Inf, N),
    d=zeros(N),
    G=ones(0, N),
    g=ones(0)) where {T}

    M = length(b)
    J = length(g)

    (M, N) == size(A) || throw(DimensionMismatch("incompatible dimension: A"))
    (J, N) == size(G) || throw(DimensionMismatch("incompatible dimension: G"))
    N == size(d, 1) || throw(DimensionMismatch("incompatible dimension: d"))
    N == size(u, 1) || throw(DimensionMismatch("incompatible dimension: u"))

    #check feasibility and redundancy of Ax=b
    rb = rank([A b])
    @assert rb == rank(A) "infeasible: Ax=b"
    #@assert M == length(getRows(A, tolN)) "redundant rows in Ax=b"   #full row rank
    @assert M == rb "redundant rows in Ax=b"       #full row rank

    @assert !any(d .== u) "downside bound == upper bound detected"
    @assert J > 0 || any(isfinite.(d)) || any(isfinite.(u)) "no inequalities and bounds"

    iu = u .< d
    if sum(iu) > 0
        @warn "swap the elements where u < d, to make sure u > d"
        t = u[iu]
        u[iu] .= d[iu]
        d[iu] .= t
    end

    LP{T}(c, A, b, G, g, d, u, N, M, J)
end


"""

        QP(V::Matrix{T}; q, u, d, G, g, A, b) where T
        QP(P::Problem{T}, L::T) where T
        QP(mu::T, P::Problem{T}) where T

Setup a quadratic programming model:

```math
    min   (1/2)z′Vz+q′z
    s.t.   Az=b ∈ R^{M}
           Gz≤g ∈ R^{J}
           d≤z≤u ∈ R^{N}
```

some variable may be free, say -Inf < zi < +Inf. No equalities if M=0. Default values: q = 0, u = +∞, d = 0, G = [], g = [], A = ones(1,N), b = [1]

See also [`Problem`](@ref), [`solveQP`](@ref)

"""
struct QP{T<:AbstractFloat}    #standard QP, or structure of QP
    V::Matrix{T}
    A::Matrix{T}
    G::Matrix{T}
    q::Vector{T}
    b::Vector{T}
    g::Vector{T}
    d::Vector{T}
    u::Vector{T}
    N::Int
    M::Int
    J::Int
end

function QP(V::Matrix{T}; N=size(V, 1), #N=convert(Int32, size(V, 1)),
    q=zeros(N),
    u=fill(Inf, N),
    d=zeros(N),
    G=ones(0, N),
    g=ones(0),
    A=ones(1, N),
    b=ones(1)) where {T}

    M = length(b)
    J = length(g)

    (N, N) == size(V) || throw(DimensionMismatch("incompatible dimension: V"))
    V = convert(Matrix{T}, (V + V') / 2)   #make sure symmetric
    @assert det(V) >= 0 "variance matrix has negative determinant"
    (M, N) == size(A) || throw(DimensionMismatch("incompatible dimension: A"))
    (J, N) == size(G) || throw(DimensionMismatch("incompatible dimension: G"))
    N == size(q, 1) || throw(DimensionMismatch("incompatible dimension: q"))
    N == size(d, 1) || throw(DimensionMismatch("incompatible dimension: d"))
    N == size(u, 1) || throw(DimensionMismatch("incompatible dimension: u"))

    #check feasibility and redundancy of Ax=b
    rb = rank([A b])
    @assert rb == rank(A) "infeasible: Ax=b"
    @assert M == rb "redundant rows in Ax=b"       #full row rank

    iu = u .< d
    if sum(iu) > 0
        @warn "swap the elements where u < d, to make sure u > d"
        t = u[iu]
        u[iu] .= d[iu]
        d[iu] .= t
    end
    #to do: J+ (num of finite d u) > 0  , when LP is introduce

    @assert !any(d .== u) "downside bound == upper bound detected"
    @assert J > 0 || any(isfinite.(d)) || any(isfinite.(u)) "no any inequalities or bounds"

    QP{T}(V, convert(Matrix{T}, copy(A)),
        convert(Matrix{T}, copy(G)),
        convert(Vector{T}, copy(q)),
        convert(Vector{T}, copy(vec(b))),
        convert(Vector{T}, copy(vec(g))),
        convert(Vector{T}, copy(d)),
        convert(Vector{T}, copy(u)), N, M, J)
end


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
    rule::Symbol    #rule for Simplex {:Dantzig, :maxImprovement}
end

Settings(; kwargs...) = Settings{Float64}(; kwargs...)

function Settings{Float64}(; maxIter=7777,
    tol=2^-26,
    tolN=2^-26,
    tolG=2^-27,
    pivot=:column, rule=:Dantzig)
    Settings{Float64}(maxIter, tol, tolN, tolG, pivot, rule)
end

function Settings{BigFloat}(; maxIter=7777,
    tol=BigFloat(2)^-76,
    tolN=BigFloat(2)^-76,
    tolG=BigFloat(2)^-77,
    pivot=:column, rule=:Dantzig)
    Settings{BigFloat}(maxIter, tol, tolN, tolG, pivot, rule)
end

