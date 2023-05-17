#Speed and Accuracy: Quadratic pragramming
#compare: SSQP, OSQP, Clarabel

using EfficientFrontier, LinearAlgebra
using TranscodingStreams, CodecXz, Serialization, Downloads


using Statistics

if length(filter((x) -> x == :uOSQP, names(Main, imported=true))) == 0
    jlFile = "./uOSQP.jl"
    if !isfile(jlFile)
        Downloads.download("https://github.com/PharosAbad/EfficientFrontier.jl/raw/main/examples/uOSQP.jl", jlFile)
    end
    include("./uOSQP.jl")
    using .uOSQP
end

if length(filter((x) -> x == :uClarabel, names(Main, imported=true))) == 0
    jlFile = "./uClarabel.jl"
    if !isfile(jlFile)
        Downloads.download("https://github.com/PharosAbad/EfficientFrontier.jl/raw/main/examples/uClarabel.jl", jlFile)
    end
    include("./uClarabel.jl")
    using .uClarabel
end


function fPortfolio(P::Problem{T}, L::T=0.0; settings=SettingsQP(P), settingsLP=SettingsLP(P)) where {T}
    Q = QP(P, L)
    solveQP(Q; settings=settings, settingsLP=settingsLP)    #z, S, iter
end

function fPortfolio(mu::T, P::Problem{T}; settings=SettingsQP(P), settingsLP=SettingsLP(P)) where {T}
    Q = QP(mu, P)
    solveQP(Q; settings=settings, settingsLP=settingsLP)    #z, S, iter
end

function testData(ds::Symbol, pd=true)
    if ds == :Ungil
        E, V = EfficientFrontier.EVdata(:Ungil, false)
        A = [1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0
            1.0 1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0]
        b = [1.0; 0.25]
        G = [0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0
            0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0]
        G[1, :] = -G[1, :]
        g = [-0.3; 0.6]
        d = vec([-0.0 -0.0 -0.0 -0.0 -0.0 -0.0 -0.0 -0.0 -0.0 -0.0 -0.1 -0.1 -0.1 -0.1])
        u = vec([0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.1 0.1 0.1 0.3 0.3 0.3 0.3])

        P = Problem(E, V, u, d, G, g, A, b)
    elseif ds == :SP500
        xzFile = joinpath(tempdir(), "sp500.jls.xz") #xzFile = "/tmp/sp500.jls.xz"
        if !isfile(xzFile)
            Downloads.download("https://github.com/PharosAbad/PharosAbad.github.io/raw/master/files/sp500.jls.xz", xzFile)
        end
        io = open(xzFile)
        io = TranscodingStream(XzDecompressor(), io)
        E = deserialize(io)
        V = deserialize(io)
        close(io)
        N = length(E)
        u = fill(3 / 32, N)

        #pd = true
        if pd
            N = 263
            ip = 1:N
            V = V[ip, ip]
            E = E[ip]
            u = fill(3 / 32, N)
        end
        println(ds, " number of stocks: ", N)
        P = Problem(E, V, u)
    else
        error("Unknown dataset")
    end
    return P
end

function saveText(fn, Sl, So, Sc, Al, Ao, Ac, Tl, To, Tc, Ol, Oo, Oc)
    redirect_stdio(stdout=fn) do
        #status
        println("\n------- Solution status ------- SSQP/OSQP/Clarabel")
        show(stdout, "text/plain", Sl)
        println("")
        show(stdout, "text/plain", So)
        println("")
        show(stdout, "text/plain", Sc)
        println("")
        #Accuracy
        println("\n------- Accuracy -------SSQP/OSQP/Clarabel   ", round.([norm(Al, Inf), norm(Ao, Inf), norm(Ac, Inf)], sigdigits=3))
        println("---- quantile 99% ----   ", round.([quantile(Al[:], 0.99), quantile(Ao[:], 0.99), quantile(Ac[:], 0.99)], sigdigits=3))
        println("------- median -------   ", round.([median(Al[:]), median(Ao[:]), median(Ac[:])], sigdigits=3))
        println("---- quantile  1% ----   ", round.([quantile(Al[:], 0.01), quantile(Ao[:], 0.01), quantile(Ac[:], 0.01)], sigdigits=3))
        show(stdout, "text/plain", round.(Al, sigdigits=3))
        println("")
        show(stdout, "text/plain", round.(Ao, sigdigits=3))
        println("")
        show(stdout, "text/plain", round.(Ac, sigdigits=3))
        println("")
        #Speed
        println("\n--- Speed (time span, smaller for faster speed) ---SSQP/OSQP/Clarabel   ", round.([norm(Tl, Inf), norm(To, Inf), norm(Tc, Inf)], sigdigits=3))
        println("---- quantile 99% ----   ", round.([quantile(Tl[:], 0.99), quantile(To[:], 0.99), quantile(Tc[:], 0.99)], sigdigits=3))
        println("------- median -------   ", round.([median(Tl[:]), median(To[:]), median(Tc[:])], sigdigits=3))
        println("---- quantile  1% ----   ", round.([quantile(Tl[:], 0.01), quantile(To[:], 0.01), quantile(Tc[:], 0.01)], sigdigits=3))
        show(stdout, "text/plain", round.(Tl, sigdigits=3))
        println("")
        show(stdout, "text/plain", round.(To, sigdigits=3))
        println("")
        show(stdout, "text/plain", round.(Tc, sigdigits=3))
        println("")
        #Objective function
        println("\n--- Objective function value (diff in sd, not variance) ---SSQP/OSQP/Clarabel   ", round.([norm(Ol, Inf), norm(Oo, Inf), norm(Oc, Inf)], sigdigits=3))
        println("---- quantile 99% ----   ", round.([quantile(Ol[:], 0.99), quantile(Oo[:], 0.99), quantile(Oc[:], 0.99)], sigdigits=3))
        println("------- median -------   ", round.([median(Ol[:]), median(Oo[:]), median(Oc[:])], sigdigits=3))
        println("---- quantile  1% ----   ", round.([quantile(Ol[:], 0.01), quantile(Oo[:], 0.01), quantile(Oc[:], 0.01)], sigdigits=3))
        show(stdout, "text/plain", round.(Ol, sigdigits=3))
        println("")
        show(stdout, "text/plain", round.(Oo, sigdigits=3))
        println("")
        show(stdout, "text/plain", round.(Oc, sigdigits=3))
        println("")
    end
end

function SpeedAccuracy(aEF, P, QPsolver, M=16)
    println("QP solver ", QPsolver, " is running")
    V = P.V
    N = length(aEF.mu) - 1
    T = zeros(N, M)     #time used, for speed
    A = zeros(N, M)     #Accuracy
    O = zeros(N, M)     #Objective function
    S = trues(N, M)     #Solution status
    for k in 1:N
        for m = 1:M
            mu = ((M + 1 - m) * aEF.mu[k] + (m - 1) * aEF.mu[k+1]) / M
            z = ePortfolio(mu, aEF)
            if QPsolver == :SSQP
                ts = @elapsed x, Sz, status = fPortfolio(mu, P)
                st = status > 0
            elseif QPsolver == :OSQP
                ts = @elapsed y = OpSpQP(mu, P)
                st = y.info.status_val == 1
                x = y.x
            elseif QPsolver == :Clarabel
                ts = @elapsed y = ClarabelQP(mu, P)
                st = Int(y.status) == 1
                x = y.x
            else
                error("Unknown QP solver")
            end
            S[k, m] = st
            T[k, m] = ts
            A[k, m] = norm(x - z, Inf)
            O[k, m] = sqrt(x' * V * x) - sqrt(z' * V * z)

        end
    end

    return T, A, O, S
end

function cmpSA(ds::Symbol, pd=true)
    P = testData(ds, pd)
    println("--- Starting EfficientFrontier ---")
    t0 = time()
    ts = @elapsed aCL = EfficientFrontier.ECL(P)
    aEF = eFrontier(aCL, P)
    t1 = time()
    println("1st run, EfficientFrontier:  ", t1 - t0, "  seconds", "\n   aCL:  ", ts, "  seconds")
    t0 = time()
    ts = @elapsed aCL = EfficientFrontier.ECL(P)
    aEF = eFrontier(aCL, P)
    t1 = time()
    println("2nd run, EfficientFrontier:  ", t1 - t0, "  seconds", "\n   aCL:  ", ts, "  seconds")

    QPsolver = :SSQP
    Tl, Al, Ol, Sl = SpeedAccuracy(aEF, P, QPsolver)
    QPsolver = :OSQP
    To, Ao, Oo, So = SpeedAccuracy(aEF, P, QPsolver)
    QPsolver = :Clarabel
    Tc, Ac, Oc, Sc = SpeedAccuracy(aEF, P, QPsolver)

    saveText("stdout.txt", Sl, So, Sc, Al, Ao, Ac, Tl, To, Tc, Ol, Oo, Oc)
    return nothing
end

function SpeedAccuracyL(aEF, P, aCL, QPsolver, M=16)
    println("QP solver ", QPsolver, " is running")
    V = P.V
    N = length(aEF.mu) - 1
    T = zeros(N, M)     #time used, for speed
    A = zeros(N, M)     #Accuracy
    O = zeros(N, M)     #Objective function
    S = trues(N, M)     #Solution status
    for k in 1:N
        i = aEF.ic[k]
        t = aCL[i]
        for m = 1:M
            #mu = ((M + 1 - m) * aEF.mu[k] + (m - 1) * aEF.mu[k+1]) / M
            L = ((M + 1 - m) * t.L1 + (m - 1) * t.L0) / M
            z = ePortfolio(P, L, aCL)
            if QPsolver == :SSQP
                ts = @elapsed x, Sz, status = fPortfolio(P, L)
                st = status > 0
            elseif QPsolver == :OSQP
                ts = @elapsed y = OpSpQP(P, L)
                st = y.info.status_val == 1
                x = y.x
            elseif QPsolver == :Clarabel
                ts = @elapsed y = ClarabelQP(P, L)
                st = Int(y.status) == 1
                x = y.x
            else
                error("Unknown QP solver")
            end
            S[k, m] = st
            T[k, m] = ts
            A[k, m] = norm(x - z, Inf)
            O[k, m] = sqrt(x' * V * x) - sqrt(z' * V * z)
        end
    end

    return T, A, O, S
end


function cmpSA_L(ds::Symbol, pd=true)
    P = testData(ds, pd)
    println("--- Starting EfficientFrontier ---")
    t0 = time()
    ts = @elapsed aCL = EfficientFrontier.ECL(P)
    aEF = eFrontier(aCL, P)
    t1 = time()
    println("1st run, EfficientFrontier:  ", t1 - t0, "  seconds", "\n   aCL:  ", ts, "  seconds")
    t0 = time()
    ts = @elapsed aCL = EfficientFrontier.ECL(P)
    aEF = eFrontier(aCL, P)
    t1 = time()
    println("2nd run, EfficientFrontier:  ", t1 - t0, "  seconds", "\n   aCL:  ", ts, "  seconds")

    QPsolver = :SSQP
    Tl, Al, Ol, Sl = SpeedAccuracyL(aEF, P, aCL, QPsolver)
    QPsolver = :OSQP
    To, Ao, Oo, So = SpeedAccuracyL(aEF, P, aCL, QPsolver)
    QPsolver = :Clarabel
    Tc, Ac, Oc, Sc = SpeedAccuracyL(aEF, P, aCL, QPsolver)

    saveText("stdoutL.txt", Sl, So, Sc, Al, Ao, Ac, Tl, To, Tc, Ol, Oo, Oc)
    return nothing
end

function cmpSnA(ds::Symbol, pd=true)
#FP(mu=mu0), Az=b contains z′E=μ, objective function L=0
#cmpSA(:Ungil)
cmpSA(ds, pd)

#FP(L=L0), , Az=b excludes z′E=μ, objective function has -L*z′E
#cmpSA_L(:Ungil)
cmpSA_L(ds, pd)
end


function main()
    cmpSnA(:Ungil)
    #cmpSnA(:SP500)
    #cmpSnA(:SP500, false)

    nothing
end

main()

nothing


