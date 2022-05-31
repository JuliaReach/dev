# https://easychair.org/publications/paper/QrVj
# Chen, X., Sankaranarayanan, S., & Abrahám, E. (2015, December).
# Flow* 1.2: More Effective to Play with Hybrid Systems. In ARCH@ CPSWeek (pp. 152-159).

using ReachabilityAnalysis, TaylorModels

# expand e^(At) into a Taylor series up to order k
function _exp_with_error(A, δ; order)
    t = Taylor1(order+1)

    n = size(A, 1)
    Iₙ = Matrix(1.0I, n, n)
    At = A*t
    W = Iₙ + At
    Q = copy(At) # stores A^i t^i / i!

    for i in 2:order
        Q = Q * At / i
        W += Q
    end

    # here Q = A^k t^k / k!
    Anorm = opnorm(A, Inf) # maximum norm
    ρ = exp(Anorm * δ)
    Mρ = IntervalArithmetic.Interval(-ρ, ρ)

    IΦ = broadcast(x -> x * Mρ, Q*At/(order+1))

    return W + IΦ
end

# recursive version of the Chen linear algorithm for linear ODEs
function _chen2015_recur(X0, A, δ, nsteps; orderT=2, orderQ=8)
    
    n = size(A, 1)
    C = Vector{Vector{TaylorN{IntervalArithmetic.Interval{Float64}}}}(undef, n)
    M = _exp_with_error(A, δ, order=orderT)

    R = ReachSet(X0, 0 .. 0)
    RTM0 = overapproximate(R, TaylorModelReachSet, orderQ=orderQ, orderT=orderT)
    X0_tm = set(RTM0)
    
    Y0 = [X.pol.coeffs[1] for X in X0_tm]
    X0 = [sum([TaylorModel1(Taylor1([TaylorN(Mi * Y0[i].coeffs)
                            for Mi in M[j,i].coeffs]), 0..0, 0..0, 0..δ) for i in 1:n]) for j in 1:n]
    RTM1 = TaylorModelReachSet(X0, 0..δ)
    
    fpv = [RTM1] 
    for i = 2:nsteps
        X1_tm = set(fpv[end])
        X1 = [sum([(X1_tm[i].pol * M[i,j]) for j in 1:n]) for i in 1:n]
        X1v = [TaylorModel1(Taylor1(Xi.coeffs), 0..0, 0..0, 0..δ) for Xi in X1]
        RTM2 = TaylorModelReachSet(X1v, (i-1)*δ..i*δ)
        push!(fpv, RTM2)
    end
    
    return Flowpipe(fpv)
end

# iterative version of the Chen linear algorithm for linear ODEs
function _chen2015_iter(X0, A, δ, nsteps; orderT=2, orderQ=8)
    
    n = size(A, 1)
    C = Vector{Vector{TaylorN{IntervalArithmetic.Interval{Float64}}}}(undef, n)
    
    R = ReachSet(X0, 0 .. 0)
    RTM0 = overapproximate(R, TaylorModelReachSet, orderQ=orderQ, orderT=orderT)
    X0_tm = set(RTM0)
    Y0 = [X.pol.coeffs[1] for X in X0_tm]
    
    M = _exp_with_error(A, δ, order=orderT)
    X0 = [sum([TaylorModel1(Taylor1([TaylorN(Mi * Y0[i].coeffs) 
                            for Mi in M[j,i].coeffs]), 0..0, 0..0, 0..δ) for i in 1:n]) for j in 1:n]
    RTM1 = TaylorModelReachSet(X0, 0..δ)
    
    fpv = [RTM1]
    for i = 2:nsteps
        M = _exp_with_error(A, i*δ, order=orderT)
        X0 = [sum([TaylorModel1(Taylor1([TaylorN(Mi * Y0[i].coeffs) 
                            for Mi in M[j,i].coeffs]), 0..0, 0..0, 0..δ) for i in 1:n]) for j in 1:n]
        RTM1 = TaylorModelReachSet(X0, (i-1)*δ..i*δ)
        push!(fpv, RTM1)
    end
    
    return Flowpipe(fpv)
    
end

# test:

X0 = BallInf(ones(2), 0.1)
A = [0 1; -1 0.]
δ = 0.01

orderT = 4
orderQ = 8;

fp_recur = _chen2015_recur(X0, A, δ, 50);

fp_iter = _chen2015_iter(X0, A, δ, 50);

prob = @ivp(x'=Ax, x(0) ∈ X0)
@time sol = solve(prob, alg=TMJets(orderQ=orderQ, orderT=orderT), tspan=0..0.5);

using Plots

plot(sol, vars=(1,2), lab="TMJets")

plot(fp_recur, vars=(0,1), lab="recursive")
plot!(fp_iter, vars=(0,1), lab="iterative")
plot!(sol, vars=(0,1), lab="TMJets")

plot(fp_recur, vars=(0,2), lab="recursive")
plot!(fp_iter, vars=(0,2), lab="iterative")
plot!(sol, vars=(0,2), lab="TMJets")

# FIXME:
#plot(fp_recur, vars=(1,2), lab="recursive")
#plot(fp_iter, vars=(1,2), lab="iterative")