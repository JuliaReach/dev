using ReachabilityAnalysis, TaylorModels

# expand e^(At) into a Taylor series
# up to order k
function _exp_with_error(A, δ; order)
    t = Taylor1(order+1)
    #t = TaylorModel1(order+1, 0, 0..1)

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

    # IΦ = Q * At / (k+1) * Mρ # FIXME
    IΦ = broadcast(x -> x * Mρ, Q*At/(order+1))

    return W + IΦ
end
