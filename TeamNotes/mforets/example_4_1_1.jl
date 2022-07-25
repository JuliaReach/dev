using ReachabilityAnalysis, Symbolics
using ReachabilityAnalysis: _new_variables,  _jacobian_expression, _jacobian_function,
                            _hessian_expression, _hessian_function, _exp, Φ₁

function f!(dx, x, p, t)
    dx[1] = -x[1] + x[1]^2
    dx
end

X0 = Interval(-1, 1)
Δt = 1.0
prob = @ivp(x' = f!(x), dim=1, x(0) ∈ X0)

var = _new_variables(1)

J_ex = _jacobian_expression(prob, var)
J_fun = _jacobian_function(prob, var)
J_fun([0.0]), J_fun([1.0])

H_ex = _hessian_expression(prob, var)
H_fun = _hessian_expression(prob, var)
H_fun[1]

xc = LazySets.center(X0)
dx = similar(xc)
zs = xc + 1/2*Δt*f!(dx, xc, 0, 0)

A = J_fun(zs)

Q = [hcat(2.0)] # H_fun
V0 = 1/2 * overapproximate(QuadraticMap(Q, X0), Zonotope)

Y = linear_map(_exp(A, Δt), X0)
Γ = Φ₁(A, Δt)
W = linear_map(Γ, V0)
out = overapproximate(minkowski_sum(Y, W), Interval) # [-0.36788, 1]

# FIXME Repeat the same as above but using X0 as a SSPZ


# Using the exact sum
# ---------------------
#
# Y_spz = convert(SparsePolynomialZonotope, Y)
# W_spz = convert(SparsePolynomialZonotope, W)
# overapproximate(exact_sum(Y_spz, W_spz), Interval)
