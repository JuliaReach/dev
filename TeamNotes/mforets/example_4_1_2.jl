using ReachabilityAnalysis, Symbolics
using ReachabilityAnalysis: _new_variables,  _jacobian_expression, _jacobian_function,
                            _hessian_expression, _hessian_function
using Plots

function vanderpol!(dx, x, p, t)
    dx[1] = x[1]*x[2]
    dx[2] = (1-x[1]^2)*x[2] - x[1]
end

X0 = Hyperrectangle(low=[-1.4, 0.6], high=[0.6, 1.4])
Δt = 1.0
prob = @ivp(x' = vanderpol!(x), dim=2, x(0) ∈ X0)

sol = solve(prob, tspan=0..3, alg=TMJets()) # to compare



# sol = solve(prob, tspan=0..3, alg=KA20()) # REF [47] in the thesis

var = _new_variables(2)

J_ex = _jacobian_expression(prob, var)
J_fun = _jacobian_function(prob, var)
J_fun([0.0, 0.0]), J_fun([1.0, 0.0])

H_ex = _hessian_expression(prob, var)
H_fun = _hessian_expression(prob, var)

# TO-DO: Repeat same steps as in example 4.1.1
# trying to reproduce the Figure 4.2.

