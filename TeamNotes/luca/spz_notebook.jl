### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ b2ba0168-16df-11ed-3969-77cc5e20a6c4
begin
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add([Pkg.PackageSpec(name="LazySets", rev="lf-spz-qm"),
			 Pkg.PackageSpec(name="ReachabilityAnalysis"),
		     Pkg.PackageSpec(name="Plots"),
		     Pkg.PackageSpec(name="Symbolics")
			])
	using LazySets, Plots, Symbolics, ReachabilityAnalysis
	using ReachabilityAnalysis: _new_variables,  _jacobian_expression, _jacobian_function, _hessian_expression, _hessian_function, _exp, Φ₁
end

# ╔═╡ 956923fc-1c1b-488a-a5ef-03f938f6f026
using LinearAlgebra

# ╔═╡ b1f14c83-8256-4f1c-96cb-7082efc0fbec
md"""
## Example 4.1.1
"""

# ╔═╡ c7bad52c-e6a8-4972-a04c-6a8cc33e5e6f
let
	function f!(dx, x, p, t)
	    dx[1] = -x[1] + x[1]^2
	    dx
	end
	X0 = convert(SparsePolynomialZonotope, Interval(-1, 1))
	Δt = 1.0
	prob = @ivp(x' = f!(x), dim=1, x(0) ∈ X0)
	var = _new_variables(1)
	
	J_fun = _jacobian_function(prob, var)
	H_fun = _hessian_function(prob, var)
	
	xc = LazySets.center(X0)
	dx = similar(xc)
	zs = xc + 1/2*Δt*f!(dx, xc, 0, 0)
	A = J_fun(zs)
	Q = [Float64.(hfun(zs)) for hfun in H_fun]
	
	QM = QuadraticMap(Q, X0)
	V0 = linear_map(0.5, overapproximate(QuadraticMap(Q, X0), SparsePolynomialZonotope))
	Y = linear_map(_exp(A, Δt), X0)
	Γ = Φ₁(A, Δt)
	W = linear_map(Γ, V0)
	out = exact_sum(Y, W) # 0.368α₁ + 0.632α₁²
end

# ╔═╡ 0ff014ac-f8d3-43d7-b55e-8184772af792
md"""
## Example 4.1.2
"""

# ╔═╡ adc65f1f-bdf0-46be-b09a-8a794df9fbe5
let
	function f!(dx, x, p, t)
		dx[1] = x[1] * x[2]
		dx[2] = (1 - x[1]^2) * x[2] - x[1]
		return dx
	end
	Δt = 1.0
	# X0 = SparsePolynomialZonotope([0.2, -0.6], [1 0;0 0.4], zeros(2, 0), [1 0;0 1])

	X0 = convert(SparsePolynomialZonotope, Hyperrectangle(low=[-1.4, 0.6], high=[0.6, 1.4]))
	prob = @ivp(x' = f!(x), dim=2, x(0) ∈ X0)
	var = _new_variables(2)
	J_fun = _jacobian_function(prob, var)
	H_fun = _hessian_function(prob, var)

	xc = LazySets.center(X0)
	dx = similar(xc)
	zs = xc + 1/2*Δt*f!(dx, xc, 0, 0)
	R = translate(X0, -zs) # translate
	A = J_fun(zs)
	Q = [Float64.(hfun(zs)) for hfun in H_fun]

	QM = QuadraticMap(Q, X0)
	V0 = linear_map(0.5, overapproximate(QuadraticMap(Q, R), SparsePolynomialZonotope))
	cnew = (f!(dx, zs, 0, 0) - A * zs) + LazySets.center(V0)
	V0 = SparsePolynomialZonotope(cnew, genmat_dep(V0), genmat_indep(V0), expmat(V0))
	
	Y = linear_map(_exp(A, Δt), X0)
	Γ = Φ₁(A, Δt)
	W = linear_map(Γ, V0)
	out_mink = minkowski(Y, W)
	out_exact = exact_sum(Y, W)
end

# ╔═╡ b39685c9-3a2c-44d5-95e4-8b57c5e8654a
md"""
## Vanderpool example (from paper)
"""

# ╔═╡ dc442c86-4f26-4967-8fc8-a401642de779


# ╔═╡ 46ac8a14-7b1a-482e-a503-bf2aa73e49b5


# ╔═╡ Cell order:
# ╠═b2ba0168-16df-11ed-3969-77cc5e20a6c4
# ╠═956923fc-1c1b-488a-a5ef-03f938f6f026
# ╟─b1f14c83-8256-4f1c-96cb-7082efc0fbec
# ╠═c7bad52c-e6a8-4972-a04c-6a8cc33e5e6f
# ╠═0ff014ac-f8d3-43d7-b55e-8184772af792
# ╠═adc65f1f-bdf0-46be-b09a-8a794df9fbe5
# ╟─b39685c9-3a2c-44d5-95e4-8b57c5e8654a
# ╠═dc442c86-4f26-4967-8fc8-a401642de779
# ╠═46ac8a14-7b1a-482e-a503-bf2aa73e49b5
