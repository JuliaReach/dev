### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 1591bb74-2e89-11ed-07f8-33695ecb3632
md"""

## Reachability algorithm (without input)

```math
\begin{align}
&\mathbf{input: } \mathcal{X}_0,t_f,\Delta t,\lambda,\nu,\rho_d,p_d\mu_d\\
&\mathbf{initialization: }t_0\leftarrow0,j\leftarrow0, \mathcal{R}^{union}\leftarrow\emptyset,\mathcal{R}^o(0)\leftarrow\mathcal{X}_0,\Psi(\tau_0)\leftarrow\mathbf{0}\\
&\mathbf{while }t<t_f\mathbf{ do}\\
&\quad x_c=center(\mathcal{R}^o(t_j)),x^*=x_c+0.5\Delta t\cdot f(x_c)\\
&\quad A = \nabla f(x^*), \mathcal{Q}=\nabla^2f(x^*)\\
&\quad\mathcal{R}^d=\mathcal{R}^o(t_j)\oplus(-x^*),\quad\mathcal{R}^d_z(t_j)=zonotope(R^d(t_j))\\
&\quad\mathcal{R}^\Delta(t_j+1)=((e^{A\Delta t}-I_n)\otimes\mathcal{R}^d(t_j))\oplus\Gamma(\Delta t)f(x^*)\\
&\quad R^{s,\Delta}=\textcolor{blue}{comb}(\mathbf{0},R^{\Delta}(t_{j+1}))\oplus(\textcolor{red}{F}\otimes R_z^d(t_j))\oplus\textcolor{red}{\hat{F}}f(x^*)\\
&\quad V(t_j) = (F(x^*)-Ax^*)\oplus\frac{1}{2}sq(Q,R^d(t_j))\\
&\quad\mathbf{repeat}\\
&\qquad \bar{\Psi}(\tau_j) = (\lambda\cdot I_n)\otimes\Psi(\tau_j)\\
&\qquad R^p(\tau_j)=\textcolor{red}{reachVarInput}(\bar{\Psi}(\tau_j), A, \Delta t, \nu)\\
&\qquad R^\Delta(\tau_j)=zonotope(R^{s,\Delta}(\tau_j)\oplus R^p(\tau_j))\\
&\qquad R^o(\tau_j)=R^o(t_j)\oplus R^\Delta(\tau_j)\\
&\qquad I = \textcolor{blue}{interval}(R^p(\tau_j))\\
&\qquad D = \textcolor{red}{bound}(\nabla^3 f(x), I)\\
&\qquad L(\tau_j) = \frac{1}{6}\textcolor{red}{poly}(D, I\oplus(-x^*))\\
&\qquad V^\Delta(\tau_j)=\frac{1}{2}\textcolor{blue}{sq}(Q, R^d_z, R^\Delta)\oplus\frac{1}{2}\textcolor{blue}{sq}(Q, R^\Delta, D^d_z)\\
&\qquad \Psi(\tau_j) = \textcolor{blue}{interval}(V(t_j)\oplus V^\Delta(\tau_j))\oplus(-f(x^*)+Ax^*)\\
&\quad \mathbf{until } \Psi(\tau_j)\subseteq\bar{\Psi}(\tau_j)\\
&\quad R^{p,\Delta}(\tau_j)=\textcolor{red}{reachVarInput}(V^\Delta(\tau_j), A, \Delta t, \nu)\\
&\quad R^o(t_{j+1})=(e^{A\Delta t}\otimes R^o(t_j))\boxplus(\Gamma(\Delta t)\otimes V(t_j))\oplus R^{p, \Delta}(\tau_j)\\
&\quad R^o(t_{j+1}) = reduce(R^o(t_{j+1}, \rho_d)\\
&\quad \mathbf{if~} \textcolor{blue}{volRatio}(R^o(t_{j+1})) > \mu_d\mathbf{~then}\\
&\qquad R^o(t_{j+1}) = \textcolor{blue}{restructure}(R^o(t_{j+1}), p_d)\\
&\quad\mathbf{end if}\\
&\quad R^{union}=R^{union}\cup R^o(\tau_j)\\
&\quad t_{j+1} = t_j + \Delta t,\quad \Psi(\tau_{j+1}) = \Psi(\tau_j),\quad j += 1\\
&\mathbf{end~while}\\
\end{align}
```


"""

# ╔═╡ 4129a82e-da07-40ae-aa57-0efbc807269a
md"""
## Missing pieces
- ``\textcolor{red}{red}``: requires more discussion
- ``\textcolor{blue}{blue}``: easy / PR open

### Easy

- ``\textcolor{blue}{interval}`` (box_overapproximation): once the PR for `support_function` is merged, this should come for free. Technically this can be already achieved by going through zonotopes, but that gives poor results
- ``\textcolor{blue}{volRatio}``: trivial using ``interval``
- ``\textcolor{blue}{restructure}``: builds on top of reduce, see prop 3.1.41
- ``\textcolor{blue}{sq}``: two-input quadratic map, see 3.2.22
- ``\textcolor{blue}{comb}``: linear combination, see 3.1.26
### More discussion needed

- ``\textcolor{red}{reachVarInput~(eq 4.14)}:\bigoplus_{k=0}^\nu\frac{\Delta t^{k+1}}{(k+1)!}(A^k\otimes V^\Delta)\oplus(\Delta t\cdot\textcolor{red}{\varepsilon}\otimes V^\Delta)``
- ``\textcolor{red}{\varepsilon~(eq 4.14)}: [-1, 1]\frac{(\Vert A\Vert_\infty\Delta t)^{\nu+1}}{(\nu+1)!}\frac{\nu+2}{\nu+2-\Vert A\Vert_\infty\Delta t}``
- ``\textcolor{red}{F~(eq 4.16)}=\bigoplus_{k=2}^\nu[(k^{-\frac{k}{k-1}}-k^{-\frac{1}{k-1}})\Delta t^k, 0]\frac{A^k}{k!}\oplus\varepsilon``
- ``\textcolor{red}{\hat{F}~(eq 4.16)}=\bigoplus_{k=2}^\nu[((k+1)^{-\frac{k+1}{k}}-(k+1)^{-\frac{1}{k}})\Delta t^{k+1}, 0]\frac{A^k}{(k+1)!}\oplus\Delta t\cdot\varepsilon``
- ``\textcolor{red}{bound}``: requires introducing RangeEnclosures as optional dependency, also `enclose` must be adapted for MIMO functions
- ``\textcolor{red}{poly}``: polynomial map, can be implemented by a series of quadratic maps (appendix A)
"""

# ╔═╡ a4dafbdd-475a-4791-a443-44349972afec


# ╔═╡ 344e4268-2915-446f-b37e-3ba37cb6dad5


# ╔═╡ 170ee9ee-d2a9-4d37-85d6-5fa5841dbaa8


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.0"
manifest_format = "2.0"
project_hash = "da39a3ee5e6b4b0d3255bfef95601890afd80709"

[deps]
"""

# ╔═╡ Cell order:
# ╠═1591bb74-2e89-11ed-07f8-33695ecb3632
# ╠═4129a82e-da07-40ae-aa57-0efbc807269a
# ╠═a4dafbdd-475a-4791-a443-44349972afec
# ╠═344e4268-2915-446f-b37e-3ba37cb6dad5
# ╠═170ee9ee-d2a9-4d37-85d6-5fa5841dbaa8
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
