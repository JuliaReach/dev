using ReachabilityAnalysis
using Plots

# Create a Hyperrectangle centered at 5 and of radius 1
H = Hyperrectangle([5.0, 5.0], [1.0, 1.0])

# Make it a reach-set the assigning time interval 0 .. 1
R = ReachSet(H, 0 .. 1)

# Convert to Taylor model representation
T = overapproximate(R, TaylorModelReachSet)

# f1(x, y) = 5.0 + 1.0 x₁ + [0, 0]
# f2(x, y) = 5.0 + 1.0 x₂ + [0, 0] 

# check test breaks
# Z = overapproximate(T, Zonotope, dom=IntervalBox(-1 .. 1)) # dimensions should match

Z = overapproximate(T, Zonotope, dom=IntervalBox(-1 .. 1, 2)) |> set
isequivalent(Z, H)

Z1 = overapproximate(T, Zonotope, dom=IntervalBox(0 .. 0.5, 0 .. 0.5)) |> set
# expected to be 1/4th but is not
H1 = overapproximate(Z, Hyperrectangle)

H3 = Hyperrectangle([5.5, 5.5], [0.5, 0.5])
isequivalent(H1, H3)

# plot incorrect evaluation
# includes zero, but evaluation is wrong
plot(H)
plot!(Z1)
plot!(H3, ls=:dash)

# doesn't include zero, can't evaluate
overapproximate(T, Zonotope, dom=IntervalBox(0.9 .. 1.0, 0.9 .. 1.0))

_overapproximate(T, Zonotope, dom=IntervalBox(0 .. 0.5, 0 .. 0.5)) |> set

using ReachabilityAnalysis: zeroBox, TaylorModelN, zeroI, fp_rpa

R = T
Δt = tspan(R)
# dimension of the reachset
n = dim(R)

# evaluate the Taylor model in time
# X_Δt is a vector of TaylorN (spatial variables) whose coefficients are intervals
X = set(R)
tdom = Δt - tstart(R)  # normalize time (to TM-internal time)
tdom = tdom ∩ domain(R)  # intersection handles round-off errors
X_Δt = evaluate(X, tdom)



# Idea (Luis): recenter the polynomials to the new domain
function _overapproximate(R::TaylorModelReachSet{N}, ::Type{<:Zonotope};
                          dom::IntervalBox=symBox(dim(R))) where {N}
    Δt = tspan(R)
    # dimension of the reachset
    n = dim(R)

    # evaluate the Taylor model in time
    # X_Δt is a vector of TaylorN (spatial variables) whose coefficients are intervals
    X = set(R)
    tdom = Δt - tstart(R)  # normalize time (to TM-internal time)
    tdom = tdom ∩ domain(R)  # intersection handles round-off errors
    X_Δt = evaluate(X, tdom)

    # builds the associated taylor model for each coordinate j = 1...n
    #  X̂ is a TaylorModelN whose coefficients are intervals
    X̂ = [TaylorModelN(X_Δt[j], zeroI, zeroBox(n), dom) for j in 1:n]

    # compute floating point rigorous polynomial approximation
    # fX̂ is a TaylorModelN whose coefficients are floats
    fX̂ = fp_rpa.(X̂)

    # LazySets can overapproximate a Taylor model with a Zonotope
    Zi = overapproximate(fX̂, Zonotope)
    return ReachSet(Zi, Δt)
end