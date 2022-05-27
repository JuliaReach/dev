using LazySets, Plots

B = rand(BallInf)
dump(B)

X = rand(Hyperrectangle)
plot(X)

dump(X)

Z = convert(Zonotope, X)
dump(Z)

P = convert(HPolytope, X)
dump(P)

constraints_list(X)
vertices_list(X)

# linear map: lazy
A = rand(2, 2)
A * X

# linear map: concrete
linear_map(A, X)

# cartesian product: lazy
I = Interval(0, 1)
X × I

# cartesian product: concrete
cartesian_product(X, I)

# they can compose:
T = BallInf(zeros(3), 100.0)
E = ((A * X) × I) ∩ T

typeof(E)

concretize(E)

plot(project(concretize(E), 1:2))

# for 3d: it is possible to use Makie.jl ....
