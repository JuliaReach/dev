Idea: Create a `JuliaReachBase.jl` package that contains common "utility" functions used in different packages, but don't naturally fit into one of them.

Methods:

-  Activate / Deactivate assertions (LazySets)
-  Macro to test if package is loaded (LazySets and ReachabilityAnalysis)
