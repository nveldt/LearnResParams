# On my laptop, run /Applications/Julia-1.0.3.app/Contents/Resources/julia/bin/julia

include("MetricOptimization-1.0.jl")

using MAT

name = "KarateA"
mat = matread(homedir()*"/data/Graphs/"*name*".mat")
A = mat["A"]

# Define parameters
GapTol = .01
ConTol = .01
filename = "testoutput"
statusFrequency = 10
gam = 10.0
maxits = 1000
lam = .2

# Testing standard lambdaCC

Output = Dykstra_lamCC_TFA(A,lam,GapTol,ConTol,gam,maxits,statusFrequency)
Xd = Output[1]
bound_d = Output[9]

include("GurobiCC-1.0.jl")

X,bound = LamCC_LP_Gurobi(A,lam,false)

@show bound, bound_d, FullTriangleCheck(X), FullTriangleCheck(Xd)



# Testing Degree-weighted lambdaCC

lam = 1/sum(nonzeros(A))

Output = Dykstra_lamCC_dw(A,lam,GapTol,ConTol,gam,maxits,statusFrequency)
Xd = Output[1]
bound_d = Output[9]

X,bound = Degree_Weighed_LamCC_LP_Gurobi(A,lam,false)

@show bound, bound_d, FullTriangleCheck(X), FullTriangleCheck(Xd)
