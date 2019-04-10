using MAT
using Random


mat = matread("ToyGraph16.mat")
A = mat["A"]
n = size(A,1)


mat = matread("UniqueOptimalClusterings_Toy16.mat")
C = mat["C"]
C = round.(Int64,C)

c = C[:,4]

# Define a few parameters
GapTol = .01
ConTol = .01
filename = "output"
statusFrequency = 10
gam = 10.0
maxits = 5000
lam = .1

include("LamCC_MinPFF_Parallel.jl")


#@time X = Parallel_Dykstra_lamCC_TFA(A,lam,GapTol,ConTol,gam,maxits,statusFrequency,filename)

Start = time()
@time lamBest, BestP, ResList, Plist = MinPFF_LamCC(A,c,false,0.01,.03,.9)
partime = time()-Start

# Serial version
include("LamCC_MinPFF-1.0.jl")

@time lamBest, BestP, ResList, Plist = MinPFF_LamCC(A,c,false,0.01,.03,.9)

@show partime
