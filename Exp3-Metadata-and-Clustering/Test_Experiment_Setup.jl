## Test the experiment setup on the simple 16 node graph


# Choose either the parallel or the serial version of the code.

# include("../src/LamCC_MinPFF_Parallel.jl")
include("../src/LamCC_MinPFF-1.0.jl")

# Load a graph
using MAT
mat = matread("../Graphs/ToyGraph16.mat")
A = mat["A"]
volA = sum(nonzeros(A))
n = size(A,1)

# The toy graph comes with several clusterings that optimize the standard
# LambdaCC objective in different parameter regimes
C = mat["C"]
C = round.(Int64,C)

# Extract a clustering and see if we can learn a good parameter for it
c = C[:,4]

# Permute the labels, so that the labeling shouldn't correspond to a
# good clustering
using Random
p = randperm(n)
cfake = c[p]

# Minimize the degree-weighted LambdaCC parameter fitness function
dweight = true
m = volA/2
error = 1e-8
right = 2/m
left = 1/(8*m)

start = time()
lamBest, BestP, ResList, Plist = MinPFF_LamCC(A,vec(round.(Int64,c)),dweight,error,left,right)
runtime = time()-start

start = time()
lamBestFake, BestPFake, ResListFake, PlistFake = MinPFF_LamCC(A,vec(round.(Int64,cfake)),dweight,error,left,right)
faketime = time()-start

println("Min P for true clustering is $BestP\n Min P for fake clustering is $BestPFake")

matwrite("Test_Experiment_Setup_Output.mat",
Dict("lamBestFake" => lamBestFake,"BestPFake" => BestPFake,
"ResListFake" => ResListFake, "PlistFake" => PlistFake, "cfake" => cfake,
"faketime" => faketime, "lamBest" => lamBest,"BestP" => BestP,
"ResList" => ResListFake, "Plist" => Plist, "c" => c,
"runtime" => runtime))
