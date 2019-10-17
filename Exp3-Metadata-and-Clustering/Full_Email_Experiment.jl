
# Choose either the parallel or the serial version of the code.

# include("../src/LamCC_MinPFF_Parallel.jl")
include("../src/LamCC_MinPFF-1.0.jl")

# Load a graph
using MAT
mat = matread("../Graphs/Email_Graph.mat")
A = mat["A"]
c = vec(round.(Int64,mat["truth"]))

# Load the pre-saved fake metadata attribute, obtained by permuting the labels
# of the metadata clustering c
mat = matread("cfake_email.mat")
cfake = vec(round.(Int64,["cfake"]))

volA = sum(nonzeros(A))
n = size(A,1)

dweight = true
m = volA/2
error = 1e-8
right = 2/m
left = 1/(8*m)

tic()
lamBest, BestP, ResList, Plist = MinPFF_LamCC(A,c,dweight,error,left,right)
runtime = toc()

tic()
lamBestFake, BestPFake, ResListFake, PlistFake = MinPFF_LamCC(A,cfake,dweight,error,left,right)
faketime = toc()

matwrite("Full_Email_Experiment_Output.mat",
Dict("lamBestFake" => lamBestFake,"BestPFake" => BestPFake,
"ResListFake" => ResListFake, "PlistFake" => PlistFake, "cfake" => cfake,
"faketime" => faketime, "lamBest" => lamBest,"BestP" => BestP,
"ResList" => ResListFake, "Plist" => Plist, "c" => c,
"runtime" => runtime))
