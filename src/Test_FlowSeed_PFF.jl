# Test the local parameter fitness function.

using MAT

include("../src/CutVol_MinPFF-1.0.jl")

datasets = ["DBLP","Amazon","LiveJournal","Orkut"]


# For a better understanding of the FlowSeed algorithm, see:
# https://arxiv.org/abs/1811.12280
#
# For this work, we are mostly just interested in using the FlowSeed routine
# because it is fast code for performing flow-based local graph clustering.

# In order to run on the other experiments, you need to download and process the
# the datasets. See the README in
# https://github.com/nveldt/FlowSeed/tree/master/Exp1-snap-communities
graph = "Amazon"
graph = "LiveJournal"
graph = "Orkut"

# DBLP is already included in the Experiment 2 folder
graph = "DBLP"

# Load the graph, and ten communities in it
mat = matread("../Exp2-snap-communities/graphs/"*graph*"-top10.mat")

# This is the adjacency matrix for DBLP
A = mat["A"]

d = sum(A,dims = 2)
volA = sum(nonzeros(A))
n = size(A,1)

# This stores the 10 largest communities, in an n x 10, 0-1 indicator matrix
comm = mat["C"]
numcom = size(comm,2)

# Conductance values for these communities
Conds = mat["Conds"]


# Select a community, and find the local resolution parameter that best fits it
# Then compare whether you get better results by minimizing conductance,
# or by minimizing an objective function with a special resolution parameter.

commID = 1

# Extract the target communities
Target = findall(x->x!=0,comm[:,commID])
TarSize = length(Target)
TargetStats = set_stats(A,vec(Target),volA)
Tcond = round(TargetStats[4],digits = 3)
Tsize = length(Target)
Tcut = TargetStats[1]
Tvol = TargetStats[2]
println("Community $commID has $Tsize nodes and a conductance of $Tcond")

## Start with 5 times the target set as a reference set
# We will search only for communities inside R:
#
#      min     cut(S) - alpha* vol(S)
#
#               s.t.  S \subseteq R

R = BFS_neighborhood(A,Target,5*TarSize)
volR = sum(d[R])

# Set a high locality parameter, so that we only search for sets that are
# subsets of R
epsilon = 10000.0
fR = volR/(volA-volR)
delta = epsilon - fR

# Report statistics for the seed set (also called the "reference set")
STATS = set_stats(A,R,volA)
condR = round(STATS[4],digits = 4)
pr,re,f1 = round.(PRF(Target,R),digits = 4)
numR = length(R)
println("Seed Set: |R| = $numR, condR = $condR, PR = $pr, RE = $re, F1 = $f1")

# Set up the penalties on excluding seed nodes.
# If we want to find the minimum conductance subset of R, then
# we should include no soft penalties  on excluding nodes from R (pR = 0),
# and also no strict penalties (i.e. RinS = 0)
RinS = zeros(numR,1)
pR = zeros(numR,1)

# Vanilla FlowSeed software returns the minimum conductance set
start = time()
S1, relcond = FlowSeed(A,R,epsilon,pR,RinS,true,true)
t1 = time() - start
cut, vol, edges, cond = set_stats(A,S1,volA)
c1 = round(cond,digits = 6)
size1 = length(S1)
prf = PRF(Target,S1)
pr1 = round(prf[1],digits = 6)
re1 = round(prf[2],digits = 6)
f11 = round(prf[3],digits = 6)
f21 = round(F2score(Target,S1),digits = 6)
println("Min Conductance: \tPR = $pr1 \t RE = $re1 \t F1 = $f11 \t F2 = $f21 \t Size = $size1 \t cond = $c1 \t Time = $t1")


## Now we learn a resolution parameter

## How well we approximately minimize the parameter fitness function
error = .05

## The minimum conductance in R defines the lower bound on meaningful
# resolution parameters to try
left = relcond

# To learn a better resolution parameter, we only need to know the cut and volume
# of the target set.
start = time()
alphaBest, BestP, Sbest, ResList, Plist = MinPFF_SimpleLocal(A,R,Tcut,Tvol,error,left,.99, 0.0, false)
t2 = time()-start
cut, vol, edges, condbest = set_stats(A,Sbest,volA)
c2 = round(condbest,digits = 6)
size2 = length(Sbest)
prf = PRF(Target,Sbest)
pr2 = round(prf[1],digits = 6)
re2 = round(prf[2],digits = 6)
f12 = round(prf[3],digits = 6)
f22 = round(F2score(Target,Sbest),digits = 6)
println("Learning alpha: \tPR = $pr2 \t RE = $re2 \t F1 = $f12 \t F2 = $f22 \t Size = $size2 \t cond = $c2 \t Time = $t2 \t alphaBest = $alphaBest")
