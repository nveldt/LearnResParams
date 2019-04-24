# Test the FlowSeed local graph clustering implementation

using MAT

include("../src/FlowSeed_MinPFF-1.0.jl")

datasets = ["DBLP","Amazon","LiveJournal","Orkut"]

# In order to run on the other experiments, you need to download and process the
# the datasets. See the README in
# https://github.com/nveldt/FlowSeed/tree/master/Exp1-snap-communities

# DBLP is included in the repo for testing
datasets = ["DBLP"]

function runall()

for graph = datasets

    outputstring = "Output/"*graph*"_5X.txt"
    outputmatrix = "Output/"*graph*"_5X.mat"

    open(outputstring,"w") do f
        write(f,"Assuming we known cut(X) and vol(X), and a superset of X of
        size 5|X|, compare how well we detect X by finding the mininmum conductance
        subset, and by learning alpha. \n")
    end

    mat = matread("graphs/"*graph*"-top10.mat")
    A = mat["A"]
    comm = mat["C"]
    Conds = mat["Conds"]
    d = sum(A,dims = 2)
    volA = sum(nonzeros(A))
    n = size(A,1)
    numcom = size(comm,2)

    mincond_sets = spzeros(n,numcom)
    mincond_stats = zeros(7,numcom)
    learnalpha_sets = spzeros(n,numcom)
    learnalpha_stats = zeros(9,numcom)

    for commID = 1:numcom

        # Extract one of the target communities
        Target = findall(x->x!=0,comm[:,commID])
        TarSize = length(Target)
        TargetStats = set_stats(A,vec(Target),volA)
        Tcond = round(TargetStats[4],digits = 3)
        Tsize = length(Target)
        println("Community $commID has $Tsize nodes and a conductance of $Tcond")

        ## Start with 5 times the target set as a reference set
        R = BFS_neighborhood(A,Target,5*TarSize)
        volR = sum(d[R])
        fR = volR/(volA-volR)
        epsilon = 10000.0
        delta = epsilon - fR
        STATS = set_stats(A,R,volA)
        condR = round(STATS[4],digits = 4)
        pr,re,f1 = round.(PRF(Target,R),digits = 4)
        numR = length(R)
        println("Seed Set: |R| = $numR, condR = $condR, PR = $pr, RE = $re, F1 = $f1")
        open(outputstring,"a") do f
              write(f,"\nCommunity $commID has $Tsize nodes and a conductance of $Tcond \n")
              write(f, "Seed Set: |R| = $numR, condR = $condR, PR = $pr, RE = $re, F1 = $f1 \n")
        end

        # Set up the penalties on excluding seed nodes, which is zero
        # in this case if we want the minimum conductance subset
        RinS = zeros(numR,1)
        pR = zeros(numR,1)

        # Run FlowSeed software to find the minimum conductance set
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

        mincond_sets[S1,commID] .= 1
        mincond_stats[:,commID] = [t1; length(S1); c1; pr1; re1; f11; f21]

        ## Now try to learn the very best alpha to use:
        error = .05
        left = relcond

        start = time()
        alphaBest, BestP, Sbest, ResList, Plist = MinPFF_SimpleLocal(A,R,Target,error,left, .99, delta, false, true)
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
        learnalpha_sets[Sbest,commID] .= 1
        learnalpha_stats[:,commID] = [t2; length(Sbest); c2; pr2; re2; f12; f22; alphaBest; BestP]

        open(outputstring,"a") do f
          write(f, "Minimum conduc: \tPR = $pr1 \t RE = $re1 \t F1 = $f11 \t F2 = $f21 \t Size = $size1 \t cond = $c1 \t Time = $t1 \n")
          write(f, "Learning alpha: \tPR = $pr2 \t RE = $re2 \t F1 = $f12 \t F2 = $f22 \t Size = $size2 \t Cond = $c2 \t Time = $t2 \t alphaBest = $alphaBest \t BestP = $BestP\n")
        end

    end

    matwrite(outputmatrix, Dict( "learnalpha_stats" => learnalpha_stats,
    "learnalpha_sets" => learnalpha_sets, "mincond_sets" => mincond_sets,
    "mincond_stats"=> mincond_stats))

end

end

runall()
