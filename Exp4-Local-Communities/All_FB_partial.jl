# /p/mnt/software/julia-0.6.0/bin/julia

include("../include/FlowSeed-1.0.jl")
include("../src/FlowSeed_MinPFF-1.0.jl")
using MAT
using MatrixNetworks
# This file loads a list of Facebook 100 network names. From that list,
# this function will run the experimental setup for networks indexed by "start"
# to "finish".
function AllFB(start::Int64,last::Int64)

f = open("Facebook_Sets.txt")

lines = readlines(f)

# Get a list of all the Facebook100 network names
FBlist = split(lines[1])

# Extract a subset of these to run experiments on
FBlist = FBlist[start:last]

# Look only at communities with at least MinSize nodes
MinSize = 20

# Using a breadth-first-search, grow a community that is "expandby" times
# the size of the original network (or at most half the graph).
expandby = 3

path_to_FB100 = homedir()*"/data/Facebook100/"
println("You will need to obtain and save the Facebook 100 networks to run these experiments.")

     for name =  FBlist

         mat = matread(path_to_FB100*name*".mat")
         A = mat["A"]

         # Metadata
         local_info = mat["local_info"]


         # Get largest connected component, to avoid 0-conductance sets.
         Acc,p = largest_component(A)
         A = Acc
         local_info = local_info[p,:]

         volA = sum(nonzeros(A))
         n = size(A,1)
         d = sum(A,dims=2)

         outputstring = "all_output/FB_metacom_"*name*"_"*string(MinSize)*"_"*string(expandby)*".txt"
         outputmat = "all_output/FB_metacom_"*name*"_"*string(MinSize)*"_"*string(expandby)*".mat"
         open(outputstring,"w") do f
              write(f,"\n Experiments for Graph: "*name*"\n")
         end

         # State the minimum community size we care about

         outputdata = zeros(1,10)

         # For each meta-data category...
         for cat = 1:6
             println("")

             open(outputstring,"a") do f
                  write(f,"\nMetadata attribute number $cat\n")
             end

             # Get all the unique labels.
             labels = unique(local_info[:,cat])

             # ignore the label 0, which means "no information"
             labels = setdiff(labels,0)

             for j = 1:length(labels)
                 l = labels[j]
                 Group = findall(x -> x == l, local_info[:,cat])
                 gs = length(Group)

                 # The group can't be too large or too small
                 if gs > MinSize && gs < n/2
                     cutG,volG,edgesG,condG = set_stats(A,Group,volA);
                     println("Group $j in category $cat has size $gs and cond = $condG")
                     open(outputstring,"a") do f
                          write(f,"\nGroup $j in category $cat has size $gs and cond = $condG\n")
                     end

                     # Now we grow this an learn the best resolution parameter alpha
                     Target = Group

                     # Get a superset of this group
                     R = BFS_neighborhood(A,Target,min(length(Target)*expandby,round(Int,n/2)))

                       X = Target
                       numR = length(R)
                       volR = sum(d[R])
                       fR = volR/(volA - volR)
                       epsilon = 10000.0
                       STATS = set_stats(A,R,volA)
                       condR = round(STATS[4],digits=4)
                       pr,re,f1 = round.(PRF(Target,R),digits=4)
                       println("Superset: |R| = $numR, condR = $condR, PR = $pr, RE = $re, F1 = $f1")
                       open(outputstring,"a") do f
                            write(f, "Superset: |R| = $numR, condR = $condR, PR = $pr, RE = $re, F1 = $f1\n")
                       end
                       # First Test: Run SimpleLocal
                       pR = zeros(numR,1)
                       RinS = zeros(numR,1)

                       start = time()
                       S1, cond1 = FlowSeed(A,R,epsilon,pR,RinS,d,true)
                       t1 = time()-start
                       c1 = round(cond1,digits=4)
                       size1 = length(S1)
                       prf = PRF(Target,S1)
                       pr1 = round(prf[1],digits=3)
                       re1 = round(prf[2],digits=3)
                       f1_mincond = round(prf[3],digits=3)
                       println("Ratio Objective: \tPR = $pr1 \t RE = $re1 \t F1 = $f1_mincond \t Size = $size1 \t Cond = $c1 \t Time = $t1")
                       open(outputstring,"a") do f
                            write(f, "Min Cond Subset: \tPR = $pr1 \t RE = $re1 \t F1 = $f1_mincond \t Size = $size1 \t Cond = $c1 \t Time = $t1 \n")
                       end
                       alphaStart = cond1


                       ## Second Approach: Find best Alpha
                       start = time()
                       alphaBest,BestP,Sbest,ResList,Plist = MinPFF_SimpleLocal(A,R,X,0.05,alphaStart,1.0,epsilon,false)
                       tLalph = time()-start

                       size1 = length(Sbest)
                       prf = PRF(Target,Sbest)
                       pr1 = round(prf[1],digits=3)
                       re1 = round(prf[2],digits=3)
                       f1_minP = round(prf[3],digits=3)
                       println("Sbest: \tPR = $pr1 \t RE = $re1 \t F1 = $f1_minP \t Size = $size1 \t Alpha = $alphaBest, \t Pbest = $BestP")
                       open(outputstring,"a") do f
                            write(f, "Learning Alpha: \tPR = $pr1 \t RE = $re1 \t F1 = $f1_minP \t Size = $size1 \t Alpha = $alphaBest \t BestP = $BestP \t time = $tLalph \n")
                       end
                       outputdata = [outputdata; cat l gs condG alphaBest BestP pr1 re1 f1_minP f1_mincond]
                 end

             end

         end

         readme = "outputdata = [category, label of group, group size, conductance of group, learned alphaBest, min of parameter fitness, precision, recall, f1_minP,f1_mincond]"
         matwrite(outputmat, Dict("outputdata" => outputdata, "readme" => readme))

     end

end

AllFB(1,100)
