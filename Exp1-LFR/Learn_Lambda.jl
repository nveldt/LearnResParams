using MAT
include("../src/LamCC_PFF_Gurobi_1.0.jl")


for mu = 20:5:50
    outputmat = "Output/Learned_SmallComm_mu_"*string(mu)*"_Gurobi.mat"
    mat1 = matread("data/Train_1_SmallComm__mu_"*string(mu)*".mat")

    A1 = mat1["A"]
    c1 = vec(round.(Int64,mat1["truth"]))

    # These parameters don't do anything if we are going to call Gurobi
    error = 0.0
    GapTol = 0.0
    ConTol = 0.0

    # Setting gam = 0.0 means the Gurobi solver is used
    gam = 0.0

    m = sum(nonzeros(A1))/2
    left = 1/(8*m)
    right = 2/(m)
    s1 = time()
    lamBest1, BestP1, ResList1, Plist1 = MinPFF_LamCC(A1,c1,true,error,left,right,GapTol,ConTol,gam)
    time1 = time()-s1

    matwrite(outputmat, Dict(
    "lamBest1" => lamBest1, "BestP1" => BestP1,"ResList1" => ResList1, "Plist1" => Plist1,"time1"=> time1))

end
