# Minimize the Parameter Fitness Function for LambdaCC
# This code uses Gurobi for solving the LambdaCC LP relaxation.
# You can actually also use it to still call the projection based method.

include("../include/MetricOptimization/MetricOptimization-1.0.jl")
include("../include/MetricOptimization/GurobiCC-1.0.jl")

# Run the parameter fitness minimizer on the LamCC objective
function MinPFF_LamCC(A::SparseMatrixCSC{Float64,Int64},c::Vector{Int64},dweight::Bool,
    error::Float64 = 1e-6, left::Float64 = 0.0, right::Float64 = 1.0,GapTol::Float64=1e-3,ConTol::Float64=1e-3,gam::Float64=10.0)

    ResList = Vector{Float64}()
    Plist = Vector{Float64}()
    current = right
    lamBest,list1,list2 = LocalBranchOne(A,c,dweight,left,right,0.0,0.0,current,error,ResList,Plist,GapTol,ConTol,gam)

    BestP = EvaluatePFF(A,c,lamBest,dweight,GapTol,ConTol,gam)

    return lamBest, BestP, ResList, Plist

end

# Evaluate the denominator of the parameter fitness function by solving a relaxed LamCC problem
function EvaluatePFF(A::SparseMatrixCSC{Float64,Int64},c::Vector{Int64},lam::Float64,dweight::Bool,GapTol::Float64=1e-3,ConTol::Float64=1e-3,gam::Float64=10.0)

    n = size(A,1)
    if dweight
        D,Wc = LamCC_degweighted(A,lam)
    else
        D,Wc = LamCC_DandW(A,lam)
    end
    Xc = zeros(n,n)
    for i = 1:n-1
        for j = i+1:n
            if c[i] != c[j]
                Xc[j,i] = 1
            end
        end
    end
    Fc = abs.(Xc-D)


    Num = LPobj(Fc,Wc)

    # Compute the denominator using either Gurobi (if gamma == 0)
    # or using DysktraCC

    if gam == 0.0

        if dweight
            X,Denom = Degree_Weighed_LamCC_LP_Gurobi(A,lam)
        else
            X,Denom = LamCC_LP_Gurobi(A,lam)
        end

    else

        # Dykstra's method in fact returns a close approximation to the lower bound,
        # but it is not exact. To better ensure at least a good bound on the
        # parameter fitness function, we adjust the output on Dykstra's method
        # using the techinques shown in
        # "A Projection Method for Metric-Constrained Optimization"
        # https://arxiv.org/abs/1806.01678
        if dweight
            Output = Dykstra_lamCC_dw(A,lam,GapTol,ConTol,gam)
            DualLB = Output[4]
            Denom = DualLB[end]/(1+1/gam)
        else
            Output = Dykstra_lamCC_TFA(A,lam,GapTol,ConTol,gam)
            DualLB = Output[4]
            Denom = DualLB[end]/(1+1/gam)
        end

    end

    return Num/Denom
end

# Just go ahead and show the full bowl-shape
function FullBowl(A::SparseMatrixCSC{Float64,Int64},c::Vector{Int64},
    lambdas::Vector{Float64},dweight::Bool,GapTol::Float64=1e-3,ConTol::Float64=1e-3,gam::Float64=10.0)

    Poutput = Vector{Float64}()
    BestP = Inf
    BestRes = lambdas[1]

    for lam = lambdas
        Pout  = EvaluatePFF(A,c,lam,dweight,GapTol,ConTol,gam)
        push!(Poutput,Pout)

        if Pout < BestP
            BestP = Pout
            BestRes = lam
        end
    end

    return BestRes, BestP, Poutput
end



# We give it a left and right bound, and an error tolerance to get below
# X is the example clustering that we are "learning" from.
#
# Every time we evaluate the parameter fitness function at a new point,
# we update ResList with the new resolution parameter and the output for it
function LocalBranchOne(A::SparseMatrixCSC{Float64,Int64},c::Vector{Int64},dweight::Bool,
    left::Float64,right::Float64, Pleft::Float64, Pright::Float64, current::Float64,
    error::Float64,ResList::Vector{Float64},Plist::Vector{Float64},GapTol::Float64=1e-3,ConTol::Float64=1e-3,gam::Float64=10.0)

    # If we weren't passed them, evaluate at the endpoints
    if Pright == 0
        Pright = EvaluatePFF(A,c,right,dweight,GapTol,ConTol,gam)
        push!(ResList,right)
        push!(Plist,Pright)
    end

    if Pleft == 0
        Pleft = EvaluatePFF(A,c,left,dweight,GapTol,ConTol,gam)
        push!(ResList,left)
        push!(Plist,Pleft)
    end


    # Get the current best resolution parameter
    if current == 0
        current = left
        if Pright < Pleft
            current = right
        end
    end

    # Base case: return the best answer so far if we're within tolerance
    if abs(right-left) < error
        return current, ResList, Plist
    end
    len = right-left
#    println("The length of the interval is $len")

    # Midpoint:
    m = (left+right)/2
    Pm = EvaluatePFF(A,c,m,dweight,GapTol,ConTol,gam)
    push!(ResList,m)
    push!(Plist,Pm)

    # Just check cases to narrow down the window in which the minimizer is located
    if Pm == Pleft && Pleft == Pright && Pright < Inf
        current = m
        return current, ResList, Plist
    elseif Pleft <= Pm && Pm < Pright
        current = left
        return LocalBranchOne(A,c,dweight,left,m,Pleft,Pm,current,error,ResList,Plist,GapTol,ConTol,gam)
    elseif Pleft > Pm && Pm >= Pright
        current = right
        return LocalBranchOne(A,c,dweight,m,right,Pm,Pright,current,error,ResList,Plist,GapTol,ConTol,gam)
    elseif Pleft > Pm && Pm < Pright
        current = m
        return LocalBranchTwo(A,c,dweight,left,m,right,Pleft,Pm,Pright,current,error,ResList,Plist,GapTol,ConTol,gam)
    else

        @show left, right, m
        @show Pleft, Pright, Pm
        println("There is an issue with the curve not following the theoretical pattern.")
        println("A potential problem is that the underlying solver may not actually be solving the relaxation perfectly.")
        println("Need to decrease the constraint tolerance and primal-dual gap tolerance for the underlying solver.")
        return current, ResList, Plist
    end

end


# Sometimes we have to check two branches.
function LocalBranchTwo(A::SparseMatrixCSC{Float64,Int64},c::Vector{Int64},dweight::Bool,
    left::Float64, m::Float64,
    right::Float64, Pleft::Float64, Pm::Float64, Pright::Float64, current::Float64,
    error::Float64,ResList::Vector{Float64},Plist::Vector{Float64},GapTol::Float64=1e-3,ConTol::Float64=1e-3,gam::Float64=10.0)

    volA = sum(nonzeros(A))
    if abs(right-left) < error
        return current, ResList, Plist
    end
    len = right-left
#    println("The length of the interval is $len")

    # Left and right midpoints
    lmid = (left+m)/2
    Plmid = EvaluatePFF(A,c,lmid,dweight,GapTol,ConTol,gam)
    push!(ResList,lmid)
    push!(Plist,Plmid)
    rmid = (right+m)/2
    Prmid = EvaluatePFF(A,c,rmid,dweight,GapTol,ConTol,gam)
    push!(ResList,rmid)
    push!(Plist,Prmid)

    if Plmid == Pm && Pm == Prmid && Pm < Inf
        current = m
        return current, ResList, Plist
    elseif Plmid == Pm && Pm < Inf
        return LocalBranchOne(A,c,dweight,lmid,m,Plmid,Pm,current,error,ResList,Plist,GapTol,ConTol,gam)
    elseif Pm == Prmid && Pm < Inf
        return LocalBranchOne(A,c,dweight,m,rmid,Pm,Prmid,current,error,ResList,Plist,GapTol,ConTol,gam)
    elseif Plmid < Pm && Pm < Prmid
        current = lmid
        return LocalBranchTwo(A,c,dweight,left,lmid,m,Pleft,Plmid,Pm,current,error,ResList,Plist,GapTol,ConTol,gam)
    elseif Plmid > Pm && Pm > Prmid
        current = rmid
        return LocalBranchTwo(A,c,dweight,m,rmid,right,Pm,Prmid,Pright,current,error,ResList,Plist,GapTol,ConTol,gam)
    else
        current = m
        return LocalBranchTwo(A,c,dweight,lmid,m,rmid,Plmid,Pm,Prmid,current,error,ResList,Plist,GapTol,ConTol,gam)
    end
end
