# Minimize the Parameter Fitness Function for Local Clustering,
# specifically in the situation where you know an estimate for cut(X) and
# vol(X), the cut and volume of the target set, but X is not known exactly.
# Knowing cut(X) and vol(X) is enough to minimize the parameter fitness function
# corresponding to the following graph clustering objective:
#
#           cut(S) - alpha* vol(S)
#
#               s.t.  S \subseteq R
#
#
# This code relies heavily on other work on localized community detection.
include("../include/FlowSeed-1.0.jl")


# Run the parameter fitness minimizer with FlowImprove
function MinPFF_SimpleLocal(A::SparseMatrixCSC{Float64,Int64},R::Vector{Int64},cutX::Float64,volX::Float64,
    error::Float64 = 0.01,left::Float64 = 0.0, right::Float64 = 1.0, delta::Float64 = 0.0, getLB::Bool = true, localFlag::Bool = true)

    # Ignore the input delta--if you are using this code, then you are
    # looking for the minimum conductance subset of the seed set

    delta = 1e8
    volA = sum(nonzeros(A))
    n = size(A,1)
    d = sum(A,dims=2)
    numR = length(R)
    volR = sum(d[R])
    fR = volR/(volA - volR)
    epsilon = fR + delta
    pR = zeros(numR,1)
    RinS = pR

    # Get the lower bound
    if getLB
        Sstar, alphaStar = FlowSeed(A,R,epsilon,pR,RinS,d,localFlag)
        # a lower bound we need is alphaStar
        left = alphaStar
    end

    # Find one-hop neighbors of R, and get the complement set
    Rn = neighborhood(A,R,1)    # get the immediate neighbors of R...
    Rn = setdiff(Rn,R)          # ...but we exclude R itself
    inRc = ones(n)
    inRc[R] .= 0
    Rc = findall(x->x!=0,inRc)            # complement of R

    ResList = Vector{Float64}()
    Plist = Vector{Float64}()
    current = right
    alphaBest,list1,list2 = LocalBranchOne(A,R,Rn,Rc,epsilon,pR,RinS,d,true,cutX,volX,left,right,0.0,0.0,current,error,ResList,Plist)

    BestP, Sbest = EvaluatePFF(A,cutX,volX,R,Rn,Rc,epsilon,pR,RinS,d,true,alphaBest)
    return alphaBest, BestP, Sbest, ResList, Plist

end

# Evaluate the denominator of the parameter fitness function by solving a flow problem
function EvaluatePFF(A::SparseMatrixCSC{Float64,Int64},cutX::Float64,volX::Float64,R::Vector{Int64},
    Rn::Vector{Int64},Rc::Vector{Int64}, epsilon::Float64,pR::Array{Float64},RinS::Array{Float64},
    d::Array{Float64},localFlag,alpha::Float64)

    volA = sum(nonzeros(A))
    n = size(A,1)
    source = zeros(n)
    sink = zeros(n)
    dr = d[R]
    drc = d[Rc]
    source[R] = alpha*(pR .+ 1).*dr
    sink[Rc] = alpha*epsilon*drc
    volR = sum(d[R])

    # The numerator  of the parameter fitness function is a linear function in alpha
    Num = cutX - alpha*volX + alpha*volR

    # Compute the min s-t cut to get the denominator
    if localFlag
        # Do it by repeatedly solving smaller problems, starting
        # by looking at the immediate neighbors Rn
        S = LocalPushRelabel(A,R,source,sink,Rn)
    else
        # Run a single min-cut computation on the whole graph
        S = NonLocalPushRelabel(A,R,source,sink)
    end

    if length(S) == 0 || length(S) == n
        return Inf, S
    else
        Denom = STcutAlpha(A,S,R,d,alpha,epsilon,volA,pR,RinS)
        ptemp = Num/Denom
        println("New evaluation: alpha = $alpha, P(alpha) = $ptemp")
        return ptemp, S
    end

end

# Just go ahead and show the full bowl-shape

function FullBowl(A::SparseMatrixCSC{Float64,Int64},R::Vector{Int64},
    epsilon::Float64,pR::Array{Float64},RinS::Array{Float64},
    d::Array{Float64},localFlag::Bool,cutX::Float64,volX::Float64,alphas::Vector{Float64})

    n = size(A,1)

    # Find one-hop neighbors of R, and get the complement set
    Rn = neighborhood(A,R,1)    # get the immediate neighbors of R...
    Rn = setdiff(Rn,R)          # ...but we exclude R itself
    inRc = ones(n)
    inRc[R] .= 0
    Rc = findall(x->x!=0,inRc)             # complement of R

    Poutput = Vector{Float64}()
    BestP = Inf
    BestRes = alphas[1]

    for alpha = alphas

        Pout, S = EvaluatePFF(A,cutX,volX,R,Rn,Rc,epsilon,pR,RinS,d,localFlag,alpha)
        push!(Poutput,Pout)

        if Pout < BestP
            BestP = Pout
            BestRes = alpha
        end
    end

    return BestRes, BestP, Poutput
end


# Starts with the parameters of FlowSeed, because we'll need to call that.
#
# Then we give it a left and right bound, and an error tolerance to get below
# X is the example clustering that we are "learning" from.
#
# Every time we evaluate the parameter fitness function at a new point,
# we update ResList with the new resolution parameter and the output for it
function LocalBranchOne(A::SparseMatrixCSC{Float64,Int64},R::Vector{Int64},
    Rn::Vector{Int64},Rc::Vector{Int64}, epsilon::Float64,pR::Array{Float64},RinS::Array{Float64},
    d::Array{Float64},localFlag::Bool, cutX::Float64,volX::Float64, left::Float64,
    right::Float64, Pleft::Float64, Pright::Float64, current::Float64,
    error::Float64,ResList::Vector{Float64},Plist::Vector{Float64})


    volA = sum(nonzeros(A))

    # If we weren't passed them, evaluate at the endpoints
    if Pleft == 0
        Pleft, S = EvaluatePFF(A,cutX,volX,R,Rn,Rc,epsilon,pR,RinS,d,localFlag,left)
        push!(ResList,left)
        push!(Plist,Pleft)
    end
    if Pright == 0
        Pright,S = EvaluatePFF(A,cutX,volX,R,Rn,Rc,epsilon,pR,RinS,d,localFlag,right)
        push!(ResList,right)
        push!(Plist,Pright)
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
    Pm,S = EvaluatePFF(A,cutX,volX,R,Rn,Rc,epsilon,pR,RinS,d,localFlag,m)
    push!(ResList,m)
    push!(Plist,Pm)

    # Just check cases to narrow down the window in which the minimizer is located
    if Pm == Pleft && Pleft == Pright && Pright < Inf
        current = m
        return current, ResList, Plist
    elseif Pleft <= Pm && Pm < Pright
        current = left
        return LocalBranchOne(A,R,Rn,Rc,epsilon,pR,RinS,d,localFlag,cutX,volX,left,m,Pleft,Pm,current,error,ResList,Plist)
    elseif Pleft > Pm && Pm >= Pright
        current = right
        return LocalBranchOne(A,R,Rn,Rc,epsilon,pR,RinS,d,localFlag,cutX,volX,m,right,Pm,Pright,current,error,ResList,Plist)
    elseif Pleft > Pm && Pm < Pright
        current = m
        return LocalBranchTwo(A,R,Rn,Rc,epsilon,pR,RinS,d,localFlag,cutX,volX,left,m,right,Pleft,Pm,Pright,current,error,ResList,Plist)
    end

end



# Sometimes we have to check two branches.
function LocalBranchTwo(A::SparseMatrixCSC{Float64,Int64},R::Vector{Int64},
    Rn::Vector{Int64},Rc::Vector{Int64}, epsilon::Float64,pR::Array{Float64},RinS::Array{Float64},
    d::Array{Float64},localFlag::Bool, cutX::Float64, volX::Float64, left::Float64, m::Float64,
    right::Float64, Pleft::Float64, Pm::Float64, Pright::Float64, current::Float64,
    error::Float64,ResList::Vector{Float64},Plist::Vector{Float64})

    volA = sum(nonzeros(A))
    if abs(right-left) < error
        return current, ResList, Plist
    end
    len = right-left
#    println("The length of the interval is $len")

    # Left and right midpoints
    lmid = (left+m)/2
    Plmid,S = EvaluatePFF(A,cutX,volX,R,Rn,Rc,epsilon,pR,RinS,d,localFlag,lmid)
    push!(ResList,lmid)
    push!(Plist,Plmid)
    rmid = (right+m)/2
    Prmid,S = EvaluatePFF(A,cutX,volX,R,Rn,Rc,epsilon,pR,RinS,d,localFlag,rmid)
    push!(ResList,rmid)
    push!(Plist,Prmid)

    if Plmid == Pm && Pm == Prmid && Pm < Inf
        current = m
        return current, ResList, Plist
    elseif Plmid == Pm && Pm < Inf
        return LocalBranchOne(A,R,Rn,Rc,epsilon,pR,RinS,d,localFlag,cutX,volX,lmid,m,Plmid,Pm,current,error,ResList,Plist)
    elseif Pm == Prmid && Pm < Inf
        return LocalBranchOne(A,R,Rn,Rc,epsilon,pR,RinS,d,localFlag,cutX,volX,m,rmid,Pm,Prmid,current,error,ResList,Plist)
    elseif Plmid < Pm && Pm < Prmid
        current = lmid
        return LocalBranchTwo(A,R,Rn,Rc,epsilon,pR,RinS,d,localFlag,cutX,volX,left,lmid,m,Pleft,Plmid,Pm,current,error,ResList,Plist)
    elseif Plmid > Pm && Pm > Prmid
        current = rmid
        return LocalBranchTwo(A,R,Rn,Rc,epsilon,pR,RinS,d,localFlag,cutX,volX,m,rmid,right,Pm,Prmid,Pright,current,error,ResList,Plist)
    else
        current = m
        return LocalBranchTwo(A,R,Rn,Rc,epsilon,pR,RinS,d,localFlag,cutX,volX,lmid,m,rmid,Plmid,Pm,Prmid,current,error,ResList,Plist)
    end
end
