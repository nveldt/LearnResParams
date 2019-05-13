include("../src/LamCC_MinPFF_Parallel.jl")
using MAT

println("Note that graph must be obtained separately from the Facebook 100 dataset")
mat = matread("../Graphs/Caltech36.mat")
A = mat["A"]
local_info = mat["local_info"]
n = size(A,1)
m = sum(nonzeros(A))/2

# We simply evaluate the denominator of the parameter fitness function at a range of values
# rather than running our algorithm over and over again for 7 attributes.

left = 1/(8*m)
right = 2/m

Lambdas = collect(range(left,stop=right,length=20))
LPbounds = zeros(length(Lambdas))
next = 1
timers = zeros(length(Lambdas))
next = 1
gam = 10.0
GapTol = 0.01
ConTol = 0.01
maxits = 10000
statusFrequency = 100
for lam = Lambdas
    global next
    start = time()
    Output = Parallel_Dykstra_lamCC_dw(A,lam,GapTol,ConTol,gam,maxits,statusFrequency)
    timers[next] = time()-start
    next+=1

    # To get a lower bound on the LP, we have to adjust a little based on
    # the value of gamma, since this is only approximately solving the LP
    #
    # For details on why this is a lower bound, see
    # https://arxiv.org/abs/1806.01678
    DualLB = Output[4]
    LPbounds[next] = DualLB[end]/(1+1/gam)
    next += 1
end

matwrite("Caltech36_LamCC_LPbounds.mat", Dict("LPbounds" => LPbounds, "Lambdas" => Lambdas,"timers" => timers))


# Get a bound for one extra small value of lambda as well
increment = Lambdas[5]-Lambdas[4]
lam = left-increment

start = time()
Output = Parallel_Dykstra_lamCC_dw(A,lam,GapTol,ConTol,gam,maxits,statusFrequency)
timer = time()-start
DualLB = Output[4]
LPbound = DualLB[end]/(1+1/gam)

matwrite("Caltech36_LamCC_smallestLam_LPbound.mat", Dict("LPbound" => LPbound, "lam" => lam,"timer" => timer))
