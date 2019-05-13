using MAT
using Random

# If a clustering has values of 0, they should be treated as individual singleton
# clusters, rather than having one large cluster with label 0, because a label
# of zero means there is not data available.
function addSingletons(c)
    cur = maximum(c) + 1;

    n = length(c)
    for i = 1:n
        if c[i] == 0
            c[i] = cur
            cur += 1
        end
    end
    return c
end


include("../src/LamCC_MinPFF-1.0.jl")
mat = matread("Caltech36.mat")
A = mat["A"]
local_info = round.(Int64,mat["local_info"])
n = size(A,1)

# Cat1 = local_info(:,1); % student-faculty flag
# Cat2 = local_info(:,2); % gender
# Cat3 = local_info(:,3); % major
# Cat4 = local_info(:,4); % 2nd major
# Cat5 = local_info(:,5); % dorm/residence
# Cat6 = local_info(:,6); % Year
# Cat7 = local_info(:,7); % high school

mat = matread("Caltech36_LamCC_LPbounds.mat")

LPbounds = mat["LPbounds"]
Lambdas = mat["Lambdas"]
timers = mat["timers"]

mat = matread("Caltech36_LamCC_smallestLam_LPbound.mat")
smallbound = mat["LPbound"]
smalllam = mat["lam"]

Lambdas = [smalllam; Lambdas]
LPbounds = [smallbound; LPbounds]

Plists = zeros(21,7)
Plists_fake = zeros(21,7)

FakeCs = zeros(n,7)
for lamInd = 1:21

    lam = Lambdas[lamInd]
    bound = LPbounds[lamInd]
    D,Wc = LamCC_degweighted(A,lam)
    for attribute = 1:7
        c = addSingletons(local_info[:,attribute])
        cfake = shuffle(c)
        FakeCs[:,attribute] = cfake
        Xc = zeros(n,n)
        Xfake = zeros(n,n)
        for i = 1:n-1
            for j = i+1:n
                if c[i] != c[j]
                    Xc[j,i] = 1
                end
                if cfake[i] != cfake[j]
                    Xfake[j,i] = 1
                end
            end
        end
        Fc = abs.(Xc-D)
        Ffake = abs.(Xfake-D)
        Num = LPobj(Fc,Wc)
        NumFake = LPobj(Ffake,Wc)
        Plists[lamInd,attribute] = Num/bound
        Plists_fake[lamInd,attribute] = NumFake/bound
    end

end

LamBest = zeros(7)
minPs = zeros(7)

fake_LamBest = zeros(7)
fake_minPs = zeros(7)
# Then for each attribute, find the minimum of the parameter fitness function,
# both for the true metadata clustering and the fake clustering (permuted labels).

for attribute = 1:7

    # Take the best lambda found...

    # ...for the real metadata clustering...
    themin, wheremin = findmin(Plists[:,attribute])
    LamBest[attribute] = Lambdas[wheremin]
    minPs[attribute] = themin

    # ...and the permuted labels.
    themin, wheremin = findmin(Plists_fake[:,attribute])
    fake_LamBest[attribute] = Lambdas[wheremin]
    fake_minPs[attribute] = themin
end

matwrite("Caltech_Best_Lams.mat", Dict("LamBest" => LamBest, "minPs" => minPs,
"FakeCs" => FakeCs, "fake_LamBest" => fake_LamBest, "fake_minPs" => fake_minPs))
