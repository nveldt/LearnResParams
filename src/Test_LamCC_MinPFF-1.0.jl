using MAT
using Random

name = "KarateA"
#name = "Netscience"
mat = matread("../Graphs/"*name*".mat")
A = mat["A"]
n = size(A,1)

# Define a few parameters
GapTol = .01
ConTol = .01
filename = "output"
statusFrequency = 10
maxits = 1000
errortol = 1e-4

# Construct an arbitrary "example" clustering
c = ones(Int64,n)
c[1:10] .= 2

## SOLVER 1: SERIAL PROJECTION METHODS

# include("LamCC_MinPFF-1.0.jl")
# Solve the LambdaCC LP relaxation once
# gam = 10.0
# lam = .1
# @time X = Dykstra_lamCC_TFA(A,lam,GapTol,ConTol,gam,maxits,statusFrequency)

## SOLVER 2: PARALLEL PROJECTION METHODS
# include("LamCC_MinPFF_Parallel.jl")
# gam = 10.0

## S0LVER 3: Gurobi software
include("LamCC_PFF_Gurobi_1.0.jl")
gam = 0.0 # setting gam = 0.0 means Gurobi is called. If gam > 0, this defaults to serial projection methods.

## Standard LambdaCC
left = 0.01
right = 0.99
lamBest, BestP, ResList, Plist = MinPFF_LamCC(A,c,false, errortol, left, right, GapTol, ConTol, gam)
println("The Standard LamCC parameter fitness function equals $BestP when evaluated at lambda = $lamBest.")

## Degree-weighted LambdaCC
m = sum(nonzeros(A))
lam = 1/m
@time X = Dykstra_lamCC_dw(A,lam,GapTol,ConTol,gam,maxits,statusFrequency)

left = 1/(4*m)
right = 2/(m)
lamBest, BestP, ResList, Plist = MinPFF_LamCC(A,c,true, errortol, left, right, GapTol, ConTol, gam)
println("The degree-weighted LamCC parameter fitness function equals $BestP when evaluated at lambda = $lamBest.")
