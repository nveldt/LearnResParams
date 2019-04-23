# README

This repo contains code for algorithms and experiments from the paper

Learning Resolution Parameters for Graph Clustering

Nate Veldt, David Gleich, Anthony Wirth
Proceedings of the 29th International Conference on World Wide Web.

arXiv preprint: [https://arxiv.org/abs/1903.05246]()

The majority of the code is in the Julia progamming language (version 1.0). Some of the code is written in Matlab. 

## Outside software

The Algorithms folder constains algorithms and code from several different sources

GenLouvain Algorithm:

    Lucas G. S. Jeub, Marya Bazzi, Inderjit S. Jutla, and Peter J. Mucha,    
    *"A generalized Louvain method for community detection implemented
    in MATLAB,"* http://netwiki.amath.unc.edu/GenLouvain (2011-2016).
    
For NMI scores: 

	https://www.mathworks.com/matlabcentral/fileexchange/29047-normalized-mutual-information

For ARI scores: 

	https://www.mathworks.com/matlabcentral/fileexchange/49908-adjusted-rand-index
	
The algorithms folder also contains code from previous projects:

* LambdaCC algorithms: [https://github.com/nveldt/LamCC]()
* FlowSeed (local graph clustering): [https://github.com/nveldt/FlowSeed]()
* Metric Optimization (solving the LambdaCC LP relaxation using projection methods): [https://github.com/nveldt/MetricOptimization]()
* Parallel Metric Optimization (solving the LP relaxation in parallel): [https://github.com/nveldt/ParallelDykstras]()


Some of the function require an academic license for Gurobi (e.g. the functions in GurobiCC-1.0.jl for solving the correlation clustering linear programming relaxation).

## Testing the algorithms

For the LambdaCC objective, several choices have to be made regarding which algorithm is used:

* Two versions of LambdaCC: Standard and Degree-weighted 
* Three solvers for the underlying LP relaxation: Gurobi optimization software (more accurate solutions, but won't work for large problems), Serial Projection methods (approximate solutions, but more memory efficient than Gurobi), Parallel Projection methods (threaded version of serial projection methods).

An example for running the parameter fitness function in different cases is given in `Test_LamCC_MinPFF-1.0.jl`

