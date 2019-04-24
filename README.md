# README

This repo contains code for algorithms and experiments from the paper

Learning Resolution Parameters for Graph Clustering

Nate Veldt, David Gleich, Anthony Wirth
Proceedings of he 29th International Conference on World Wide Web.

arXiv preprint: [https://arxiv.org/abs/1903.05246]()

The majority of the code is in the Julia progamming language (version 1.0). Some of the code is written in Matlab. 

## Outside software

Software from outside projects is included in the `include` folder.

* GenLouvain Algorithm:

    Lucas G. S. Jeub, Marya Bazzi, Inderjit S. Jutla, and Peter J. Mucha,    
    *"A generalized Louvain method for community detection implemented
    in MATLAB,"* [http://netwiki.amath.unc.edu/GenLouvain]() (2011-2016).
    
* For NMI scores: 

	[https://www.mathworks.com/matlabcentral/fileexchange/29047-normalized-mutual-information]()

* For ARI scores: 

	[https://www.mathworks.com/matlabcentral/fileexchange/49908-adjusted-rand-index]()
	

* LambdaCC algorithms: [https://github.com/nveldt/LamCC]()
* FlowSeed (local graph clustering): [https://github.com/nveldt/FlowSeed]()
* Metric Optimization (solving the LambdaCC LP relaxation using projection methods): [https://github.com/nveldt/MetricOptimization]()
* Parallel Metric Optimization (solving the LP relaxation in parallel): [https://github.com/nveldt/ParallelDykstras]()


Some routines require an academic license for Gurobi (e.g. the functions in GurobiCC-1.0.jl for solving the correlation clustering linear programming relaxation).

## Testing the algorithms

For the LambdaCC objective, there are several versions that can be run depending on a few choices:

* Two versions of LambdaCC: Standard and Degree-weighted 
* Three solvers for the underlying LP relaxation: Gurobi optimization software (more accurate solutions, but won't work for large problems), Serial Projection methods (approximate solutions, but more memory efficient than Gurobi), Parallel Projection methods (threaded version of serial projection methods).

An example for running the parameter fitness function in different cases is given in `Test_LamCC_MinPFF-1.0.jl`

To see how to minimize the local parameter fitness function in `Test_FlowSeed_PFF.jl`.

## Experiments

Code and instructions for reproducing each experiment from the paper are included in separate folders.
