## README

This is the readme file for Experiment 1 on LFR graphs from the paper "Learning Resolution Parameters for Graph Clustering"


## File Explanations

`generate_LFR.m`: Matlab interface for generating a single LFR graph

`Generate_data.m`: Generates all the data for the LFR experiments, for each mixing parameter we generate 1 for training and 5 for testing.

`Learn_Lambda.jl`: Learns Lambda on the training graph for each value of the mixing parameter. Calls Gurobi optimization software. 

`Run_GenLouvain.m`: Loads the lambda value learned from the training LFR graph, and runs GenLouvain with the learned lambda and the lambda corresponding to modularity. Stores output ARI and NMI scores. 

`ARI_Plot.m`: Plots the ARI scores obtained by running GenLouvain with the two different choices of lambda

`LFR_visualize.m`: Gives a way to visualize a network's community structure
