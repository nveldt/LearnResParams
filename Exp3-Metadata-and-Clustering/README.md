## README

In this set of experiments we measure the relationship between community structure and metadata attributes in two real-world networks: Email, and Caltech36.

### Datasets

1. Email [1][2]: dataset downloaded from the SNAP large networks repository: [https://snap.stanford.edu/data/email-Eu-core.html](). The original graph is directed. We make edges undirected and extract the largest connected component. Metadata corresponds to departmental affilition of researchers (nodes) at a European university. We include the undirected largest component and corresponding metadata labels in the data folder.
2. Caltech36: subset of the Facebook network at the California Institute of Technology, from the Facebook 100 dataset [3].


[1] Hao Yin, Austin R. Benson, Jure Leskovec, and David F. Gleich. "Local Higher-order Graph Clustering." In Proceedings of the 23rd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining. 2017.

[2] J. Leskovec, J. Kleinberg and C. Faloutsos. Graph Evolution: Densification and Shrinking Diameters. ACM Transactions on Knowledge Discovery from Data (ACM TKDD), 1(1), 2007.

[3] Traud, A. L.; Mucha, P. J. & Porter, M. A. Social structure of Facebook networks Physica A: Statistical Mechanics and its Applications, 2012, 391, 4165-4180

## Experiments

Test out the experimental setup on a small graph using `Test_Experimen_Setup.jl`. 

Run the full experiment on the email graph with `Full_Email_Experiment.jl`.

For Caltech, we start by solving the LP relaxation of LambdaCC for a range of values of lambda.

`Caltech_LamCCLP_bounds.jl` shows how we obtained lower bounds.

`Caltech_Best_Lams.jl` shows how we computed the minimum of the parameter fitness function for each metadata attribute. This will save the results in `Caltech_Best_Lams.mat`.