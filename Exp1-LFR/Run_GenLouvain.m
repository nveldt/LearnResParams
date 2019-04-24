% Load the learned lambdas for the small community training sets,
% and then run louvain with standard and learned parameter settings

addpath('../src/')
addpath('../include/GenLouvain-2.1/')

for muval = 20:5:50
    
    fprintf('\nNew Mu Val\n')
    load(strcat('Output/Learned_SmallComm_mu_',num2str(muval),'_Gurobi.mat'))
    
    [~,loc1] = min(Plist1);
    lam1 = ResList1(loc1);

    % Run GenLouvain k times and take the best output with respect to the
    % LambdaCC objective function score
    k = 5;
    Modularity = zeros(5,2);
    Learned = zeros(5,2);
    for testcase = 1:5
        load(strcat('data/Test_',num2str(testcase),'_SmallComm_','_mu_',num2str(muval),'.mat'))
        fprintf('\nNew Test Case\n')
        c = truth;
        m = nnz(A)/2;
        lam = 1/(2*m);
        [c_mod,modscore] = many_lam_louvain(A,k,lam);
        fprintf('Modularity: %f %f %f \n',NMI(c_mod,c),ARI(c_mod,c),modscore)
        Modularity(testcase,1:2) = [NMI(c_mod,c),ARI(c_mod,c)];
        
        [c_mod,modscore] = many_lam_louvain(A,k,lam1);
        fprintf('Learn Lam1: %f %f %f \n',NMI(c_mod,c),ARI(c_mod,c),modscore)
        Learned(testcase,1:2) = [NMI(c_mod,c),ARI(c_mod,c)];
        
        
    end
    info = 'column 1 is NMI scores, column 2 is ARI scores';
    save(strcat('Output/Scores_ari_nmi_SmallComm_mu_',num2str(muval),'_Gurobi.mat'),'Modularity','Learned','info')
    
end