addpath('../src/')
addpath('../include/GenLouvain-2.1/')

fprintf('Small Com: \n')


for mu = 20:5:50
        
    load(strcat('Output/Learned_SmallComm_mu_',num2str(mu),'_Gurobi.mat'))    
    fprintf('mu = %d, Time to solve using Gurobi solver = %f \n',mu,time1)

end

