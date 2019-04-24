%% Generate 1 training LFR graphs and 5 test examples

n = 200;
k = 10;     % average degree
maxk = 20;  % maximum degree
minc = 5;  % minimum number of nodes in the community
maxc = 20;  % max community size

for mu = .2:.05:.5

    for j = 1
        [A, truth,EdgeList] = generate_LFR(n,k,maxk,mu,minc,maxc);
        A = sparse(A);
        save(strcat('data/Train_',num2str(j),'_SmallComm_','_mu_',num2str(mu*100),'.mat'),'k','maxk','mu','minc','maxc','A','truth')
    end

    % 5 testing examples
    for j = 1:5
        [A, truth,EdgeList] = generate_LFR(n,k,maxk,mu,minc,maxc);
        A = sparse(A);
        save(strcat('data/Test_',num2str(j),'_SmallComm_','_mu_',num2str(mu*100),'.mat'),'k','maxk','mu','minc','maxc','A','truth')
    end

end