function [cBest,BestLCC,LCCS] = many_lam_louvain(A,k,lam,cinit)
% Run the gen-louvain algorithm over and over with different settings
% and take the best outcome
%addpath('GenLouvain-2.1/')

n = size(A,1);
m = nnz(A)/2;
if nargin < 4
    cinit = (1:n)';
end
if nargin < 3
    lam = 1/(2*m);
end

w = sum(A,1);
B = @(i) A(:,i) - lam*w'*w(i);
%B = full(A - w'*w*lam);
limit = 100000;

% Run the deterministic version
cBest = iterated_genlouvain(B,limit,0,0,0,cinit);
BestLCC = lamCCobj(A,lam,cBest);

LCCS = [];
for times = 1:k
    c = iterated_genlouvain(B,limit,0,1,1,cinit);
    lcc = lamCCobj(A,lam,c);
    LCCS = [LCCS; lcc];
    if lcc < BestLCC
        BestLCC = lcc;
        cBest = c;
    end
end

BestLCC = BestLCC*lam;
LCCS = LCCS*lam;


end