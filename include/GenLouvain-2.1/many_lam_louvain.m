function [cBest,BestMod,Mods] = many_lam_louvain(A,k,lam,cinit)
% Run the gen-louvain algorithm over and over with different settings
% and take the best outcome

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
limit = 1000000;

% Run the deterministic version
[cBest,BestMod] = iterated_genlouvain(B,limit,0,0,0,cinit);

Mods = [];
for times = 1:k
    [c,mod] = iterated_genlouvain(B,limit,0,1,1,cinit);
    Mods = [Mods; mod];
    if mod > BestMod
        BestMod = mod;
        cBest = c;
    end
end


end