function [A, truth,EdgeList] = generate_LFR(n,k,maxk,mu,minc,maxc,path_to_binary_networks)
% This calls the LFR benchmark generator to generate a graph with n nodes,
% average degree k, max degree maxk, mixing parameter mu, smallest
% community of size minc, largest community of size maxc

if nargin < 7
    path_to_binary_networks = '../include';
end
command = strcat(path_to_binary_networks,'/binary_networks/./benchmark -N ',{' '},num2str(n),{' '},'-k',{' '}, ...
    num2str(k), {' '},'-maxk',{' '},num2str(maxk), {' '}, '-mu',{' '},num2str(mu),{' '},...
    '-minc',{' '},num2str(minc),{' '},'-maxc',{' '},num2str(maxc));

system(char(command))

[A,truth,EdgeList] = readinLFR('network.dat','community.dat');

!rm network.dat community.dat statistics.dat time_seed.dat
end

function [A,truth,EdgeList] = readinLFR(netfile,comfile)
% [A,truth] = readinLFR(filename): reads in an LFR network and outputs the
% adjacency matrix corresponding to that network and the ground truth
% community vector

net = fopen(netfile,'r');

com = fopen(comfile);
ComList = fscanf(com,'%d %d',[2, Inf])';
n = ComList(end,1);
truth = ComList(:,2);
fclose(com);

A = zeros(n);
net = fopen(netfile);
EdgeList = fscanf(net,'%d %d',[2, Inf])';
fclose(net);

for i = 1:size(EdgeList,1)
    A(EdgeList(i,1),EdgeList(i,2)) = 1;
end

end