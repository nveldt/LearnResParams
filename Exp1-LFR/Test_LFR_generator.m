
path_to_binary_networks = '../include';
n = 500;
[A, truth,EdgeList] = generate_LFR(n,15,50,.1,20,50,path_to_binary_networks);

A = sparse(A);

% Before sorting by ground truth communities
figure(1)
spy(A)

[~,p] = sort(truth);

% After sorting by ground truth communities
figure(2)
spy(A(p,p))