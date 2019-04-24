%% This code should build a graph from an MRI image
% Developed by David Gleich with Nate Veldt 2016-01-22
%
% Includes slight updates made by Nate Veldt to apply to more types of MRIs
% than just the brain MRI scans from the MICCAI 2012 challenge 2018-11-06
function [ei,ej,evi,evd,N,Nmap] = mri2graph(X,r)
% X is the MRI image
% n1-by-n2-by-n3
% N will be the node id of each voxel
% Nmap is a map from node ids to voxels, stored as columns
% Nmap(:,nodeid) is the associated voxel
% r is the nearest neighbor distance

dims = size(X);
nnodes = prod(dims);
edges_per_node = (2*r+1)^3 - 1; % we are going manhattan distance k
nedges = nnodes*edges_per_node; % this is an upper-bound, will only use some
ei = ones(nedges,1);
ej = ones(nedges,1);
evi = zeros(nedges,1); 
evd = zeros(nedges,1); 
% note that we always call sparse(ei,ej,ev) and any unused edges will
% drop out because they have value 0.

% build the node ids, this will 
N = zeros(dims);
N(:) = 1:nnodes; % this will assign each node a linear index
Nmap = zeros(3,nnodes); % this will map back, will build below.

% Find the maximum value in X, used for normalization
M = double(max(max(max(X))));

fprintf('graph_from_mri: building edge list ...\n');

% let's build the node map too
nzi = 1;

t0 = tic;
tt = tic;
nextedgecheck = 50000;
nexttime = 5;

for i=1:dims(1)
    for j=1:dims(2)
        for k=1:dims(3) % oof, that's a lotta loops!

            % Get our node id
            src = N(i,j,k);
            Nmap(:,src) = [i,j,k];
            
            for oi=-r:r
                for oj=-r:r
                    for ok=-r:r % even more loops!
                        
                        u = i+oi;
                        v = j+oj;
                        w = k+ok;
                        
                        if u <= 0 || u > dims(1)
                            continue
                        end
                        
                        if v <= 0 || v > dims(2)
                            continue
                        end
                        
                        if w <= 0 || w > dims(3)
                            continue
                        end
                        
                        % Insert an edge!
                        dst = N(u,v,w);
                        
                        ei(nzi) = src;
                        ej(nzi) = dst;
                        
                        xi = min(sqrt(double(X(i,j,k)))/sqrt(M),1);
                        xj = min(sqrt(double(X(u,v,w)))/sqrt(M),1);
                        dist = norm([oi,oj,ok])^2;
                        
                        evi(nzi) = xi-xj;
                        evd(nzi) = dist;
                        %ev(nzi) = exp(-double((X(i,j,k) - X(u,v,w)))^2/sigma1)*exp(-norm([u,v,w])^2/sigma2);
                        %ev(nzi) = exp(-(log(1.+double(X(i,j,k))) - log(1.+double(X(u,v,w))))^2/sigma1^2)*exp(-norm([oi,oj,ok])^2/sigma2^2);
                        
                        % next index!
                        nzi = nzi+1;
                        
                        
                    end
                end
            end
         
            if nzi >= nextedgecheck
                if toc(tt) >= nexttime
                    % update nexttime and print status message
                    nexttime = min(nexttime*2,3600);
                    tt = tic;
                    
                    edges_per_sec = nzi/toc(t0); % edges per second
                    est_time_completion = (nedges-nzi)*1/edges_per_sec;
                    
                    fprintf('graph_from_mri: %10i / %10i edges done (%5.1f%%) %8.1f sec (%8.1f sec left)\n', ...
                        nzi, nedges, 100*nzi/nedges, toc(t0), est_time_completion);
                end
                % always check next 50k edges
                nextedgecheck = nextedgecheck + 50000;
            end

        end
    end
end

