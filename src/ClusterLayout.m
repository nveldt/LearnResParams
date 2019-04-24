% Create a simple ring of circles graph visualization based on a clustering
function xy = ClusterLayout(c)


n = numel(c);

% Get a circle for the cluster centers
numclus = max(c);

centers = zeros(numclus,2);
for i = 1:numclus
    centers(i,1) = sin(2*pi/numclus*(i-1));
    centers(i,2) = cos(2*pi/numclus*(i-1));
end

% This will control how large the radius is allowed to be for the smaller
% circles
dist = norm(centers(1,:) - centers(2,:));
%scatter(centers(:,1),centers(:,2))

xy = zeros(n,2);
r = dist*.35;
for cluster = 1:numclus
    clusInds = find(c == cluster);
    nc = numel(clusInds);
    xy_clus = zeros(nc,2);
    
    for k = 1:nc
       xy_clus(k,1) = centers(cluster,1) + r*sin(2*pi/nc*(k-1));
       xy_clus(k,2) = centers(cluster,2) + r*cos(2*pi/nc*(k-1));  
    end
    xy(clusInds,:) = xy_clus;
    
end  




end