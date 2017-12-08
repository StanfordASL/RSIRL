function [means, sigmas] = segment_actions(u_expert, K, nb_clusters)
% takes u_expert, a set of demonstrations of dimensions 2 \times K
% put each action in dimension 1 \times 2*K
% then, sliding window with increment equal to 1 to get all possible
% actions
% run K-Means over that set of actions to extract nb_clusters centroids
% also compute empirical covariance matrix of each gaussian distribution
% (useless if we do not use GMM)

d = 2;
numExps = length(u_expert);

U = zeros(d, numExps*K);
for i=1:numExps
    u_i = u_expert{i};
    U(:, (i-1)*K+1:i*K) = [u_i(1,:); u_i(2,:)];
end

actions = zeros(numExps*K-K+1, 2*K);
for i=1:(numExps*K-K+1)
    actions(i,:) = [U(1,i:i+K-1), U(2,i:i+K-1)];
end

% get mean actions
[idx, means] = kmeans(actions, nb_clusters);

% compute empirical covariance matrices
size_clusters = zeros(nb_clusters, 1);
sigmas = zeros(2*K, 2*K, nb_clusters);
%{
for i=1:size(actions,1)
    size_clusters(idx(i)) = size_clusters(idx(i)) + 1;
    sigmas(:,:,idx(i)) = sigmas(:,:,idx(i)) + (actions(i,:)-means(idx(i)))'*(actions(i,:)-means(idx(i)));
end

for i=1:nb_clusters
    sigmas(:,:,i) = sigmas(:,:,i) / size_clusters(i);
end
%}

end









