function [actions, opt_action] = map_actions(u_expert, means)
% maps to an action from continuous space to the closest mean action, using
% Euclidean distance

nb_actions = length(u_expert);
actions = {};
opt_action = zeros(nb_actions, 1);

K = size(u_expert{1},2);
nb_means = size(means, 1);
for i=1:nb_actions
    u_i = u_expert{i};
    u_i = [u_i(1,:), u_i(2,:)];
    dist = zeros(nb_means, 1);
    for j=1:nb_means
        dist(j,1) = norm(u_i - means(j,:));
    end
    [~, idx] = min(dist);
    opt_action(i,1) = idx;
    actions{i} = [means(idx,1:K); means(idx,K+1:2*K)];
end

end


