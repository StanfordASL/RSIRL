function [u, cost, tau, idx] = compute_H(x, k, dynamics, disturbances, counter, V, u_disc, w_c, par)

%% get control space depending on the value of k
action_space = u_disc{k+1};

%% for each control input, compute the value function
n_A = length(action_space);
tau = zeros(n_A, 1);

for j = 1:n_A
    u = action_space{j};
    tau(j) = compute_tau(x, u, k, dynamics, disturbances, counter, V, u_disc, w_c, par);
end

%% get best control input
[val,idx] = min(tau);

%% store values
u = action_space{idx};
cost = val;

end

