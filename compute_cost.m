function [c, features] = compute_cost(state_fwd, dynamics, u, w_c, K)

%compute features
features = compute_features(state_fwd, dynamics, u, K);

%compute sum
c = w_c'*features;

end