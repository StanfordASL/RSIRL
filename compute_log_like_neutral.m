function l_array = compute_log_like_neutral(Idx_do,p_true,w_c_0,dynamics, disturbances, data, u_disc, beta, par)
%% parameters
n = par(1); L = par(3); N = par(4); K = par(5);
x_expert = data.x;
u_expert = data.u;
a_expert = data.opt_action;
dist_last = data.dist_last;
counters = data.counter;
% numExps = length(u_expert);

%% Outliers

numExps = numel(Idx_do);

%% initialize
w_c = w_c_0;
l_array = zeros(numExps,1);

%compute gradient of log-likelihood
parfor j=1:numExps
    i = Idx_do(j);
    x_obs = x_expert(:,(i-1)*K+1);
    counter = counters(i); % indicates if robot on left or right lane
    a_obs = a_expert(i,1);
    
    [~, tau, ~, ~] = est_H_neutral(x_obs, 0, dynamics, disturbances, dist_last(i),...
        counter, p_true, u_disc, w_c, beta, 0, par);
    
    B_weights = exp(-beta*tau)./sum(exp(-beta*tau));
    
    l_array(j) = log( B_weights(a_obs) ); % log likelihood of i-th observation
end

end



