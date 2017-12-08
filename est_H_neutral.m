function [cost, tau, D_tau_wc, D_cost_wc] = est_H_neutral(x, k, dynamics, disturbances, dist_last,...
                        counter, p_true, u_disc, w_c, beta, grad_compute, par)

%% 
action_space = u_disc{k+1};

n_A = length(action_space);
tau = zeros(n_A, 1 );
D_tau_wc = zeros(length(w_c),n_A);

for j = 1:n_A
    u = action_space{j};
    %disp(['action'])
    %disp(u)
    [tau(j), D_tau_wc(:,j)] = est_tau_neutral(x, u, k, dynamics, disturbances,dist_last, counter, p_true, u_disc, w_c, beta, grad_compute, par);
end 

%% apply softmax

tau_s = beta*tau;

sigma_H = exp(-tau_s)./sum(exp(-tau_s));
cost = sigma_H'*tau;

D_cost_wc = zeros(length(w_c),1);

if (grad_compute)
    D_cost_wc = D_tau_wc*sigma_H;
end

end
