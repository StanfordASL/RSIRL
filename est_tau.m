function [tau, D_tau_w, D_tau_wc] = est_tau(x, u, k, M, dynamics, disturbances,dist_last, counter, P_w, V_w, u_disc, w_c, beta, grad_compute, par, LP_sol, LP_prob)
%% parameters 
n = par(1); L = par(3); N = par(4); K = par(5);
%% compute next trajectory for each disturbance realization
x_fwd = zeros(n, K+1, L);

disturbances.last = dist_last;
%j_repeat = dist_last;
for j=1:L
    robot_action = get_robot_action(disturbances, j, K, counter);
    x_fwd(:,:,j) = compute_trajectory(x, u, dynamics, robot_action, K);
end
%x_fwd(:,:,L) = x_fwd(:,:,j_repeat);

%% compute cost vector

%last element is the "repeated cost"
costs = zeros(L,1);
features = zeros(length(w_c),L);

is_term = 0;
if (k == N-1) %terminal node
    is_term = 1;
    for j = 1:L
        [costs(j), features(:,j)] = compute_cost(x_fwd(:,:,j), dynamics, u, w_c, K);
    end    
    %costs(L) = costs(j_repeat);
    %features(:,L) = features(:,j_repeat);
else
    if (grad_compute)
        D_cost_w = zeros(M,L);
        D_cost_wc = zeros(length(w_c),L);
    end
    
    for j = 1:L
        counter_ = counter;
        if (j == 4) % if robot turns
            counter_ = -counter; % change sign of counter
        end
        %recursive call
        x_j = x_fwd(:,:,j);
        disturbances.last = j;
        [cost_to_go, ~, ~, D_cost_w(:,j),~,D_cost_wc(:,j)] = est_H(x_j(:,end), k+1, M, dynamics, disturbances, j,...
                    counter_, P_w, V_w, u_disc, w_c, beta,grad_compute, par, LP_sol, LP_prob);
        %net cost                            
        [costs(j), features(:,j)] = compute_cost(x_j, dynamics, u, w_c, K);
        costs(j) = costs(j) + cost_to_go;
    end   
    %{
    costs(L) = costs(j_repeat);
    features(:,L) = features(:,j_repeat);
    D_cost_w(:,L) = D_cost_w(:,j_repeat);
    D_cost_wc(:,L) = D_cost_wc(:,j_repeat);    
    %}
end

%% now evaluate risk-sensitive cost
b_w = P_w.b;
[tau,idx] = max(V_w*costs);

%% compute gradients 
D_tau_w = zeros(M,1);
D_tau_wc = zeros(length(w_c),1);

eps = 1e-3;
I_M = eye(M);

if (grad_compute) 
    %derivative due to change in constraints
    [tau_pert,~,~] = compute_LP_yalmip(LP_prob,costs,kron(ones(1,M),b_w)+eps*I_M,M,L);
    %[tau_pert, ~, ~] = compute_LP_tomlab(LP_prob,costs,kron(ones(1,M),b_w)+eps*I_M,M,L,LP_sol);
    D_tau_w = (tau_pert - tau)/(-eps);
    
    q_star = V_w(idx,:)';

    D_tau_wc = features*q_star;
    
    %derivative due to change in continuation cost
    if (~is_term)
        D_tau_w = D_tau_w + D_cost_w*q_star;
        D_tau_wc = D_tau_wc + D_cost_wc*q_star;
    end
end

end
