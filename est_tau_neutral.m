function [tau, D_tau_wc] = est_tau_neutral(x, u, k, dynamics, disturbances,dist_last, counter, p_true, u_disc, w_c, beta, grad_compute, par)
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
costs = zeros(L,1);
features = zeros(length(w_c),L);

is_term = 0;
if (k == N-1) %terminal node
    is_term = 1;
    for j = 1:L
        %disp(['relative distance'])
        %disp(x_fwd(6,:,j)-x_fwd(1,:,j));
        [costs(j), features(:,j)] = compute_cost(x_fwd(:,:,j), dynamics, u, w_c, K);
        %disp(['features'])
        %disp(features(:,j))
    end
else
    if (grad_compute)
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
        [cost_to_go, ~, ~,D_cost_wc(:,j)] = est_H_neutral(x_j(:,end), k+1, dynamics, disturbances, j,...
                        counter_, p_true, u_disc, w_c, beta, grad_compute, par);
                
        %net cost                            
        [costs(j), features(:,j)] = compute_cost(x_j, dynamics, u, w_c, K);
        costs(j) = costs(j) + cost_to_go;
    end                      
end

%% now evaluate risk-sensitive cost

tau = p_true'*costs;

%% compute gradients 

D_tau_wc = zeros(length(w_c),1);

if (grad_compute)
    D_tau_wc = features*p_true;
    
    %derivative due to change in continuation cost
    if (~is_term)
        D_tau_wc = D_tau_wc + D_cost_wc*p_true;
    end
end

end
