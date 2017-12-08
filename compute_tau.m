function tau = compute_tau(x, u, k, dynamics, disturbances, counter, V, u_disc, w_c, par)
% par = [n, m, L, N, K]
n = par(1); L = par(3); N = par(4); K = par(5);
%% compute next trajectory for each disturbance realization
x_fwd = zeros(n, K+1, L);

for j=1:L
    robot_action = get_robot_action(disturbances, j, K);
    x_fwd(:,:,j) = compute_trajectory(x, u, dynamics, robot_action, K);
end

%% compute cost vector
costs = zeros(L,1);
if (k == N-1) %terminal node
    for j = 1:L
        [costs(j),~] = compute_cost(x_fwd(:,:,j), dynamics, u, w_c, K);
    end    
else
    for j = 1:L
        %recursive call
        x_j = x_fwd(:,:,j);
        disturbances.last = j;
        counter_ = counter;
        if (j == 4)
            counter_ = counter_ - 1;
        elseif (j == 5)
            counter_ = counter_ + 1;
        end
        [~, cost_to_go, ~, ~] = compute_H(x_j(:,end), k+1, dynamics, disturbances, counter_, V, u_disc, w_c, par);
        %net cost
        costs(j) = compute_cost(x_j, dynamics, u, w_c, K) + cost_to_go;
    end
end

%% now evaluate risk-sensitive cost
% disp(costs);
tau = max(V*costs);

end


