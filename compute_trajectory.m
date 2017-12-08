function state_fwd = compute_trajectory(curr_state, u, dynamics, robot_action, K)

n = length(curr_state);
state_fwd = zeros(n, K+1);
state_fwd(:,1) = curr_state;

for k=1:K
    u_k = u(:,k);
    w_k = robot_action(:,k);
    state_fwd(6:10,k+1) = dynamics.A*state_fwd(6:10,k) + dynamics.B*w_k;
    state_fwd(1:5,k+1) = dynamics.f(state_fwd(1:5,k)) + dynamics.g(state_fwd(1:5,k))*u_k;
end
end

