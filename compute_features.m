function features = compute_features(state_fwd, dynamics, u, K)
% state_fwd = [x_f, y_f, theta_f, v_f, delta_f, x_l, vx_l, y_l, vy_l, ay_l]
% code for human non linear dynamics and double/triple integrator for robot

x_rel = state_fwd(6,2:end) - state_fwd(1,2:end);
y_rel = state_fwd(8,2:end) - state_fwd(2,2:end);
y_f = state_fwd(2,2:end);
v_rel = state_fwd(4,2:end).*cos(state_fwd(3,2:end)) - state_fwd(7,2:end);
u_a = u(1,:);
u_diff = u_a(2:end) - u_a(1:end-1);

r_1 = 0.05;
r_2 = 1.;
r_3 = 0.1;
r_4 = 1.0;
r_5 = 0.1;
r_6 = 0.5;
features = zeros(6,K);
features(1,:) =  (x_rel > 2.5).*(log(1+exp(r_1*(x_rel-2.5))) - log(2) ); % not too far from leader
features(2,:) = (x_rel < 2.5).*( log(1+exp(-r_2*(x_rel-2.5))) -log(2) ); % leader should remain first
% features(2,:) = 
features(3,:) = log(1+exp(r_3*sqrt(v_rel.^2) )) - log(2); % relative velocity should be low
features(4,:) = r_4 * abs([0,u_diff]); %penalize high differentials in inputs
features(5,:) = log(1+exp(r_5*abs(y_rel)))-log(2);
features(6,:) = (y_f > 2).*( log(1+exp(r_6*(y_f-2.))) - log(2) ) + (y_f < -2).*( log(1+exp(-r_6*(y_f+2.))) - log(2) );

features = sum(features,2);
end



