ss = 6; % number of disturbances fully discarded by pruning
%% Load data
load(['Data/human_robot_data_', name_str,'_train.mat']);
load(['Data/robot_data_10Hz_', name_str,'_train.mat']);

% prune and subsample (x,y)-position of robot
x_enforced = state_robot(1,:);
y_enforced = state_robot(3,:);

x_l = x_l + 140.;
y_l = -(y_l - 3.);
x_f = x_f + 140;
y_f = -(y_f - 3.);
vy_f = -vy_f;
ay_f = -ay_f;

tt = length(x_l);
tt = tt - mod(tt, 6);
nn = tt/6 - 5*K;
nn = nn - mod(nn, K);

prune = 1+ss*K + overlap;
x_l_d = zeros(1, nn-prune+1);
y_l_d = zeros(1, nn-prune+1);

for i=1:nn-prune+1
    x_l_d(1,i) = x_l( (i+prune-1-k_delay-1)*6+1) ;
    y_l_d(1,i) = y_l( (i+prune-1-k_delay-1)*6+1) ;
end

if (check_offset)
    plot(prune:nn, y_enforced(prune:nn), 'r');
    hold on;
    plot(prune:nn, y_l_d, 'b');
    save(['Data/Inference/Means/k_delay_train_',name_str,'.mat'], 'k_delay');
    keyboard;
end


%% Prune and subsample actions of human
% define parameters of model: time step and length l
l = 3.476;
dt = dt(2);
% define theta, delta, theta_dot, delta_dot, u_a, u_s
nt = length(x_f);
theta_f = atan(vy_f./vx_f);
v_f = sqrt(vx_f.^2+vy_f.^2);
u_a = ax_f.*cos(theta_f) + ay_f.*sin(theta_f);
theta_f_dot = (ay_f.*cos(theta_f) - ax_f.*sin(theta_f) ) ./ v_f;
delta_f = atan(-l*theta_f_dot ./ v_f);
theta_f_ddot = zeros(1,nt);
theta_f_ddot(1:nt-1) = (theta_f_dot(1,2:nt) - theta_f_dot(1,1:(nt-1))) ./ dt;
theta_f_ddot(1,nt) = theta_f_ddot(nt-1);
delta_f_dot = (- l./(1+ (l^2*theta_f_dot.^2 ./ v_f.^2)) ).*( (theta_f_ddot.*v_f - theta_f_dot.*u_a)./(v_f.^2) ); 
%delta_f_dot = -(l./(v_f.*(sec(delta_f)).^2)).*(theta_f_ddot + (u_a/l).*tan(delta_f));
u_s = delta_f_dot;
% subsample actions from 60Hz to 10Hz
u_s_d = zeros(1, nn-prune+1);
u_a_d = zeros(1, nn-prune+1);
for i=1:nn-prune+1
    u_s_d(1,i) = u_s( (i+prune-1-k_delay-1)*6+1) ;
    u_a_d(1,i) = u_a( (i+prune-1-k_delay-1)*6+1) ;
end

% subsample follower's state from 60Hz to 10Hz
for i=1:nn-prune+1
    x_f_d(1,i) = x_f((i+prune-1-k_delay-1)*6+1);
    y_f_d(1,i) = y_f((i+prune-1-k_delay-1)*6+1);
    theta_f_d(1,i) = theta_f((i+prune-1-k_delay-1)*6+1);
    v_f_d(1,i) = v_f((i+prune-1-k_delay-1)*6+1);
    delta_f_d(1,i) = delta_f((i+prune-1-k_delay-1)*6+1);
end

%% Store subsampled data in required format for inference
state_robot_p = state_robot(:,prune:nn);
x_expert = [x_f_d; y_f_d; theta_f_d; v_f_d; delta_f_d; state_robot_p];
u_expert = [u_a_d; u_s_d];
numExps = (nn - ss*K)/K-5; % total number of experiments
L = 4; % number of disturbances
n = size(x_expert,1);
m = size(u_expert,1);
par = [n, m, L, N, K];
% define disturbances
disturbances.overlap = overlap;
disturbances.list{1} = zeros(m,K);

accelerate = zeros(1,K);
for k=1:K
    accelerate(1,k) = k*(K+1-k);
end
accelerate = 4*(1/K^2)*accelerate;
acc_factor = 3;
dec_factor = 5.5;
disturbances.list{2} = acc_factor*[accelerate; zeros(1,K)];
disturbances.list{3} = -dec_factor*[accelerate; zeros(1,K)];

turn_right = [44.117 25.21 9.211 -3.878 -14.05 -21.33 -25.69 -27.14 -25.69 -21.33 -14.05 -3.878 9.211 25.21 44.117];
disturbances.list{4} = [zeros(1,K); -turn_right]; % turn left 
disturbances.list{5} = [zeros(1,K); turn_right]; % turn right
disturbances.last = js_noise(ss+1) +1 ;
dist_last = js_noise(ss+1:ss+numExps) + 1;
for i=1:length(dist_last)
    if dist_last(i) == 5
        dist_last(i) = 4;
    end
end
counters = counters(ss+2:ss+numExps+1);








    





