
%% Load results
load(['Data/Inference/results_neutral_N_',num2str(N),'_',name_str,'.mat'], 'p_true', 'infres', 'par', 'disturbances', 'dynamics', 'beta');
w_c = infres.wc_f;
n = par(1); L = par(3); K = par(5);
par(4) = N;
%% compute boltzmann weights of each action at beginning of each trajectory
B_weights_n = zeros(length(u_disc{1}), numExps);
parfor i=1:numExps
    x_i = x_expert(:,(i-1)*K+1);
    [~, tau_i, ~, ~] = est_H_neutral(x_i, 0, dynamics, disturbances, dist_last(i),...
        counters(i), p_true, u_disc, w_c, beta, 0, par);
    B_weights_n(:,i) = exp(-beta*tau_i)./sum(exp(-beta*tau_i));
end

%% compute prediction errors
ac_space = u_disc{1};
err_x_n = zeros(1,numExps);
err_x_inf_n = zeros(1,numExps);
err_soft_x_n = zeros(1,numExps);
err_y_n = zeros(1,numExps);
err_y_inf_n = zeros(1,numExps);
err_soft_y_n = zeros(1,numExps);
x_pred_n = zeros(1,K+1,numExps);
v_pred_n = zeros(1,K+1,numExps);
err_vx_n = zeros(1,numExps);
err_vy_n = zeros(1,numExps);
y_pred_n = zeros(1,K+1,numExps);
theta_pred_n = zeros(1,K+1,numExps);
for i=1:numExps
    x_i = x_expert(1:5,(i-1)*K+1:i*K+1); % human states only over i-th trajectory 
    x_init = x_i(:,1); % state of human at beginning of i-th trajectory
    B_i = B_weights_n(:,i); % probabilities of playing each action at that state
    nb_ac = length(ac_space);
    for j=1:nb_ac
        u_j = ac_space{j};
        x_i_pred = zeros(5, K+1);
        x_i_pred(:,1) = x_init;
        for k=1:K
            x_i_pred(:,k+1) = dynamics.f(x_i_pred(:,k)) + dynamics.g(x_i_pred(:,k))*u_j(:,k);
        end
        err_x_n(i) = err_x_n(i) + B_i(j)*norm(x_i_pred(1,:)-x_i(1,:));
        err_x_inf_n(i) = err_x_inf_n(i) + B_i(j)*max(abs((x_i_pred(1,:)-x_i(1,:))));
        err_y_n(i) = err_y_n(i) + B_i(j)*norm(x_i_pred(2,:)-x_i(2,:));
        err_y_inf_n(i) = err_y_inf_n(i) + B_i(j)*max(abs((x_i_pred(2,:)-x_i(2,:))));
        err_vx_n(i) = err_vx_n(i) + B_i(j)*norm(x_i_pred(4,:).*cos(x_i_pred(3,:))-x_i(4,:).*cos(x_i(3,:)));
        err_vy_n(i) = err_vy_n(i) + B_i(j)*norm(x_i_pred(4,:).*sin(x_i_pred(3,:))-x_i(4,:).*sin(x_i(3,:)));
    end
    [~, idx_i] = max(B_i);
    x_soft_i = zeros(5,K+1);
    x_soft_i(:,1) = x_init;
    u_soft_i = ac_space{idx_i};
    for k=1:K
        x_soft_i(:,k+1) = dynamics.f(x_soft_i(:,k)) + dynamics.g(x_soft_i(:,k))*u_soft_i(:,k);
    end
    x_pred_n(1,:,i) = x_soft_i(1,:);
    y_pred_n(1,:,i) = x_soft_i(2,:);
    v_pred_n(1,:,i) = x_soft_i(4,:);
    theta_pred_n(1,:,i) = x_soft_i(3,:);
    err_soft_x_n(i) = norm(x_soft_i(1,:) - x_i(1,:));
    err_soft_y_n(i) = norm(x_soft_i(2,:) - x_i(2,:));
end
numExps_neutral = numExps;

save(['Data/Inference/errs_neutral_N_',num2str(N),'_', name_str,'.mat'], 'B_weights_n','K','numExps_neutral', 'err_soft_x_n', 'err_soft_y_n', 'err_x_n', 'err_y_n', 'err_vx_n', 'err_vy_n', 'x_pred_n', 'y_pred_n', 'theta_pred_n', 'v_pred_n');






