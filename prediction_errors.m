
%% Load results
load(['Data/Inference/results_N_',num2str(N),'_',name_str,'.mat'], 'infres', 'par', 'disturbances', 'dynamics', 'beta');
P_w = infres.P_w; V_w = P_w.V;
w_c = infres.wc_f;
n = par(1); L = par(3); K = par(5); 
parameterized_polytope;

%% compute boltzmann weights of each action at beginning of each trajectory
B_weights = zeros(length(u_disc{1}), numExps);
parfor i=1:numExps
    x_i = x_expert(:,(i-1)*K+1);
    [~, tau_i, ~, ~, ~, ~] = est_H(x_i, 0, M, dynamics, disturbances, dist_last(i),...
        counters(i), P_w, V_w, u_disc, w_c, beta, 0, par, LP_sol, LP_prob);
    tau_s = beta*tau_i;
    B_weights(:,i) = exp(-beta*tau_s)./sum(exp(-beta*tau_s));
end

%% compute prediction errors
ac_space = u_disc{1};
err_x = zeros(1,numExps);
err_x_inf = zeros(1,numExps);
err_soft_x = zeros(1,numExps);
err_y = zeros(1,numExps);
err_y_inf = zeros(1,numExps);
err_soft_y = zeros(1,numExps);
x_pred = zeros(1, K+1, numExps);
v_pred = zeros(1, K+1, numExps);
err_vx = zeros(1,numExps);
err_vy = zeros(1,numExps);
theta_pred = zeros(1, K+1, numExps);
y_pred = zeros(1, K+1, numExps);
for i=1:numExps
    x_i = x_expert(1:5,(i-1)*K+1:i*K+1); % human states only over i-th trajectory
    u_i = u_expert(1:2, (i-1)*K+1:i*K+1);
    x_init = x_i(:,1); % state of human at beginning of i-th trajectory
    B_i = B_weights(:,i); % probabilities of playing each action at that state
    nb_ac = length(ac_space);
    %disp(['true trajectory'])
    %disp(x_i)
    for j=1:nb_ac
        u_j = ac_space{j};
        %disp(['action'])
        %disp(u_j)
        x_i_pred = zeros(5, K+1);
        x_i_h = zeros(5,K+1);
        x_i_pred(:,1) = x_init;
        x_i_h(:,1) = x_init;
        for k=1:K
            x_i_pred(:,k+1) = dynamics.f(x_i_pred(:,k)) + dynamics.g(x_i_pred(:,k))*u_j(:,k);
            x_i_h(:,k+1) = dynamics.f(x_i_h(:,k)) + dynamics.g(x_i_h(:,k))*u_i(:,k);
        end
        %disp(['predicted trajectory'])
        %disp(x_i_pred)
        %disp(['imitated trajectory'])
        %disp(x_i_h)
        err_x(i) = err_x(i) + B_i(j)*norm(x_i_pred(1,:)-x_i(1,:));
        err_x_inf(i) = err_x_inf(i) + B_i(j)*max(abs((x_i_pred(1,:)-x_i(1,:))));
        err_y(i) = err_y(i) + B_i(j)*norm(x_i_pred(2,:)-x_i(2,:));
        err_y_inf(i) = err_y_inf(i) + B_i(j)*max(abs((x_i_pred(2,:)-x_i(2,:))));
        err_vx(i) = err_vx(i) + B_i(j)*norm(x_i_pred(4,:).*cos(x_i_pred(3,:))-x_i(4,:).*cos(x_i(3,:)));
        err_vy(i) = err_vy(i) + B_i(j)*norm(x_i_pred(4,:).*sin(x_i_pred(3,:))-x_i(4,:).*sin(x_i(3,:)));
    end
    [~, idx_i] = max(B_i);
    x_soft_i = zeros(5,K+1);
    x_soft_i(:,1) = x_init;
    u_soft_i = ac_space{idx_i};
    %disp('most plausible action')
    %disp(u_soft_i)
    for k=1:K
        x_soft_i(:,k+1) = dynamics.f(x_soft_i(:,k)) + dynamics.g(x_soft_i(:,k))*u_soft_i(:,k);
    end
    x_pred(1,:,i) = x_soft_i(1,:);
    y_pred(1,:,i) = x_soft_i(2,:);
    v_pred(1,:,i) = x_soft_i(4,:);
    theta_pred(1,:,i) = x_soft_i(3,:);
    err_soft_x(i) = norm(x_soft_i(1,:) - x_i(1,:));
    err_soft_y(i) = norm(x_soft_i(2,:) - x_i(2,:));
end

save(['Data/Inference/errs_N_',num2str(N),'_', name_str,'.mat'], 'B_weights', 'K','numExps', 'err_soft_x', 'err_soft_y', 'err_x', 'err_y', 'err_vx', 'err_vy', 'x_pred', 'y_pred', 'v_pred', 'theta_pred');











