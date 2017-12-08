%load results_bis_N_1.mat ;
% tt = 1:520;

load(['Data/Inference/results_N_',num2str(N),'_',name_str,'_','s','.mat'], 'u_disc');

%% Plot x relative over total time
cmap_data = jet;
indx_exp = 35;
x_init = x_expert(1:5, 1+(indx_exp-1)*K);
t_k = (indx_exp-1)*K+1:indx_exp*K+1;
[weights, index] = sort(B_weights(:,indx_exp), 'descend');
weights_sum = cumsum(weights);
index_stop = min(find(weights_sum > 0.95));
figure()
subplot(2,1,1)
for j=1:index_stop
    x_pred_j = zeros(5,K+1);
    x_pred_j(:,1) = x_init;
    B_j = B_weights(index(j),indx_exp);
    u_j = u_disc{1}{index(j)};
    for k=1:K
        x_pred_j(:,k+1) = dynamics.f(x_pred_j(:,k)) + dynamics.g(x_pred_j(:,k))*u_j(:,k);
    end
    
%     ind_j = ceil((64/exp(1))*exp(B_j));
%     ind_j = ceil(64*B_j);
    plot(1:16, x_expert(6,t_k)-x_pred_j(1,:), 'Color', cmap_data(ind_j,:), 'LineWidth', 3);
    hold on;
end
plot(1:16, x_expert(6,t_k)-x_expert(1,t_k), 'r', 'LineWidth', 1);

%%
load(['Data/Inference/results_neutral_N_',num2str(N),'_',name_str,'_','s','.mat'], 'u_disc');

[weights, index] = sort(B_weights_n(:,indx_exp), 'descend');
weights_sum = cumsum(weights);
index_stop = min(find(weights_sum > 0.95));
subplot(2,1,2)
for j=1:index_stop
    x_pred_j = zeros(5,K+1);
    x_pred_j(:,1) = x_init;
    B_j = B_weights_n(index(j),indx_exp);
    u_j = u_disc{1}{index(j)};
    for k=1:K
        x_pred_j(:,k+1) = dynamics.f(x_pred_j(:,k)) + dynamics.g(x_pred_j(:,k))*u_j(:,k);
    end
    ind_j = ceil((64/exp(1))*exp(B_j));
%     ind_j = ceil(64*B_j);
    plot(1:16, x_expert(6,t_k)-x_pred_j(1,:), 'Color', cmap_data(ind_j,:), 'LineWidth', 3);
    hold on;
end
plot(1:16, x_expert(6,t_k)-x_expert(1,t_k), 'r', 'LineWidth', 1);














