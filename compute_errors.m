%% load data
clear; clc;close all;

name_str = 'p2'; 

K = 15; overlap = 8;

%%%%%%%%%%%%%% Parameters for processing and learning %%%%%%%%%%%%%%%%%
N = 1;
k_delay = 3; % offset (tune it manually)
check_offset = false;
if (~check_offset)
    load(['Data/Inference/Means/k_delay_test_',name_str,'.mat'], 'k_delay');
end
process_test_data;

load(['Data/Inference/Means/means_', num2str(N),'_',name_str,'.mat']);

%% Performance of risk-sensitive model
prediction_errors;

%% Performance of risk-neutral model
prediction_errors_neutral;

%% Process and save results
tt = 1:520; %% total time
x_rel_actual = x_expert(6,tt) - x_expert(1,tt);
x_rel = zeros(1,K+1,numExps);
x_rel_n = zeros(1,K+1,numExps);
y_rel = zeros(1,K+1,numExps);
y_rel_n = zeros(1,K+1,numExps);
vx_rel = zeros(1,K+1,numExps);
vx_rel_n = zeros(1,K+1,numExps);
for i=1:numExps
    t_i = (i-1)*K+1:i*K+1;
    x_rel(1,:,i) = x_expert(6,t_i) - x_pred(1,:,i);
    x_rel_n(1,:,i) = x_expert(6,t_i) - x_pred_n(1,:,i);
    y_rel(1,:,i) = x_expert(8,t_i) - y_pred(1,:,i);
    y_rel_n(1,:,i) = x_expert(8,t_i) - y_pred_n(1,:,i);
    vx_rel(1,:,i) = x_expert(7,t_i) - v_pred(1,:,i).*cos(theta_pred(1,:,i));
    vx_rel_n(1,:,i) = x_expert(7,t_i) - v_pred_n(1,:,i).*cos(theta_pred_n(1,:,i));
end

save(['Data/Inference/Performance/perf_N_',num2str(N),'_',name_str,'.mat']) ;




