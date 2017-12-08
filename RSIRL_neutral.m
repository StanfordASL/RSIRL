%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Risk-sensitive Inverse Reinforcement Learning %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Choose participant to process data from
clear;clc;close all;

name_str = 'p2';

l = 3.476; dt = 0.1; K = 15;
overlap = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%   Parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 2;

nb_clusters_1 = 15;
nb_clusters_2 = 5;

beta = 3.0; % inverse temperature of Boltzmann
T = 10; % number of time steps in gradient ascent

check_offset = false;
load(['Data/Inference/Means/k_delay_train_',name_str,'.mat'], 'k_delay');
k_delay
process_training_data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load means; computed in RS-IRL

load(['Data/Inference/Means/means_', num2str(N),'_',name_str,'.mat']);

%% Load dynamics

define_dynamics; 

%% Create data structure for gradient ascent
data.x = x_expert;
data.counter = counters;
data.dist_last = dist_last;
u = {};
for i=1:numExps
    u{i} = u_expert(:,(i-1)*K+1:i*K);
end
[data.u, data.opt_action] = map_actions(u, means_1);

%% Parameters for gradient ascent

alpha = 0.7;
wc_initial = [0.028368186502538
   0.367823198560694
   0.297345595794308
   0.264149922027673
   0.025538473699593
   0.016774623415196];
% wc_initial = [1;5;2;3;1;2]; wc_initial = wc_initial/sum(wc_initial);

p_true = [0.3; 0.3; 0.3; 0.1];

%% Outliers

%Exps = [1,3:28,30:43];
Exps = 1:numExps;

%% Learn cost weights by maximizing log-likelihood using full-batch gradient ascent
delta_l = 1;
delta_l_eps = 0.01;
wc_f = wc_initial;

while (delta_l > delta_l_eps) 
    res = subgradientdescent_neutral(Exps,wc_initial, p_true, dynamics, disturbances, data, u_disc, beta, alpha, T, par);
    
    %net change (use for convergence test)
    delta_l = res.log_likelihood(end)-res.log_likelihood(1);
    fprintf('Delta(log_like): %f\n',delta_l);
    disp('**********************************************');
    
    if (delta_l > 0) %successful sweep
        %update record
        l_array = res.log_like_array(:,end);
        l_w = res.log_likelihood;
        wc_f = res.wc_s(:,end-1); % cost weights
        
        %update starting point
        wc_initial = wc_f;
    end
end

%% Save results
%keyboard;
infres.l_w = l_w;
infres.wc_f = wc_f;

par_ascent.Exps = Exps;
par_ascent.wc_initial = wc_initial;
par_ascent.T = T;
par_ascent.beta = beta; 
par_ascent.alpha = alpha;
save(['Data/Inference/results_neutral_N_',num2str(N),'_',name_str,'.mat'], 'p_true', 'beta', 'infres', 'disturbances', 'dynamics', 'data', 'x_expert', 'u_disc', 'par', 'par_ascent');

